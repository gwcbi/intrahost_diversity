#! /bin/bash

dest=$1
btargs=$2
r1=$3
r2=$4
[[ -n ${5+x} ]] && rU=$5 || unset rU

echo "Destination.............$dest"
echo "Bowtie2 args............$btargs"
echo "Read1...................$r1"
echo "Read2...................$r2"
[[ -n ${rU+x} ]] && echo "ReadU...................$rU"

rtmp=$(mktemp -d  -t tmp-align_and_fix-XXXXXX)
echo "Temporary dir...........$rtmp"

# Copy and index initial reference
cp $dest/consensus.fasta $rtmp/consensus.fasta
samtools faidx $rtmp/consensus.fasta
picard CreateSequenceDictionary R=$rtmp/consensus.fasta O=$rtmp/consensus.dict

module load bowtie2
bowtie2-build $rtmp/consensus.fasta $rtmp/consensus

# Set standard args
[[ -n ${samp+x} ]] && rgid=$samp || rgid=$(dirname $dest)
stdargs="--rg-id $rgid --rg SM:$rgid --rg LB:1 --rg PU:1 --rg PL:illumina"

# Set -U argument
[[ -n ${rU+x} ]] && uarg=" -U $rU " || uarg=""

echo "Aligning with bowtie2"
cmd="bowtie2 -p $NCPU $stdargs $btargs -x $rtmp/consensus -1 $r1 -2 $r2 $uarg -S $rtmp/a1.sam"
echo -e "Command:\n$cmd"
$cmd

samtools view -uS $rtmp/a1.sam | samtools sort - $rtmp/a1
samtools index $rtmp/a1.bam
picard MarkDuplicates REMOVE_DUPLICATES=false CREATE_INDEX=true \
    M=$rtmp/rmdup.metrics.txt \
    I=$rtmp/a1.bam \
    O=$rtmp/rmdup.bam

samtools index $rtmp/rmdup.bam

java -Xmx50g -jar $GATK -T RealignerTargetCreator \
    -I $rtmp/rmdup.bam \
    -R $rtmp/consensus.fasta \
    -o $rtmp/tmp.intervals

java -Xmx50g -jar $GATK -T IndelRealigner \
    -maxReads 1000000 -dt NONE \
    -I $rtmp/rmdup.bam \
    -R $rtmp/consensus.fasta \
    -targetIntervals $rtmp/tmp.intervals \
    -o $rtmp/all.bam

samtools index $rtmp/all.bam

# Copy BAM with all reads
cp $rtmp/all.bam $dest && samtools index $dest/all.bam

### Create filtered alignment ############################################################

# Exclude: UNMAP,MUNMAP,SECONDARY,QCFAIL,DUP,SUPPLEMENTARY
samtools view -b -F 3852 $rtmp/all.bam > $rtmp/filtered.bam
samtools index $rtmp/filtered.bam

### Call polymorphisms using freebayes ###################################################
module load freebayes
module load bcftools/1.3

fbopts=" --min-alternate-fraction 0.01 --pooled-continuous --standard-filters --ploidy 1 --haplotype-length 0 "

### Parallel freebayes
# module load bamtools
# samtools faidx $rtmp/consensus.fasta
# bamtools coverage -in $rtmp/filtered.bam | coverage_to_regions.py $rtmp/consensus.fasta.fai $NCPU > $rtmp/regions.bed
# freebayes-parallel $rtmp/regions.bed $NCPU $fbopts -f $rtmp/consensus.fasta $rtmp/filtered.bam > $rtmp/init_fb.vcf

echo "Calling variants with freebayes"
freebayes $fbopts -f $rtmp/consensus.fasta $rtmp/filtered.bam | \
    bcftools view -Oz | \
        bcftools filter -m '+' -Oz -e "QA > 4000" -s 'HQ' | \
        bcftools filter -m '+' -Oz -e "AO/DP > 0.50" -s 'GT50' | \
        bcftools filter -m '+' -Oz -e "AO/DP <= 0.20 & AO/DP >= 0.05" -s 'LT20' | \
        bcftools filter -m '+' -Oz -e "AO/DP < 0.05" -s 'LT05' > $rtmp/filtered_fb.vcf.gz

### Call polymorphisms using UG ##########################################################
echo "Calling variants with UnifiedGenotyper"
java -Xmx50g -jar $GATK -T UnifiedGenotyper \
    --num_threads $NCPU \
    -glm BOTH -dt NONE  \
    -ploidy 4 -maxAltAlleles 4 \
    -A AlleleBalance --min_base_quality_score 15 \
    -I $rtmp/filtered.bam \
    -R $rtmp/consensus.fasta \
    -o $rtmp/filtered_ug.vcf.gz

### Call polymorphisms using HC #####################################################
echo "Calling variants with HaplotypeCaller"
java -Xmx50g -jar $GATK -T HaplotypeCaller \
    -dt NONE  \
    -ploidy 4 -maxAltAlleles 4 \
    -A AlleleBalance --min_base_quality_score 15 \
    -I $rtmp/filtered.bam \
    -R $rtmp/consensus.fasta \
    -o $rtmp/filtered_hc.vcf.gz

# Copy filtered data
cp $rtmp/filtered.bam $dest && samtools index $dest/filtered.bam
cp $rtmp/filtered_fb.vcf.gz $dest && bcftools index $dest/filtered_fb.vcf.gz
cp $rtmp/filtered_ug.vcf.gz $dest && bcftools index $dest/filtered_ug.vcf.gz
cp $rtmp/filtered_hc.vcf.gz $dest && bcftools index $dest/filtered_hc.vcf.gz

if [[ -n ${rU+x} ]]; then
    # Exclude: PAIRED,UNMAP,MUNMAP,SECONDARY,QCFAIL,DUP,SUPPLEMENTARY
    samtools view -b -F 3853 $rtmp/filtered.bam > $rtmp/final.bam
    samtools index $rtmp/final.bam
else
    # Include: PAIRED Exclude: UNMAP,MUNMAP,SECONDARY,QCFAIL,DUP,SUPPLEMENTARY
    # samtools view -b -f 1 -F 3852 $rtmp/all.bam > $rtmp/paired.bam
    # Include:PAIRED,PROPER_PAIR Exclude: UNMAP,MUNMAP,SECONDARY,QCFAIL,DUP,SUPPLEMENTARY
    samtools view -b -f 3 -F 3852 $rtmp/filtered.bam > $rtmp/final.bam
    samtools index $rtmp/final.bam
fi

### Call final polymorphisms using freebayes #############################################
echo "Calling final variants with freebayes"
freebayes $fbopts -f $rtmp/consensus.fasta $rtmp/final.bam | \
    bcftools view -Oz | \
        bcftools filter -m '+' -Oz -e "QA > 4000" -s 'HQ' | \
        bcftools filter -m '+' -Oz -e "AO/DP > 0.50" -s 'GT50' | \
        bcftools filter -m '+' -Oz -e "AO/DP <= 0.20 & AO/DP >= 0.05" -s 'LT20' | \
        bcftools filter -m '+' -Oz -e "AO/DP < 0.05" -s 'LT05' > $rtmp/final_fb.vcf.gz

### Call final polymorphisms using UG #####################################################
echo "Calling final variants with UnifiedGenotyper"
java -Xmx50g -jar $GATK -T UnifiedGenotyper \
    --num_threads $NCPU \
    -glm BOTH -dt NONE  \
    -ploidy 4 -maxAltAlleles 4 \
    -A AlleleBalance --min_base_quality_score 15 \
    -I $rtmp/final.bam \
    -R $rtmp/consensus.fasta \
    -o $rtmp/final_ug.vcf.gz

### Call final polymorphisms using HC #####################################################
echo "Calling final variants with HaplotypeCaller"
java -Xmx50g -jar $GATK -T HaplotypeCaller \
    -dt NONE  \
    -ploidy 4 -maxAltAlleles 4 \
    -A AlleleBalance --min_base_quality_score 15 \
    -I $rtmp/final.bam \
    -R $rtmp/consensus.fasta \
    -o $rtmp/final_hc.vcf.gz

# Copy filtered data
cp $rtmp/final.bam $dest && samtools index $dest/final.bam
cp $rtmp/final_fb.vcf.gz $dest && bcftools index $dest/final_fb.vcf.gz
cp $rtmp/final_ug.vcf.gz $dest && bcftools index $dest/final_ug.vcf.gz
cp $rtmp/final_hc.vcf.gz $dest && bcftools index $dest/final_hc.vcf.gz