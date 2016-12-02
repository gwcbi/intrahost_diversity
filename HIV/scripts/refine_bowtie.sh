#! /bin/bash

dest=$1
btargs=$2
mincov=$3
r1=$4
r2=$5
[[ -n ${6+x} ]] && rU=$6 || unset rU

echo "Destination.............$dest"
echo "Bowtie2 args............$btargs"
echo "Min coverage............$mincov"
echo "Read1...................$r1"
echo "Read2...................$r2"
[[ -n ${rU+x} ]] && echo "ReadU...................$rU"

rtmp=$(mktemp -d  -t tmp-refine-XXXXXX)
echo "Temporary dir...........$rtmp"

# Copy and index initial reference
cp $dest/initial.fasta $rtmp/initial.fasta
samtools faidx $rtmp/initial.fasta
picard CreateSequenceDictionary R=$rtmp/initial.fasta O=$rtmp/initial.dict

module load bowtie2
bowtie2-build $rtmp/initial.fasta $rtmp/initial

# Set standard args
[[ -n ${samp+x} ]] && rgid=$samp || rgid=$(dirname $dest)
stdargs="--no-unal --rg-id $rgid --rg SM:$rgid --rg LB:1 --rg PU:1 --rg PL:illumina"

# Set -U argument
[[ -n ${rU+x} ]] && uarg=" -U $rU " || uarg=""

echo "Aligning with bowtie2"
cmd="bowtie2 -p $NCPU $stdargs $btargs -x $rtmp/initial -1 $r1 -2 $r2 $uarg -S $rtmp/aligned.sam"
echo -e "Command:\n$cmd"
$cmd

samtools view -uS $rtmp/aligned.sam | samtools sort - $rtmp/aligned
samtools index $rtmp/aligned.bam
picard MarkDuplicates REMOVE_DUPLICATES=true CREATE_INDEX=true \
    M=$rtmp/rmdup.metrics.txt \
    I=$rtmp/aligned.bam \
    O=$rtmp/rmdup.bam

java -Xmx50g -jar $GATK -T RealignerTargetCreator \
    -I $rtmp/rmdup.bam \
    -R $rtmp/initial.fasta \
    -o $rtmp/tmp.intervals

java -Xmx50g -jar $GATK -T IndelRealigner \
    -maxReads 1000000 -dt NONE \
    -I $rtmp/rmdup.bam \
    -R $rtmp/initial.fasta \
    -targetIntervals $rtmp/tmp.intervals \
    -o $rtmp/realign.bam

java -Xmx50g -jar $GATK -T UnifiedGenotyper \
    --num_threads $NCPU \
    -out_mode EMIT_ALL_SITES \
    -glm BOTH --baq OFF --useOriginalQualities -dt NONE  \
    -stand_call_conf 0 -stand_emit_conf 0 \
    -A AlleleBalance --min_base_quality_score 15 -ploidy 4 \
    -I $rtmp/realign.bam \
    -R $rtmp/initial.fasta \
    -o $rtmp/tmp.vcf.gz

mc=$(python ../scripts/calc_cov.py $rtmp/tmp.vcf.gz $mincov)
echo "Min coverage............$mc"

assembly.py vcf_to_fasta \
    --min_coverage $mc \
    $rtmp/tmp.vcf.gz \
    $rtmp/refined.fasta

cp $rtmp/tmp.vcf.gz $dest/initial.vcf.gz
cp $rtmp/refined.fasta $dest/refined.fasta
samtools faidx $dest/refined.fasta
picard CreateSequenceDictionary R=$dest/refined.fasta O=$dest/refined.dict

# Remove temporary directory
# rm -rf $rtmp
