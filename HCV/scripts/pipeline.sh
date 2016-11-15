#! /bin/bash

module load picard
module load flash
module load samtools

module load viral-ngs
module load picard
module load bless
module load blast+

# How many CPUs to use. If we are on a login node we do not want to use all CPUs.
[[ $(hostname) == login* ]] && NCPU=8 || NCPU=$(nproc)
# Temporary directory to use
[[ -d /scratch ]] && CURTMP="/scratch" || CURTMP="$TMPDIR"

rroot=$(git rev-parse --show-toplevel)

##########################################################################################
# Step 0: Create unaligned BAM file for raw reads
##########################################################################################
mkdir -p $samp/00_raw

picard FastqToSam ALLOW_AND_IGNORE_EMPTY_LINES=true \
    SM=$samp \
    F1=$samp/00_raw/original_1.fastq \
    F2=$samp/00_raw/original_2.fastq \
    O=$samp/00_raw/original.bam

##########################################################################################
# Step 1: Trim reads. Not implemented for now.
##########################################################################################
mkdir -p $samp/01_trimmed

##########################################################################################
# Step 2: Join reads with FLASh.
##########################################################################################
mkdir -p $samp/02_flash

flash -M 150 -d $samp/02_flash -o flash \
    $samp/00_raw/original_1.fastq \
    $samp/00_raw/original_2.fastq 2>&1 | tee $samp/02_flash/flash.log

# Make BAM for extended fragments
picard FastqToSam ALLOW_AND_IGNORE_EMPTY_LINES=true \
    SM=$samp \
    F1=$samp/02_flash/flash.extendedFrags.fastq \
    O=$samp/02_flash/flash.extendedFrags.bam

# Make BAM for reads that were not combined
picard FastqToSam ALLOW_AND_IGNORE_EMPTY_LINES=true \
    SM=$samp \
    F1=$samp/02_flash/flash.notCombined_1.fastq \
    F2=$samp/02_flash/flash.notCombined_2.fastq \
    O=$samp/02_flash/flash.notCombined.bam

# Merge the two BAM files
samtools merge $samp/02_flash/flash.bam $samp/02_flash/flash.extendedFrags.bam $samp/02_flash/flash.notCombined.bam

##########################################################################################
# Step 3: Error correction using BLESS2
##########################################################################################
mkdir -p $samp/03_bless

$(which bless)  -kmerlength 31 \
    -read1 $samp/00_raw/original_1.fastq \
    -read2 $samp/00_raw/original_2.fastq  \
    -prefix $samp/03_bless/bless

# Make BAM with corrected reads
picard FastqToSam ALLOW_AND_IGNORE_EMPTY_LINES=true \
    SM=$samp \
    F1=$samp/03_bless/bless.1.corrected.fastq \
    F2=$samp/03_bless/bless.2.corrected.fastq \
    O=$samp/03_bless/bless.bam

##########################################################################################
# Step 4: Denovo assembly using Trinity
##########################################################################################
mkdir -p $samp/04_assembly

ADAPTERFILE=$rroot/adapters/adapters.fasta

assembly.py assemble_trinity --threads $NCPU --JVMmemory 50g --tmp_dir $CURTMP \
    --n_reads 1000000 \
   $samp/00_raw/original.bam \
   $ADAPTERFILE \
   $samp/04_assembly/contigs.fa

##########################################################################################
# Step 5: Assign contigs to subtypes
##########################################################################################
mkdir -p $samp/05_subtype

python $rroot/HCV/scripts/assign_contigs.py \
    $samp/04_assembly/contigs.fa \
    $rroot/HCV/references/HCV_subtype_refs \
    $samp/05_subtype

mkdir -p $samp/05b_subtype
read_utils.py novoalign --JVMmemory 50g --tmp_dir $CURTMP \
    $samp/00_raw/original.bam \
    $rroot/HCV/references/HCV_subtype_consensus.fasta \
    $samp/05b_subtype/aligned_consensus.bam

samtools index $samp/05b_subtype/aligned_consensus.bam    
samtools idxstats $samp/05b_subtype/aligned_consensus.bam 

##########################################################################################
# Step 6: Scaffold contigs
##########################################################################################
mkdir -p $samp/06_scaffold

# If both amplicons are present we expect the scaffold to be 5957 bp (H77 coordinates)
# with 784 bp gap between amplicons.
# This is 63.5% of the reference length, with 13% ambiguous bases
# If only amplicon 2 is present, we expect the scaffold to be 3119 bp with few ambiguous
# bases.
# This is 33.3% of the reference length.

# Order, orient, and impute each scaffold
for f in $samp/05_subtype/*.fasta; do
    subt=$(basename ${f%%.*})
    if [[ ! $subt == "unassigned" ]]; then
        echo -e "\nOrder and orient subtype $subt\n"
        ref=$(ls ../references/subtypes/$subt*.fasta)
        assembly.py order_and_orient --tmp_dir $CURTMP \
            $f \
            $ref \
            $samp/06_scaffold/$subt.scaffold.fasta
        
        assembly.py impute_from_reference --tmp_dir $CURTMP \
            --newName "$subt.$samp" \
            --minLengthFraction 0 --minUnambig 0 \
            $samp/06_scaffold/$subt.scaffold.fasta \
            $ref \
            $samp/06_scaffold/$subt.imputed.fasta
    fi
done

# Combine all scaffolds into one reference
cat $samp/06_scaffold/*.imputed.fasta > $samp/06_scaffold/imputed.fasta

# Index reference
read_utils.py novoindex $samp/06_scaffold/imputed.fasta
read_utils.py index_fasta_samtools $samp/06_scaffold/imputed.fasta
read_utils.py index_fasta_picard $samp/06_scaffold/imputed.fasta

##########################################################################################
# Step 7: Refine assembly 1
##########################################################################################
mkdir -p $samp/07_refine1

# read_utils.py align_and_count_hits  --threads $NCPU --JVMmemory 50g --tmp_dir $CURTMP
assembly.py refine_assembly --threads $NCPU --JVMmemory 50g --tmp_dir $CURTMP \
    --min_coverage 2 \
    --novo_params "-r Random -l 30 -g 40 -x 20 -t 502" \
    --outBam $samp/07_refine1/refine1.bam \
    --outVcf $samp/07_refine1/refine1.vcf.gz \
    $samp/06_scaffold/imputed.fasta \
    $samp/00_raw/original.bam \
    $samp/07_refine1/refine1.fasta

samtools index $samp/07_refine1/refine1.bam
samtools idxstats $samp/07_refine1/refine1.bam

##########################################################################################
# Step 8: Refine assembly 2
##########################################################################################
mkdir -p $samp/08_refine2

assembly.py refine_assembly --threads $NCPU --JVMmemory 50g --tmp_dir $CURTMP \
    --min_coverage 30 \
    --novo_params "-r Random -l 40 -g 40 -x 20 -t 100" \
    --outBam $samp/08_refine2/refine2.bam \
    --outVcf $samp/08_refine2/refine2.vcf.gz \
    $samp/07_refine1/refine1.fasta \
    $samp/00_raw/original.bam \
    $samp/08_refine2/refine2.fasta

##########################################################################################
# Step 9:
##########################################################################################
mkdir -p $samp/09_fixed

# Copy the refined assembly
cp $samp/08_refine2/refine2.fasta $samp/09_fixed/consensus.fasta

# Index reference
read_utils.py novoindex $samp/09_fixed/consensus.fasta
read_utils.py index_fasta_samtools $samp/09_fixed/consensus.fasta
read_utils.py index_fasta_picard $samp/09_fixed/consensus.fasta

read_utils.py align_and_fix --threads $NCPU --JVMmemory 50g --tmp_dir $CURTMP \
    --outBamAll $samp/09_fixed/all_raw.bam \
    --outBamFiltered $samp/09_fixed/aligned_raw.bam \
    --novoalign_options "-r Random -l 40 -g 40 -x 20 -t 100 -k" \
    $samp/00_raw/original.bam \
    $samp/09_fixed/consensus.fasta

# FLASh joined reads
read_utils.py align_and_fix --threads $NCPU --JVMmemory 50g --tmp_dir $CURTMP \
    --outBamAll $samp/09_fixed/all_flash.bam \
    --outBamFiltered $samp/09_fixed/aligned_flash.bam \
    --novoalign_options "-r Random -l 40 -g 40 -x 20 -t 100 -k" \
    $samp/02_flash/flash.bam \
    $samp/09_fixed/consensus.fasta

# BLESS corrected reads
read_utils.py align_and_fix --threads $NCPU --JVMmemory 50g --tmp_dir $CURTMP \
    --outBamAll $samp/09_fixed/all_bless.bam \
    --outBamFiltered $samp/09_fixed/aligned_bless.bam \
    --novoalign_options "-r Random -l 40 -g 40 -x 20 -t 100 -k" \
    $samp/03_bless/bless.bam \
    $samp/09_fixed/consensus.fasta

##########################################################################################
# Step 10: Call SNP
##########################################################################################





##########################################################################################

mkdir -p $samp/10_test
assembly.py refine_assembly --threads $NCPU --JVMmemory 50g --tmp_dir $CURTMP \
    --min_coverage 30 \
    --novo_params "-r Random -l 40 -g 40 -x 20 -t 100" \
    --outBam $samp/10_test/refine_ec.bam \
    --outVcf $samp/10_test/refine_ec.vcf.gz \
    $samp/08_refine2/refine2.fasta \
    $samp/03_bless/bless.bam \
    $samp/10_test/refine_ec.fasta

mkdir -p $samp/11_test
assembly.py refine_assembly --threads $NCPU --JVMmemory 50g --tmp_dir $CURTMP \
    --min_coverage 30 \
    --novo_params "-r Random -l 40 -g 40 -x 20 -t 100" \
    --outBam $samp/11_test/refine_mg.bam \
    --outVcf $samp/11_test/refine_mg.vcf.gz \
    $samp/08_refine2/refine2.fasta \
    $samp/02_flash/flash.bam \
    $samp/11_test/refine_mg.fasta












