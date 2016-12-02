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
export TMPDIR=$CURTMP

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

# Subsample original reads
# samtools view -bs 123.01 $samp/00_raw/original.bam > $samp/00_raw/onepct1.bam
# samtools view -bs 987.01 $samp/00_raw/original.bam > $samp/00_raw/onepct2.bam
# samtools view -bs 123.10 $samp/00_raw/original.bam > $samp/00_raw/tenpct1.bam
# samtools view -bs 987.10 $samp/00_raw/original.bam > $samp/00_raw/tenpct2.bam

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

# Subsample flash reads
# samtools view -bs 123.01 $samp/02_flash/flash.bam > $samp/02_flash/flash.onepct1.bam
# samtools view -bs 123.10 $samp/02_flash/flash.bam > $samp/02_flash/flash.tenpct1.bam

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

# Subsample bless reads
# samtools view -bs 123.01 $samp/03_bless/bless.bam > $samp/03_bless/bless.onepct1.bam
# samtools view -bs 123.10 $samp/03_bless/bless.bam > $samp/03_bless/bless.tenpct1.bam    

##########################################################################################
# Step 4: Denovo assembly using Trinity
##########################################################################################
ADAPTERFILE=$rroot/adapters/adapters.fasta

# Assemble with the original reads
mkdir -p $samp/04_assembly
assembly.py assemble_trinity --threads $NCPU --JVMmemory 50g --tmp_dir $CURTMP \
    --n_reads 1000000 \
   $samp/00_raw/original.bam \
   $ADAPTERFILE \
   $samp/04_assembly/contigs.fa

python $rroot/HIV/scripts/assembly_stats.py < $samp/04_assembly/contigs.fa | tee $samp/04_assembly/summary.txt

# Assemble with error-corrected reads
mkdir -p $samp/04c_assembly
assembly.py assemble_trinity --threads $NCPU --JVMmemory 50g --tmp_dir $CURTMP \
    --n_reads 1000000 \
   $samp/03_bless/bless.bam \
   $ADAPTERFILE \
   $samp/04c_assembly/contigs.fa

python $rroot/HIV/scripts/assembly_stats.py < $samp/04c_assembly/contigs.fa | tee $samp/04c_assembly/summary.txt

# Assemble with extended reads
mkdir -p $samp/04x_assembly
source activate viral-ngs
trintmp=$(mktemp -d  -t tmp-assembly-assemble_trinity-XXXXXX)
Trinity --CPU $NCPU --bflyHeapSpace 50G \
    --min_contig_length 300 \
    --seqType fq \
    --single $PWD/$samp/02_flash/flash.extendedFrags.fastq \
    --output $trintmp

cp $trintmp/Trinity.fasta $samp/04x_assembly/contigs.fa
rm -rf $trintmp

source deactivate

python $rroot/HIV/scripts/assembly_stats.py < $samp/04x_assembly/contigs.fa | tee $samp/04x_assembly/summary.txt

##########################################################################################
# Step 5: Assign contigs to subtypes
##########################################################################################
mkdir -p $samp/05_subtype
python $rroot/HIV/scripts/assign_contigs.py \
    --cutoff 400 \
    $samp/04_assembly/contigs.fa \
    $rroot/HIV/references/HIV_subtype_refs \
    $samp/05_subtype | tee $samp/05_subtype/summary.txt

mkdir -p $samp/05c_subtype
python $rroot/HIV/scripts/assign_contigs.py \
    --cutoff 400 \
    $samp/04c_assembly/contigs.fa \
    $rroot/HIV/references/HIV_subtype_refs \
    $samp/05c_subtype | tee $samp/05c_subtype/summary.txt

mkdir -p $samp/05x_subtype
python $rroot/HIV/scripts/assign_contigs.py \
    --cutoff 400 \
    $samp/04x_assembly/contigs.fa \
    $rroot/HIV/references/HIV_subtype_refs \
    $samp/05x_subtype | tee $samp/05x_subtype/summary.txt


##########################################################################################
# Step 6: Scaffold contigs
##########################################################################################

mkdir -p $samp/06_scaffold
cat $samp/05_subtype/subtypes.config | while read l; do
    # Subtype
    subt=$(cut -f1 <<<"$l")
    # Selected isolate for scaffolding
    isol=$(cut -f2 <<<"$l")
    echo -e "\nOrder and orient subtype $subt"
    echo -e "Isolate $isol\n"
    ref=../references/subtypes/$isol.fasta
    assembly.py order_and_orient --tmp_dir $CURTMP \
            $samp/05_subtype/$subt.fasta \
            $ref \
            $samp/06_scaffold/$subt.scaffold.fasta
    if [[ -s $samp/06_scaffold/$subt.scaffold.fasta ]]; then
        assembly.py impute_from_reference --tmp_dir $CURTMP \
            --newName "$subt.$samp" \
            --minLengthFraction 0 --minUnambig 0 \
            $samp/06_scaffold/$subt.scaffold.fasta \
            $ref \
            $samp/06_scaffold/$subt.imputed.fasta
    fi
done

# Combine all scaffolds into one reference
cat $samp/06_scaffold/*.scaffold.fasta > $samp/06_scaffold/scaffold.fasta
python $rroot/HIV/scripts/scaffold_stats.py < $samp/06_scaffold/scaffold.fasta | tee $samp/06_scaffold/summary.txt

# Combine all imputed scaffolds
cat $samp/06_scaffold/*.imputed.fasta > $samp/06_scaffold/imputed.fasta

# Index reference
read_utils.py novoindex $samp/06_scaffold/imputed.fasta  
read_utils.py index_fasta_samtools $samp/06_scaffold/imputed.fasta  
read_utils.py index_fasta_picard $samp/06_scaffold/imputed.fasta

##########################################################################################

mkdir -p $samp/06c_scaffold
cat $samp/05c_subtype/subtypes.config | while read l; do
    # Subtype
    subt=$(cut -f1 <<<"$l")
    # Selected isolate for scaffolding
    isol=$(cut -f2 <<<"$l")
    echo -e "\nOrder and orient subtype $subt"
    echo -e "Isolate $isol\n"
    ref=../references/subtypes/$isol.fasta
    assembly.py order_and_orient --tmp_dir $CURTMP \
            $samp/05c_subtype/$subt.fasta \
            $ref \
            $samp/06c_scaffold/$subt.scaffold.fasta
    if [[ -s $samp/06c_scaffold/$subt.scaffold.fasta ]]; then
        assembly.py impute_from_reference --tmp_dir $CURTMP \
            --newName "$subt.$samp" \
            --minLengthFraction 0 --minUnambig 0 \
            $samp/06c_scaffold/$subt.scaffold.fasta \
            $ref \
            $samp/06c_scaffold/$subt.imputed.fasta
    fi
done

# Combine all scaffolds into one reference
cat $samp/06c_scaffold/*.scaffold.fasta > $samp/06c_scaffold/scaffold.fasta
python $rroot/HIV/scripts/scaffold_stats.py < $samp/06c_scaffold/scaffold.fasta | tee $samp/06c_scaffold/summary.txt

# Combine all imputed scaffolds
cat $samp/06c_scaffold/*.imputed.fasta > $samp/06c_scaffold/imputed.fasta

# Index reference
read_utils.py novoindex $samp/06c_scaffold/imputed.fasta  
read_utils.py index_fasta_samtools $samp/06c_scaffold/imputed.fasta  
read_utils.py index_fasta_picard $samp/06c_scaffold/imputed.fasta

##########################################################################################

mkdir -p $samp/06x_scaffold
cat $samp/05x_subtype/subtypes.config | while read l; do
    # Subtype
    subt=$(cut -f1 <<<"$l")
    # Selected isolate for scaffolding
    isol=$(cut -f2 <<<"$l")
    echo -e "\nOrder and orient subtype $subt"
    echo -e "Isolate $isol\n"
    ref=../references/subtypes/$isol.fasta
    assembly.py order_and_orient --tmp_dir $CURTMP \
            $samp/05x_subtype/$subt.fasta \
            $ref \
            $samp/06x_scaffold/$subt.scaffold.fasta
    if [[ -s $samp/06x_scaffold/$subt.scaffold.fasta ]]; then
        assembly.py impute_from_reference --tmp_dir $CURTMP \
            --newName "$subt.$samp" \
            --minLengthFraction 0 --minUnambig 0 \
            $samp/06x_scaffold/$subt.scaffold.fasta \
            $ref \
            $samp/06x_scaffold/$subt.imputed.fasta
    fi
done

# Combine all scaffolds into one reference
cat $samp/06x_scaffold/*.scaffold.fasta > $samp/06x_scaffold/scaffold.fasta
python $rroot/HIV/scripts/scaffold_stats.py < $samp/06x_scaffold/scaffold.fasta | tee $samp/06x_scaffold/summary.txt

# Combine all imputed scaffolds
cat $samp/06x_scaffold/*.imputed.fasta > $samp/06x_scaffold/imputed.fasta

# Index reference
read_utils.py novoindex $samp/06x_scaffold/imputed.fasta  
read_utils.py index_fasta_samtools $samp/06x_scaffold/imputed.fasta  
read_utils.py index_fasta_picard $samp/06x_scaffold/imputed.fasta


##########################################################################################
# Steps 7 & 8: Refine assembly
##########################################################################################

### Raw data #############################################################################
mkdir -p $samp/07_refine1

cp $samp/06_scaffold/imputed.fasta $samp/07_refine1/initial.fasta
. $rroot/HIV/scripts/refine_bowtie.sh \
    $samp/07_refine1 \
    "--very-fast-local" \
    5 \
    $samp/00_raw/original_1.fastq \
    $samp/00_raw/original_2.fastq

mkdir -p $samp/08_refine2

cp $samp/07_refine1/refined.fasta $samp/08_refine2/initial.fasta
. $rroot/HIV/scripts/refine_bowtie.sh \
    $samp/08_refine2 \
    "--very-sensitive-local -N 1 " \
    0.01 \
    $samp/00_raw/original_1.fastq \
    $samp/00_raw/original_2.fastq

python $rroot/HIV/scripts/scaffold_stats.py < $samp/07_refine1/refined.fasta | tee $samp/07_refine1/summary.txt
python $rroot/HIV/scripts/scaffold_stats.py < $samp/08_refine2/refined.fasta | tee $samp/08_refine2/summary.txt

### Error corrected ######################################################################
mkdir -p $samp/07c_refine1

cp $samp/06c_scaffold/imputed.fasta $samp/07c_refine1/initial.fasta
. $rroot/HIV/scripts/refine_bowtie.sh \
    $samp/07c_refine1 \
    "--very-fast-local" \
    5 \
    $samp/03_bless/bless.1.corrected.fastq \
    $samp/03_bless/bless.2.corrected.fastq

mkdir -p $samp/08c_refine2

cp $samp/07c_refine1/refined.fasta $samp/08c_refine2/initial.fasta
. $rroot/HIV/scripts/refine_bowtie.sh \
    $samp/08c_refine2 \
    "--very-sensitive-local -N 1 " \
    0.01 \
    $samp/03_bless/bless.1.corrected.fastq \
    $samp/03_bless/bless.2.corrected.fastq

python $rroot/HIV/scripts/scaffold_stats.py < $samp/07c_refine1/refined.fasta | tee $samp/07c_refine1/summary.txt
python $rroot/HIV/scripts/scaffold_stats.py < $samp/08c_refine2/refined.fasta | tee $samp/08c_refine2/summary.txt

### Extended #############################################################################
mkdir -p $samp/07x_refine1

cp $samp/06x_scaffold/imputed.fasta $samp/07x_refine1/initial.fasta
. $rroot/HIV/scripts/refine_bowtie.sh \
    $samp/07x_refine1 \
    "--very-fast-local" \
    5 \
    $samp/02_flash/flash.notCombined_1.fastq \
    $samp/02_flash/flash.notCombined_2.fastq \
    $samp/02_flash/flash.extendedFrags.fastq

mkdir -p $samp/08x_refine2

cp $samp/07x_refine1/refined.fasta $samp/08x_refine2/initial.fasta
. $rroot/HIV/scripts/refine_bowtie.sh \
    $samp/08x_refine2 \
    "--very-sensitive-local -N 1 " \
    0.01 \
    $samp/02_flash/flash.notCombined_1.fastq \
    $samp/02_flash/flash.notCombined_2.fastq \
    $samp/02_flash/flash.extendedFrags.fastq

python $rroot/HIV/scripts/scaffold_stats.py < $samp/07x_refine1/refined.fasta | tee $samp/07x_refine1/summary.txt
python $rroot/HIV/scripts/scaffold_stats.py < $samp/08x_refine2/refined.fasta | tee $samp/08x_refine2/summary.txt

##########################################################################################
# Step 9: Fix consensus
##########################################################################################
mkdir -p $samp/09_fixed
cp $samp/08_refine2/refined.fasta $samp/09_fixed/consensus.fasta
. $rroot/HIV/scripts/align_and_fix.sh \
    $samp/09_fixed \
    "--very-sensitive-local -N 1 " \
    $samp/00_raw/original_1.fastq \
    $samp/00_raw/original_2.fastq

mkdir -p $samp/09c_fixed
cp $samp/08c_refine2/refined.fasta $samp/09c_fixed/consensus.fasta
. $rroot/HIV/scripts/align_and_fix.sh \
    $samp/09c_fixed \
    "--very-sensitive-local -N 1 " \
    $samp/03_bless/bless.1.corrected.fastq \
    $samp/03_bless/bless.2.corrected.fastq

mkdir -p $samp/09x_fixed
cp $samp/08x_refine2/refined.fasta $samp/09x_fixed/consensus.fasta
. $rroot/HIV/scripts/align_and_fix.sh \
    $samp/09x_fixed \
    "--very-sensitive-local -N 1 " \
    $samp/02_flash/flash.notCombined_1.fastq \
    $samp/02_flash/flash.notCombined_2.fastq \
    $samp/02_flash/flash.extendedFrags.fastq



