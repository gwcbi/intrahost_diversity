# HCV

The HCV data is located here: `/import/cbi/Projects/UMD_HCV/fastq`

## Pilot study

Selected 10 samples:

### Sample 01

1a, only 2nd amplicon

ID: `P153SEQ0008_P189M0015_S15`

### Sample 02

1b, both amplicons present

ID: `P153SEQ0010_P189M0048_S8`

### Sample 03 

1b, only 2nd amplicon
Low diversity (54 vars)

ID: `P153SEQ0059_P189M0115_S22`

gunzip -c fastq/

### Sample 04

4m, both amplicons

ID: `P153SEQ0122_P189M0171_S22`

### Sample 05

1a, both amplicons
59% alignment rate to all, 33% alignment rate to selected reference.
In subtyping, 15820 map to 1a and 14,516 map to next best subtype.

ID: `P153SEQ0062_P189M0089_S1`

### Sample 06

1a, both amplicons
66% alignment rate to all, 30% alignment rate to selected reference.
In subtyping, 27,065 map to 1b and 10,777 map to next best subtype.

ID: `P153SEQ0119_P189M0199_S3`

### Sample 07

4a, both amplicons
10,710 map to next best subtype

ID: `P153SEQ0111_P189M0169_S27`

### Sample 08

1b, both amplicons
Very deep, over 1M reads

ID: `P153SEQ0079_P189M0149_S21`

### Sample 09

1a, both amplicons
Low diversity (8 vars)

ID: `P153SEQ0119_P189M0237_S40`

### Sample 10

1b, both amplicons

ID: `P153SEQ0119_P189M0223_S26`

## Setup pilot directory

```bash
while read l; do 
    samp=$(awk '{print $1}' <<<$l)
    id=$(awk '{print $2}' <<<$l)
    mkdir -p pilot/$samp
    gunzip -c $(ls fastq/$id*R1*.fastq.gz) > pilot/$samp/original_1.fastq
    gunzip -c $(ls fastq/$id*R2*.fastq.gz) > pilot/$samp/original_2.fastq
done < pilot/key.txt
```