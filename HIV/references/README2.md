# HIV References

## 1. Download HIV references sequences from LANL

The [Los Alamos HIV Database](https://hiv.lanl.gov/) is used to obtain reference sequences.
LANL has [curated alignments](https://www.hiv.lanl.gov/content/sequence/NEWALIGN/align.html)
available for download at the above link. We will download complete **genomes** for the 
**subtype references** from **2010** and use this to identify the subtypes of our samples.
We will also download complete **genomes** for **consensus/ancestral** sequences from **2002**.

The fields of interest are *Alignment type*, *Region*, *Subtype*, and *Year*.

![Screenshot of LANL download page](../../docs/lanl_hiv_download.png)

For the subtype reference, select *Subtype reference*, *GENOME*,
*M group without recombinants*, and *2010*.
Save the downloaded file as `HIV1_REF_2010_genome_DNA.fasta` in the `LANL` directory.

For the consensus reference, select *Consensus/Ancestral*, *GENOME*,
*M group without recombinants*, and *2002*.
Save the downloaded file as `HIV1_CON_2002_genome_DNA.fasta` in the `LANL` directory.


## 2. Extract references

The strategy we will use for identifying the subtype of each sample is to compare
sequences from the sample to a reference database containing a few representatives from
each subtype. The downloaded file contains approximately 4 representative sequences for 
each subtype.

The following script reformats the subtype reference names and sequences:

```python
from Bio import SeqIO
import re

"""
Load all sequences in FASTA file into a dictionary with new ID as the key
"""
seqs = {}
for s in SeqIO.parse('LANL/HIV1_REF_2010_genome_DNA.fasta','fasta'):
    subt = s.id.split('.')[1]
    acc = s.id.split('.')[-1]
    name = s.id.split('.')[-2]
    newid = 'HIV_%s.%s.%s' % (subt, acc, name)
    seqs[newid] = s

with open('HIV_subtype_refs.fasta', 'w') as outh:
    for rn in sorted(seqs.keys()):
        seqstr = str(seqs[rn].seq)
        seqstr = seqstr.replace('-', '').upper()
        seqstr = re.sub('\?','N', seqstr)
        print >>outh, '>%s' % rn
        for i in range(0,len(seqstr),100):
            print >>outh, seqstr[i:i+100]
        with open('subtypes/%s.fasta' % rn, 'w') as outs:
            print >>outs, '>%s' % rn
            for i in range(0,len(seqstr),100):
                print >>outs, seqstr[i:i+100]
``` 

We will also create references for the consensus sequences.

```python
from Bio import SeqIO
import re

"""
Load all sequences in FASTA file into a dictionary with the sequence ID as the key
"""
seqs = {}
for s in SeqIO.parse('LANL/HIV1_CON_2002_genome_DNA.fasta','fasta'):
    seqs[s.id] = s

"""
Extract the subtype consensus sequences to HIV_subtype_consensus.fasta
"""
connames = ['CONSENSUS_%s' % st for st in ['A1','A2','B','C','D','F1','G','H','O']]
with open('HIV_subtype_consensus.fasta', 'w') as outh:
    for rn in connames:
        seqstr = str(seqs[rn].seq)
        seqstr = seqstr.replace('-', '').upper()
        seqstr = re.sub('\?','N', seqstr)
        id = 'HIV_%s.con' % rn.split('_')[1]
        print >>outh, '>%s' % id
        for i in range(0,len(seqstr),100):
            print >>outh, seqstr[i:i+100]

```

There should be two fasta files in the current directory: 
`HIV_subtype_consensus.fasta`, and `HIV_subtype_refs.fasta`.
There should also be one `*.fasta` file for
each subtype in the `subtypes` directory.

## 3. Index references

Some software programs require special "index" files in order to use sequence files.
Here we are going to prebuild indexes used by different programs. This
builds the `novoalign`, `samtools` and `picard` indexes.

```bash
module load viral-ngs
for f in *.fasta; do
    read_utils.py novoindex $f
    read_utils.py index_fasta_samtools $f
    read_utils.py index_fasta_picard $f
done

for f in subtypes/*.fasta; do
    read_utils.py novoindex $f
    read_utils.py index_fasta_samtools $f
    read_utils.py index_fasta_picard $f
done
```

Here we build the `blast+` indexes.

```bash
module load blast+
for f in *.fasta; do
    makeblastdb -in $f -dbtype nucl -out ${f%.*}
done

for f in subtypes/*.fasta; do
    makeblastdb -in $f -dbtype nucl -out ${f%.*}
done
```

Finally, we build `bowtie2` indexes:

```bash
module load bowtie2
for f in *.fasta; do
    bowtie2-build $f ${f%.*}
done

for f in subtypes/*.fasta; do
    bowtie2-build $f ${f%.*}
done
```
