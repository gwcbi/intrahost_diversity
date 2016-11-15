# HIV References
Create ids.txt file


Want to put subtype consensus sequences in fasta file --
	HIV1_CON_2002_genome_DNA
	Has consensus sequences for each subtype 


HIV consensus data:
Alignment type: Consensus/Ancestral
Year: 2002 
Organism: HIV-1/SIVcpz 
DNA/Protein: DNA 
Region: genome 
Subtype: ALL 
Format: FASTA 
Alignment ID : 102CG1
Number of sequences: 20


HIV subtype data:
Alignment type: Subtype reference
Year: 2010 
Organism: HIV-1/SIVcpz 
DNA/Protein: DNA 
Region: genome 
Subtype: ALL 
Format: FASTA 
Alignment ID : 110RG1
Number of sequences: 170






Download HIV references sequences from LANL
Go to: Alignments, Curated Alignments
Whole genome reference
Select: Consensus/Ancestral, HIV-1/SIVcpz, GENOME, All, DNA, 2002, fasta
Save downloaded file as HIV1_CON_2002_genome_DNA.fasta
Go to: Alignments, Curated Alignments
Subtype reference
Select: Subtype reference, HIV-1/SIVcpz, GENOME, All, DNA, 2010, fasta
Save downloaded file as HIV_SUB_2010_genome_DNA.fasta
Extract subset of references


















from Bio import SeqIO
import re




"""
Load all sequences in FASTA file into a dictionary with the sequence ID as the key
"""
seqs = {}
for s in SeqIO.parse('LANL/HIV1_CON_2002_Genome_DNA.fasta','fasta'):
    seqs[s.id] = s
"""

Extract the subtype consensus sequences to HIV_subtype_consensus.fasta
***********in the HIV consensus fasta file, there are consensus sequences for each subtype******

"""
connames = [k for k in seqs.keys() if k.split('.')[1].startswith('CON')]
connames.sort(key=lambda x:x.split('.')[1].split('(')[0])
with open('HIV_subtype_consensus.fasta', 'w') as outh:
    for rn in connames:
        seqstr = str(seqs[rn].seq)
        seqstr = seqstr.replace('-', '').upper()
        seqstr = re.sub('\?','N', seqstr)
        id = 'HIV_%s.con' % rn.split('.')[1].split('(')[0].split('_')[1]
        print >>outh, '>%s' % id
        for i in range(0,len(seqstr),100):
            print >>outh, seqstr[i:i+100]




"""
Extract one sequence from each subtype (from HIVids.txt) to HIV_subtype_refs.fasta.
Also, write each sequence within its own file in the "subtypes" directory
"""
refnames = [l.strip() for l in open('HIVids.txt','r')]
refnames.sort(key=lambda x: x.split('.')[1].lower())
with open('HIV_subtype_refs.fasta', 'w') as outh:
    for rn in refnames:
        seqstr = str(seqs[rn].seq)
        seqstr = seqstr.replace('-', '').upper()
        seqstr = re.sub('\?','N', seqstr)
        id = 'HIV_%s' % rn.split('.')[1].lower()
        id = '%s.%s' % (id, rn.split('.')[-1])
        print >>outh, '>%s' % id
        for i in range(0,len(seqstr),100):
            print >>outh, seqstr[i:i+100]        
        with open('subtypes/%s.fasta' % id, 'w') as outs:
            print >>outs, '>%s' % id
            for i in range(0,len(seqstr),100):
                print >>outs, seqstr[i:i+100]



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
module load blast+
for f in *.fasta; do
    makeblastdb -in $f -dbtype nucl -out ${f%.*}
done

for f in subtypes/*.fasta; do
    makeblastdb -in $f -dbtype nucl -out ${f%.*}
done


module load bowtie2
for f in *.fasta; do
    bowtie2-build $f ${f%.*}
done

for f in subtypes/*.fasta; do
    bowtie2-build $f ${f%.*}
done








TXT
Save as HIVids.txt


Ref.A1.AU.03.PS1044_Day0.DQ676872
Ref.A2.CD.97.97CDKTB48.AF286238
Ref.B.FR.83.HXB2_LAI_IIIB_BRU.K03455
Ref.C.BR.92.BR025_d.U52953
Ref.D.CD.83.ELI.K03454
Ref.F1.BE.93.VI850.AF077336
Ref.G.BE.96.DRCBL.AF084936
Ref.H.BE.93.VI991.AF190127
