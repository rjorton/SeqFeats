### SeqFeats

CURRENTLY UNDERGOING TESTING - NOT FINISHED

SeqFeats is python program to create various SEQuence FEATureS from FASTA sequence files. 

This is replacing [VirusFeatures](https://github.com/rjorton/), and being expanded to included various different way to measure codon bias, codon pair bias, and dinucleotide bias, in both viral and host sequences. SeqFeats will also incorporate the original codon pair deoptimization and dinucleotide deoptimization routines.

The first script seqfeats_cai.py is to create Codon Adaptation Index (CAI) for individual sequences.

---
    python seqfeats_cai.py fasta_file allow_Ns
---

An example using the provided human_cds_new_sel.fa file to create human_cds_new_sel_cai.txt:

---
    python seqfeats_cai.py human_cds_new_sel.fa
    
    python seqfeats_cai.py human_cds_new_sel.fa y
---

The script takes:

Argument 1: A FASTA sequence file. This must be of in frame coding sequences - it can NOT strip out ORFs from genomes for you.

Argument 2 (optional):  allow_Ns - allow sequences with Ns in them - if allow, codons with Ns will be discounted but the rest of the codons will be analysed for the sequence. Default is TRUE. Allowable options are (case insensitive): true/yes/y/1 or false/no/n/0

It outputs the following messages to screen:

---
    python seqfeats_cai.py human_cds_new_sel.fa 
    
    seqfeats_cai started...
    fasta_calculate_cai
    Input FASTA file = human_cds_new_sel.fa 
    Output CU table file = human_cds_new_sel_cu_table.txt
    Output CAI data file = human_cds_new_sel_cai.txt
    Sequences with Ns allowed [1 = True, 0 = False] = 1
    Calculated ref codon table: processed sequences = 23312, total sequences = 23312
    Outputted CU table
    Calculated CAI for each sequence: processed sequences = 23312, total sequences = 23312
    ...finished seqfeats_cai
---

File 1 _cu_table.txt: codon usgae table of all the sequences

File 2 _cai.txt: the CAI score for each sequence

I would recommend [COUSIN](http://cousin.ird.fr) (COdon Usage Similarity INdex) for calculating CAI and a wealth of other codon bias mesures. 
The seqfeats_cai.py script was created as COUSIN failed to process sequences will an N in them. It should give the same results as the CAI_59 metric in COUSIN.


seqfeats_isg: observed/expected and raw dinucleotide frequencies of host sequences

seqfeats_cpb: nucleotide, dinculeotide frequencies, and codon, aa, and codon-pairs biases of host sequences