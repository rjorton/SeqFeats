### SeqFeats

SeqFeats is python program to create various SEQuence FEATureS from FASTA sequence files. 

This is replacing [VirusFeatures](https://github.com/rjorton/), and being expanded to included various different way to measure codon bias, codon pair bias, and dinucleotide bias, in both viral and host sequences. SeqFeats will also incorporate the original codon pair deoptimization and dinucleotide deoptimization routines.

The first script seqfeats_cai.py is to create Codon Adaptation Index (CAI) for individual sequences

---
    python seqfeats_cai.py fasta_file allow_Ns
---

An example using the provided human_cds_new_sel.fa file to create human_cds_new_sel_cai.txt:

---
    python seqfeats_cai.py human_cds_new_sel.fa y
---

The script takes:

Argument 1: A FASTA sequence file. This must be of in frame coding sequences - it can NOT strip out ORFs from genomes for you.

Argument 2:  allow_Ns - allow sequences with Ns in them - if allow, codons with Ns will be discounted but the rest of the codons will be analysed for the sequence. Default is TRUE. Allowable options are (case insensitive): true/yes/y/1 or false/no/n/0

I would recommend [COUSIN](http://cousin.ird.fr) (COdon Usage Similarity INdex) for calculating CAI and a wealth of other codon bias mesures. 
The seqfeats_cai.py script was created as COUSIN failed to process sequences will an N in them. It should give the same results as the CAI_59 metric in COUSIN.
