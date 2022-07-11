# Compute PN/PS expected ratio

Here is stored a tool in C++ that computes the estimated transitions from one nucleotide to another. To do so you have to:

 - Estimate the introns.
 - Overlap introns with a vcf/gvcf file.
 - Obtain de transition matrix by identifying the SNPs in the vcf/gvcf file in the intronic regions. The matrix has to be ordered like in the following example:
 
 ```
         A       C       G       T
A       0   46247  101347   74995
C   55021       0   35614  117275
G  116998   35607       0   55342
T   75316  101451   46050       0
 ```
 
 In the example matrix the number of SNPs that change from an adenosine to a cytosine (A -> C) will be stored at cell (1,2) (in the example the value is 55021).
 
 - Extract the CDSs of the sequences that you want to estimate in fasta format.
 - With the transition matrix and the CDS run the program computePNPS_v2_10_03_2020.
 
To run the software on the example data type:

```
computePNPS_v2_10_03_2020 Transition_matrix_introns.tsv CDS_reference_few.fna > outFile.txt
```
The program outputs the transition matrix, the probability matrix and the espected transition of each CDS. The columns in the last section are: 

Protein_name, PN, PS, PN/PS, PN/(PN+PS)*length_of_the_sequence, PS/(PN+PS)*length_of_the_sequence
