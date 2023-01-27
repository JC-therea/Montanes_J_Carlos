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
 
 In the example matrix the number of SNPs that change from an adenosine to a cytosine (A -> C) will be stored at cell (1,2) (in the example the value is 46247).
 
 - Extract the CDSs of the sequences that you want to estimate in fasta format.
 - With the transition matrix and the CDS run the program computePNPS_v2_10_03_2020.
 
To run the software on the example data type:

```
computePNPS_v2_10_03_2020 Transition_matrix_introns.tsv CDS_reference_few.fna > outFile.txt
```
The program outputs the transition matrix, the probability matrix and the espected transition of each CDS. The columns in the last section are: 

Protein_name, PN, PS, PN/PS, PN/(PN+PS)*length_of_the_sequence, PS/(PN+PS)*length_of_the_sequence

# Estimate amino acid expected substitution

Here I provide the 3 programs that we also used to estimate the expected amino acid substitution. The three scripts are in Compute_expected_amino_acid_substitution.

## tsvExtractorFromMSAfiles.py

Program that use several files to extract the alignment betwen 1 or more species and produce a tsv file from them. The program requieres a config file that is created the first time that you run the program without any existing config file.
The config file (argument -c) requires the following:
- Nodes from which you want to obtain the alignments.
- Species.
- Directory with all the MSA alignments.
- A summary file that has the information about all the transcripts, their Orthogroup, their ages and the last duplications that the transcripts had.

Additionally there is 2 more arguments requires by the program:
- Mode to work with, can be used for alignments among different species (0), duplications of the same species (1), or rooted alignments with the reference species in the first position of the config file (2) (argument -m).
- path off the output/s (-o).

Usage:

```
CFILE="Compute_expected_amino_acid_substitution/Example/configFile.txt"

echo "Nodes=N1,Conserved" > $CFILE
echo "Species=Drosophila_melanogaster,Drosophila_simulans" >> $CFILE
echo "MSAdir=Compute_expected_amino_acid_substitution/Example/MultipleSequenceAlignments/" >> $CFILE
echo "SummaryTable=Compute_expected_amino_acid_substitution/Example/Dmelanogaster_summary.tsv" >> $CFILE

python Compute_expected_amino_acid_substitution/tsvExtractorFromMSAfiles.py -c $CFILE -m 0 -o "Compute_expected_amino_acid_substitution/Example/"
```

## ExtractCodonFrequency.py

Program that outputs the codon frequency of all the transcripts in a given multifasta file. Penultimate column is the absolut frequency of all codons and last column is the frequency of each codon (extracted from the absolute frequencies). Usage:

```
CONSFA="Compute_expected_amino_acid_substitution/Example/Dmel_CDSs_conserved.fa"
CONSTSV="Compute_expected_amino_acid_substitution/Example/Dmel_CDSs_conserved.tsv"
RELCONSTSV="Compute_expected_amino_acid_substitution/Example/Dmel_CDSs_conserved_relative.tsv"
N1FA="Compute_expected_amino_acid_substitution/Example/Dmel_CDSs_N1.fa"
N1TSV="Compute_expected_amino_acid_substitution/Example/Dmel_CDSs_N1.tsv"
RELN1TSV="Compute_expected_amino_acid_substitution/Example/Dmel_CDSs_N1_relative.tsv"

python Compute_expected_amino_acid_substitution/ExtractCodonFrequency.py $CONSFA $CONSTSV
python Compute_expected_amino_acid_substitution/ExtractCodonFrequency.py $N1FA $N1TSV

# This part is to extract the last column of the files
LASTC=$(awk '{print NF}' $CONSTSV | sort -nu | tail -n 1)
cut -f 1,$LASTC $CONSTSV | tail -n +2 > $RELCONSTSV

LASTC=$(awk '{print NF}' $N1TSV | sort -nu | tail -n 1)
cut -f 1,$LASTC $N1TSV | tail -n +2 > $RELN1TSV
```
## expected_aa_changes_cons.py

Program that creates a tsv file with observed and expected frequencies from the provided observed files.
Requires pandas and argparse packages.
Also requires :
- A table of frequencies with all the coodons in the subset of proteins that you are searching for (-c). Files are the output of ExtractCodonFrequency.py.
- A substitution matrix of nucleotides (-ns). 
- Mode that you want to work (with or without direction of amino acid changes) (-m 'Dir' or 'noDir')
- Any number of inputs files with the aminoacid changing frequencies (-i). Files are the output of tsvExtractorFromMSAfiles.py.
- Directory to output the files (-o) (default = current directory).

Usage:

```
RELCONSTSV="Compute_expected_amino_acid_substitution/Example/Dmel_CDSs_conserved_relative.tsv"
RELN1TSV="Compute_expected_amino_acid_substitution/Example/Dmel_CDSs_N1_relative.tsv"
TRANSMATRIX="Compute_expected_PNPS/Example/Transition_matrix_introns_full_info.tsv"
CONSALN="Compute_expected_amino_acid_substitution/Example/Conserved_Drosophila_melanogaster_Drosophila_simulans.tsv"
N1ALN="Compute_expected_amino_acid_substitution/Example/N1_Drosophila_melanogaster_Drosophila_simulans.tsv"

mkdir -p "Compute_expected_amino_acid_substitution/Example/N1"
mkdir -p "Compute_expected_amino_acid_substitution/Example/Cons"

python Compute_expected_amino_acid_substitution/expected_aa_changes_cons.py -c $RELCONSTSV -ns $TRANSMATRIX -m noDir -i $CONSALN -o "Compute_expected_amino_acid_substitution/Example/Cons/"
python Compute_expected_amino_acid_substitution/expected_aa_changes_cons.py -c $RELN1TSV -ns $TRANSMATRIX -m noDir -i $N1ALN -o "Compute_expected_amino_acid_substitution/Example/N1/"
```


