# Microbial Genomes and Metagenomes

## Setup: Reserve a computer from Google Cloud

To carry out our bioinformatics work, we'll order up a virtual machine
(VM) from Google Cloud. Once you're signed up, go to the console. In
the "Quick Access" area, click on the box for "Compute Engine." The
screen should now say "Compute Engine" in the upper left. Next to
"Compute Engine," it should say "VM instances." This is your dashboad
for managing virtual machines.

We need to create a new virtual machine instance to do our work, so
click on "CREATE INSTANCE." We need to change a few settings for our
bioinformatics work. Here's a summary of the important stuff:

* Change the name to something you like.
* The default machine configuration, E2, is fine.
* Under "Machine type," go to the "CUSTOM" tab and crank up the
  memory all the way to 16 GB.
* The default operating system and version, Debian GNU/Linux 12 (bookworm),
  is fine for this workshop.
* Under "Firewall," check the boxes for "Allow HTTP traffic" and
  "Allow HTTPS traffic."

The machine will take a minute to start up. Once it does, you'll see a
button that says "SSH" under the "Connect" column on the right. When
you click on this button, your web browser will launch a command-line
terminal window. We are now ready to get started with bioinformatics.

## Setup: Install Conda

Conda is our system for installing third-party software needed for the
pipeline. Conda allows us to install whatever we need inside our home
directory, and places the software we need into separate environments
so the different software components don't clash with each other.

```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
# License agreement: hit return to see the license, hit "q" to exit the
#   license, then type "yes" to accept the terms.
# Miniconda3 will be installed into this location: hit return to confirm the
#   location.
# Do you wish to update your shell profile to automatically initialize conda?
#   Type "yes" and hit return.
```

At this point you'll be instructed to close and re-open the terminal
window. **Please do this.**

When the window re-opens, you should see "(base)" on the left-hand
side of your prompt. This indicates that you're in the base
environment for Conda, and that your installation was most likely
successful.

Clean up.

```bash
rm Miniconda3-latest-Linux-x86_64.sh
```

## 1. Microbial sequence reads (Kyle)

Assign a sequence read to each participant in the workshop.

```
>read_1  VH01989:2:22257VMNX:2:1101:26106:4636
CTTGTACACACCGCCCGTCACACCATGGGAGTGGGTTGCAAAAGAAGTAGGTAGCTTAACCTTCGGGAGGGCGCTTACCACTTTGTGATTCATGACTGGGGTGAAGTCGTAACAAGGTAACCGTAGGGGAACCTGCGGTTGGATCACCTCC
>read_2  VH01989:2:22257VMNX:2:1101:53631:7672
CCGCTGCACAACCCGGCTCACCTGATCGGTATCGAAGAAGCTCTGAAATCTTTCCCACAGCTGAAAGACAAAAACGTTGCTGTATTTGACACCGCGTTCCACCAGACTATGCCGGAAGAGTCTTACCTCTACGCCCTGCCGTACAACCTGT
>read_3  VH01989:2:22257VMNX:2:1101:27856:10742
GTGATCATCTGGTCGCTGGGGAATGAGTCAGGCCACGGCGCTAATCACGACGCACTCTATCGCTGGATCAAATCTGTCGATCCTTCCCGCCCGGTACAGTATGAAGGCGGCGGAGCCGACACCTCCGCAACCGATATTATTTGCCCGATGT
>read_4  VH01989:2:22257VMNX:2:1101:20241:30364
GTTTCAGCCTTAGTCATTATCGACTTTTGTTCGAGTGGAGTCCGCCGTGTCACTTTCGCTTTGGCAGCAGTGTCTTGCCCGATTGCAGGATGAGTTACCAGCCACAGAATTCAGTATGTGGATACGCCCATTGCAGGCGGAACTGAGCGAT
>read_5  VH01989:2:22257VMNX:2:1101:65344:36984
GTATTTGAACTGATCAAAGTGCCGGTTAAATCGCATACCAATCACGGGCGCAATGCGGAATATTTTGCCTGGGTGCAAAAACATTTACGTGAACACCCCGTCGATAGAGTCGTTGGATTTAATAAAATGCCGGGGCTGGACGTTTATTATG
>read_6  VH01989:2:22257VMNX:2:1101:41266:44325
GCGTTGATGCGCGTCAGGATCAGACTGACATTGAGATGTTCGAACTGCTGGAGCCAATTGCTGACGGTTTCCGTAACTATCGCGCTCGTCTGGACGTTTCCACCACCGAGTCACTGTTGATTGATAAAGCACAGCAACTGACGCTGACCGC
>read_7  VH01989:2:22257VMNX:2:1102:9712:42850
GCCTGCGCGAACATCCCCGCGTCTTTGACTCCCTGTTTTGTTCGATGGTCGCTGCCGGAGAAAAATCCGGACATCTCGACGTGGTGCTCAATCGCCTGGCGGATTACACCGAACAGCGGCAGCGTCTGAAATCACGCCTGCTGCAGGCCAT
>read_8  VH01989:2:22257VMNX:2:1102:69768:46932
GTGTGAACGTGCTGTCGAAAAACGATTCGATGAAGATTCAGATTGGTGCCAATGATAACCAGACGATCAGCATTGGCTTGCAACAAATCGACAGTACCACTTTGAATCTGAAAGGATTTACCGTGTCCGGCATGGCGGATTTCAGCGCGGC
>read_9  VH01989:2:22257VMNX:2:1103:52293:55388
AAGTATTGGAGCAGGTACAACCTATGTTTATGTTAATCTCGATCCTGTAATACAACCGGGCCAGAATCTGGTTGTAGACTTGTCTCAGCATATAAGTTGCTGGAATGATTACGGCGGCTGGTACGACACTGATCATATAAACCTGGTACAA
>read_10 VH01989:2:22257VMNX:2:1304:24271:31376
ATATTGAACCGTATAAGCTGGGGAATATTGCATGAGTATATTTATCTCATGGCTTGTTCTGATTATTTCGGTGGTCTGCGCCATTGGGATTATGCAAATTATTCATTCAGTAAAAAAGATTGAACGCTTTTTCACTGGCGAATAACAGCGC
>read_11 VH01989:2:22257VMNX:2:1611:32074:6489
CCGGAACGCGCCTCCACTTTCTTCCCGAGCCCGGATGGTGGAATCGGTAGACACAAGGGATTTAAAATCCCTCGGCGTTCGCGCTGTGCGGGTTCAAGTCCCGCTCCGGGTACCATGGGAAAGATAAGAATAAAATCAAAGCAATAAGCCT
```

[BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi) the read at the NCBI
website.

With audience participation, we’ll discuss the genes associated with
each read and learn about their role in bacteria.

## 2. Assemble a microbial genome (Ahmed)

We'll download our data files from Zenodo, a site for depositing
scientific data.

Download sequence data for a microbial genome.

```bash
wget 'https://zenodo.org/records/14229558/files/marc.2536_1.fastq.gz?download=1' -O marc.2536_1.fastq.gz
wget 'https://zenodo.org/records/14229558/files/marc.2536_2.fastq.gz?download=1' -O marc.2536_2.fastq.gz
```

Check out a FASTQ file.

```bash
zless marc.2536_1.fastq.gz
# Hit q to quit
```

Assemble the genome with Spades.

```bash
conda install -c bioconda -c conda-forge spades
```

```bash
spades.py --isolate -1 marc.2536_1.fastq.gz -2 marc.2536_2.fastq.gz -o spades_marc
```

Tour of output files.

Annotate with Prokka.

```bash
conda install -c bioconda -c conda-forge prokka
```

```bash
prokka spades_marc/contigs.fasta --outdir prokka_marc --prefix marc.2536
```

Tour of output files. Discuss genome size, GC content, total number of
genes, rRNA genes, etc.

Find the gene in the genome that contains the example read.

```bash
nano example_read.fasta
# Paste in your seqeunce read
# Control-x to exit
```

```bash
makeblastdb -dbtype nucl -in prokka_marc/marc.2536.ffn
```

```bash
blastn -query example_read.fasta -db prokka_marc/marc.2536.ffn -outfmt 7
```

Save the nucleotide and protein sequence of the matching gene.

```bash
less prokka_marc/marc.2536.ffn
# Copy the nucleotide sequence to the clipboard
```

```bash
nano example_gene.ffn
# Paste in the nucleotide sequence
# Control-x to exit
```

```bash
less prokka_marc/marc.2536.faa
# Copy the protein sequence to the clipboard
```

```bash
nano example_gene.faa
# Paste in the protein sequence
# Control-x to exit
```

## 3. Metagenome analysis (Kyle)

Download a human-filtered FASTQ file from infant feces collected hours
after birth.

```bash
wget 'https://zenodo.org/records/14229558/files/s188.STL.V01_1.4d_R1.fastq.gz?download=1' -O s188.STL.V01_1.4d_R1.fastq.gz
wget 'https://zenodo.org/records/14229558/files/s188.STL.V01_1.4d_R2.fastq.gz?download=1' -O s188.STL.V01_1.4d_R2.fastq.gz
```

Assemble with MEGAHIT.

```bash
conda install -c bioconda -c conda-forge megahit
```

```bash
megahit -1 s188.STL.V01_1.4d_R1.fastq.gz -2 s188.STL.V01_1.4d_R2.fastq.gz -o megahit_s188
```

Annotate with Prokka.

```bash
prokka megahit_s188/final.contigs.fa --outdir prokka_s188 --prefix s188.STL.V01
```

```bash
less prokka_s188/s188.STL.V01.tsv
```

Characterize taxonomy with Sourmash.

```bash
conda install -c bioconda -c conda-forge sourmash
```

```bash
wget 'https://osf.io/4f8n3/download' -O genbank-k31.lca.json.gz
```

```bash
sourmash sketch dna -p scaled=10000,k=31,abund s188.STL.V01_1.4d_R1.fastq.gz --name-from-first
```

```bash
sourmash gather -k 31 s188.STL.V01_1.4d_R1.fastq.gz.sig genbank-k31.lca.json.gz
```

This sample contains a representative of the same species we assembled
in the morning session. Find matches for the example genes in the
isolated organism.

Search using the nucleotide sequence.

```bash
makeblastdb -dbtype nucl -in prokka_s188/s188.STL.V01.ffn
```

```bash
blastn -query example_gene.ffn -db prokka_s188/s188.STL.V01.ffn -outfmt 7
```

Search using the protein sequence.

```bash
makeblastdb -dbtype prot -in prokka_s188/s188.STL.V01.faa
```

```bash
blastp -query example_gene.faa -db prokka_s188/s188.STL.V01.faa -outfmt 7
```

## 4. Exploring MicrobiomeDB (Dan)

Download a table of taxonomic abundances from the infant study
[here](igram_birth_1m.biom).

Import the BIOM file into MicrobiomeDB.

Compare sourmash results with abundances in BIOM file. Compare birth
and one month time points.
