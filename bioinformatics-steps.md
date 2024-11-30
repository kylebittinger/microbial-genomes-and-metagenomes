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
rm Miniconda3-py311_23.10.0-1-Linux-x86_64.sh
```

## 1. Microbial sequence reads (Kyle)

Each student will be assigned a sequence read. We’ll go over what a
sequence read looks like and how the data was generated. Students will
take their sequence read and BLAST it at NCBI. All the sequence reads
will be from different genes within an E. coli genome, but the
students won’t know that yet.

With audience participation, we’ll discuss each gene and learn about
its role.  I plan to have housekeeping genes, something from the lac
operon, something related to oxygen stress, and one unknown
hypothetical protein.

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
zless marc.bacteremia.2536_1.fastq.gz
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
prokka spades_marc/contigs.fasta --outdir prokka_marc
```

Tour of output files. Discuss genome size, GC content, total number of
genes, rRNA genes, etc.

Find the genes in the E. coli genome that contain the student-assigned
reads.

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
prokka megahit_s188/final.contigs.fa --outdir prokka_s188
```

This sample will contain E. coli, but not the same E. coli that was
assembled in isolation.

Find homologs of genes that were in the isolated E. coli genome and
discuss the accessory genome of E. coli.

## 4. Exploring MicrobiomeDB (Dan)

Download a table of taxonomic abundances from the infant study and
import them into MicrobiomeDB.

Students will see how the infant gut microbiome compares to adult gut
microbiome, other human body sites, and the universe of microbiomes in
the database. It would be great if we could circle back to E. coli
somehow at the very end, maybe by summarizing E. coli abundance across
human body sites and other sample types.
