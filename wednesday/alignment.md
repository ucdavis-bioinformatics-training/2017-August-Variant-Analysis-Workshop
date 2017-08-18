Alignment of High-Throughput Sequencing Data using BWA
-------------------------------------------------------

In this section, we will use BWA (Burrows-Wheeler Aligner) to align all of our trimmed reads to a genome. We will align both the paired-end and single-end reads to the genome separately, and then merged the two alignments together. We will do this for all the samples.

**1\.** First go back to your variant_example directory and create a new directory called '02-Alignment':

    cd ~/variant_example
    mkdir 02-Alignment
    cd 02-Alignment

Link to the files we will be using for alignment:

    ln -s ../01-Trimming/*.sickle.fastq .
    ln -s ../01-Trimming/*.singles.fastq .

-----

**2\.** Now, in order to align any data, we need a reference to align against. We have reduced the size of our dataset (by selecting only one chromosome) so that the steps will occur relatively quickly. We have also reduced the genome down to just one chromosome, so that the alignment steps will happen quickly. First, create a directory for the reference and then download the reference:

    cd ~/variant_example
    mkdir ref
    cd ref
    wget https://ucdavis-bioinformatics-training.github.io/2017-August-Variant-Analysis-Workshop/wednesday/chr18.fa

First, we must index the reference to be able to align the data. Load the BWA module and look at the options:

    module load bwa
    bwa
    bwa index

Index the reference:

    bwa index chr18.fa 

This will produce 5 files in the reference directory that BWA will use during the alignment phase.

-----

**3\.** 
