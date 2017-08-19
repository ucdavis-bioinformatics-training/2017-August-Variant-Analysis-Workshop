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

Take a look at it:

    less chr18.fa

Press 'q' to quit out of less. We must index the reference to be able to align the data. Load the BWA module and look at the options:

    module load bwa
    bwa
    bwa index

Index the reference:

    bwa index chr18.fa 

This will produce 5 files in the reference directory that BWA will use during the alignment phase.

-----

**3\.** Now, go back to your alignment directory and list the files:

    cd ../02-Alignment
    ls -l

You should see 5 sets of files, one for each sample. Each set should contain two paired-end reads (R1 & R2), and a single-end reads file. So, let's start by running bwa on the first sample. We will be using the 'bwa mem' subcommand with our files. Take a look at the options:

    bwa mem
    
Note that the Usage shows that we need to give bwa a location for the 'idxbase', which is the path to the reference. Now, we will align the two paired-end files and redirect the alignment output (in SAM format) to a file. We will use 4 threads (processors) and add read group (i.e sample ID) information to the alignment:

    bwa mem -t 4 -R "@RG\tID:A8100\tSM:A8100" ../ref/chr18.fa A8100.chr18.R1.sickle.fastq A8100.chr18.R2.sickle.fastq > A8100.chr18.paired.sam

This step will take about 10 minutes to run.

-----

**4\.** Next, we will align the single-end file for the same sample:

    bwa mem -t 4 ../ref/chr18.fa A8100.chr18.singles.fastq > A8100.chr18.singles.sam

Then, we need to convert the sam files into bam files for downstream processing. We will use a tool called 'samtools' to do this. Load the samtools module and take a look at the various subcommands and options:

    module load samtools
    samtools
    samtools view

We will use 'samtools view' to convert the sam files into a bam files (binary sam)... using 4 threads and the '-b' flag to output bam format:

    samtools view -@ 4 -b A8100.chr18.paired.sam > A8100.chr18.paired.bam

Do the same for the single-end alignment file:

    samtools view -@ 4 -b A8100.chr18.singles.sam > A8100.chr18.singles.bam

-----

**5\.** Now we sort the files by alignment position so that they are easy to merge later on... with 4 threads and specifying the output file:

    samtools sort -@ 4 -o A8100.chr18.paired.sorted.bam A8100.chr18.paired.bam

Do the same for the single-end file:

    samtools sort -@ 4 -o A8100.chr18.singles.sorted.bam A8100.chr18.singles.bam

Finally, we merge the two files into one alignment file:

    samtools merge -@ 4 A8100.chr18.all.bam A8100.chr18.singles.sorted.bam A8100.chr18.paired.sorted.bam

And index the final alignment file. This will allow downstream programs to easily read data from the bam file:

    samtools index A8100.chr18.all.bam

-----

**6\.** In the next step, we will use another Slurm script to run all the alignment commands on the rest of the samples. First download the script:

    wget https://ucdavis-bioinformatics-training.github.io/2017-August-Variant-Analysis-Workshop/wednesday/bwa.sh

Take a look at it:

    cat bwa.sh
