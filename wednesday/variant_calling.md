Variant Calling using Freebayes and Delly
==========================================

Now we will call variants (SNPs and indels) using two different programs, 'freebayes' and 'delly'. We will use the output from the alignment step as input into these programs. In the end you will get a VCF (Variant Call Format) file with genotype information for every variant across all the samples.

![fc04](fc04.png)

---

**1\.** First, create a new directory for calling variants:

    cd ~/variant_example
    mkdir 03-freebayes
    cd 03-freebayes

Now, let's link in the relevant files from our alignments, along with their indices:

    ln -s ../02-Alignment/*.all.bam .
    ln -s ../02-Alignment/*.all.bam.bai .

---

**2\.** First we will use software called 'freebayes' to find SNPs and short indels. Load the module and take a look at the help text:

    module load freebayes
    freebayes -h
