Effect Prediction using SnpEff
===============================

In this section we will be using a software called 'snpeff' to do effect prediction of our variants. We will build an effect prediction database using our reference and annotation and then use that database to run effect prediction. This will give us a VCF file with an extra "ANN" field per variant, which will give us the effect of that variant.

![fc05](fc05.png)

---

**1\.** As before, let's create a directory for this section:

    cd ~/variant_example
    mkdir 05-snpeff
    cd 05-snpeff

And let's link in our variant file from our freebayes run:

    ln -s ../03-freebayes/all.chr18.vcf

---

**2\.** Normally we would use the snpeff module to run snpeff, but in this case, since we want to create our own database, we need to download our own instance of snpeff. Go to the [SnpEff sourceforge page](http://snpeff.sourceforge.net/SnpEff_manual.html) to look at the documentation. Use 'wget' to get the source code to your 05-snpeff directory:

    wget https://downloads.sourceforge.net/project/snpeff/snpEff_latest_core.zip

Unzip the file:

    unzip snpEff_latest_core.zip
    
Look at the directories created and then go into the 'snpEff' directory:

    ls -ltrh
    cd snpEff 

---

**3\.** Now go back to the [snpeff documentation](http://snpeff.sourceforge.net/SnpEff_manual.html) and click on the "Building Databases" section. We will be following directions from this section. First, we need to add lines to the snpeff configuration file. Use nano to open the file:

    nano snpEff.config

And then add the following lines to the file in the "Non-Standard Databases" section:

> # chr18
> chr18.genome : chr18

Save and Exit.

---

**4\.** Next we need to make directories where the database will reside:

    mkdir data
    mkdir data/chr18 
    cd data/chr18

Now copy over the reference file, but call it "sequences.fa":

    cp ../ref/chr18.fa sequences.fa
    cd data/chr18

Now copy over the reference file, but call it "sequences.fa":

    cp ../ref/chr18.fa sequences.fa
    cd data/chr18

Now, link to the reference file, but call it "sequences.fa":

    ln -s ../../../../ref/chr18.fa sequences.fa

And get the annotation file:

    wget https://raw.githubusercontent.com/ucdavis-bioinformatics-training/2017-August-Variant-Analysis-Workshop/master/thursday/chr18.gtf

Change its name to 'genes.gtf':

    mv chr18.gtf genes.gtf

---

**5\.** We are ready to create the database. Go back to the snpEff directory and run the 'build' command for snpeff:

    cd ../..
    java -jar snpEff.jar build -gtf22 -v chr18

Now we're ready to do effect prediction.

---

**6\.** Go back to your '05-snpeff' directory and look at the options for the effect prediction command for snpeff:

    cd ..
    java -jar snpEff/snpEff.jar eff -h

Now run the prediction:

    java -jar snpEff/snpEff.jar eff chr18 all.chr18.vcf > snpeff.chr18.vcf

This will take about 8 minutes to run.

---

**7\.** Take a look at the output file:

    less snpeff.chr18.vcf

Also, download (to your laptop) and take a look at the 'snpEff_summary.html' file. The VCF file is the same file as the input, except every variant has an 'ANN' field added to it. Look at the [snpeff docs](http://snpeff.sourceforge.net/SnpEff_manual.html#input) and [this detailed pdf](http://snpeff.sourceforge.net/VCFannotationformat_v1.0.pdf) to get information about the format of the 'ANN' field.

