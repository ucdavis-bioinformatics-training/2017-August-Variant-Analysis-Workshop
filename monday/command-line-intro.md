Intro to the Command Line
==========================

Joe Fass

jnfass@ucdavis.edu


Download the slides [here](CLintro.pdf).


Getting There
--------------------------------

Secure SHell ... ssh.

    ssh [class#@]cabernet.genomecenter.ucdavis.edu
    # DON't WORK ON THE HEAD NODE!
    srun -t 1440 -c 4 -n 1 --mem 8000 --reservation workshop --pty /bin/bash

You'll see, eventually:

    srun: job ####### has been allocated resources

Multiple ways to clear the clutter, and have an empty screen:

    <enter> <enter> <enter>
    <ctrl-l>
    <ctrl-k?>
    clear

And if you're really done working on the command line:

    exit  # kills current shell!
    # Note that any text *following* a "#" symbol is ignored.
    # Note that, in this case, the first 'exit' will exit the srun job
    # ... you'd need another 'exit' to log out of cabernet

Command Line Basics
--------------------

    pwd  # present working directory ... where am I?
    ls   # list files here ... you should see nothing since your 'class##' homes are empty
    ls /tmp/  # list files somewhere else
    sleep 1000  # wait for 1000 seconds!
    <ctrl-c>  # shows as '^C'; exits command entry without executing
    # In some cases, a different key sequence is required
    python  # enter Python interactive session
    <ctrl-d>  # escape from some repl's ... "Read, Execute, & Print Loop"
    R  # enter R interactive session
    <ctrl-d>
    # Try this fun command:
    yes  # that's a lot of yes
    <ctrl-c>
    yes | more  # "pipe," or direct, the yes output to the 'more' command
    <spacebar>
    <spacebar>  # next page
    <q>  # quits from more, less, 'man' pages, etc.
    yes | less  # pipe to 'less' paginator, instead of 'more'
    <spacebar>
    <arrow keys, pgup, pgdn>  # forward or back through file
    g *or* G  # beginning or end of file
    /es  # '/' enters search mode, "es" is pattern sought for (could be any string)
    /s  
    /no
    # After searching for a pattern, using the '/' key, you can use:
    n *or* N  # next or previous pattern match


You've Got Options
-------------------

    ls -R /software
    # ack! too much going to the screen!
    <ctrl-c>
    ls -R /software/scythe  # lists directories and files *recursively*
    # how do I know which options can do what?
    man ls
    # navigate like in "less" (up/down,pgup/dn,g,G,/pattern,n,N,q)
    # try the following:
    ls -l
    ls -a
    ls -l -a
    ls -la  # option "smushing"
    ls -ltrha
    ls -ltrha --color  # single letter (smushed) vs word options
    # what if I want to use same options repeatedly? and be lazy?
    alias  # lists current *command aliases* or *shortcuts*
    alias ll='ls -ltrhaF --color'
    ll
    alias l='ls --color'


Getting Around
----------------

    cd  # no effect? that's because it sends you home (to ~)
    cd /  # go to root of tree's root system
    cd home  # go to where everyone's homes are
    pwd
    cd class10  # use your actual home, not class10
    pwd
    cd /
    pwd
    cd ~  # a shortcut to home, from anywhere
    pwd
    cd .  # "." always means *this* directory
    pwd
    cd ..  # ".." always means *one directory up*
    pwd

Don't get confused between the "." directory name and filenames that start with the "." character. The latter are just valid filenames or directory names, but are usually hidden (use ls's "-a" option to see hidden files).


Absolute and Relative Paths
----------------------------

    cd ~  # let's start at home
    cd ../../home/class10/  # *relative* (start here, take steps)
    pwd
    cd /home/class10/  # *absolute* (start at root, take steps)
    pwd


Tab Completion - A Real Tendon-Saver
-------------------------------------

Using tab-completion will literally save your life. Hours of it. A single <tab> auto-completes file or directory names when there's only one name that could be completed correctly. If multiple files could satisfy the tab-completion, then nothing will happen after the first <tab>. In this case, press <tab> a second time to list all the possible completing names. Note that if you've already made a mistake, and no files will ever be completed correctly, then <tab>'s will do nothing.

    touch one seven september  # create three empty files using "touch" command
    cat o<tab><no enter>  # will complete to "one"
    <enter>
    cat s<tab><no enter>  # nothing!
    <tab><no enter>  # should list seven and september
    ev<tab><no enter>  # should complete to seven
    <enter>  # runs "cat seven" command
    # we often literally autocomplete commands letter by letter
    # comes in handy if we don't exactly remember what name we want
    ls /hom<tab>j<tab><tab>f<tab>  # completes to my home directory, /home/jfass/


Create and Destroy
-------------------

    cd  # another way to return to your home
    mkdir temp
    cd temp/
    echo "Hello, world!" > first.txt
    file first.txt  # tells us what it is
    cat first.txt
    cd ../
    rmdir temp  # shouldn't work!
    rm temp/first.txt  # clear directory first
    rmdir temp


Piping and Redirection
-----------------------

    mkdir CLB
    cd CLB/
    echo "first" > test.txt
    cat test.txt
    echo "second" > test.txt
    cat test.txt
    echo "third" >> test.txt
    cat test.txt

The ">" character redirects output of a command that would normally go to the screen instead into a specified file. ">" replaces, ">>" appends.

    cut -c 1-3 test.txt  # cuts three character columns from file
    cat test.txt | cut -c 1-3  # same thing, piping output of one command into input of another
    cat test.txt | cut -c 1-3 | sort -r
    cat test.txt | cut -c 1-3 | sort -r | grep s
    # pipes cat to cut to sort (-r means reverse order sort, grep searches for pattern matches)

This is a great way to build up a set of operations while inspecting the output of each step in turn.


History Repeats Itself
-----------------------

    <up>  # last command
    <up>  # next-to-last command
    <down>  # last command, again
    <down>  # current command, empty or otherwise
    # for a more global view of your history:
    history
    history | head
    history | tail
    history | tail -n 30
    history | less
    cat test.txt | cut -c 1-3 | sort -r | grep -s > reallyImportantResult.txt
    history | tail
    !560  # re-executes 560th command (yours will have different numbers)
    <ctrl-r>fir  # should find most recent command with "fir" string ... 'echo "first" > test.txt'?
    <enter>  # to run command
    <ctrl-c>  # get out of recursive search
    <ctr-r>  # repeat <ctrl-r> to find successively older previous string matches


Editing Yourself
-----------------

    <up><up>  # go to some previous command
    <ctrl-a>  # go to the beginning of the line!
    <ctrl-e>  # go to the end of the line!
    <left><right>  # move to a single word (surrounded by whitespace)
    <ctrl-k>  # delete from here to end of line
    <ctrl-w>  # delete from here to beginning of preceeding word
    blah blah blah<ctrl-w><ctrl-w>  # leaves you with only one "blah"


Compression
-------------

    gzip test.txt
    file test.txt
    gunzip -c test.txt | more
    ll
    gunzip test.txt
    bzip2 test.txt; bunzip2 test.txt.bz2


Archives
--------

Tape archives, or .tar files, that are further compressed are called "tarballs." Let's grab one.

    wget ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com/PhiX/Illumina/RTA/PhiX_Illumina_RTA.tar.gz
    # ... or ...
    curl ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com/PhiX/Illumina/RTA/PhiX_Illumina_RTA.tar.gz > PhiX.tgz
    # .tar.gz and .tgz are *commonly used* extensions for compressed tar files
    tar -xzvf PhiX_Illumina_RTA.tar.gz  # or whatever you called it, if you used curl
    # -x = extract, -z = use gzip/unzip, -v = verbose (show each file in archive), -f = use a file, not a tape drive(!)
    

Forced Removal
---------------

This gets a heading all its own. Because when you're on the command line, there's no "Recycle Bin". 

    rm -rf PhiX  # be sure ... there's no going back!
    # -r = recursively remove sub-directories, -f means *force* (auto-"yes")
    # We actually want to use those directories, so un-archive them again!


BASH Wildcard Characters and Find
----------------------------------

    ls ?hiX/Illumina  # list files in Illumina sub-directory of any directory ending in "hiX"
    ls PhiX/Illumina/RTA/Sequence/*/*.fa  # list all .fa files a few directories down
    # So, "?" fills in for zero or one character, "*" fills in for zero or more characters
    find . -name "*.fa"
    find . -name "*.f?"  # how is this different from the previous command?


Manipulation of a FASTA File
-----------------------------

We just found the phiX-174 genome, so let's copy it to our current directory so we can play with it:

    cp ./PhiX/Illumina/RTA/Sequence/WholeGenomeFasta/genome.fa phix.fa
    # Note how we copied the 'genome.fa' file to a different name: 'phix.fa'
    wc -l phix.fa

We can use the 'grep' command to search for matches to patterns (more flexibly than by using less's '/' key). 'grep' comes from "**g**lobally search for a **r**egular **e**xpression and **p**rint." 

    grep -c ">" phix.fa  # only one FASTA sequence entry, since only one header line (">gi|somethingsomething...")
    cat phix.fa  # this may not be useful for anything larger than a virus!
    # let's look at start codon and 2 following:
    grep "ATG......" phix.fa  # '.' characters are the single-character wildcards for grep
    # use the '-o' option to **o*nly print the pattern matches, one per line
    grep -o "ATG......" phix.fa
    # use the 'cut' command with '-c' to select characters 4-6, the second codon
    grep -o "ATG......" phix.fa | cut -c4-6
    # 'sort' the second codon sequences (default order is same as ASCII table; see 'man ascii')
    grep -o "ATG......" phix.fa | cut -c4-6 | sort
    # combine successive identical sequences, but count them ('-c' option)
    grep -o "ATG......" phix.fa | cut -c4-6 | sort | uniq -c
    # and finally sort using only the 1st "field" as a key ('-k1,1'), in reverse numeric order ('-rn')
    grep -o "ATG......" phix.fa | cut -c4-6 | sort | uniq -c | sort -rn -k1,1
    # ... which gives us the most common codons first

This may not be a particularly useful thing to do with a genomic FASTA file, but it illustrates the process by which one can build up a string of operations, using pipes, in order to ask quantitative questions about sequence content. More generally than that, this process allows one to ask questions about files and file contents and the operating system, and verify at each step that the process so far is working as expected. The command line is, in this sense, really a modular workflow management system.


Symbolic Links
---------------

Since copying or even moving large files (like sequence data) around your filesystem may be impractical, we can use links to reference "distant" files without duplicating the data in the files. Symbolic links are disposable pointers that refer to other files, but behave like the referenced files in commands.

    ln -s PhiX/Illumina/RTA/Sequence/WholeGenomeFasta/genome.fa .
    ll  # notice the symbolic link pointing at its target
    grep -c ">" genome.fa


Bioinformatics, At Last! ... ?
-------------------------------

OK, let's try to do some sequence alignment (similar to a BLAST search).

    bwa mem genome.fa phix.fa > aln.sam
    # We get the following:
    # [E::bwa_idx_load] fail to locate the index files

What went wrong? We redirected output to the 'aln.sam' file using the '>' character. But there's nothing in the output file:

    cat aln.sam
    # <nothing>


STDOUT & STDERR
-----------------

Programs can write to two separate output streams, "standard out" (STDOUT), and "standard error" (STDERR). The former is generally for direct output of a program, while the latter is supposed to be used for reporting problems. I've seen some bioinformatics tools use STDERR to report summary statistics about the output, but this is probably bad practice. Default behavior in a lot of cases is to dump both STDOUT and STDERR to the screen, unless you specify otherwise. In order to nail down what goes where, and record it for posterity:

    bwa mem genome.fa phix.fa 1> aln.sam 2> aln.err
    # the 1st output, STDOUT, goes to 'aln.sam'
    # the 2nd output, STDERR, goes to 'aln.err'
    cat aln.sam
    # <nothing>
    cat aln.err
    # [E::bwa_idx_load] fail to locate the index files

That didn't really help us with our problem, but at least you can direct output where you want it now. Saving STDOUT is pretty routine (you want your results, yes?), but remember that explicitly saving STDERR is important in a computing cluster, since you may not directly see the "screen" when you're running jobs. And even if you do, it could be a mix of STDERR from multiple jobs (if it's dumped to screen), so explicitly naming each STDERR file lets you diagnose what went wrong with which set of data, after everything's run or failed to run.

For the curious, our problem was that we didn't index the "target" or "reference" of our sequence alignment. Try this:

    bwa index genome.fa
    bwa mem genome.fa phix.fa 1> aln.sam 2> aln.err
    cat aln.err
    # see? summary of the alignment process / results
    cat aln.sam

The resulting SAM file shows a perfect match between the query sequence ('phix.fa') and the reference ('genome.fa') ... which is good, because they're the same file (see the symbolic links section above)! Beyond that, don't worry about the SAM format; we'll get into that tomorrow.


Loops
-------

Loops are useful for quickly telling the shell to perform one operation after another, in series. For example:

    for i in {1..21}; do echo $i >> a; done
    cat a
    # <1 through 21 on separate lines>

The general form is:

    for name in {list}; do
        commands
    done

The list can be a sequence of numbers or letters, or a group of files specified with wildcard characters:

    # imagine you have 20 sequence files, in a 'fastqs' directory:
    bwa index reference.fa
    for sample in fastqs/\*.fastq; do
        bwa mem reference.fa $sample 1> $sample.sam 2> $sample.err
    done
    # this would produce, for example, ./fastqs/sample1.fastq.sam and ./fastqs/sample1.fastq.err, etc.

Sometimes a "while" loop is more convenient than a "for" loop ... if you don't readily know how many iterations of the loop you want:

    while {condition}; do
        commands
    done

Or, imagining a file that contains the filenames (one per line) of samples' sequence data:

    cat file-of-filenames.txt | while read sample; do
        bwa mem reference.fa $sample 1> $sample.sam 2> $sample.err
    done


Paste Command
--------------

The paste command is useful in creating tables.

    echo "WT1" >> b
    echo "WT2" >> b
    echo "control1" >> b
    echo "control2" >> b
    # now we can number our four samples to conveniently refer to them in order
    for i in {1..4}; do echo $i >> a; done  # semicolons terminate separate lines of multi-line commands like loops
    paste a b > c
    cat c


Running in the Background
--------------------------

Sometimes it's useful to continue working on the command line, even after you've executed a command that's going to take a while to finish. Normally this command would occupy the shell, and prevent you from typing in commands and receiving results. But we can "put jobs in the background" so that they don't occupy your shell directly:

    sleep 1000000
    <control-z>  # shows as '^Z'
    # [1]+  Stopped                 sleep 1000000
    bg
    # [1]+ sleep 1000000 &

'^Z' first suspends the sleep command (or whatever command is running when you hit <control-z>). Then, 'bg' resumes running that command *in the background*, so that it doesn't occupy the terminal. The output of the 'bg' command tells you that you have one command running in the background. You could start more, suspend them, then resume them in the background, and query what background jobs are running or are suspended, not running:

    jobs
    # [1]+  Running                 sleep 1000000 &
    
We can also start a job in the background in one step, without having to suspend then resume it, using the '&' character at the end of the command:

    sleep 5000000 &

If we want to delete these jobs for any reason, we can kill them using the numbering that 'jobs' reveals:

    jobs
    # [1]-  Running                 sleep 1000000 &
    # [2]+  Running                 sleep 5000000 &
    kill %1
    jobs
    # [1]-  Terminated              sleep 1000000
    # [2]+  Running                 sleep 5000000 &
    kill %2
    jobs
    # [2]+  Terminated              sleep 5000000

Finally, the 'nohup' command (from "no hangup"!) makes jobs extra resistant to lost connections or terminal problems. In other words, even jobs running in the background can be terminated if one's shell dies. 'nohup' separates the running job from the shell, so it'll keep running until it dies or is actively killed by you.

    nohup sleep 1000000 &
    # [1] 34993
    # class##@c4-0:~/CLB$ nohup: ignoring input and appending output to ‘nohup.out’
    jobs
    # [1]+  Running                 nohup sleep 1000000 &
    # output is dumped into the 'nohup.out' file unless specifically redirected in your command
    kill %1


Table of Processes
-------------------

The 'top' command prints a self-updating table of running processes and system stats. Use 'q' to exit top, 'z' to toggle better color contrast, 'M' to sort by memory use, 'P' to sort by processor use, and 'c' to toggle display of the full commands. Hit '1' to toggle display of all processors, and hit 'u' followed by typing in a username in order to only show processes (jobs) owned by that user.

![top-example](top.png)

Shell Scripts, File Permissions
--------------------------------

Often it's useful to define a whole string of commands to run on some input, so that (1) you can be sure you're running the same commands on all data, and (2) you don't have to type the same commands in over and over! Let's use the 'nano' text editor program that's pretty reliably installed on most linux systems (MacOS??).

    nano test.sh
    # nano now occupies the whole screen; see commands at the bottom
    # type/paste in the following ...
    # (note that '#!' is an interpreted command to the shell, not a comment)
    # (also note that you may have to add and delete spaces and <enter>s to get the spacing right, though that's not critical)
    #!/bin/bash
    grep -o . $1 | \
        sort | \
        uniq -c | \
        sort -rn -k1,1
    <control-o><control-x>  # to write **o**ut to test.sh, and then e**x**it nano

Note that '$1' means 'the value of the 1st argument to the shell script' ... in other words, the text that follows the shell script name when we run it (see below).

Though there are ways to run the commands in test.sh right now, it's generally useful to give yourself (and others) "execute" permissions for test.sh, really making it a shell script. Note the characters in the first (left-most) field of the file listing:

    ll test.sh
    # -rw-r--r--  1 jfass biocore   79 Aug 19 15:05 test.sh

The first '-' becomes a 'd' if the "file" is actually a directory. The next three characters represent **r**ead, **w**rite, and e**x**ecute permissions for the file owner (you), followed by three characters for users in the owner's group, followed by three characters for all other users. Run the 'chmod' command to change permissions for the 'test.sh' file, adding execute permissions ('+x') for the user (you) and your group ('ug'):

    chmod ug+x test.sh
    ll test.sh
    # -rwxr-xr-- 1 jfass biocore 79 Aug 19 15:05 test.sh*

OK! So let's run this script, feeding it the phiX genome. When we put the genome file 1st after the name of the script, this filename becomes variable '1', which the script can access by specifying '$1'.

    ./test.sh genome.fa
    #   1686 T
    #   1292 A
    #   1253 G
    #   1155 C
    #      1 x
    #      1 p
    #      1 i
    #      1 h
    #      1 >

The script's grep command splits out every character in the file on a separate line, then sorts them so it can count the occurrences of every unique character and show the most frequent characters first ... a quick and dirty way to get at GC content.


Installing (Simple) Software
-------------------------------

Let's install a straightforward tool ... another read aligner by the author of BWA that's intended for longer, lower accuracy reads (such as from the PacBio or Oxford Nanopore sequencers): minimap2.

    cd  # returns you to your home directory, since we've been working in the 'CLB' directory
    mkdir tools
    cd tools/
    git clone https://github.com/lh3/minimap2.git
    cd minimap2/
    make
    ./minimap2
    cd

Now you can run this tool (an executable file) by fully specifying the path to your 'tools' directory from wherever you're running the tool. In addition, you could copy, move, or put symbolic links to the tool in a common place, like /usr/bin/, which should be in everybody's path.


