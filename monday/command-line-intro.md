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
    srun -t 1440 -c 4 -n 1 --mem 8000 --reservation workshop --pty /bin/bish

You'll see, eventually:

    srun: job ####### has been allocated resources

Multiple ways to clear the clutter:

    <enter> <enter> <enter>
    <ctrl-l>
    <ctrl-k?>
    clear
    exit  # kills current shell!
    # Note that any text *following* a "#" symbol is ignored.


Command Line Basics
--------------------

    pwd  # present working directory ... where am I?
    ls   # list files here
    ls /tmp/  # list files somewhere else
    sleep 1000  # wait for 1000 seconds!
    <ctrl-c>  # shows as ^C
    command fubar <ctrl-c>  # exits command entry without executing
    python  # enter Python interactive session
    <ctrl-d>  # escape from some repl's
    R  # enter R interactive session
    <ctrl-d>
    yes  # that's a lot of yes
    <ctrl-c>
    yes | more  # "pipe" yes' output to the "more" command
    <spacebar>
    <spacebar>  # next page
    <q>  # quits from more, less, "man" pages, etc.
    yes | less  # pipe to "less" paginator, instead of "more"
    <spacebar>
    <arrow keys, pgup, pgdn>  # forward or back through file
    g *or* G  # beginning or end of file
    /es  # "/" enters search mode, "es" is pattern sought for
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


Wildcard Characters and Find
-----------------------------

    ls ?hiX/Illumina  # list files in Illumina sub-directory of any directory ending in "hiX"
    ls PhiX/Illumina/RTA/Sequence/*/*.fa  # list all .fa files a few directories down
    # So, "?" fills in for zero or one character, "*" fills in for zero or more characters
    find . -name "*.fa"
    find . -name "*.f?"  # how is this different from the previous command?


![name of image](URL or folder/filename.png)



