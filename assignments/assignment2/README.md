## Assignment 2: Genome Assembly
Assignment Date: Monday, September 12, 2022 <br>
Due Date: Monday, September 19, 2022 @ 11:59pm <br>

### Assignment Overview

In this assignment, you are given a set of unassembled reads from a mysterious pathogen that contains a secret message encoded someplace in the genome. The secret message will be recognizable as a novel insertion of sequence not found in the reference. Your task is to assess the quality of the reads, assemble the genome, identify, and decode the secret message. If all goes well the secret message should decode into a recognizable 
english text, otherwise doublecheck your coordinates and try again. As a reminder, any questions about the assignment 
should be posted to [Piazza](https://piazza.com/class/l7dg3c82ftw1d/).

For this assignment, we recommend you install and run the tools using [bioconda](https://www.nature.com/articles/s41592-018-0046-7). There are some tips below in the Resources section. Alternatively, you can try running the tools using [Docker](https://www.docker.com/). Docker is a powerful containerization tool to make software easier to distribute. This will become even more important later in the semester when you are running your assignments in the cloud. For this assignment, you can use this docker instance that has these tools preinstalled. ***Note: this does not work on new Apple M1 chips :(:*** 
[https://github.com/mschatz/wga-essentials](https://github.com/mschatz/wga-essentials)


#### Question 1. Coverage Analysis [20 pts]

Download the reads and reference genome from: [https://github.com/schatzlab/appliedgenomics2022/blob/main/assignments/assignment2/asm.tgz?raw=true](https://github.com/schatzlab/appliedgenomics2022/blob/main/assignments/assignment2/asm.tgz?raw=true)

Note we have provided both paired-end and mate-pairs reads (see included README for details). Make sure to look at all of the reads for the coverage analysis and kmer analysis, as well as in the assembly.

- Question 1a. How long is the reference genome? [Hint: Try `samtools faidx`]
- Question 1b. How many reads are provided and how long are they? Make sure to measure each file separately [Hint: Try `FastQC`]
- Question 1c. How much coverage do you expect to have? [Hint: A little arthmetic]
- Question 1d. Plot the average quality value across the length of the reads [Hint: Screenshot from `FastQC`]

#### Question 2. Kmer Analysis [20 pts]

Use `Jellyfish` to count the 21-mers in the reads data. Make sure to use the "-C" flag to count cannonical kmers, 
otherwise your analysis will not correctly account for the fact that your reads come from either strand of DNA.

- Question 2a. How many kmers occur exactly 50 times? [Hint: try `jellyfish histo`]
- Question 2b. What are the top 10 most frequently occurring kmers [Hint: try `jellyfish dump` along with `sort` and `head`]
- Question 2c. What is the estimated genome size based on the kmer frequencies? [Hint: upload the jellyfish histogram to [GenomeScope2](http://qb.cshl.edu/genomescope/genomescope2.0/) and report the min "Genome Haploid Length" in the "Results" section]
- Question 2d. How well does the GenomeScope genome size estimate compare to the reference genome? [Hint: In a sentence or two]

#### Question 3. De novo assembly [20 pts]

Assemble the reads using `Spades`. Spades will *not* run on Windows, you must use a linux/mac environment. 

- Question 3a. How many contigs were produced? [Hint: try `grep -c '>' contigs.fasta`]
- Question 3b. What is the total length of the contigs? [Hint: try `samtools faidx`, plus a short script/excel]
- Question 3c. What is the size of your large contig? [Hint: check `samtools faidx` plus `sort -n`]
- Question 3d. What is the contig N50 size? [Hint: Write a short script, or use excel]

#### Question 4. Whole Genome Alignment [15 pts]

- Question 4a. What is the average identity of your assembly compared to the reference? [Hint: try `dnadiff`]
- Question 4b. What is the length of the longest alignment [Hint: try `nucmer` and `show-coords`]
- Question 4c. How many insertions and deletions are in the assembly? [Hint: try `dnadiff`]

#### Question 5. Decoding the insertion [20 pts]
- Question 5a. What is the position of the insertion on the reference? [Hint: try `show-coords`]
- Question 5b. How long is the novel insertion? [Hint: try `show-coords`]
- Question 5c. What is the DNA sequence of the encoded message? [Hint: try `samtools faidx` to extract the insertion]
- Question 5d. What is the secret message? [Hint: run `dna-encode.pl -d` to decode the string from 5c]


### Packaging

The solutions to the above questions should be submitted as a single PDF document that includes your name, email address, and 
all relevant figures (as needed). Submit your solutions by uploading the PDF to [GradeScope](https://www.gradescope.com/courses/431817), and remember to select where in your submission each question/subquestion is. The Entry Code is: J37JKW. 

If you submit after this time, you will use your late days. Remember, you are only allowed 4 late days for the entire semester!



### Resources


#### [Bioconda](https://bioconda.github.io/) - Package manager for bioinformatics software

On linux or mac I *highly* recommend that you use bioconda to install the packages rather than installing from source. 

The easiest way to install conda is with [Miniconda](https://docs.conda.io/en/latest/miniconda.html). For M1 macs, I recommend using the x86 installation in emulation mode since  M1/arm support is still limited. I also recommend using [mamba](https://github.com/mamba-org/mamba) instead of the default `conda` command for installing new packages:

```
## Replace Linux-x86_64 with the version you downloaded!
$ chmod +x Miniconda2-latest-Linux-x86_64.sh
$ ./Miniconda2-latest-Linux-x86_64.sh

## After conda is installed add some default channels
$ conda config --add channels bioconda
$ conda config --add channels conda-forge
$ conda install -c conda-forge mamba
```

Once bioconda is configured, all of the tools needed for this assignment can be installed using:

```
$ conda install fastqc jellyfish spades mummer samtools
```


#### [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) - Raw read quality assessment

```
$ fastqc /path/to/reads.fq
```

If you have problems, make sure java is installed (`sudo apt-get install default-jre`)


#### [Jellyfish](http://www.genome.umd.edu/jellyfish.html) - Fast Kmer Counting

When counting kmers, make sure to count "canonical" kmers on both strands (-C):

```
$ jellyfish count -m 21 -C -s 1000000 /path/to/reads*.fq
$ jellyfish histo mer_counts.jf > reads.histo
```

#### [GenomeScope](http://www.genomescope.org/) - Analyze Kmer Profile to determine genome size and other properties

GenomeScope is a web-based tool so there is nothing to install. Hooray! Just make sure to use the `-C` when running jellyfish count so that the reads are correctly processed.

####  [Spades](http://cab.spbu.ru/software/spades/) - Short Read Assembler. 

Normally spades would try several values of k and merge the results together, but here we will force it to just use k=31 to save time. The assembly should take a few minutes.

```
$ spades.py --pe1-1 frag180.1.fq --pe1-2 frag180.2.fq --mp1-1 jump2k.1.fq --mp1-2 jump2k.2.fq -o asm -t 4 -k 31
```


***Note: On mac you may need to run spades like this with the -m flag:***
```
$ spades.py -m 1024 --pe1-1 frag180.1.fq --pe1-2 frag180.2.fq --mp1-1 jump2k.1.fq --mp1-2 jump2k.2.fq -o asm -t 4 -k 31
```


#### [MUMmer](http://mummer.sourceforge.net/) - Whole Genome Alignment

```
$ dnadiff /path/to/ref.fa /path/to/qry.fa
$ nucmer /path/to/ref.fa /path/to/qry.fa
$ show-coords out.delta
```

**WARNING: nucmer and related tools do not like it if/when you have spaces or special characters ('@') in the path to the binaries***

If you have problems with mummer not running perl, reinstall it like this:

```
$ mamba install mummer=3.23=h6de7cb9_11
```


#### [SAMTools](http://www.htslib.org/) - Extract part of a genome sequence using 'samtools faidx' (this will extract from contig_id bases 1234 through 5678)

```
$ ./samtools faidx /path/to/genome.fa contig_id:1234-5678
```



