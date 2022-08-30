## Assignment 1: Chromosome Structures
Assignment Date: Wednesday, August 31, 2022 <br>
Due Date: Wednesday, Sept. 7, 2022 @ 11:59pm <br>

### Assignment Overview

In this assignment you will profile the overall structure of the genomes of several important species and then study human chromosome 22 in more detail.
As a reminder, any questions about the assignment should be posted to [Piazza](https://piazza.com/class/l7dg3c82ftw1d/).

### Question 1: Chromosome structures [10 pts]

Download the chomosome size files for the following genomes (Note these have been preprocessed to only include main chromosomes):

1. [Arabidopsis thaliana (TAIR10)](https://schatz-lab.org/appliedgenomics2020/assignments/assignment1/TAIR10.chrom.sizes) - An important plant model species [[info]](https://en.wikipedia.org/wiki/Arabidopsis_thaliana)
2. [Tomato (Solanum lycopersicum v4.00)](https://schatz-lab.org/appliedgenomics2020/assignments/assignment1/tomato.chrom.sizes) - One of the most important food crops [[info]](https://en.wikipedia.org/wiki/Tomato)
3. [E. coli (Escherichia coli K12)](https://schatz-lab.org/appliedgenomics2020/assignments/assignment1/ecoli.chrom.sizes) - One of the most commonly studied bacteria [[info]](https://en.wikipedia.org/wiki/Escherichia_coli)
4. [Fruit Fly (Drosophila melanogaster, dm6)](https://schatz-lab.org/appliedgenomics2020/assignments/assignment1/dm6.chrom.sizes) - One of the most important model species for genetics [[info]](https://en.wikipedia.org/wiki/Drosophila_melanogaster)
5. [Human (hg38)](https://schatz-lab.org/appliedgenomics2020/assignments/assignment1/hg38.chrom.sizes) - us :) [[info]](https://en.wikipedia.org/wiki/Homo_sapiens)
6. [Wheat (Triticum aestivum, IWGSC)](https://schatz-lab.org/appliedgenomics2020/assignments/assignment1/wheat.chrom.sizes) - The food crop which takes up the largest land area [[info]](https://en.wikipedia.org/wiki/Wheat)
7. [Worm (Caenorhabditis elegans, ce10)](https://schatz-lab.org/appliedgenomics2020/assignments/assignment1/ce10.chrom.sizes) - One of the most important animal model species [[info]](https://en.wikipedia.org/wiki/Caenorhabditis_elegans)
8. [Yeast (Saccharomyces cerevisiae, sacCer3)](https://schatz-lab.org/appliedgenomics2020/assignments/assignment1/yeast.chrom.sizes) - an important eukaryotic model species, also good for bread and beer [[info]](https://en.wikipedia.org/wiki/Saccharomyces_cerevisiae)

Using these files, make a table with the following information per species:

- Question 1.1. Total genome size

- Question 1.2. Number of chromosomes

- Question 1.3. Largest chromosome size and name

- Question 1.4. Smallest chromosome size and name

- Question 1.5. Mean chromosome length


### Question 2. Coverage simulator [20 pts]

- Question 2.1. How many 100bp reads are needed to sequence a 1Mbp genome to 5x coverage?

- Question 2.2. In the language of your choice, simulate sequencing 5x coverage of a 1Mbp genome with 100bp reads and plot the histogram of coverage. Note you do not need to actually output the sequences of the reads, you can just randomly sample positions in the genome and record the coverage. You do not need to consider the strand of each read. The start position of each read should have a uniform random probabilty at each possible starting position (1 through 999,901). You can record the coverage in an array of 1M positions. Overlay the histogram with a Poisson distribution with lambda=5

- Question 2.3. Using the histogram from Q2.2, how much of the genome has not been sequenced (has 0x coverage)? How well does this match Poisson expectations?

- Question 2.4. Now repeat the analysis with 15x coverage: 1. simulate the appropriate number of reads, 2. make a histogram, 3. overlay a Poisson distribution with lambda=15, 4. compute the number of bases with 0x coverage, and 5. evaluate how well it matches the Poisson expectation.



### Question 3: Sequence content [20 pts]

Download the human chromosome 22 from here: [https://schatz-lab.org/appliedgenomics2022/assignments/assignment1/chr22.fa.gz](https://schatz-lab.org/appliedgenomics2022/assignments/assignment1/chr22.fa.gz)

- Question 3.1. How many As, Cs, Gs, Ts and Ns are found in the entire chromosome?

- Question 3.2. Divide the chromosome into 100bp non-overlapping windows/bins. How many bins contain at least one N?

- Question 3.3. Make a histogram of the number of bins of a given %GC: x-axis = %GC, y-axis = # bins with this %GC. You should skip any bins that contain Ns.

- Question 3.4. Recall that Illumina sequencing performs poorly when the %GC is <= 30% or >= 65%. Based on what you saw in Q3.3, how many bins (that do not contain any Ns) do we expect to sequence poorly?

- Question 3.5. Assuming the rest of the human genome has the same sequence composition as chromosome 22, how many bases of the human genome do we expect to sequence poorly? Hint: extrapolate the results from 3.4 to the size of the whole genome computed in 1.1. 



### Hints

- Many of the questions can be addressed with standard command line [tools](http://lh3lh3.users.sourceforge.net/biounix.shtml) such as `grep`, `wc`, `awk`, `sort`, `fold`, etc
- You may wish to try out [`datamash`](https://www.gnu.org/software/datamash/)
- You may find [`samtools`](http://www.htslib.org/) and especially `samtools faidx` helpful for indexing the fasta files
- Plotting can be done in any language; R or Python are recommended; Excel is okay but ugly :-P
- We highly recommend you use [Jupyter notebooks](https://jupyter.org/) that you can then "print" to a PDF.
- Be sure to clearly mark each question and subquestion

### Packaging

The solutions to the above questions should be submitted as a single PDF document that includes your name, email address, and 
all relevant figures (as needed). Make sure to clearly label each of the subproblems and give the exact commands and/or code snippets you used for solving the question. Submit your solutions by uploading the PDF to [GradeScope](https://www.gradescope.com/courses/431817), and remember to select where in your submission each question/subquestion is. The Entry Code is: J37JKW. 

If you submit after this time, you will start to use up your late days. Remember, you are only allowed 96 hours (4 days) for the entire semester!



