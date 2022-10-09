## Assignment 4: Variant Analysis and Bedtools
Assignment Date: Monday, October 10, 2022 <br>
Due Date: Monday, Oct 17, 2022 @ 11:59pm <br>

### Assignment Overview

In this assignment you will consider the fundamental properties for variant calling, and get some experience with small variant analysis.

In this assignment you will explore WDLs as a workflow language to orchestrate a variety of read mapping tasks. You can execute your WDLs using [miniwdl](https://github.com/chanzuckerberg/miniwdl) after running `pip install miniwdl`. You may also find [learn-wdl](https://github.com/openwdl/learn-wdl) to be very helpful.

If you are using a Mac, make sure to disable [gRPC FUSE for file sharing](https://github.com/chanzuckerberg/miniwdl/issues/145). There is a currently a very weird bug in Docker desktop 4.12.0 (docker engine 20.10.17) on Apple M1 where if you try to disable "Use gRPC FUSE for file sharing" via the GUI it will turn itself back on when you try to activate it. Docker is aware of the problem and are working on a fix. In the meantime there is a workaround available: [https://github.com/docker/for-mac/issues/6467](https://github.com/docker/for-mac/issues/6467)

The bioinformatics tools you will need for the assignment (bowtie, samtools, etc) are bundled into a [Docker](https://www.docker.com) container so you wont need to install other software packages. This will get you ready to run your code in the cloud for future assignments. Instructions for running the container are available here: [https://github.com/mschatz/wga-essentials](https://github.com/mschatz/wga-essentials). This is known to work on new macs (M1) and older macs (intel chip). It should also run on Linux and Windows but let us know if you have any issues.

As a reminder, any questions about the assignment should be posted to [Piazza](https://piazza.com/class/l7dg3c82ftw1d/).


### Question 1. Binomial Distribution [10 pts]

- 1a. For coverage n = 10 to 200, calculate the maximum number of minor allele reads (round down) that would make your one-sided binomial test reject the null hypothesis p=0.5 at 0.05 significance. Plot coverage on the x-axis and (number of reads)/(coverage) on the y-axis. Note this is the minimum number of reads that are necessary to believe we might have a heterozygous variant on the second haplotype rather than just mere sequencing error.

- 1b. What asymptote does the plot seem to approach? Why is this?


### Question 2. Hello World on AnVIL [5 pts]

- Using AnVIL, run the [Hello World WDL](https://github.com/openwdl/learn-wdl/blob/master/1_script_examples/1_hello_worlds/3_hello_input_task/input-task.wdl) to print out "Hello, YOUR_NAME". Please submit your code as well as a screenshot of the window in which you run the code. Note for this analysis you should sign up for [$300 worth of free credits](https://support.terra.bio/hc/en-us/articles/360046295092). 



### Question 3. Identifying a Specific Variant [25 pts]

*Three friends, Ezra, Sabine, and Zeb go out to lunch together and, in discussing the menu, all realize that they all have a specific culinary preference. They believe the cause of this preference to be genetic. Can you identify what SNP causes this preference and what the preference is?*

Download the read set from here: [https://github.com/schatzlab/appliedgenomics2022/tree/main/assignments/assignment4/input_files](https://github.com/schatzlab/appliedgenomics2022/tree/main/assignments/assignment4/input_files)

For this question, you may find this tutorial helpful: [http://clavius.bc.edu/~erik/CSHL-advanced-sequencing/freebayes-tutorial.html](http://clavius.bc.edu/~erik/CSHL-advanced-sequencing/freebayes-tutorial.html)

To answer the following questions, use a WDL file to run the required tools. You should use one task per computational question. Please include your WDL file in your submission. You can use your local computer or AnVIL to run the WDL.

- 3a. You are able to narrow down the SNP to a specific region in chromosome 11. Using bowtie2, align their FASTQ files to the provided reference file. How many reads does each file have and how many are successfully mapped? [Hint: Use `samtools flagstat`.] 

- 3b. How many SNPs and indels does each file have? [Hint: Sort the SAM file first. Then, call variants with `freebayes`. Summarize using `bcftools stats`.]

- 3c. You now know which SNPs and indels each friend has. However, you want to know which variant they all share. How many variants are shared between all three friends? [Hint: Use `bcftools isec`.]

- 3d. Between the variants shared between all three friends, which is likeliest to cause a phenotype of interest? [Hint: You can search for variants at a certain chromosome and position at https://www.ncbi.nlm.nih.gov/snp/advanced/. Remember, the position in the intersected VCF is the position within the region we're looking at, so you will have to find the starting location of the region!]

- 3e. What is the phenotype? [Hint: Search the name of the gene associated with the variant you found in 3d.]



#### Question 4. De novo mutation analysis [20 pts]

For this question, we will be focusing on the de novo variants identified in this paper:<br>
[http://www.nature.com/articles/npjgenmed201627](http://www.nature.com/articles/npjgenmed201627)

Download the de novo variant positions from here (Supplementary Table S4):<br>
[https://github.com/schatzlab/appliedgenomics2022/blob/main/assignments/assignment4/41525_2016_BFnpjgenmed201627_MOESM431_ESM.xlsx?raw=true](https://github.com/schatzlab/appliedgenomics2022/blob/main/assignments/assignment4/41525_2016_BFnpjgenmed201627_MOESM431_ESM.xlsx?raw=true)

Download the gene annotation of the human genome here: <br>
[ftp://ftp.ensembl.org/pub/release-87/gff3/homo_sapiens/Homo_sapiens.GRCh38.87.gff3.gz](ftp://ftp.ensembl.org/pub/release-87/gff3/homo_sapiens/Homo_sapiens.GRCh38.87.gff3.gz)

Download the annotation of regulatory variants from here:<br>
[ftp://ftp.ensembl.org/pub/release-87/regulation/homo_sapiens/homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20161111.gff.gz](ftp://ftp.ensembl.org/pub/release-87/regulation/homo_sapiens/homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20161111.gff.gz)

**NOTE** The variants are reported using version 37 of the reference genome, but the annotation is for version 38. Fortunately, you can 'lift-over' the variants to the coordinates on the new reference genome using several avaible tools. I recommmend the [UCSC liftover tool](https://genome.ucsc.edu/cgi-bin/hgLiftOver) that can do this in batch by converting the variants into BED format. Note, some variants may not successfully lift over, especially if they become repetitive and/or missing in the new reference, so please make a note of how many variants fail liftover.

- 4a. How much of the genome is annotated as a gene?

- 4b. What is the sequence of the shortest gene on chromosome 22? [Hint: `bedtools getfasta`]

- 4c. How much of the genome is an annotated regulatory sequence (any type)? [Hint `bedtools merge`]

- 4d. How much of the genome is neither gene nor regulatory sequences? [Hint: `bedtools merge` + `bedtools subtract`]

- 4e. Using the [UCSC liftover tool](https://genome.ucsc.edu/cgi-bin/hgLiftOver), how many of the variants can be successfully lifted over to hg38?

- 4f. How many variants are in genes? [Hint: convert xlsx to BED, then `bedtools`]

- 4g. How many variants are in *any* annotated regulatory regions? [Hint: `bedtools`]

- 4h. What type of annotated regulatory region has the most variants? [Hint: `bedtools`]

- 4i. Is this a statistically significant number of variants (P-value < 0.05)? [Hint: If you don't want to calculate this analytically, you can do an experiment. Try simulating the same number of variants as the original file 100 times, and see how many fall into this regulatory type. If at least this many variants fall into this feature type more than 5% of the trials, this is not statistically significant]


### Packaging

The solutions to the above questions should be submitted as a single PDF document that includes your name, email address, and all relevant figures (as needed). Submit your solutions by uploading the PDF to [GradeScope](https://www.gradescope.com/courses/431817), and remember to select where in your submission each question/subquestion is. The Entry Code is: J37JKW. 

If you submit after this time, you will use your late days. Remember, you are only allowed 4 late days for the entire semester!



## Resources

#### [FreeBayes](https://github.com/ekg/freebayes) - Small variant identification

```
conda install freebayes
freebayes -f chr22.fa sample.bam > sample.vcf
```

#### [bcftools](https://samtools.github.io/bcftools/bcftools.html) - VCF summary

```
conda install bcftools
bcftools stats sample.vcf > stats.txt
```

#### [BEDTools](http://bedtools.readthedocs.io/en/latest/) - Genome arithmetic

Get bedtools from conda, if you haven't already.

```
conda install bedtools
```
