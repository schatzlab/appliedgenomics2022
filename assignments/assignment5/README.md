## Assignment 5: RNA-seq
Assignment Date: Monday October 24, 2022 <br>
Due Date: Monday October 31 2022 @ 11:59pm <br>

### Assignment Overview

In this assignment you will explore a couple of aspects of RNA-seq (with a small introduction to clustering). For this assignment, you will have to generate some visualizations - we recommend R or Python, but use a language you are comfortable with! 

 **Make sure to show your work/code in your writeup!**

As a reminder, any questions about the assignment should be posted to [Piazza](https://piazza.com/class/l7dg3c82ftw1d/).


#### Question 1. Time Series [20 pts]

[This file](http://schatz-lab.org/teaching/exercises/rnaseq/rnaseq.1.expression/expression.txt) contains normalized expression values for 100 genes over 10 time points. Most genes have a stable background expression level, but some special genes show increased
expression over the timecourse and some show decreased expression.

- Question 1a. Cluster the genes using an algorithm of your choice. Which genes show increasing expression and which genes show decreasing expression, and how did you determine this? What is the background expression level (numerical value) and how did you determine this? [Hint: K-means and hierarchical clustering are common clustering algorithms you could try.]

- Question 1b. Calculate the first two principal components of the expression matrix. Show the plot and color the points based on their cluster from part (a). Does the PC1 axis, PC2 axis, neither, or both correspond to the clustering?

- Question 1c. Create a heatmap of the expression matrix. Order the genes by cluster, but keep the time points in numerical order.

- Question 1d. Visualize the expression data using t-SNE.

- Question 1e. Using the same data, visualize the expression data using UMAP.

- Question 1f. In a few sentences, compare and contrast the (1) heatmap, (2) PCA, (3) t-SNE and (4) UMAP results. Be sure to comment on understandability, relative positioning of clusters, runtime, and any other significant factors that you see.


#### Question 2. Sampling Simulation [10 pts]

A typical human cell has ~250,000 transcripts, and a typical bulk RNA-seq experiment may involve millions of cells. Consequently in an RNAseq experiment you may start with trillions of RNA molecules, although your sequencer will only give a few tens of millions of reads. 
Therefore your RNAseq experiment will be a small sampling of the full composition. We hope the sequences will be a representative sample of the total population, but if your sample is very unlucky or biased it may not represent the true distribution. We will explore this concept by sampling a small subset of transcripts (500 to 50000) out of a much larger set (1M) so that you can evaluate this bias.

In [data1.txt](data1.txt) with 100,000 lines we provide an abstraction of RNA-seq data where normalization has been performed and the number of times a gene name occurs corresponds to the number of transcripts in the sample.

- Question 2a. Randomly sample 500 rows. Do this simulation 10 times and record the relative abundance of each of the 15 genes. Make a scatterplot the mean vs. variance of each gene (x-axis=mean of gene_i, y-axis=variance of gene_i)

- Question 2b. Do the same sampling experiment but sample 5000 rows each time. Again plot the mean vs. variance.

- Question 2c. Do the same sampling experiment but sample 50000 rows each time. Again plot the mean vs. variance.

- Question 2d. Is the variance greater in (a), (b) or (c)? What is the relationship between mean abundance and variance? Why?


#### Question 3. Differential Expression [20 pts]

- Question 3a. Using the file from question 2 (data1.txt) along with [data2.txt](data2.txt), randomly sample 5000 rows from each file. Sample 3 times for each file (this emulates making experimental replicates) and conduct a paired t-test for differential expression of each of the 15 genes. Which genes are significantly differentially expressed at the 0.05 level and what is their mean fold change?

- Question 3b. Make a volano plot of the data from part a: x-axis=log2(fold change of the mean expression of gene_i); y-axis=-log_10(p_value comparing the expression of gene_i). Label all of the genes that show a statistically siginificant change

- Question 3c. Now sample 5000 rows 10 times from each file, equivalent to making more replicates. Which genes are now significant at the 0.05 level and what is their mean fold change?

- Question 3d. Make a volcano plot using the results from part c (label any statistically significant genes)

- Question 3e. Perform the simulations from parts a/c but sample 50000 rows each time from each file. Which genes are significant and what is their mean fold change? 

- Question 3f. Make a volcano plot from 5e (label any statistically significant genes)

- Question 3g. Now examine the complete files: compare the fold change in the complete files vs the different subsamples making sure to address replicates and the size of the random sample. 


### Packaging

The solutions to the above questions should be submitted as a single PDF document that includes your name, email address, and all relevant figures (as needed). Submit your solutions by uploading the PDF to [GradeScope](https://www.gradescope.com/courses/431817), and remember to select where in your submission each question/subquestion is. The Entry Code is: J37JKW. 

If you submit after this time, you will use your late days. Remember, you are only allowed 4 late days for the entire semester!

