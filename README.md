# Differential expression of small RNAs in humans and chimpanzees
## Description
This project, analysis differential gene expression between two groups (humans and chimpanzees). The procedure for this analysis is done in four steps: 
 - Preprocessing of raw reads from fastq files
 - Run a Snakemake pipeline
 - Run gene_name.py script
 - Differential expression analysis (DEA) in R using DESeq2 package 


# Installation
The Snakemake pipeline (Snakefile), gene_name.py and DEA.R are ready for download.

## Programs used in the analysis:

 - Python v.3.6.6
 - snakemake v.5.2.4
 - R v.4.0.0
 
Programs run in snakemake:
 - cutadapt v.2.9
 - fastqc v.0.11.3
 - MultiQC v.1.2
 - STAR v.2.6.0c
 - SAMtools v.1.6
 - deepTools v.2.5.4
 - Subread v.1.6.3

Programs run in R:
 - Bioconductor v.3.11
 - DESeq2 v.1.28.1

# Usage

## Preprocessing of raw reads:
Preprocessing of fastq raw reads is done by trimming adapters and bases from each end of the reads. How this is done is usually described by the manufacturer of the sequencing kit.

The cutadapt program was used for trimming:

```
#trim 4 bases of each end
 ls *.fastq | sed s/.fastq// | while read file; do cutadapt -u 4 -u -4 -o $file.trimmed.fastq $file.fastq; echo $file processed; done
```
```
#trim adapter
ls *.trimmed.fastq | sed s/.trimmed.fastq// | while read file; do cutadapt -a what to trim -o $file.trimmed_1.fastq --minimum-length length $file.trimmed.fastq; echo $file processed; done
```
After trimming, the trimmed_1.fastq files are ready to be run with the Snakemake pipeline. 

## Snakemake pipeline (Snakefile):

The snakemake pipeline takes trimmed fastq files and performs the following work:

 - Quality control: FastQC and MultiQC.
 - Align reads to a reference genome (hg38) with STAR aligner.
 - Sort and index bam files.
 - Converts bam files to BigWig files.
 - Read quantification with featureCounts.

The Snakefile is designed to run with modules. If the user does not use modules; make sure that the programs are executable in the working directory, and erase lines starting with "module load" in the "shell:" section of each rule in the Snakefile.  

Directory structure to run snakemake pipeline; the Snakefile is suppose to be in the same directory as the data and bin directories.

**Example:**

**&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;projectdir/**

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Snakefile

**&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;config.yaml**

**&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;GTF annotation file**

**&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;data/**

**&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;bin/**

A file named config.yaml is created by the user, and the file has to include a list called "samples" with the names of the fastq sample files used, as shown below.

    #config.yaml
    
    samples: ["sample_1.fastq", "sample_2.fastq", sample_3.fastq"]

Snakefile is executed in the project directory:

    #run snakefile
    snakemake



The output files from the Snakemake pipeline, are a single .csv file and multiple .bw files. The .csv file is used directly with DEA.R for the differential expression analysis. Bw files are used with the UCSC genome browser.


##  Gene_info script
**gene_info.py** will extract, the gene id, type and name from the GTF annotation file. Three output files are generated, gene_id.txt, gene_type.txt and gene_name.txt. These files are used as input files for the **DEA.R** script.

    #run the gene_info.py
    python3 gene_info.py -i file.gtf -o gene_id.txt -t gene_type.txt -n gene_name.txt

## Differential expression analysis with DESeq2:
The **DEA.R** script performs differential expression analysis and visualizes the results.

The .csv file from snakemake, along with the three files from **gene_info.py**, are input files for the **DEA.R** script in R.

In the **DEA.R** script the user has to: 
 - Set a directory from where the input files are stored.
 - Install packages used in the script.
 - Edit some parts of the script to fit the users analysis:
	 - Sections marked with \*
		 - Change column names in the "Rename columns" section to match sample names.
		 - Add or remove columns based on sample numbers in the "Delete columns" section.
		 - Edit the "Metadata" to fit sample conditions
		 - In the "Set factor and run DESeq" section change factor conditions.
		 - The "Pool together small RNAs types" section can be skipped.
