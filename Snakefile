# define config file, that contains samples list with sample names.
configfile: "config.yaml"

rule all:
	input:
		"results/QC/trimmed/multiQC/multiqc.html",
		expand("results/Sorted_bam/default_map/{sample}_sorted.bam.bai", sample=config["samples"]),
		expand("results/Sorted_bam/unique_map/{sample}_sorted.bam.bai", sample=config["samples"]),
		expand("results/Sorted_bam/default_map/{sample}_sorted.bw", sample=config["samples"]),
		expand("results/Sorted_bam/unique_map/{sample}_sorted.bw", sample=config["samples"]),
		"results/Counts/trimmed/counts_default.csv",
		"results/Counts/trimmed_unique/counts_unique.csv"

rule fastqc:
	# Run FastQC on a FASTQ file.
	input:
		"data/trimmed/{sample}.fastq"
	output:
		"results/QC/trimmed/fastQC/{sample}_fastqc.html",
		"results/QC/trimmed/fastQC/{sample}_fastqc.zip"
	priority: 12
	shell:
	# Run fastQC
		"""
		module load fastqc/0.11.3
		fastqc {input} -t 6 -o results/QC/trimmed/fastQC/
		module purge
		"""

rule multiqc:
	# Run MultiQC on the fastqc.zip files
	input:
		expand("results/QC/trimmed/fastQC/{sample}_fastqc.zip", sample=config["samples"])
	output:
		html = "results/QC/trimmed/multiQC/multiqc.html",
		stats = "results/QC/trimmed/multiQC/multiqc_general_stats.txt"
	priority: 11
	log:
		"results/logs/trimmed/multiQC/multiqc.log"
	shell:
	# Run MulitQC
		"""
		module load icc/2017.4.196-GCC-6.4.0-2.28
		module load impi/2017.3.196
		module load ifort/2017.4.196-GCC-6.4.0-2.28
		module load impi/2017.3.196
		module load MultiQC/1.2-Python-2.7.14

		multiqc -n multiqc.html {input} 2> {log}
		mv multiqc.html {output.html}
		mv multiqc_data/multiqc_general_stats.txt {output.stats}
		mv multiqc_data results/QC/trimmed/multiQC
		module purge
		"""

rule mapping:
	# Align reads to reference genome with default parameters
	input:
		sample="data/trimmed/{sample}.fastq",
		genome="/projects/fs1/common/genome/lunarc/indicies/star/human/hg38",
		GTF="data/gencode.v30.annotation.gtf"

	output:
		"results/Mapping/STAR/trimmed/{sample}_Aligned.out.bam",
		"results/logs/STAR/{sample}_Log.final.out"
	priority: 10
	shell:
	# Run STAR aligner
		"""
		module load GCC/5.4.0-2.26
		module load OpenMPI/1.10.3
		module load STAR/2.6.0c

		STAR --runThreadN 10 \
		--readFilesIn {input.sample} \
		--genomeDir {input.genome} \
		--outSAMtype BAM Unsorted \
		--outSAMunmapped Within \
		--sjdbGTFfile {input.GTF} \
		--outFileNamePrefix {wildcards.sample}_

		mv {wildcards.sample}_Aligned.out.bam {output[0]}
		mv {wildcards.sample}_Log.final.out {output[1]}
		mv {wildcards.sample}_Log.out {wildcards.sample}_Log.progress.out {wildcards.sample}_SJ.out.tab results/logs/STAR
		mv {wildcards.sample}__STARgenome results/Mapping/STAR/trimmed

		module purge
		"""
rule sort_bam:
	# Sort default bam files
	input:
		"results/Mapping/STAR/trimmed/{sample}_Aligned.out.bam"
	output:
		"results/Sorted_bam/default_map/{sample}_sorted.bam"
	priority: 9
	shell:
	# Run samtools to sort
		"""
		module load icc/2017.4.196-GCC-6.4.0-2.28  impi/2017.3.196
		module load ifort/2017.4.196-GCC-6.4.0-2.28  impi/2017.3.196
		module load SAMtools/1.6

		samtools sort {input} -o {output}

		module purge
		"""

rule index_bam:
	# Index sorted bam files
	input:
		"results/Sorted_bam/default_map/{sample}_sorted.bam"
	output:
		"results/Sorted_bam/default_map/{sample}_sorted.bam.bai"
	priority: 8
	shell:
	# Run samtools to index files
		"""
		module load icc/2017.4.196-GCC-6.4.0-2.28  impi/2017.3.196
		module load ifort/2017.4.196-GCC-6.4.0-2.28  impi/2017.3.196
		module load SAMtools/1.6

		samtools index {input}

		module purge

		"""

rule bigwig:
	# Convert sorted bam files to bigwig files
	input:
		"results/Sorted_bam/default_map/{sample}_sorted.bam"
	output:
		"results/Sorted_bam/default_map/{sample}_sorted.bw"
	priority: 7

	shell:
	# Run bamCoverage
		"""
		module load GCC/7.3.0-2.30  OpenMPI/3.1.1
		module load deepTools/2.5.4-Python-3.6.6

		bamCoverage -b {input} -o {output}

		module purge

		"""

rule unique_mapping:
	# Align reads to reference genome with unique parameters
	input:
		sample="data/trimmed/{sample}.fastq",
		genome="/projects/fs1/common/genome/lunarc/indicies/star/human/hg38",
		GTF="data/gencode.v30.annotation.gtf"

	output:
		"results/Mapping/STAR_unique/trimmed/{sample}_Aligned.out.bam",
		"results/logs/STAR_unique/{sample}_Log.final.out"
	priority: 6

	shell:
	# Run STAR aligner
		"""
		module load GCC/5.4.0-2.26
		module load OpenMPI/1.10.3
		module load STAR/2.6.0c

		STAR --runThreadN 10 \
		--readFilesIn {input.sample} \
		--genomeDir {input.genome} \
		--outSAMtype BAM Unsorted \
		--outSAMunmapped Within \
		--sjdbGTFfile {input.GTF} \
		--outFilterMultimapNmax 1 \
		--outFileNamePrefix {wildcards.sample}_

		mv {wildcards.sample}_Aligned.out.bam {output[0]}
		mv {wildcards.sample}_Log.final.out {output[1]}
		mv {wildcards.sample}_Log.out {wildcards.sample}_Log.progress.out {wildcards.sample}_SJ.out.tab results/logs/STAR_unique
		mv {wildcards.sample}__STARgenome results/Mapping/STAR_unique/trimmed

		module purge
		"""
rule sort_bam_unique:
	# Sort unique bam files
	input:
		"results/Mapping/STAR_unique/trimmed/{sample}_Aligned.out.bam"
	output:
		"results/Sorted_bam/unique_map/{sample}_sorted.bam"
	priority: 5
	shell:
	# Run samtools to sort files
		"""
		module load icc/2017.4.196-GCC-6.4.0-2.28  impi/2017.3.196
		module load ifort/2017.4.196-GCC-6.4.0-2.28  impi/2017.3.196
		module load SAMtools/1.6

		samtools sort {input} -o {output}

		module purge
		"""
rule index_bam_unique:
	# Index sorted unique files
	input:
		"results/Sorted_bam/unique_map/{sample}_sorted.bam"
	output:
		"results/Sorted_bam/unique_map/{sample}_sorted.bam.bai"
	priority: 4
	shell:
	# Run samtools to index files
		"""
		module load icc/2017.4.196-GCC-6.4.0-2.28  impi/2017.3.196
		module load ifort/2017.4.196-GCC-6.4.0-2.28  impi/2017.3.196
		module load SAMtools/1.6

		samtools index {input}


		module purge

		"""

rule bigwig_unique:
	# Convert unique bam files to bigwig files
	input:
		"results/Sorted_bam/unique_map/{sample}_sorted.bam"
	output:
		"results/Sorted_bam/unique_map/{sample}_sorted.bw"
	priority: 3

	shell:
	# Run bamCoverage to convert to bigwig files
		"""
		module load GCC/7.3.0-2.30  OpenMPI/3.1.1
		module load deepTools/2.5.4-Python-3.6.6

		bamCoverage -b {input} -o {output}

		module purge

		"""

rule featureCounts:
	# Quantify default bam files
	input:
		GTF="data/gencode.v30.annotation.gtf",
		bam=expand("results/Mapping/STAR/trimmed/{sample}_Aligned.out.bam", sample=config["samples"])
	output:
		"results/Counts/trimmed/counts_default.csv"
	priority: 2

	shell:
	# Run featureCounts
		"""
		module load GCC/7.3.0-2.30
		module load OpenMPI/3.1.1
		module load Subread/1.6.3

		featureCounts -T 10 -s 2 -a {input.GTF} -o {output} {input.bam}


		module purge
		"""

rule unique_featureCounts:
	# Quantify unique bam files
	input:
		GTF="data/gencode.v30.annotation.gtf",
		bam=expand("results/Mapping/STAR_unique/trimmed/{sample}_Aligned.out.bam", sample=config["samples"])
	output:
		"results/Counts/trimmed_unique/counts_unique.csv"
	priority: 1

	shell:
	# Run featureCounts
		"""
		module load GCC/7.3.0-2.30
		module load OpenMPI/3.1.1
		module load Subread/1.6.3

		featureCounts -T 10 -s 2 -a {input.GTF} -o {output} {input.bam}


		module purge
		"""
