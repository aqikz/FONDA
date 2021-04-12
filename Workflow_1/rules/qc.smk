# Check quality of fastq-files before trimming
rule fastqc_before_trim:
    input:
        config["sample_dir"] + "{sample}.fastq.gz"
    output:
        html="output/stats/fastqc_original/{sample}.html",
        zip="output/stats/fastqc_original/{sample}_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params: ""
    log:
        "logs/fastqc_original/{sample}.log"
    wrapper:
        "0.60.1/bio/fastqc"


# Perform trimming on fastq-files
rule trimmomatic_pe:
    input:
        r1= config["sample_dir"] + "{sample}_1.fastq.gz",
        r2= config["sample_dir"] + "{sample}_2.fastq.gz"
    output:
        r1="output/trimmed/{sample}_1.fastq.gz",
        r2="output/trimmed/{sample}_2.fastq.gz",
        
        # reads where trimming entirely removed the mate
        r1_unpaired="output/trimmed/{sample}_1.unpaired.fastq.gz",
        r2_unpaired="output/trimmed/{sample}_2.unpaired.fastq.gz"
    log:
        "logs/trimmomatic/{sample}.log"
    params:
        # list of trimmers (see manual)
        trimmer=["TRAILING:3"],
        # optional parameters
        extra="",
        compression_level="-9"
    threads:
        32
    wrapper:
        "0.60.1/bio/trimmomatic/pe"


# Check quality of fastq-files after trimming
rule fastqc_after_trim:
    input:
        "output/trimmed/{sample}.fastq.gz",
    output:
        html="output/stats/fastqc_trimmed/{sample}.html",
        zip="output/stats/fastqc_trimmed/{sample}_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params: ""
    log:
        "logs/fastqc_trimmed/{sample}.log"
    wrapper:
        "0.60.1/bio/fastqc"


# Qualimap2 -> get a statistic for the mapping results
rule qualimap_statistics:
	input:
		"output/mapped/{sample}.sorted.bam"
	output:
		"output/stats/qualimap2/{sample}/qualimapReport.html"
	params:
		out_dir = "output/stats/qualimap2/{sample}"
	conda:
		"../envs/env.yaml"
	shell:
		"qualimap bamqc -bam {input} -outdir {params.out_dir}"

# Bind all fastqc & Qualimap statistics together
rule multiqc:
	input:
		expand("output/stats/qualimap2/{sample}/qualimapReport.html", sample=list(samples.index)),
		expand("output/stats/fastqc_original/{sample}_1.html", sample=list(samples.index)),
		expand("output/stats/fastqc_original/{sample}_2.html", sample=list(samples.index)),
		expand("output/stats/fastqc_trimmed/{sample}_1.html", sample=list(samples.index)),
		expand("output/stats/fastqc_trimmed/{sample}_2.html", sample=list(samples.index))
	output:
		"output/stats/multiqc.html"
	params:
		""  # Optional: extra parameters for multiqc.
	log:
		"logs/multiqc.log"
	wrapper:
		"0.60.1/bio/multiqc"