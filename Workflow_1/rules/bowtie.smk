# Build bowtie-index for reference file
rule bowtie_build_index:
	input:
		reference_file
	output:
		expand("{index_name}.{index_ext}", index_name=index, index_ext=INDEX_EXT)
	conda:
		"../envs/env.yaml"
	threads:	4
	shell:
		"bowtie2-build {input} {index}"

		 
# Perform mapping of paired ends 
rule bowtie2_pairedEnd_mapping:
    input:
        sample=["output/trimmed/{sample}_1.fastq.gz", "output/trimmed/{sample}_2.fastq.gz"]
    output:
        "output/mapped/{sample}.bam"
    log:
        "logs/bowtie2/{sample}.log"
    params:
        index=expand("{index_files}", index_files=config["index"]),  # prefix of reference genome index (built with bowtie2-build)
        extra=expand("-p {p_param} --np {np_param}", p_param=config["bowtie_parameters"]["p_param"],
         np_param=config["bowtie_parameters"]["np_param"]) 	# optional parameters
    threads: 8  # Use at least two threads
    wrapper:
        "0.60.1/bio/bowtie2/align"