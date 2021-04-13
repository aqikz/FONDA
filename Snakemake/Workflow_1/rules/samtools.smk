# Converts sam-files to bam-files
rule samtools_view:
    input:
        "output/mapped/{sample}.sam"
    output:
        "output/mapped/{sample}.bam"
    params:
        "-b" # optional params string
    wrapper:
        "0.60.1/bio/samtools/view"

# Sort bam-files
rule samtools_sort:
    input:
        "output/mapped/{sample}.bam"
    output:
        "output/mapped/{sample}.sorted.bam"
    params:
        "-m 2G"
    threads:  # Samtools takes additional threads through its option -@
        8     # This value - 1 will be sent to -@.
    wrapper:
        "0.60.1/bio/samtools/sort"
        	
        	
# Create samtools index
rule samtools_index:
    input:
         "output/mapped/{sample}.sorted.bam"
    output:
        "output/mapped/{sample}.sorted.bam.bai"
    params:
        "" # optional params string
    wrapper:
        "0.60.1/bio/samtools/index"


# Get mapping stats
rule samtools_idxstats:
    input:
        bam="output/mapped/{sample}.sorted.bam",
        idx="output/mapped/{sample}.sorted.bam.bai"
    output:
        "output/mapped/{sample}.bam.idxstats"
    log:
        "logs/samtools/idxstats/{sample}.log"
    wrapper:
        "0.60.1/bio/samtools/idxstats"