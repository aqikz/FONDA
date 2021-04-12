rule ivar_create_consensus_sequence:
    input:
        "output/mapped/{sample}.sorted.bam"
    output:
        "output/ivar/{sample}.fa",
        "output/ivar/{sample}.qual.txt"
    threads:
        2
    conda:
        "../envs/ivar_env.yaml"
    params:
        output_dir = "output/ivar/{sample}"
    shell:
        "samtools mpileup -d 1000 -A -Q 0 {input}| ivar consensus -p {params.output_dir} -q 20 -t 0"
        