# Create multiple alignment
rule augur_multiple_alignment:
    input:
        #sequences = rules.filter.output.sequences,
        sequences = expand("output/ivar/{sample}.fa", sample=list(samples.index)),
        reference = reference_file
    output:
        alignment = "output/augur/aligned.fasta"
    conda:
        "../envs/augur_env.yaml"
    shell:
        """
        augur align \
            --sequences {input.sequences} \
            --reference-sequence {input.reference} \
            --output {output.alignment}
        """

# Create phylogenetic tree -> tree is stored in Newick format with Branch lengths == nucleotide divergence
rule augur_create_tree:
    input:
        alignment = rules.augur_multiple_alignment.output.alignment
    output:
        tree = "output/augur/tree_raw.nwk"
    conda:
        "../envs/augur_env.yaml"
    shell:
        """
        augur tree \
            --alignment {input.alignment} \
            --output {output.tree}
        """


# 
rule newick_tree_to_png:
    input:
        tree = rules.augur_create_tree.output.tree
    output:
        image = "output/augur/phylogenetic_tree.png"
    conda:
        "../envs/augur_env.yaml"
    script:
        "../scripts/tree_transforming.py"
        

rule compute_entropies_and_plotting:
    input:
        multiple_alignment = rules.augur_multiple_alignment.output.alignment
    output:
        entropy_file = "output/augur/entropies.txt",
        plotting = "output/augur/entropies.png"
    params:
        window_size = config["window_size_x"]
    conda:
        "../envs/augur_env.yaml"
    script:
        "../scripts/compute_entropy.py"
    
    
