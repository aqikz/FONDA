

Input: Samples with CSV file


Two external python scripts: used for computing entropy and tree transformation.

"rules/bowtie.smk"		# Creates bowtie-index, 4. runs mapping
"rules/qc.smk"			# 1. FastQC before trim 2. Trim 3. FastQC after trim, 8. Qualimap for mappings, 9. MultiQC
"rules/samtools.smk"		# 5. sort bam-files, 6. Samtools-index, 7. Samtools-mapping-stats
"rules/ivar.smk"		# Genome assembly
"rules/augur.smk"		# Multiple alignment, Phylo Tree, Entropies

Note: still contains errors, no successful run over tiny data so far.
