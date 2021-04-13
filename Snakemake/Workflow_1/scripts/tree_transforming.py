import sys
from ete3 import Tree

tree_file = snakemake.input[0]
png_output = snakemake.output[0]

t = Tree(tree_file)
t.render(png_output, w=183, units="mm")
