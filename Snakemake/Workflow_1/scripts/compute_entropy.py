from Bio import SeqIO
from math import log2
import matplotlib.pyplot as plt
import numpy as np

# Compute Shannon entropies for every position in multiple alignment
def compute_entropies(output_file):
    fasta_sequences = SeqIO.parse(open(multiple_alignment_file), 'fasta')

    # Collect sequences
    sequences_array = []
    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        sequences_array.append(sequence)
            
    # Same length for all sequences?
    general_length = len(sequences_array[0])
    for seq in sequences_array:
        assert(len(seq) == general_length)
            
            
    # Calculate entropies
    entropy_array = [0] * general_length
    nr_of_sequences = len(sequences_array)
    for i in range(0, general_length):
        count_A = 0.0
        count_C = 0.0
        count_G = 0.0
        count_T = 0.0
        count_N = 0.0
        count_gap = 0.0
        
        # Count appearances
        for seq in sequences_array:
            if seq[i] == "A":
                count_A += 1
            elif seq[i] == "C":
                count_C += 1
            elif seq[i] == "G":
                count_G += 1
            elif seq[i] == "T":
                count_T += 1
            elif seq[i] == "N":
                count_N += 1
            elif seq[i] == "-":
                count_gap += 1
                
        count_array = [count_A, count_C, count_G, count_T, count_N, count_gap]
        count_prod_array = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        
        # calculate partial entropy for every base
        for x in range(0, 6):
            if (count_array[x] > 0):
                product = count_array[x]/nr_of_sequences * log2(count_array[x]/nr_of_sequences)
                count_prod_array[x] = product
        
        # Sum up to calculate total entropy for x-position
        current_entropy = 0.0
        for prod in count_prod_array:
            current_entropy = current_entropy - prod    # minus to turn it positive
        entropy_array[i] = current_entropy
        
    # Output
    with open(output_file, 'w') as f:
        for item in entropy_array:
            f.write("%s\n" % item)
    
    return entropy_array


# Plot entropies
def plot_entropies(entropy_array, output_png_file, window_size):
    fig = plt.figure()
    
    index = np.arange(len(entropy_array))
    x_label = "Alignment-Position regarding window size of " + str(window_size) + "."
    
    plt.bar(index, entropy_array)
    plt.xlabel(x_label, fontsize=10)
    plt.ylabel('Shannon-Entropy-Value', fontsize=10)
    plt.title('Variability along sequenced genome')
    plt.savefig(output_png_file)
    #plt.show()
    

# Collect all params
multiple_alignment_file = snakemake.input[0]
output_txt_file = snakemake.output[0]
output_png_file = snakemake.output[1]
window_size = snakemake.params[0]

# Compute entropies
entropy_array = compute_entropies(output_txt_file)

# Split array according to set X and add entropies up
split_array = np.array_split(entropy_array, (len(entropy_array)/window_size)+1)
updated_entries = []
for array in split_array:
    sum = 0.0
    for x in array:
        sum = sum + x
    updated_entries.append(sum)


# Plot entropies for chosen X
plot_entropies(updated_entries, output_png_file, window_size)
