import argparse
from itertools import combinations
import pyfaidx
import itertools
import math
import numpy as np
import pandas as pd
from scipy.special import comb


def parse_args():
    parser = argparse.ArgumentParser(
        description="""Evaluate entropy and overlap complexity of seeds""")
    # parser.add_argument("-g", "--genomes", type=str, help="Paths to the genomes, comma delimited", required=True)
    #parser.add_argument("-S", "--sseeds", type=str, help="Set of seeds to evaluate", required=True)
    parser.add_argument("-e", "--entropy-bits", type=int, required=True,
                        help="Specify number of bits per seed for entropy calculation")
    parser.add_argument("-f", "--file", type=str, required=True,
                        help="inputfile")
    args = parser.parse_args()
    return args

def calculate_entropy(seed, s_size):
    bit_counts = {}
    for i in itertools.product(['0', '1'], repeat=s_size):
        key = ''.join(i)
        bit_counts[key] = 0
    for i in range(len(seed) - s_size + 1):
        bit = seed[i:i+s_size]
        bit_counts[bit] += 1
    # print(bit_counts)
    # Calculate Shannon entropy
    entropy = 0
    for bit, count in bit_counts.items():
        if count == 0:
            continue
        num_substrings = len(seed) - s_size + 1
        entropy -= float(count)/float(num_substrings) * math.log(float(count)/float(num_substrings), 2)
    return entropy

def nCr(n, r):
    return int(math.factorial(n)/math.factorial(r)/math.factorial(n-r))

def overlap_complexity(seeds):
    # For now assume 2 seeds
    complexity = 0
    for seed1, seed2 in combinations(seeds, 2):
        for i in range(1, len(seed1)):
            sigma1 = get_sigma(seed1[:i], seed2[len(seed2)-i:])
            sigma2 = get_sigma(seed1[len(seed1)-i:], seed2[:i])
            complexity += 2 ** sigma1 + 2 ** sigma2
        final_sigma = get_sigma(seed1, seed2)
        complexity += 2 ** final_sigma
    return complexity

def get_sigma(s1, s2):
    sigma = 0
    for i, j in zip(s1, s2):
        if i == '1' and j == '1':
            sigma += 1
    return sigma

def main():
    args = parse_args()
    fh = open(args.file, "r")
    for line in fh:
        seedSet = line.strip().split(" ")
        print(line.strip() + "\t"+ str(overlap_complexity(seedSet)))
#     seedSet = args.sseeds.split(" ")
#     for seed in seedSet:
#         print(calculate_entropy(seed, args.entropy_bits))
#     print(overlap_complexity(seedSet))
#     calculate_entropy_vect = np.vectorize(calculate_entropy, excluded=['s_size'])
#     entropies_1 = calculate_entropy_vect(seeds, 2)
#     entropies_2 = calculate_entropy_vect(seeds, 3)
#     entropies_3 = calculate_entropy_vect(seeds, 4)
#     data = pd.DataFrame({'2bit': entropies_1, '3bit': entropies_2, '4bit': entropies_3},
#                  index=seeds,
#                  columns=['2bit', '3bit', '4bit'])
#     data.to_csv(args.output, sep='\t')
#     print(data.describe())


if __name__ == '__main__':
    main()
