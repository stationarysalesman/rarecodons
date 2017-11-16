import numpy as np
import matplotlib.pyplot as plt
from Bio import SeqIO

fname = 'ecoli-nt.fasta'
numCodons = 59 # number of degenerate codons
def getSequences():
    lst = []
    for record in SeqIO.parse(fname, 'fasta'):
        lst.append(record.seq)
    return lst


def countCodons(sequences, maxLength, numSequences):
    """Count the times a codon occurs in a transcript
    at all positions."""
    
    # Keep track of which codons (first axis) occur at which locations (second axis)
    n = np.zeros((59, numSequences)) 


def MLE():
    """Estimate the parameters of codon bias in E. coli based 
    on Krumpp et. al. exponential model."""

    # Get metadata about the sequences
    sequences = getSequences() 
    lengths = map(lambda x: len(x), sequences)
    maxLength = max(lengths)
    numSequences = len(sequences)

    # Find out the codon distribution across the genes
    n = countCodons()
