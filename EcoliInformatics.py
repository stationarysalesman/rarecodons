import numpy as np
import matplotlib.pyplot as plt
import csv
from Bio import SeqIO

fname = 'ecoli-nt.fasta'
numCodons = 64 # number of degenerate codons

def generateCodonMap():
    """Generate a mapping of codons onto integers."""
    codonMap = dict()
    f = open('codonfreqs.csv', 'rB')
    csvfile = csv.DictReader(f, delimiter=',')
    for i, item in enumerate(csvfile):
        codonMap[str(item['Codon'])] = i 
    return codonMap


def getSequences():
    lst = []
    for record in SeqIO.parse(fname, 'fasta'):
        lst.append(record.seq.transcribe())
    return lst


def countCodonsInGene(gene, maxLength, codonMap):
    """Count the codons in a gene based on their positions."""
    n = np.zeros((64, maxLength))
    for i in range(0, len(gene)-2, 3):
       cdn = str(gene[i:i+3])
       n[codonMap[cdn]][i/3] += 1
    return n

def countCodons(sequences, maxLength, codonMap):
    """Count the times a codon occurs in a transcript
    at all positions."""
    
    # Keep track of which codons (first axis) occur 
    # at which locations (second axis)
    n = np.zeros((64, maxLength)) 
    for seq in sequences:
        n += countCodonsInGene(seq, maxLength, codonMap) 

    return n


def MLE():
    """Estimate the parameters of codon bias in E. coli based 
    on Krumpp et. al. exponential model."""

    # Get metadata about the sequences
    codonMap = generateCodonMap()
    sequences = getSequences() 
    lengths = map(lambda x: len(x), sequences)
    maxLength = max(lengths)
    numSequences = len(sequences)

    # Find out the codon distribution across the genes
    n = countCodons(sequences, maxLength, codonMap)
    for codon in n:
        print codon[:20]


MLE()

