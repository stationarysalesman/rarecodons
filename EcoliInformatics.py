import numpy as np                  
import matplotlib.pyplot as plt     # Plotting, if you want it
import csv                          # Read CSVs
from scipy.special import comb      # N choose K
import scipy.optimize               # minimize
from Bio import SeqIO               # Biopython sequence I/O

fname = 'ecoli-nt.fasta'
numCodons = 64 # number of degenerate codons, change later(?)
                # Clearly wrong, but might be interesting to look at
numAA = 21      # include stop codon

aaMap = {'F': 0, 'L': 1, 'I': 2, 'M': 3, 'V': 4, 'S': 5, 'P': 6, 'T': 7, 'A': 8, 'Y': 9, \
            'STOP': 10, 'H': 11, 'Q': 12, 'N': 13, 'K': 14, 'D': 15, 'E': 16, 'C': 17, 'R': 18, \
            'G': 19, 'W': 20}

def generateCodonMap():
    """Generate a mapping of codons/amino acids onto integers."""
    codonMap = dict()
    f = open('codonfreqs.csv', 'rB')
    csvfile = csv.DictReader(f, delimiter=',')
    for i, item in enumerate(csvfile):
        codonMap[str(item['Codon'])] = (i, item['Amino Acid'])
    return codonMap


def getSequences():
    lst = []
    for record in SeqIO.parse(fname, 'fasta'):
        lst.append(record.seq.transcribe())
    return lst


def countsInGene(gene, maxLength, codonMap):
    """Count the codons and amino acids in a gene based on their positions."""
    n = np.zeros((numCodons, maxLength))
    N = np.zeros((numAA, maxLength))
    for i in range(0, len(gene)-2, 3):
        cdn = str(gene[i:i+3])
        n[codonMap[cdn][0]][i/3] += 1
        N[aaMap[codonMap[cdn][1]]][i/3] += 1
    return n,N


def counts(sequences, maxLength, codonMap):
    """Count the times a codon occurs in a transcript
    at all positions, and times an amino acid occurs."""
    
    # Keep track of which codons/amino acid (first axis) occur 
    # at which locations (second axis)
    n = np.zeros((numCodons, maxLength)) 
    N = np.zeros((numAA, maxLength))
    for seq in sequences:
        n_add, N_add = countsInGene(seq, maxLength, codonMap) 
        n += n_add
        N += N_add 
    return n,N


def exponentialModel(thetas, k):
    """Exponential model for codon distribution.
    Calculates the probabilty of seeing codon i at 
    position k. See Krumpp et. al. materials and methods."""

    return thetas[0] * np.exp(-k / thetas[1]) + thetas[2] 

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
    n,N = counts(sequences, maxLength, codonMap)

    # Now we have the numbers. It's time to estimate the parameters
    # for each distribution separately. We'll define a new function
    # for each codon, allowing us to write a function whose only 
    # varying parameters are the thetas (in the main text, these are 
    # a, tau, and c). Then, we can run an optimization algorithm on 
    # the log-likelihood function.
    for codon, v in codonMap.items():
        i, aa = v
        n_i = n[i] 
        N_j = N[aaMap[aa]]
        def L(thetas):
            total = 0.0
            for k in range(maxLength):
                binomialCoeff = comb(N_j[k], n_i[k])
                modelVal = exponentialModel(thetas, k)
                factor1 = np.power(modelVal, n_i[k])
                factor2 = np.power(1-modelVal, N_j[k]-n_i[k]) 
                result = binomialCoeff * factor1 * factor2
                total += np.log(result)
            return -total

        # Run 5 times to help avoid getting stuck at local maxima 
        for run in range(5): 
            thetas = [np.random.random()*5, \
                      np.random.random()*5, \
                      np.random.random()*5]

            result = scipy.optimize.minimize(L, thetas)
            print "Results for codon " + str(codon)
            print str(result['x']) 

MLE()

