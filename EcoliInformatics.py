import numpy as np                  
import matplotlib.pyplot as plt     # Plotting, if you want it
import csv                          # Read CSVs
import scipy.optimize               # minimize
import scipy.stats                  # binomial pmf
import signal, os                   # signaling
from Bio import SeqIO               # Biopython sequence I/O

class TooLongException(Exception):
    """Raise this when things take too long. Impatience is a virtue."""

    def __init__(self):
        pass

def handler(signum, frame):
    raise TooLongException

def generateCodonMap():
    """Generate a mapping of codons onto amino acids."""
    codonMap = dict()
    f = open('codonfreqs.csv', 'rB')
    csvfile = csv.DictReader(f, delimiter=',')
    for i, item in enumerate(csvfile):
        codonMap[str(item['Codon'])] = item['Amino Acid']
    return codonMap


def getSequences():
    lst = []
    for record in SeqIO.parse('ecoli-nt.fasta', 'fasta'):
        lst.append(record.seq.transcribe())
    return lst


def updateCountsInGene(gene, maxLength, codonMap, codonCounts, aaCounts):
    """Update the counts for codons and amino acids in a gene based on their positions."""

    for i in range(0, min(len(gene)-2, maxLength-2), 3):
        cdn = str(gene[i:i+3])
        codonCounts[cdn][i/3] += 1
        aaCounts[codonMap[cdn]][i/3] += 1

    return

def counts(sequences, maxLength, codonMap):
    """Count the times a codon occurs in a transcript
    at all positions, and times an amino acid occurs."""
    
    # Set up the dictionary that will keep track of the counts for 
    # each codon
    codonCounts = dict()
    aaCounts = dict()
    for cdn in codonMap:
        codonCounts[cdn] = np.zeros(maxLength, dtype=float)
    aaCounts = dict()
    for aa in set(codonMap.values()):
        aaCounts[aa] = np.zeros(maxLength, dtype=float)

    # Keep track of which codons/amino acid (first axis) occur 
    # at which locations (second axis)
    for seq in sequences:
        updateCountsInGene(seq, maxLength, codonMap, codonCounts, aaCounts) 

    return codonCounts, aaCounts 


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
    lengths = map(lambda x: len(x)/3, sequences)
#    maxLength = max(lengths)
    maxLength = 700*3
    numSequences = len(sequences)

    # Find out the codon distribution across the genes
    codonCounts, aaCounts = counts(sequences, maxLength, codonMap)

    # Now we have the numbers. It's time to estimate the parameters
    # for each distribution separately. We'll define a new function
    # for each codon, allowing us to write a function whose only 
    # varying parameters are the thetas (in the main text, these are 
    # a, tau, and c). Then, we can run an optimization algorithm on 
    # the log-likelihood function.
    print 'Estimating...' 
    for codon, aa in codonMap.items():
        n_i = codonCounts[codon] 
        N_j = aaCounts[aa] 
        def L(thetas):
            total = 0.0
            for k in range(maxLength):
                total += scipy.stats.binom.logpmf(n_i[k], N_j[k], exponentialModel(thetas, k))
            return -1 * total

        # Run 5 times to help avoid getting stuck at local maxima 
        for run in range(5): 
            thetas = [np.random.uniform(-1.0, 1.0), \
                      np.random.uniform(0.0, 1.0), \
                      np.random.uniform(0.0, 1.0)]
            signal.signal(signal.SIGALRM, handler)
#            signal.alarm(10)
            try:
                result = scipy.optimize.fmin(L, thetas)
                print "Results for codon " + str(codon)
                print result 
                
                # Plot the model against the actual data
                data_xs = np.array(range(maxLength), dtype=float)
                data_ys = np.zeros(data_xs.size, dtype=float)
                for k in range(maxLength):
                    data_ys[k] += n_i[k]
                data_ys /= N_j
                plt.plot(data_xs[1:], data_ys[1:], color='black')
                model_ys = map(lambda x: exponentialModel(thetas, x), data_xs)
                plt.plot(data_xs[1:], model_ys[1:], color='red')
                plt.show()

            except TooLongException:
                print "aw jeez Rick it took too long"


MLE()

