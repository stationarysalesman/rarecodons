import numpy as np                  
import matplotlib.pyplot as plt     # Plotting, if you want it
import csv                          # Read CSVs
import scipy.optimize               # minimize
import scipy.stats                  # binomial pmf
import signal, os                   # signaling, pid
import time                         # for seeding PRNG
from Bio import SeqIO               # Biopython sequence I/O
from multiprocessing import Pool

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


def parallelCodonMLE(l):
    return codonMLE(l[0], l[1], l[2], l[3], l[4], l[5])


def codonMLE(codon, aa, n_i, N_j, maxLength, num):
    """Estimate model parameters for a single codon distribution."""
   
    def L(thetas):
        total = 0.0
        for k in range(maxLength):
            total += scipy.stats.binom.logpmf(n_i[k], N_j[k], exponentialModel(thetas, k))
        return -1 * total

    # Run 5 times to help avoid getting stuck at local maxima 
    print "Performing MLE for codon " + str(codon) + "(" + str(aa) + ")..."
    thetas_lst = []
    fname_lst = []
    
    # compensate for workers receiving jobs at same time
    s = hash(int(os.getpid()+time.time()+maxLength))
    np.random.seed(s)
    thetas = [np.random.uniform(0.0, 1.0), \
              np.random.uniform(0.0, 1.0), \
              np.random.uniform(0.0, 1.0)]

    
    result = scipy.optimize.fmin(L, thetas)
    print "Optimization complete. Result:"
    print result 
    xs = np.array(range(maxLength), dtype=float)
    ys = np.array(map(lambda x: exponentialModel(result, x), xs))
    plt.plot(xs[1:700], ys[1:700], color='red')
    plt.ylabel('P({}|{})'.format(codon, aa))
    plt.xlabel('Distance from start codon')
    plt.ylim((0.0, 1.0)) 
    plt.xlim((0.0, 700))
    data_ys = [sum(n_i[i:i+8])/sum(N_j[i:i+8]) for i in range(0, ys.size-8, 8)] 
    plt.plot(xs[:696:8], data_ys[:87], color='black')
    fname = 'data/{}{}{}.png'.format(codon, aa, num)
    plt.savefig(fname)
    plt.clf()
    return (result, fname) 


    
def test():
    codonMap = generateCodonMap()
    sequences = getSequences() 
    lengths = map(lambda x: len(x)/3, sequences)
    maxLength = 700*3
    numSequences = len(sequences)
    codonCounts, aaCounts = counts(sequences, maxLength, codonMap)
    n_i = codonCounts['CCG']
    N_j = aaCounts[codonMap['CCG']]
    itemList = [] 
    for run in range(5):
        thetas, fname = codonMLE('CCG', 'P', n_i, N_j, maxLength, run)
        itemList.append((thetas, fname))
    bestThetas, fname, minKey = calcBestParams(n_i,N_j,itemList)
    print str(bestThetas)
    print fname
    print str(minKey)

def calcBestParams(n_i, N_j, itemList):
    """Use RMSE to determine the parameters that fit best to the 
    conditional probability P(codon | amino acid)."""

    scores = dict()
    for item in itemList: 
        thetas = item[0]
        fname = item[1]
        xs = np.array(range(700), dtype=float)
        ys = np.array(map(lambda x: exponentialModel(thetas, x), xs))
        xs_smooth = xs[:692:8]
        ys_smooth = ys[:692:8]
        cprob_ys_smooth = np.array([sum(n_i[i:i+8]) / sum(N_j[i:i+8])  \
                                    for i in range(0, ys.size-8, 8)])
        rmse = np.sqrt(np.divide(np.sum(np.square(ys_smooth - cprob_ys_smooth)), ys.size))
        scores[rmse] = (thetas, fname)
    minKey = sorted(scores)[0]
    bestThetas, fname = scores[minKey]
    return bestThetas, fname, minKey
     

def analyze():
    """Analyze data stored in a csv file."""
    
    # Get metadata about the sequences
    codonMap = generateCodonMap()
    sequences = getSequences() 
    lengths = map(lambda x: len(x)/3, sequences)
    maxLength = 700*3
    numSequences = len(sequences)

    # Find out the codon distribution across the genes
    codonCounts, aaCounts = counts(sequences, maxLength, codonMap)

    # Get parameters from a csv file
    fname = 'parameters.csv'
    f = open(fname, 'rB')
    csvReader = csv.DictReader(f, delimiter=',')
    parameterMap = dict() 
    for row in csvReader:
        thetas = [row['Theta1'], row['Theta2'], row['Theta3']]
        thetaFloats = map(lambda x: float(x), thetas)
        try:
            parameterMap[row['Codon']].append((thetaFloats, row['filename']))
        except KeyError:
            parameterMap[row['Codon']] = [(thetaFloats, row['filename'])] 

    of = open('bestparameters.csv', 'wB')
    csvWriter = csv.writer(of, delimiter=',')
    csvWriter.writerow(['Codon', 'Amino Acid', 'Theta1', 'Theta2', 'Theta3', 'RMSE', 'filename'])
    for codon, itemList in parameterMap.items():
        n_i = codonCounts[codon]
        N_j = aaCounts[codonMap[codon]] 
        bestThetas, fname, rmse = calcBestParams(n_i, N_j, itemList)
        csvWriter.writerow([codon, codonMap[codon], bestThetas[0], bestThetas[1], bestThetas[2], 
                            rmse, fname])

def MLE():
    """Estimate the parameters of codon bias in E. coli based 
    on Krumpp et. al. exponential model."""

    # Get metadata about the sequences
    codonMap = generateCodonMap()
    sequences = getSequences() 
    lengths = map(lambda x: len(x)/3, sequences)
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
    fieldnames = ['Codon', 'Amino Acid', 'Theta1', 'Theta2', 'Theta3', 'filename']
    outfile = open('parameters.csv', 'wB')
    csvWriter = csv.writer(outfile, delimiter=',')
    csvWriter.writerow(fieldnames)
    argsList = []
    for codon, aa in codonMap.items():
        n_i = codonCounts[codon] 
        N_j = aaCounts[aa]
        for num in range(5): 
            argsList.append([codon, aa, n_i, N_j, maxLength,num])

    p = Pool(16)
    resultList = p.map(parallelCodonMLE, argsList)    
    for i,item in enumerate(argsList):
        csvWriter.writerow([item[0], item[1], resultList[i][0][0], resultList[i][0][1], 
                            resultList[i][0][2], resultList[i][1]])       

 
np.random.seed(789345678)
MLE()

#test() 
#analyze()
