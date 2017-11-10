import csv
from Bio import SeqIO
import numpy as np


# Read the csv
codonMap = dict()
with open('codonfreqs.csv', 'rb') as csvfile:
    reader = csv.DictReader(csvfile, delimiter=',')
    for row in reader:
        cdn = row['Codon']
        ratio = row['Ratio']
        freq = row['Freq']
        codonMap[cdn] = {'Ratio': ratio, 'Frequency':freq}

# Generate the weight vector
w = np.zeros(64)
codonList = list()
for i,k in enumerate(sorted(codonMap)):
    v = codonMap[k]
    print v
    w[i] = v['Frequency'] 
    codonList.append(k)



def noiseEstimator(seq, weight, P, a, d):
    """Estimate the noisy translational burst size from 
    an e coli gene rare codon bias.
    @param seq: sequence object
    @param weight: weight vector for codons
    @param P: mean protein value
    @a,d: free parameters"""
    counts = np.zeros(64) 
    for i,cdn in enumerate(codonList): # sorted codons
        counts[i] = seq.count(cdn)
    mag = counts.size
    counts /= float(mag)
    numerator = np.vdot(counts, w)
    fcw = d * (numerator / mag)

    return a * np.power(P, fcw)

print codonList
# Read in the FASTA files for the genes
with open('ecoli-nt.fasta', 'rU') as handle:
    for record in SeqIO.parse(handle, 'fasta'):
        seq = record.seq
        print "Noise estimator for " + str(record.id) + ":"
        print noiseEstimator(seq,w,100,0.3,0.7)
 

       
