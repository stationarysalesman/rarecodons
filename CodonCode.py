import csv
from Bio import SeqIO
import numpy as np
import re

# Read the csv
codonMap = dict()
with open('codonfreqs.csv', 'rb') as csvfile:
    reader = csv.DictReader(csvfile, delimiter=',')
    for row in reader:
        cdn = row['Codon']
        ratio = row['Ratio']
        freq = row['Freq']
        codonMap[cdn] = {'Ratio': float(ratio), 'Frequency':float(freq)}

# Generate the weight vector
w = np.zeros(64)
codonList = list()
for i,k in enumerate(sorted(codonMap)):
    v = codonMap[k]
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


def sequenceInfo(seq):
    """Get some information about a sequence."""
    rare = 0
    common = 0
    counts = np.zeros(64)

    for i,cdn in enumerate(codonList): # sorted codons
        # Simple count
        count = seq.count(cdn)
        if codonMap[cdn]['Frequency'] < 0.5:
            rare += count
        else:
            common += count

        # Scoring
        counts[i] = count

    rcp = rare/float(rare+common)
    total = rare+common
    score = np.vdot(counts, w)
    return [rare, common, total, rcp, score]
    

    # Read in the FASTA files for the genes
writingDict = dict()
with open('ecoli-nt.fasta', 'rU') as handle:
    for record in SeqIO.parse(handle, 'fasta'):
        seq = record.seq.transcribe()
        writingDict[record.description] = sequenceInfo(seq) 

# Read in noise data
noiseMap = dict()
noisefile = open('ecoli-noise-data.csv', 'rB')
noisecsvfile = csv.DictReader(noisefile, delimiter=',')
for row in noisecsvfile:
    k = row['SYMBOL'].lower()
    p = row['<P> (1)']
    kp = row['kp (sec-1)']
    noiseMap[k[1:-1]] = [p, kp]
    

with open('ecoli-naive-data.csv', 'wb') as csvfile:
    ecoliwriter = csv.writer(csvfile, delimiter=',')
    ecoliwriter.writerow(['Gene', 'Rare Codons', 'Common Codons', 'Total Codons', 'RC%', 'Score', '<P>', 'kp'])
    for k,v in writingDict.items():
        r, c, t, rcp, score = v
        
        # Get gene symbol
        l = re.split(' ', k)
        sym = l[1]
        try:
            p, kp = noiseMap[sym.lower()]
        except KeyError:
            continue
        
        ecoliwriter.writerow([k, r, c, t, rcp, score, p, kp])
       
