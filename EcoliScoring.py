# Use exponential model to score a gene in E. coli for rarity
import csv
from Bio import SeqIO
import numpy as np
import re
from EcoliInformatics import exponentialModel, linearModel


def generateCodonMap():
    # Read the csv
    codonMap = dict()
    with open('bestparameters.csv', 'rb') as csvfile:
        reader = csv.DictReader(csvfile, delimiter=',')
        for row in reader:
            codon = row['Codon']
            model = row['Model']
            if model == 'exponential':
                thetas = [float(row['Theta1']), float(row['Theta2']), float(row['Theta3'])]
            elif model == 'linear':
                thetas = [float(row['Theta1']), float(row['Theta2'])]
            else:
                continue
            codonMap[codon] = (model, thetas)
    return codonMap



def naiveIndex(seq, codonMap):
    """A naive scoring with no clustering."""
    score = 0.0
    for i in range(3, len(seq)-3, 3):
        try:
            model, thetas = codonMap[seq[i:i+3]]
        except KeyError:
            print 'No corresponding model found for codon ' + seq[i:i+3] + '.'
            continue

        pos = i / 3 # codon position, NOT nucleotide position
        if model == 'exponential':
            p = exponentialModel(thetas, pos)
        elif model == 'linear':
            p = linearModel(thetas, pos)
        else:
            print 'error. bad things.'
            print codonMap[seq[i:i+3]]
            return

        if p < .15:
            score += p
    return score


def rarityIndex(seq, codonMap):
    """Gives a measure of how rare the codons of a gene are."""
    score = 1.0
    scoreVector = np.zeros(len(seq)/3)
    for i in range(3, len(seq)-3, 3):
        try:
            model, thetas = codonMap[seq[i:i+3]]
        except KeyError:
            print 'No corresponding model found for codon ' + seq[i:i+3] + '.'
            continue

        pos = i / 3 # codon position, NOT nucleotide position
        if model == 'exponential':
            p = exponentialModel(thetas, pos)
        elif model == 'linear':
            p = linearModel(thetas, pos)
        else:
            print 'error. bad things.'
            print codonMap[seq[i:i+3]]
            return

        if p < .15:
            scoreVector[pos] = p 
     
    thresh = 5 
    i = 0
    while i < scoreVector.size:
        if scoreVector[i] != 0:
            clusterScore = 0.0
            k = i
            j = i+1
            once = False 
            while j < scoreVector.size and j - k < thresh:
                if scoreVector[j] != 0:
                    clusterScore += scoreVector[j] * 100
                    once = True
                    k = j
                j += 1
            if once:
                clusterScore += scoreVector[i] * 100
                score += clusterScore
            i = j
        else:
            i += 1
    
    
    return float(score) / (len(seq)/3)


def main():
    # Get map of codons to estimated parameters
    codonMap = generateCodonMap()

    # Read in the FASTA files for the genes
    writingDict = dict()
    with open('ecoli-nt.fasta', 'rU') as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            seq = record.seq.transcribe()
            if len(seq) % 3 != 0:
                print "ERROR! bad data :("
                continue
            else:
                writingDict[record.description] = (len(seq)/3, rarityIndex(seq, codonMap)) 

    # Read in noise data
    noiseMap = dict()
    noisefile = open('ecoli-noise-data.csv', 'rB')
    noisecsvfile = csv.DictReader(noisefile, delimiter=',')
    for row in noisecsvfile:
        k = row['SYMBOL'].lower()
        p = row['<P> (1)']
        kp = row['kp (sec-1)']
        b = row['burst rate'] 
        noiseMap[k[1:-1]] = [p, kp, b]
        

    with open('ecoli-data.csv', 'wb') as csvfile:
        ecoliwriter = csv.writer(csvfile, delimiter=',')
        ecoliwriter.writerow(['Gene', 'Total Codons', 'Rarity Index', '<P>', 'burst rate', 'kp'])
        for k,v in writingDict.items():
            totalCodons, score = v
            
            # Get gene symbol
            l = re.split(' ', k)
            sym = l[1]
            try:
                p, b, kp = noiseMap[sym.lower()]
            except KeyError:
                continue
            
            ecoliwriter.writerow([k, totalCodons, score, p, b, kp])
          

if __name__ == '__main__':
    main()

