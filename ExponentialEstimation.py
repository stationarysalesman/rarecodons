# Use exponential model to score a gene in E. coli for rarity
import csv
from Bio import SeqIO
import numpy as np
import re
from EcoliInformatics import exponentialModel


def generateCodonMap():
    # Read the csv
    codonMap = dict()
    with open('bestparameters.csv', 'rb') as csvfile:
        reader = csv.DictReader(csvfile, delimiter=',')
        for row in reader:
            codon = row['Codon']
            thetas = [float(row['Theta1']), float(row['Theta2']), float(row['Theta3'])]
            codonMap[codon] = thetas
    return codonMap


def rarityIndex(seq, codonMap):
    """Gives a measure of how rare the codons of a gene are."""
    score = 0.0
    bitmap = np.zeros(len(seq)/3)
    for i in range(3, len(seq), 3):
        if len(seq[i:i+3]) < 3:
            print seq.back_transcribe()
        thetas = codonMap[seq[i:i+3]]
        pos = i / 3 # codon position, NOT nucleotide position
        p = exponentialModel(thetas, pos)
        if p < .15:
            bitmap[pos] = 1
    indices = [i for i in range(bitmap.size) if bitmap[i] == 1]
    diff = np.diff(indices)
    score = np.sum(np.array([i for i in diff if i < 5]))
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
                writingDict[record.description] = ('nan', 'nan')
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
        noiseMap[k[1:-1]] = [p, kp]
        

    with open('ecoli-data.csv', 'wb') as csvfile:
        ecoliwriter = csv.writer(csvfile, delimiter=',')
        ecoliwriter.writerow(['Gene', 'Total Codons', 'Rarity Index', '<P>', 'kp'])
        for k,v in writingDict.items():
            totalCodons, score = v
            
            # Get gene symbol
            l = re.split(' ', k)
            sym = l[1]
            try:
                p, kp = noiseMap[sym.lower()]
            except KeyError:
                continue
            
            ecoliwriter.writerow([k, totalCodons, score, p, kp])
          

if __name__ == '__main__':
    main()

