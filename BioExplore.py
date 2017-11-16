from Bio import SeqIO
import numpy as np
import re
import csv

def getCodonMap():
    codonMap = dict()
    with open('codonfreqs.csv', 'rb') as csvfile:
        reader = csv.DictReader(csvfile, delimiter=',')
        for row in reader:
            cdn = row['Codon']
            ratio = row['Ratio']
            freq = row['Freq']
            codonMap[cdn] = {'Ratio': float(ratio), 'Frequency':float(freq)}
    return codonMap

def getSequenceMap():
    seqMap = dict()
    with open('ecoli-nt.fasta', 'rU') as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            seqMap[record.description] = record.seq
    return seqMap



