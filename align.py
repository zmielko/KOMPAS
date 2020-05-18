#!/usr/bin/env python
#TF KOMPAS: kmer alignment and orientation
#Author: Zachery Mielko
#Version: 5/17/2020

import argparse
programDescription = 'K.O.M.P.A.S: 1Model Alignment'
parser = argparse.ArgumentParser(description=programDescription)
req = parser.add_argument_group('required arguments')
req.add_argument('-a','--alignModel',required = True, type=str, help='PWM file')
req.add_argument('-k','--kmerFile',required = True, type=str, help='Non-gapped kmer file')
parser.add_argument('-o','--output', type=str, default='' ,help='Save file, stdout by default')
parser.add_argument('-p','--palindrome', action = 'store_true', help='Indicate if the TF is palindromic')
parser.add_argument('-m','--meme', action = 'store_true', help='Input is MEME formated probability matrix')
parser.add_argument('-log', type=str, help='Save location for log')
args = parser.parse_args()

# Parameters
PWMFile = args.alignModel
kmerFile = args.kmerFile
output = args.output
isPalindrome = args.palindrome
isMEME = args.meme
logFile = args.log

# Imports
import pandas as pd
import numpy as np
import os
import sys

# Sequence Manipulation Functions
def revComp(sequence):
    """
    Input: String of DNA sequences in any case
    Output: Reverse Complement in upper case
    """
    # Define dictionary, make input string upper case
    rcDict = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'}
    seqUp = sequence.upper()
    rcSeq = ''
    for letter in seqUp:
        if letter not in rcDict:
            raise ValueError('Error: nucleotide not A, C, G, T, N found in string')
        else:
            rcSeq = rcSeq + rcDict[letter]
    return(rcSeq[::-1])


# Output
if output == '':
    output = sys.stdout
        
# Get kmer data and length
kmerData = pd.read_csv(kmerFile, delim_whitespace=True)
kmerF = kmerData.iloc[:,0]
k = len(kmerF[0])

# Read PWM
if isMEME:
    f=open(PWMFile, "r")
    for num, line in enumerate(f):
        if line.strip().startswith('letter-probability'):
            info = line.strip()
            skip = num + 1
            numRows = int(info.split(' ')[5])
            break
    f.close()
    PWM = pd.read_csv(PWMFile, delim_whitespace=True, header = None, skiprows = (skip), nrows = numRows)
    PWMout = PWM.copy()
    PWMout = PWMout.rename(columns={0:'A:', 1:'C:', 2:'G:', 3:'T:'}).T.reset_index()
    PWM = PWM.T
    PWM = PWM.to_numpy()
else:
    PWM = pd.read_csv(PWMFile, delim_whitespace=True, header = None)
    PWM = PWM.sort_values(by=0)
    PWMout = PWM.copy() # Copy to send to output
    PWM = PWM.to_numpy()
    PWM = np.delete(PWM,0,1).astype('float')  

# Pad the array
eqiProb = np.array(([0.25],[0.25],[0.25],[0.25]))
padding = np.repeat(eqiProb, k, axis=1)
PWM = np.concatenate((padding,PWM,padding), axis=1)
# Convert to dictionary
PWMDict = {'A':PWM[0],'C':PWM[1],'G':PWM[2],'T':PWM[3]}
# Alignment for one kmer
def alignment(kmer):
    """
    Input: kmer (str)
    Output: maximum score (float), kPositions (list of integers)
    Takes a kmer and scans it along a PWM, returns all kPosition with the max score
    """
    best_score = 0
    kposition = 0
    scoreList = []
    for PWMPos, position in enumerate(range(0, PWM.shape[1] - k +1)):
        score = 1
        for letterPos, letter in enumerate(kmer):
            score = score * PWMDict[letter][letterPos + PWMPos]
        scoreList.append(score)
    maxScore = max(scoreList)
    kPosition = [idx for idx, score in enumerate(scoreList) if score == maxScore]
    return((maxScore,kPosition))

# Align for all kmers
oriented_kmer = []
kposition = []
orientation = []
Escore = []
Zscore = []
modelScore = []
classified = []
for kmerScore in zip(kmerF,kmerData.iloc[:,2],kmerData.iloc[:,4]):
    fwd = alignment(kmerScore[0])
    rev = alignment(revComp(kmerScore[0]))
    if isPalindrome == True: # Don't classify by orientation of model
        oriented_kmer.append(kmerScore[0])
        kposition.append(fwd[1])
        orientation.append("F")
        Escore.append(kmerScore[1])
        Zscore.append(kmerScore[2])
        modelScore.append(fwd[0])
        classified.append(1)
        oriented_kmer.append(revComp(kmerScore[0]))
        kposition.append(rev[1])
        orientation.append("R")
        Escore.append(kmerScore[1])
        Zscore.append(kmerScore[2])
        modelScore.append(rev[0])
        classified.append(1)
    elif isPalindrome == False: # Classify Orientation
        oriented_kmer.append(kmerScore[0])
        kposition.append(fwd[1])
        orientation.append("F")
        Escore.append(kmerScore[1])
        Zscore.append(kmerScore[2])
        modelScore.append(fwd[0])
        if fwd[0] >= rev[0]:
            classified.append(1)
        else:
            classified.append(0)
        oriented_kmer.append(revComp(kmerScore[0]))
        kposition.append(rev[1])
        orientation.append("R")
        Escore.append(kmerScore[1])
        Zscore.append(kmerScore[2])
        modelScore.append(rev[0])
        if fwd[0] <= rev[0]:
            classified.append(1)
        else:
            classified.append(0)
        
results = pd.DataFrame({'kmer':oriented_kmer,'kposition':kposition,'orientation':orientation,'classified':classified,'Escore':Escore,'Zscore':Zscore})

# Write header
if output == sys.stdout:
    f = sys.stdout
else:
    f = open(output, "w")
    
if isPalindrome:
    f.write("#Palindrome Alignment\n")
else:
    f.write("#NonPalindrome Alignment\n")
f.write('#PWM Model\n')
if output != sys.stdout:
    f.close()
    
# Write PWM
PWMout.to_csv(output, sep = '\t', header = None, index = False, mode = 'a')

# Write second header
if output != sys.stdout:
    f = open(output, "a")
    f.write("#Aligned kmers\n")
    f.close()
else:
    f.write("#Aligned kmers\n")
    
#Write aligned model
results.to_csv(output, sep = '\t', index = False,mode='a')

# If log file is True
if logFile:
    nAmb = sum(results['kposition'].apply(lambda x: len(x) != 1))
    Amb = results[results['kposition'].apply(lambda x: len(x) != 1)]['kposition']
    Amb = pd.DataFrame(Amb.value_counts())
    Amb = Amb.reset_index().rename(columns={'index':'kPositions', 'kposition':'Count'})
    # Write to Log File
    f = open(log, "w")
    f.write("##### Parameters #####\n")
    f.write(f"PWM file: {PWMFile}\n")
    f.write(f"kmerFile: {kmerFile}\n")
    f.write(f"Palindrome Mode: {isPalindrome}\n")
    f.write(f"Model kPositions: {(k, PWM.shape[1] - (k))}\n")
    f.write(f"Output File: {output}\n")
    f.write("##### Multi-position kmers #####\n")
    f.write(f"Number of kmers: {nAmb}\n")
    f.write(f"Unique positions and number: \n\n")
    f.close()
    Amb.to_csv(log, sep = '\t', mode = "a", index = False)
    
