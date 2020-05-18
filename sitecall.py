#!/usr/bin/env python

# TF KOMPAS: Site Caller
# Author: Zachery Mielko
# Version: 5/18/2020

import argparse
programDescription = 'Calls TFBS from bed/genome or fasta files'
parser = argparse.ArgumentParser(description=programDescription,add_help=False)
req = parser.add_argument_group('parameter arguments')
req.add_argument('-k', '--kmerFile',type = argparse.FileType('r'), default = '-', help='aligned kmer file to use, stdin if not specified')
req.add_argument('-c','--core', type=int,required = True, nargs=2, help='Start and end positions of the core, relative to the model (2 values)')
req.add_argument('-o','--output', type=str,required = True, help='Output file')
# Bed or fasta
search = parser.add_argument_group('search space arguments [-b/-g or -f]')
search.add_argument('-b','--bed', type=str, help='Coordinates to search (.bed)')
search.add_argument('-g','--genome', type=str,help='Genome file (.fasta/.fa)')
search.add_argument('-f','--fasta', type=str, help='Fasta file (.fasta/.fa)')
# Optional
opt = parser.add_argument_group('optional arguments')
opt.add_argument('-gc','--gcParse',action='store_true', help='If using a fasta input with genomic coordinates, returns bed file of coordinates')
opt.add_argument('-t','--threshold', type=float,default = 0.4, help='Escore threshold for calling (float, -0.5 to 0.5)')
opt.add_argument('-l', '--log',type=str, help='Generate a log file')
opt.add_argument('-rM','--rankModel' ,action='store_true', help='add alignment model score')
opt.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS, help='Show this help message and exit.')
                    
args = parser.parse_args()
if args.bed is None and args.genome is None and args.fasta is None:
    parser.error("KOMPAS requires either (-b and -g) or (-f)")
if (args.bed and args.genome is None) or (args.genome and args.bed is None):
    parser.error("-b and -g need to be specified together")
if (args.fasta and args.genome) or (args.fasta and args.bed):
    parser.error("-f argument cannot be used with -b or -g")
if (args.gcParse and args.genome) or (args.gcParse and args.bed):
    parser.error("-gc argument cannot be used with -b or -g")
    
output = args.output
kmerFile = args.kmerFile
isFasta = args.fasta
if args.bed:
    peakFile = args.bed
    genomeFile = args.genome
elif args.fasta:
    gcParse = args.gcParse
    fastaFile = args.fasta 
userCore = args.core
threshold = args.threshold
logFile = args.log
rM = args.rankModel

# Imports
import pandas as pd
import numpy as np
import os
import sys
from io import StringIO
import itertools
from pyfaidx import Fasta

# Store input data and parse
inputData = kmerFile.read()
kmerFile.close()

# Parse palindrome setting
f = StringIO(inputData)
palSetting = f.readline().strip()
if not palSetting.startswith('#Palindrome') and not palSetting.startswith('#Non'):
    raise ValueError('Invalid input, file needs to be the output from a KOMPAS alignment tool')
if palSetting.startswith('#Palindrome'):
    isPalindrome = True
else:
    isPalindrome = False
f.close()

# Read in kmer data
def convList(string):
    """
    Input: String of a list of integers
    Output: List of integers
    """
    toList = string[1:-1].split(', ')
    toInts = [int(i) for i in toList]
    return(toInts)
kmer = pd.read_csv(StringIO(inputData), sep = '\t',converters={'kposition': lambda x: convList(x)}, skiprows = 7)
k = len(kmer['kmer'][0])

# Read in model length
PWM = pd.read_csv(StringIO(inputData), delim_whitespace=True, skiprows = 2, nrows = 4, header = None)
mStart = k
mLength = (len(PWM.columns) -1)
mEnd = mStart + mLength

# Check if given positions are valid
if userCore[0] > mLength or userCore[1] > mLength:
    raise ValueError('Core position(s) greater than model length')
if userCore[1] <= userCore[0]:
    raise ValueError('Core start must be greater than core end')

# Define core in kPosition
core = []
core.append(mStart + (userCore[0] -1))
core.append(mStart + (userCore[1]))
coreLen = core[1] - core[0]
centerPos = mStart

# Find the kPositions required, any would be sufficient to call
if k > coreLen: 
    searchEnd = core[1]
    checkK = 0
    ReqKpos = set() #
    while checkK != core[0]:
        checkK = searchEnd - k
        if checkK <= core[0]:
            ReqKpos.add(checkK)
            searchEnd = searchEnd + 1
            
# Or find the group of all kPositions that are needed, all or none
else:
    searchStart = core[0]
    checkK = 0
    ReqKpos = set()
    while searchStart + k <= core[1]:
        ReqKpos.add(searchStart)
        searchStart = searchStart + 1
        
# Determine flanks of ReqKPos for threshold score reporting
ScoredKpos = ReqKpos.copy()
if k >= coreLen:
    ScoredKpos.add(min(ReqKpos) - 1)
    ScoredKpos.add(max(ReqKpos) + 1)  
    
# Filter dataframe for kmers with ScoredKpos
def filtPos(inputList, ScoredKpos = ScoredKpos):
    """
    Input: List of integers (kPosition)
    Output: Intersected list of ScoredKPos
    """
    result = []
    for position in inputList:
        if position in ScoredKpos:
            result.append(position)
    return(result)

thrKmers = kmer.copy(deep = True)
thrKmers = thrKmers.query("classified == 1") # Use classified kmers
thrKmers['kposition'] = thrKmers['kposition'].apply(lambda x: filtPos(x)) # Only ScoredKpos
thrKmers = thrKmers[thrKmers['kposition'].apply(lambda x: len(x) != 0)] # Remove empty lists
thrKmers = thrKmers[thrKmers['Escore'] >= threshold]
kDict = dict(zip(thrKmers['kmer'],zip(thrKmers['kposition'],thrKmers['Escore'])))

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

# Scoring and Calling Functions
def kmerMatch(seq):
    """
    Input: Sequence to search for kPos and Escore
    Output[0] = Lists of all possible kPosition lists
    Output[1] = Lists of all possible Escore lists
    """
    kPosLists = []
    escoreLists = []
    def recursiveKmerMatch(seq, crntKposList,crntEList, seqPos):
        """
        Input: DNA sequence, list of previous kPositions, Escores, and
        current position in the sequence
        Output: Appends kPositions and Escores to initiallized empty list
        in kmerMatch
        """
        for i in range(len(seq) - k + 1 - seqPos): #range of current fork
            iPos = i + seqPos # start from the current fork
            window = seq[iPos:iPos+k]
            if window in kDict:
                if len(kDict[window][0]) == 1: # If only 1 kPos
                    crntKposList.append(kDict[window][0][0])
                    crntEList.append(kDict[window][1])
                else:
                    for j in kDict[window][0]:
                        frkdKposList = crntKposList.copy()
                        frkdEList = crntEList.copy()
                        frkdKposList.append(j)
                        frkdEList.append(kDict[window][1])
                        recursiveKmerMatch(seq, frkdKposList,frkdEList, iPos+1)
                    return None
            else:
                crntKposList.append(0)
                crntEList.append(-0.5)
        kPosLists.append(crntKposList)
        escoreLists.append(crntEList)
    recursiveKmerMatch(seq,[],[], 0)
    return((kPosLists, escoreLists))

def findConsArrays(matchOutput):
    """
    Input: Output from kmerMatch
    Output: Sequence positions, kPositions
    """
    consArrays = []
    for kpos, kscore in zip(matchOutput[0], matchOutput[1]):
        kpos = np.array(kpos)
        kscore = np.array(kscore)
        if k >= coreLen:
            position = list(filter(lambda x: len(x) != 1,np.split(np.r_[:len(kpos)], np.where(np.diff(kpos) != 1)[0]+1)))
            kpos = list(filter(lambda x: len(x) != 1,np.split(kpos, np.where(np.diff(kpos) != 1)[0]+1)))
        elif k < coreLen:
            reqLen = len(ReqKpos)
            position = list(filter(lambda x: len(x) == reqLen,np.split(np.r_[:len(kpos)], np.where(np.diff(kpos) != 1)[0]+1)))
            kpos = list(filter(lambda x: len(x) == reqLen,np.split(kpos, np.where(np.diff(kpos) != 1)[0]+1)))
        kScore = []
        for pos in position:
            kScore.append(kscore[pos])
        consArrays.append(zip(position, kpos, kScore))
    return(consArrays)

def findCenter(zipList, orient, seqLen):
    """
    Given a zip of match position, kposition, and kscore
    Returns the center sites
    """
    centerSites = set()
    for zipMatch in zipList:
        for pos, kpos, kScore in zipMatch:
            centerSite = (centerPos - kpos[0]) + pos[0]
            if orient == 'rc':
                centerSite = (seqLen - centerSite) -1
            centerSites.add(centerSite)
    return(sorted(list(centerSites)))

def siteCallWrap(seq, orient, seqLen):
    """
    Wraps kmerMatch, findConsArrays, and findCenter together
    Input: sequence, orientation, sequence length
    Output: list of called sites
    """
    return(findCenter(findConsArrays(kmerMatch(seq)), orient, seqLen))

# read model information if need be
if rM:
    PWM = pd.read_csv(StringIO(inputData), delim_whitespace=True, skiprows = 2, nrows = 4, header = None)
    PWM = PWM.sort_values(by=0)
    PWM = PWM.to_numpy()
    PWM = np.delete(PWM,0,1).astype('float') 
    PWMDict = {'A':PWM[0],'C':PWM[1],'G':PWM[2],'T':PWM[3]}
    def PWMScore(kmer):
        """
        Input: Sequence with length == model length
        Output: score
        """
        score = 1
        for letterPos, letter in enumerate(kmer):
            score = score * PWMDict[letter][letterPos]
        return(score) 

def parseID(inputID):
    """
    Takes the concatinated names given in fasta outputs from bedtool's getfasta 
    and turns them into bed compatible columns
    Input: ID from fasta file generated from a bed file
    Output: chromosome and start position
    """
    cvp = inputID.split(':')
    pos = cvp[1].split('-')
    return((cvp[0], int(pos[0]), int(pos[1])))

# File reading and parsing functions, wrap calling functions
def writeCalls(siteList,ID,orient,o,isFasta,seqLen, seq):
    """
    Outputs the results from genCallingFasta to a file
    Filters out sites where:
    1. The model flanking the core is beyond the search space
    2. There are ambiguous bases 'N' within the model range
    """
    for siteNum, center in enumerate(siteList):
        # Get ID and position from fasta or bed/genome files
        if isFasta and gcParse == False:
            chrom = ID
            IDout = ID
            if orient == '+':
                start = center
            else:
                start = (seqLen - (center+mLength))
            site = seq[start:start+mLength]
            if orient == '-':
                        site = revComp(site)
            viableStart = 0
            viableEnd = seqLen
        elif isFasta and gcParse == True:
            parseList = parseID(ID)
            IDout = ID
            chrom = parseList[0]
            if orient == '+':
                start = parseList[1] + center
            else:
                start = parseList[1] + (seqLen - (center+mLength))
            site = seq[start - parseList[1]:start+mLength - parseList[1]]
            if orient == '-':
                        site = revComp(site)
            viableStart = parseList[1]
            viableEnd = parseList[1] + seqLen
        else:
            parseList = ID
            chrom = parseList[0]
            IDout = f'{ID[0]}:{ID[1]}-{ID[2]}'
            if orient == '+':
                start = parseList[1] + center
            else:
                start = parseList[1] + (seqLen - (center+mLength))
            site = seq[start - parseList[1]:start+mLength - parseList[1]]
            if orient == '-':
                        site = revComp(site)
            viableStart = parseList[1]
            viableEnd = parseList[1] + seqLen            
        #If viable, write output
        if 'N' not in site and start >= viableStart and start+mLength <= viableEnd:
            # If palindrome, write '.' instead of +
            if isPalindrome:
                orientOut = '.'
            else:
                orientOut = orient
            # If adding PWM score, generate score
            if rM:
                score = ''
                score = PWMScore(site)
                o.write(f"{chrom}\t{start}\t{start + mLength}\t{IDout}_{orientOut}_{siteNum}\t1000\t{orientOut}\t{score}\n") 
            else:
                o.write(f"{chrom}\t{start}\t{start + mLength}\t{IDout}_{orientOut}_{siteNum}\t1000\t{orientOut}\n")

            
# Overwrite output name if need be
o = open(output, "w")  
o.close()
# Start appending to new output name
o = open(output, "a")  

# Find the minimum length of a sequence to scan
if mLength > k+1:
    minLength = mLength
else:
    minLength = k + 1

if isFasta:
    fastaRead = open(fastaFile, 'r')
    for header,seq in itertools.zip_longest(fastaRead, fastaRead):
        if header.startswith('>') and seq:
            ID = header.rstrip()[1:]
            seq = seq.rstrip()
            seq = seq.upper()
            seqLen = len(seq)
            # If possible, look at forward sequences always. If not a palindrome, look at revComp
            if seqLen >= minLength:
                sitesFwd = siteCallWrap(seq, '+', seqLen)
                if len(sitesFwd) != 0:
                    writeCalls(sitesFwd, ID, '+', o, isFasta, seqLen=seqLen, seq=seq)
                if isPalindrome == False:
                    sitesRC = siteCallWrap(revComp(seq), '-', seqLen)
                    if len(sitesRC) != 0:
                        writeCalls(sitesRC, ID, '-', o, isFasta, seqLen=seqLen, seq=seq) 
    fastaRead.close()
else:
    # pyfaidx indexing of genome file
    genome = Fasta(genomeFile)
    peakFileIO = open(peakFile, 'r')
    for line in peakFileIO:
        col = line.split('\t')
        ID= [col[0], int(col[1]), int(col[2])]
        seq = str(genome[ID[0]][ID[1]:ID[2]])
        seq = seq.upper()
        seqLen = len(seq)
        # If possible, look at forward sequences always. If not a palindrome, look at revComp
        if seqLen >= minLength:
            sitesFwd = siteCallWrap(seq, '+', seqLen)
            if len(sitesFwd) != 0:
                writeCalls(sitesFwd, ID, '+', o, isFasta, seqLen=seqLen, seq=seq)
            if isPalindrome == False:
                sitesRC = siteCallWrap(revComp(seq), '-', seqLen)
                if len(sitesRC) != 0:
                    writeCalls(sitesRC, ID, '-', o, isFasta, seqLen=seqLen, seq=seq) 
    peakFileIO.close()
o.close()    
        

# Log File
if logFile:
    f = open(logFile, "a")
    f.write("##### Parameters ##### \n")
    if isFasta == True:
        f.write(f"Fasta file: {fastaFile}"+ "\n")
    else:
        f.write(f"Peak file: {peakFile}"+ "\n")
        f.write(f"Genome file: {genomeFile}"+ "\n")
    f.write(f"User defined core: {userCore}"+ "\n")
    f.write(f"kPosition core: {core}"+ "\n")
    f.write(f"Model Position: {mStart} {mEnd}")
    f.write(f"Threshold: {threshold}"+ "\n")
    f.write(f"Output file: {output}"+ "\n")
    if rP:
        f.write(f"Sites scored by alignment models")
    else:
        f.write(f"Default output")
    f.close()   

