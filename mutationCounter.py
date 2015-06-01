#!/usr/bin/env python

#############################################################################
#   mutationCounter.py
#   2015 James A. Stapleton, Justin R. Klesmith
#
#   This program takes short reads from shotgun sequencing of mutant
#       libraries and creates FASTQ files compatible with ENRICH.
#       
#
#############################################################################

import argparse
import time
import subprocess
import itertools
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio.Emboss.Applications import WaterCommandline

def main(infile_F, infile_R):

    COMPLEMENT_DICT = {'A': 'T', 'G': 'C', 'T': 'A', 'C': 'G', 'N': 'N'}
    fakeFASTQ = ''
    with open('wt.fasta', 'rU') as wildtype:
        wildtype = wildtype.read()
    if wildtype[0] == ">":
        wildtype = wildtype.split('\n',1)[1]
    wt = [''.join(line) for line in wildtype.split('\n')]
    wt = ''.join(wt)

    #take trimmed paired-end FASTQ files as input
    #run FLASH to combine overlapping read pairs
    # or remove this line and add to shell script
    subprocess.call(["flash", "-M", "140", "-t", "1", infile_F, infile_R])

    # run water on merged reads
    water_cline = WaterCommandline(asequence="wt.fasta", bsequence="out.extendedFrags.fastq", gapopen=10, gapextend=0.5, outfile="waterA.txt")
    stdout, stderr = water_cline()

    water_cline = WaterCommandline(asequence="wt.fasta", bsequence="out.extendedFrags.fastq", gapopen=10, gapextend=0.5, outfile="waterB.txt")
    stdout, stderr = water_cline()
    with open('waterA.txt', 'rU') as waterA:
        with open('waterB.txt', 'rU') as waterB:
            for lineA, lineB in itertools.izip(waterA, waterB):
                if lineA.split()[1] == 'Identity:':
                    
    indexList = indexFinder('water.txt')

    # iterate over merged reads from FLASH
    with open("out.extendedFrags.fastq", 'rU') as merged:
        f_iter = FastqGeneralIterator(merged)
        for (title, seq, qual), (index1, index2) in itertools.izip(f_iter, indexList):
            fakeFASTQ.append(title + '\n')
            fakeSeq = buildFakeSeq(seq, 0, wt, index1, index2, 0, 0)
            fakeFASTQ.append(fakeSeq + '\n')
            fakeFASTQ.append('+\n')
            fakeQual = ''
            fakeQual = ['' + 'A' for ch in fakeSeq]
            fakeQual = ''.join(fakeQual)
            fakeFASTQ.append(fakeQual)

    # run water on unmerged reads
    water_cline = WaterCommandline(asequence="wt.fasta", bsequence="out.notCombined_1.fastq", gapopen=10, gapextend=0.5, outfile="water_F.txt")
    stdout, stderr = water_cline()
    indexList1 = indexFinder('water_F.txt')
    water_cline = WaterCommandline(asequence="wt.fasta", bsequence="out.notCombined_2.fastq", gapopen=10, gapextend=0.5, outfile="water_R.txt")
    stdout, stderr = water_cline()
    indexList2 = indexFinder('water_R.txt')

    # iterate over unmerged reads from FLASH
    with open("out.notCombined_1.fastq", 'rU') as unmerged_F:
        with open("out.notCombined_2.fastq", 'rU') as unmerged_R:
            f_iter = FastqGeneralIterator(unmerged_F)
            r_iter = FastqGeneralIterator(unmerged_R)
            for (title, seq, qual), (title_R, seq_R, qual_R), (index1, index2), (index3, index4) in itertools.izip(f_iter, r_iter, indexList1, indexList2):
                fakeFASTQ.append(title + '\n')
                seq_R_rc = revcomp(seq_R)
                fakeSeq = buildFakeSeq(seq, seq_R_rc, wt, index1, index2, index3, index4)
                fakeFASTQ.append(fakeSeq + '\n')
                fakeFASTQ.append('+\n')
                fakeQual = ''
                fakeQual = ['' + 'A' for ch in fakeSeq]
                fakeQual = ''.join(fakeQual)
                fakeFASTQ.append(fakeQual)
    return 0


######## Function definitions ##############

def timer(readCount, start_time):
    # Print time elapsed every 100000 reads processed
    readCount += 1
    if readCount % 100000 == 0:
        print readCount
        print time.time() - start_time, "seconds"
        start_time = time.time()
    return 0


def revcomp(seq, COMPLEMENT_DICT):
    rc = []
    rc.append(''.join([COMPLEMENT_DICT[base] for base in seq])[::-1])
    return rc


def buildFakeSeq(seq_F, seq_R_rc, wt, index1, index2, index3, index4):
    '''Construct a FASTQ compatible with Enrich'''
    if seq_R_rc:
        fakeRead = ''.join(wt[:index1], seq_F, wt[index2:index3], seq_R_rc, wt[index4:])
    else:
        fakeRead = ''.join(wt[:index1], seq_F, wt[index2:])
    return fakeRead


def indexFinder(infile):
    ''' Searches output file from water for alignment position indexes '''
    with open(infile, 'rU') as waterdata:
        indexList = []
        for line in waterdata:
            if line[:2] == 'wt':
                indexList.append(line.split()[1], line.split()[3])
    return indexList

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('infile_F')
    parser.add_argument('infile_R')
    args = parser.parse_args()
    main(args.infile_F, args.infile_R)
