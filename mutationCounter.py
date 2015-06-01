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
import collections
import itertools
import copy
from Bio.SeqIO.QualityIO import FastqGeneralIterator


def main(infileF, infileR):

    COMPLEMENT_DICT = {'A': 'T', 'G': 'C', 'T': 'A', 'C': 'G', 'N': 'N'}
    
    #take trimmed paired-end FASTQ files as input
    #run FLASH to combine overlapping read pairs
    # or remove this line here and add to shell script
    subprocess.call(["flash", "-M", "140", "-t", "1", infile_F, infile_R])

    # iterate over merged reads from FLASH
    with open("out.extendedFrags.fastq", 'rU') as merged:
        for title, seq, qual in FastqGeneralIterator(merged):
            buildFakeFASTQ(seq, 0, wt)

    # iterate over unmerged reads from FLASH
    with open("out.notCombined_1.fastq", 'rU') as unmerged_F:
        with open("out.notCombined_2.fastq", 'rU') as unmerged_R:
            f_iter = FastqGeneralIterator(unmerged_F)
            r_iter = FastqGeneralIterator(unmerged_R)
            for (title, seq, qual), (title_R, seq_R, qual_R) in itertools.izip(f_iter, r_iter):
                fakeFASTQ.append(title + '\n')
                seq_R_rc = revcomp(seq_R)
                fakeSeq = buildFakeSeq(seq_F, seq_R_rc, wt)
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
    rc = 
    return rc


def buildFakeSeq(seq_F, seq_R_rc, wtseq):
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
                indexList.append(line.split()[1])
                indexList.append(line.split()[3])
    return indexList

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('infile_F')
    parser.add_argument('infile_R')
    main(args.infile_F, args.infile_R)
