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
import os
import itertools
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio.Emboss.Applications import WaterCommandline

def main(infile_F, infile_R):

    alignCheck = 0
    COMPLEMENT_DICT = {'A': 'T', 'G': 'C', 'T': 'A', 'C': 'G', 'N': 'N'}
    fakeFASTQ = ''
    notAligned = 0
    wrongLength = 0

    # open file containing wild-type sequence and pull it out as a string
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

    # run water on each merged read to align to wt
    with open("fakeFASTQ.fastq", "w") as fakeFASTQ:
        with open("out.extendedFrags.fastq", "rU") as merged:
            for (title, seq, qual) in FastqGeneralIterator(merged):
                # trim sequences with N's
                seq = Ntest(seq)
                # create a file to write each read to as fasta
                with open("read.fasta", "w") as readsfile:
                    readsfile.write(">read\n")
                    readsfile.write(seq)
                # generate water command line and call it
                if os.path.isfile('water.txt'):
                    os.remove('water.txt')
                water_cline = WaterCommandline(asequence="wt.fasta", bsequence="read.fasta", gapopen=10, gapextend=0.5, outfile="water.txt", aformat='markx10')
                stdout, stderr = water_cline()
                # Check whether the read was in the right orientation
                #   If the Identity score of the alignment is low, 
                #   take the revcomp and try aligning again
                alignCheck = 0 
                indexList = []
                with open("water.txt", "rU") as waterfile:
                    for line in waterfile:
                        if len(line.split()) > 1:
                            if line.split()[1] == 'Identity:':
                                identity = line.split()[3]
                                identity = identity[1:4]
                                if float(identity) > 90:
                                    alignCheck = 1
                if not alignCheck:
                    with open("read.fasta", "w") as readfile:
                        readfile.write(">read\n")
                        seq = revcomp(seq, COMPLEMENT_DICT)
                        readfile.write(seq)
                    water_cline = WaterCommandline(asequence="wt.fasta", bsequence="read.fasta", gapopen=10, gapextend=0.5, outfile="water.txt", aformat='markx10')
                    stdout, stderr = water_cline()
                    with open("water.txt", "rU") as waterfile:
                        for line in waterfile:
                            if len(line.split()) > 1:
                                if line.split()[1] == 'Identity:':
                                    identity = line.split()[3]
                                    identity = identity[1:4]
                                    if float(identity) > 90:
                                        alignCheck = 1
                if not alignCheck:
                    notAligned += 1
                    continue
                else:
                    # find wt positions where the read aligns
                    index1, index2 = indexFinder('water.txt')
                    fakeSeq = buildFakeSeq(str(seq), 0, wt, index1, index2, 0, 0)
                    if len(fakeSeq) != len(wt):
                        wrongLength += 1
                        continue
                    fakeFASTQ.write('@' + title + '\n')
                    fakeFASTQ.write(fakeSeq + '\n')
                    fakeFASTQ.write('+\n')
                    fakeQual = ''.join(['A' for ch in fakeSeq])
                    fakeFASTQ.write(fakeQual + '\n')

    print notAligned, wrongLength
                        

#
#    # run water on unmerged reads
#    water_cline = WaterCommandline(asequence="wt.fasta", bsequence="out.notCombined_1.fastq", gapopen=10, gapextend=0.5, outfile="water_F.txt")
#    stdout, stderr = water_cline()
#    indexList1 = indexFinder('water_F.txt')
#    water_cline = WaterCommandline(asequence="wt.fasta", bsequence="out.notCombined_2.fastq", gapopen=10, gapextend=0.5, outfile="water_R.txt")
#    stdout, stderr = water_cline()
#    indexList2 = indexFinder('water_R.txt')
#
#    # iterate over unmerged reads from FLASH
#    with open("out.notCombined_1.fastq", 'rU') as unmerged_F:
#        with open("out.notCombined_2.fastq", 'rU') as unmerged_R:
#            f_iter = FastqGeneralIterator(unmerged_F)
#            r_iter = FastqGeneralIterator(unmerged_R)
#            for (title, seq, qual), (title_R, seq_R, qual_R), (index1, index2), (index3, index4) in itertools.izip(f_iter, r_iter, indexList1, indexList2):
#                fakeFASTQ.append(title + '\n')
#                seq_R_rc = revcomp(seq_R)
#                fakeSeq = buildFakeSeq(seq, seq_R_rc, wt, index1, index2, index3, index4)
#                fakeFASTQ.append(fakeSeq + '\n')
#                fakeFASTQ.append('+\n')
#                fakeQual = ''
#                fakeQual = ['' + 'A' for ch in fakeSeq]
#                fakeQual = ''.join(fakeQual)
#                fakeFASTQ.append(fakeQual)
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
    rc = ''.join([COMPLEMENT_DICT[base] for base in seq])[::-1]
    return rc


def buildFakeSeq(seq_F, seq_R_rc, wt, index1, index2, index3, index4):
    '''Construct a FASTQ compatible with Enrich'''
    if seq_R_rc:
        fakeRead = wt[:index1-1] + seq_F + wt[index2:index3] + seq_R_rc + wt[index4:]
    else:
        fakeRead = wt[:index1-1] + seq_F + wt[index2:]
    return fakeRead 


def indexFinder(infile):
    ''' Searches output file from water for alignment position indexes '''
    with open(infile, 'rU') as waterdata:
        for line in waterdata:
            if len(line.split()) > 1:
                if line.split()[1] == 'al_start:':
                    start = int(line.split()[2])
                if line.split()[1] == 'al_stop:':
                    stop = int(line.split()[2])
                    break
    return start, stop


def Ntest(seq):
    "trim sequences with N's"
    if seq[0] == 'N':
        seq = seq[1:]
    Ntest = 0
    for i, ch in enumerate(seq):
        if ch == 'N':
            Ntest = 1
            break
    if Ntest == 1:
        seq = seq[:i-1]
    return seq


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('infile_F')
    parser.add_argument('infile_R')
    args = parser.parse_args()
    main(args.infile_F, args.infile_R)
