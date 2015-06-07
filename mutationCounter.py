#!/usr/bin/env python

#############################################################################
#   mutationCounter.py
#   2015 James A. Stapleton, Justin R. Klesmith
#
#   This program takes short reads from shotgun sequencing of mutant
#       libraries and creates FASTQ files compatible with ENRICH.
#
#   The program fills in wild-type sequence around the short reads
#       to create full-length sequences, and puts in a made-up
#       quality score.
#
#   Overlapping read pairs are merged by FLASH.
#
#############################################################################


import argparse
import subprocess
import os
import itertools
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio.Emboss.Applications import WaterCommandline


def main(forward_paired, forward_unpaired, reverse_paired, reverse_unpaired):

    fakeFASTQ = ''
    notAligned = 0
    wrongLength = 0

    # open file containing wild-type sequence and pull it out as a string
    wt = wtParser()

    # take trimmed paired-end FASTQ files as input
    # run FLASH to combine overlapping read pairs
    # or remove this line and add to shell script
    subprocess.call(["flash", "-M", "140", "-t", "1", forward_paired, reverse_paired])

    # merged read pairs
    notAligned = align_and_index('out.extendedFrags.fastq', notAligned)
    with open("fakeFASTQ.fastq", "w") as fakeFASTQ:
        with open('indexes.txt', 'rU') as indexes:
            with open('out.extendedFrags.fastq', 'rU') as merged:
                f_iter = FastqGeneralIterator(merged)
                for (title, seq, qual), indexline in itertools.izip(f_iter, indexes):
                    index1, index2, rc_flag = indexline.split()
           #         print title, seq, qual, index1, index2, rc_flag
                    if index1 and index2:
                        if int(rc_flag):
                            seq = revcomp(seq)
                        fakeSeq = buildFakeSeq(seq, 0, wt, index1, index2, 0, 0)
                        if len(fakeSeq) != len(wt):
                            wrongLength += 1
           #                 print fakeSeq
           #                 print rc_flag, seq, index1, index2
                            continue
                        fakeFASTQwriter(fakeSeq, title, fakeFASTQ)

    notAligned = align_and_index(forward_unpaired, notAligned)
    with open("fakeFASTQ.fastq", "a") as fakeFASTQ:
        with open('indexes.txt', 'rU') as indexes:
            with open(forward_unpaired, "rU") as merged:
                f_iter = FastqGeneralIterator(merged)
                for (title, seq, qual), indexline in itertools.izip(f_iter, indexes):
                    index1, index2, rc_flag = indexline.split()
                    if index1 and index2:
                        if int(rc_flag):
                            seq = revcomp(seq)
                        fakeSeq = buildFakeSeq(seq, 0, wt, index1, index2, 0, 0)
                        if len(fakeSeq) != len(wt):
                            wrongLength += 1
                            continue
                        fakeFASTQwriter(fakeSeq, title, fakeFASTQ)

    notAligned = align_and_index(reverse_unpaired, notAligned)
    with open("fakeFASTQ.fastq", "a") as fakeFASTQ:
        with open('indexes.txt', 'rU') as indexes:
            with open(reverse_unpaired, "rU") as merged:
                f_iter = FastqGeneralIterator(merged)
                for (title, seq, qual), indexline in itertools.izip(f_iter, indexes):
                    index1, index2, rc_flag = indexline.split()
                    if index1 and index2:
                        if int(rc_flag):
                            seq = revcomp(seq)
                        fakeSeq = buildFakeSeq(seq, 0, wt, index1, index2, 0, 0)
                        if len(fakeSeq) != len(wt):
                            wrongLength += 1
                            continue
                        fakeFASTQwriter(fakeSeq, title, fakeFASTQ)


#
#
#        # unmerged (non-overlapping) read pairs
#        with open("out.notCombined_1.fastq", 'rU') as unmerged_F:
#            with open("out.notCombined_2.fastq", 'rU') as unmerged_R:
#                f_iter = FastqGeneralIterator(unmerged_F)
#                r_iter = FastqGeneralIterator(unmerged_R)
#                for (title, seq, qual), (title_R, seq_R, qual_R) in itertools.izip(f_iter, r_iter):
#                    index1, index2, notAligned, seq = align_and_index(seq, notAligned)
#                    if index1 and index2:
#                        index3, index4, notAligned, seq_R = align_and_index(seq_R, notAligned)
#                        if index3 and index4:
#                            fakeSeq = buildFakeSeq(seq, seq_R, wt, index1, index2, index3, index4)
#                            if len(fakeSeq) != len(wt):
#                                wrongLength += 1
#    #                            print fakeSeq
#     #                           print index1, index2, index3, index4
#      #                          print seq, seq_R
#                                continue
#                            fakeFASTQwriter(fakeSeq, title, fakeFASTQ)

    print notAligned, wrongLength

    return 0


######## Function definitions ##############


def revcomp(seq):
    """Returns the reverse complement of a DNA sequence."""
    COMPLEMENT_DICT = {'A': 'T', 'G': 'C', 'T': 'A', 'C': 'G', 'N': 'N'}
    rc = ''.join([COMPLEMENT_DICT[base] for base in seq])[::-1]
    return rc


def buildFakeSeq(seq_F, seq_R, wt, index1, index2, index3, index4):
    """Builds a fake full-length DNA sequence line consisting of one
    merged read or two short reads filled in with wild-type sequence.
    """
    index1 = int(index1)
    index2 = int(index2)
    index3 = int(index3)
    index4 = int(index4)
    if seq_R:
        diff = 0
        if index1 < index3:
            if index2 > index3 - 1:
                diff = index2 - index3 + 1
                index2 = index3 - 1
            fakeRead = wt[:index1 - 1] + seq_F + wt[index2:index3 - 1] + seq_R[diff:] + wt[index4:]
        else:
            if index4 > index1 - 1:
                diff = index4 - index1 + 1
                index4 = index1 -1
            fakeRead = wt[:index3 - 1] + seq_R + wt[index4:index1 - 1] + seq_F[diff:] + wt[index2:]
    else:
        fakeRead = wt[:index1-1] + seq_F + wt[index2:]
    return fakeRead.upper()


def index_finder(line):
    """Searches the water output line
    for alignment position indexes.
    """
    index = 0
    if len(line.split()) > 1:
        if line.split()[1] == 'al_start:':
            index = int(line.split()[2])
        elif line.split()[1] == 'al_stop:':
            index = int(line.split()[2])
    return index


def Ntest(seq):
    """Trims sequences with N's.
    Removes N from the first position,
    truncates the sequence at subsequent N's.
    """
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


def runWater(fastq, out):
    """Removes existing water.txt file,
    generates a water command line, and runs water.
    """
    if os.path.isfile(out):
        os.remove(out)
    water_cline = WaterCommandline(asequence="wt.fasta", bsequence=fastq, gapopen=10, gapextend=0.5, outfile=out, aformat='markx10')
    stdout, stderr = water_cline()
    return 0


def wtParser():
    """Takes a wild-type DNA sequence in FASTA format
    and reads it into a string.
    """
    with open('wt.fasta', 'rU') as wildtype:
        wildtype = wildtype.read()
    if wildtype[0] == ">":
        wildtype = wildtype.split('\n', 1)[1]
    wt = ''.join([line.strip() for line in wildtype.split('\n')])
    return wt


def identity_finder(line):
    identity = 0
    if len(line.split()) > 1:
        if line.split()[1] == 'Identity:':
            identity = line.split()[3]
            identity = identity[1:4]
    return identity


def align_and_index(fastq, notAligned):
    """Runs a pipeline to align a sequence (merged or unmerged
    sequencing reads) to a wild-type reference with the EMBOSS
    water local alignment program, align the reverse complement
    if the first alignment was poor, and parse and return the
    wt positions where the alignment begins and ends.
    """
    # generate water command line and call it
    runWater(fastq, 'water_fwd.txt')
    # reverse-complement the reads in fastq
    with open('fastq_rc.fastq', 'w') as reverse:
        with open(fastq, 'rU') as forward:
            next_line_is_seq = 0
            for line in forward:
                if next_line_is_seq:
                    line_rc = revcomp(line.strip())
                    reverse.write(line_rc + '\n')
                    next_line_is_seq = 0
                elif line[0] == '@':
                    next_line_is_seq = 1
                    reverse.write(line)
                else:
                    reverse.write(line)
    # run water on the reverse complements
    runWater('fastq_rc.fastq', 'water_rc.txt')
    # Write only the index
    #  and identity lines to new files
    with open('water_fwd.txt', 'rU') as forward:
        with open('water_fwd_indexes.txt', 'w') as forward_index_lines:
            for line in forward:
                if identity_finder(line) or index_finder(line):
                    forward_index_lines.write(line)
    with open('water_rc.txt', 'rU') as forward:
        with open('water_rc_indexes.txt', 'w') as forward_index_lines:
            for line in forward:
                if identity_finder(line) or index_finder(line):
                    forward_index_lines.write(line)
    # Check whether the read was in the right orientation:
    # Iterate over the water outfiles and pick the best match
    # Write the alignment start and stop of the best matches
    with open('water_fwd_indexes.txt', 'rU') as forward:
        with open('water_rc_indexes.txt', 'rU') as reverse:
            with open('indexes.txt', 'w') as outfile:
                find_index_F = 0
                find_index_R = 0
                index1 = 0
                index2 = 0
                for line_F, line_R in itertools.izip(forward, reverse):
                    if not find_index_F and not find_index_R:
                        identity_F = identity_finder(line_F)
                        identity_R = identity_finder(line_R)
                        if float(identity_F) > 90:
                            find_index_F = 1
                            rev_flag = 0
                        elif float(identity_R) > 90:
                            find_index_R = 1
                            rev_flag = 1
                        elif identity_F and identity_R:
                            outfile.write('0 0 0\n')
                            notAligned += 1
                    elif find_index_F:
                        if not index1 and not index2:
                            index1 = index_finder(line_F)
                        elif index1:
                            index2 = index_finder(line_F)
                            outfile.write(str(index1) + ' ' + str(index2) + ' ' + str(rev_flag) + '\n')
                            find_index_F = 0
                            index1 = 0
                            index2 = 0
                    elif find_index_R:
                        if not index1 and not index2:
                            index1 = index_finder(line_R)
                        elif index1:
                            index2 = index_finder(line_R)
                            outfile.write(str(index1) + ' ' + str(index2) + ' ' + str(rev_flag) + '\n')
                            find_index_R = 0
                            index1 = 0
                            index2 = 0
    return notAligned


def fakeFASTQwriter(fakeSeq, title, handle):
    """Writes the four lines of a fake FASTQ."""
    handle.write('@' + title + '\n')
    handle.write(fakeSeq + '\n')
    handle.write('+\n')
    fakeQual = ''.join(['A' for ch in fakeSeq])
    handle.write(fakeQual + '\n')
    return 0


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('forward_paired')
    parser.add_argument('forward_unpaired')
    parser.add_argument('reverse_paired')
    parser.add_argument('reverse_unpaired')
    args = parser.parse_args()
    main(args.forward_paired, args.forward_unpaired, args.reverse_paired, args.reverse_unpaired)
