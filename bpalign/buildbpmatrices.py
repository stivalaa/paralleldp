#!/usr/bin/env python
###############################################################################
#
# buildbpmatrices.py - run RNAfold to build basepair matrices for test input
#
# File:    buildbpmatrices.py
# Author:  Alex Stivala
# Created: June 2009
#
# $Id: buildbpmatrices.py 2612 2009-07-06 00:54:51Z astivala $
#
# The input data is the set of sequence pairs (data set 2) from the
# BRAliBase II data set. See 
#
# Gardner PP, Wilm A & Washietl S (2005)
# A benchmark of multiple sequence alignment
# programs upon structural RNAs. Nucleic Acids Research. 33(8):2433-2439
#
# http://www.binf.ku.dk/~pgardner/bralibase/bralibase2.html
#
#
# Uses BioPython (developed with BioPython version 1.43) for reading Fasta
# records. http://www.biopython.org
# Uses the readrnafold.py module to read RNA sequences from 
# Vienna RNA package RNAFold -p output.
#
# RNAfold must be in PATH
#
# The output is the _ss.ps and _dp.ps files from RNAfold -p in cwd,
# and list of "seq1name seq2name" pairs to  file pairs.list
# where the names are
# the two sequence idneitifers in each input FASTA file.
#
# Note that not all the sequence in the files are unique, i.e. some sequences
# appear in more than one file, but each file has a unique pair of sequences.
# So we check if we have already RNAfolded each sequence and don't 
# do it twice (though doesn't matter much if we did, would j ust overwrite
# old one)
#
###############################################################################

import sys,os,glob
import warnings # so we can suppress the annoying tempnam 'security' warning

from readrnafold import *
from Bio import Fasta

BRALIBASE2 = "/local/cluster/users/astivala/bralibase2/data-set2/unaligned"

listfilename = "pairs.list" # output list of identifier pairs here

warnings.filterwarnings('ignore', 'tempnam', RuntimeWarning) 

fasta_files = glob.glob(os.path.join(BRALIBASE2, '*.fasta'))

seqid_dict = {} # dict of {seqid: True} to mark that identifier as done

parser = Fasta.RecordParser()
listfile = open(listfilename, 'w')
for fasta_file in fasta_files:
    iterator = Fasta.Iterator(open(fasta_file), parser)
    for i in range(2): # only 2 sequences (pairwise)
        # RNAfold requires FASTA input with NO line breaks in 
        # sequence (no matter how long) so we parse input FASTA and write
        # with no line break to temp files
        rec = iterator.next()
        name = rec.title
        seq = rec.sequence
        tmpfile = os.tempnam(None, "rna")
        newrec = Fasta.Record(colwidth = 99999) 
        newrec.title = name
        newrec.sequence = seq
        fd = open(tmpfile, 'w')
        fd.write(str(newrec))
        fd.close()
        listfile.write(name)
        if i == 0:
            listfile.write(" ")
        else:
            listfile.write("\n")
        if not seqid_dict.has_key(name):
            os.system("RNAfold -p < " + tmpfile)
        seqid_dict[name] = True
        os.unlink(tmpfile)
listfile.close()

