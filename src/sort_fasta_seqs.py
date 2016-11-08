#!/usr/bin/python
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import Restriction
from Bio.Restriction import *
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna, generic_protein
from Bio.SeqUtils import MeltingTemp as mt
from subprocess import call
from StringIO import StringIO
import os
import sys
import copy

elements_to_seq = {}
seq_names = []

fasta_file = sys.argv[1]

for record in SeqIO.parse(fasta_file, "fasta"):
    seq_name = record.id
    seq_names.append(seq_name)
    seq = record.seq
    # add backbone sequences to dictionary with all genetic elements
    elements_to_seq[seq_name] = seq

seq_names_sorted = sorted(seq_names)

for seq_name in seq_names_sorted:
    print ">%s" % seq_name
    print elements_to_seq[seq_name]
