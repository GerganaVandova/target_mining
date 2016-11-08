#!/usr/bin/python
import sys
import os
from string import *
from Bio import SeqIO
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq

PRIMER_TM = 60  # desired melting temperature for all primers generated
PRIMER_NUM = 117  # primer number to begin at


comp = {'a': 't', 't': 'a', 'c': 'g', 'g': 'c'}


def revcomplement(seq):
    return "".join([comp[base] for base in reversed(seq)])


def get_primer(seq, direction, name):
    # Tm_NN: Calculation based on nearest neighbor thermodynamics. Several
    # tables for DNA/DNA, DNA/RNA and RNA/RNA hybridizations are included.
    # Correction for mismatches, dangling ends, salt concentration and other
    # additives are available.
    # Tm_staluc is the 'old' NN calculation and is kept for compatibility.
    # It is, however, recommended to use Tm_NN instead, since Tm_staluc may be
    # depreceated in the future. Also, Tm_NN has much more options. Using
    # Tm_staluc and Tm_NN with default parameters gives (essentially) the same results.

    global PRIMER_NUM
    global PRIMER_TM
    PRIMER_LENGTH = 15  # min primer lenght
    if direction == "fwd":
        while mt.Tm_staluc(seq[0:PRIMER_LENGTH]) <= PRIMER_TM and PRIMER_LENGTH <= 65:
            PRIMER_LENGTH += 1
        primer_seq = seq[0:PRIMER_LENGTH]
        primer_tm = mt.Tm_staluc(primer_seq)
    elif direction == "rev":
        while mt.Tm_staluc(seq[-PRIMER_LENGTH:]) <= PRIMER_TM and PRIMER_LENGTH <= 65:
            PRIMER_LENGTH += 1
        primer_seq = revcomplement(seq[-PRIMER_LENGTH:]).lower()
        primer_tm = mt.Tm_staluc(primer_seq)
    primer_seq = str(primer_seq)
    primer_name = "{}_{}_{}".format(PRIMER_NUM, name, direction)
    primer = list([primer_name, primer_seq, primer_tm, PRIMER_LENGTH])
    PRIMER_NUM += 1
    return primer

fasta_file = sys.argv[1]
primers = []

for record in SeqIO.parse(fasta_file, "fasta"):
    primer_name = record.id
    seq = record.seq
    seq = seq.lower()
    primer_fwd = get_primer(seq, "fwd", primer_name)
    primer_rev = get_primer(seq, "rev", primer_name)
    print("%s_%0.1fTm_%sb\t%s" % (primer_fwd[0], primer_fwd[2], primer_fwd[3], primer_fwd[1]))
    print("%s_%0.1fTm_%sb\t%s" % (primer_rev[0], primer_rev[2], primer_rev[3], primer_rev[1]))
    # print primer_fwd
    # print primer_rev
    primers.append(primer_fwd)
    primers.append(primer_rev)


# Write primers in a file. Thus far it doesn't work, >> output in a file
with open((fasta_file + "_primers_generated.txt"), "w") as f:
    for primer in primers:
        f.write("%s_%0.1f_Tm_%sb\t%s\n" % (primer[0], primer[2], primer[3], primer[1]))

# #
# # Primer 38_CP12600-Fr4_3_rev_60.6_Tm_62b shrinked to 60b
# mystring1 = 'aacattctttctgattcaatatgtatatctccttcttatacttaactaatatactaagat'  # Tm = 59
# print('%0.2f' % mt.Tm_staluc(mystring1))
#
# # Primer 39_CP12600-Fr4_4_fwd_60.3_Tm_65b srinked to 60b
# mystring2 = 'atcttagtatattagttaagtataagaaggagatatacatattgaatcagaaagaatgtt'  # Tm = 59
# print('%0.2f' % mt.Tm_staluc(mystring2))
