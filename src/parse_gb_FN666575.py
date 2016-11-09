#!/usr/bin/env python

# Convert gb file to fasta file, write dna sequence, as not all genbank files
# have predicted protein sequences.

from Bio import SeqIO
import glob

ff = open("FN666575.fasta", "w")

for record in SeqIO.parse(open("dna.fasta", "rU"), "fasta"):
    gbid = record.id
    seq = record.seq
    if gbid == "FN666575":
            print len(seq)
            print seq[3768956:3770566]
            ff.write(">%s" % gbid)
            ff.write("\n")
            ff.write(str(seq))

ff.close()
