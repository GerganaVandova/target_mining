#!/usr/bin/env python

# Extract taxa info for each genbank id

from Bio import SeqIO
import glob

taxafile = "taxa.txt"
t = open(taxafile, "w")

gbfilenames = glob.glob("gbdir/*.gb")
# gbfilenames = glob.glob("gbdir/NZ_KE354369.gb")
for gbfilename in gbfilenames:
    # print gbfilename
    f = open(gbfilename, 'r')
    data = f.read()
    for record in SeqIO.parse(open(gbfilename, "rU"), "genbank"):
        gbid = record.id.split(".")[0]
        seq = record.seq
        t.write("%s\t%s" % (gbid, str(record.annotations["taxonomy"])))
        t.write("\n")
        # print gbid
        # print record.annotations["taxonomy"]
