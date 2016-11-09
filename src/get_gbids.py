#!/usr/bin/env python
from Bio import SeqIO

filename = "/mnt/gnpn/gnpn/projects/orphanpks/TargetMining/Blast/blast_results_seqs/blast_results.KS.fasta.cleanName.cdhit.99"
outfile = "gbids.txt"

for record in SeqIO.parse(filename, "fasta"):
    recordid = record.id
    # gbid = recordid.rsplit('_', 1)[0]
    gbid = recordid.split('.')[0]
    print gbid
    f = open(outfile, 'a')
    f.write(gbid)
    f.write("\n")

f.close()
