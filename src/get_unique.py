#!/usr/bin/python
import sys

# Append taxa to out file. To execute script:
# ./get_unique.py ../data/out.targets.7.filtered.10000.all.pks.nrps ../data/out.targets.7.filtered.10000.all.pks.nrps.unique

# 17889	mfM6_squalene_synthase_squalestatin	NZ_LIQS01000151	1	1-17890	t1pks	AOK10_RS14390_MFS
# 1	1399	39.0	547.0	466	0.07	9e-09	4819	6045	3420	3	2	2	0	['Bacteria', 'Actinobacteria', 'Streptomycetales', 'Streptomycetaceae', 'Streptomyces']
input_filename = sys.argv[1]
input_file = open(input_filename).readlines()
data = set()

for line in input_file:
    line = line.strip()
    features = line.split("\t")
    # print features
    cluster_length, target, cluster_name, cluster_num, coord, cluster_type, gene = features[:7]
    # print cluster_length, target, cluster_name, cluster_num, coord, cluster_type, gene
    identity, evalue, _, _, dist, kscount, _, _, _, taxa = features[12:]
    # print identity, evalue, dist, kscount, taxa
    data.add((target, cluster_name, cluster_num, cluster_type, gene, identity, evalue, kscount, taxa))

# print '\n'.join([str(i) for i in data])

f = sys.argv[2]
outf = open(f, "w")

for i in data:
    outf.write("\t". join(list(i)))
    outf.write("\n")

outf.close()
