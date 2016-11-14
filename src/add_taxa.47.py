#!/usr/bin/python
import sys

# Append taxa to out file. To execute script:
# ./add_taxa.py /mnt/gnpn/gnpn/projects/orphanpks/TargetMining/Genbank/taxa.txt ../data/out.targets.619.evalue.e-8.filtered.10000.t1pks.transatpks.unique.1901.domaincounts ../data/out.targets.619.evalue.e-8.filtered.10000.t1pks.transatpks.unique.1901.domaincounts.taxa

id_to_features = {}
# /mnt/gnpn/gnpn/projects/orphanpks/TargetMining/Genbank/taxa.txt
# NZ_JYJF01001045	['Bacteria', 'Actinobacteria', 'Pseudonocardiales', 'Pseudonocardiaceae', 'Saccharothrix']
taxa_filename = sys.argv[1]
taxa_file = open(taxa_filename).readlines()
for line in taxa_file:
    line = line.strip()
    features = line.split("\t")
    # print features
    gb_id, taxa = features
    id_to_features[gb_id] = (taxa)
    # print id_to_features[gb_id]


# 100088	sp|Q47146|FADE_ECOLI	NZ_KK073768	7	3530291-3630379	nrps-t1pks-fabH	ctg1_orf06070_-	3579692	3580876	78.0	814.0	394	0.10	2e-03573855	3575132	4560	4	5	0	0
input_filename = sys.argv[2]
input_file = open(input_filename).readlines()

output_file = sys.argv[3]
outf = open(output_file, "w")

for line in input_file:
    line = line.strip()
    features = line.split("\t")
    # cluster_length, target, cluster_name, cluster_num = features[:4] # for e coli targets
    target, cluster_name, cluster_num = features[:3] # for 47 targets
    feats = str('\t'.join(features))
    if cluster_name not in id_to_features.keys():
        outf.write("%s" % feats)
        outf.write("\n")
        print "No taxa", cluster_name
    if cluster_name in id_to_features.keys():
        # taxa = id_to_features[cluster_name]
        taxa = str(''.join(id_to_features[cluster_name]))
        # print taxa
        # sys.exit(0)
        outf.write("%s\t%s" % (feats, taxa))
        outf.write("\n")
        # print feats, "\t", taxa

outf.close()
