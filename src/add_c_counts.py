#!/usr/bin/python
import sys

#To execute script: ./add_ks_counts.py all.cluster.counts.txt out.targets.7.filtered.10000.t1pks.transatpks out.targets.7.filtered.10000.t1pks.transatpks.kscounts

id_to_features = {}
# all.cluster.counts.txt
# AASO01003060.1	2	1	1	1
ks_count_filename = sys.argv[1]  # all.cluster.counts.txt
ks_count_file = open(ks_count_filename).readlines()
for line in ks_count_file:
    line = line.strip()
    features = line.split("\t")
    print features
    full_id, ks, kr, dh, er, c = features
    id_to_features[full_id] = (c)
    # print full_id, ks, kr, dh, er
    # print id_to_features[full_id]


#target	NZ_LBDA01000059	2	fabH-t1pks-t2fas	VT52_RS15400_acetylornithine	0.23	9e-25	2	taxa
# 100088	sp|Q47146|FADE_ECOLI	NZ_KK073768	7	3530291-3630379	nrps-t1pks-fabH	ctg1_orf06070_-	3579692	3580876	78.0	814.0	394	0.10	2e-03573855	3575132	4560
input_filename = sys.argv[2]
input_file = open(input_filename).readlines()

output_file = sys.argv[3]
outf = open(output_file, "w")

for line in input_file:
    line = line.strip()
    features = line.split("\t")
    target, cluster_name, cluster_num = features[:3]
    full_id = cluster_name + "." + cluster_num
    if full_id in id_to_features.keys():
        # print full_id, id_to_features[full_id], features
        cdomain_count = id_to_features[full_id]
        feats = str('\t'.join(features))
        outf.write("%s\t%s" % (feats, cdomain_count))
        outf.write("\n")
        print full_id, "\t", feats, cdomain_count

        # print full_id, "\t", '\t'.join(id_to_features[full_id]), '\t'.join(features)
        # print target, cluster_name, cluster_num, full_id
outf.close()
