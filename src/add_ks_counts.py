#!/usr/bin/python
import sys

id_to_features = {}
# all.cluster.counts.txt
# AASO01003060.1	2	1	1	1
ks_count_filename = sys.argv[1]
ks_count_file = open(ks_count_filename).readlines()
for line in ks_count_file:
    line = line.strip()
    features = line.split("\t")
    # print features
    full_id, ks, kr, dh, er = features
    id_to_features[full_id] = (ks, kr, dh, er)
    # print full_id, ks, kr, dh, er
    # print id_to_features[full_id]

# LpxC,NZ_LGRC01000036,1,22846-70432,t1pks,ACU09_RS39030_UDP-3-O-[3-hydroxymyristoyl],50429,51352,95,305,307,0.31,8.00E-48424332,5.00E-39,2.00E-1680E-169
input_filename = sys.argv[2]
input_file = open(input_filename).readlines()
for line in input_file:
    line = line.strip()
    features = line.split("\t")
    target, cluster_name, cluster_num = features[:3]
    full_id = cluster_name + "." + cluster_num
    if full_id in id_to_features.keys():
        # print full_id, id_to_features[full_id], features
        print full_id, "\t", '\t'.join(id_to_features[full_id]), '\t'.join(features)
    # print target, cluster_name, cluster_num, full_id
