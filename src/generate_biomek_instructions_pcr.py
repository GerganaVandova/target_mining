#!/usr/bin/python

# Howto execute script:
# ./generate_biomek_instructions_pcr.py GV_01_Twist_plate.txt GV_P01_IDT_primers.txt >GV_01_Biomek_instructions.txt

import sys

twist_order_file = sys.argv[1]
primers_file = sys.argv[2]
primers_plate = "GV_P01_IDT_primers"
dest_plate = "GV_P02_PCR"


primer_to_pos = {}
primer_to_fragment = {}
primer_to_fwdrev = {}
primer_to_fullname = {}
ff = open(primers_file).readlines()
for line in ff:
    line = line.strip()
    # print line
    well, primer_name, seq = line.split("\t")
    s = primer_name.rsplit('_', 4)[0]  # 84_Kosan_DEBS3_6
    num = int(s.split("_")[0])  # 84
    s1 = s.split("_")[1:]  # ['Kosan', 'DEBS3', '6']
    s2 = "_".join(s1)  # Kosan_DEBS3_6
    # print s, num, s1, s2
    primer_to_pos[num] = well
    primer_to_fragment[num] = s2
    primer_to_fullname[num] = primer_name

for i in xrange(len(primer_to_fragment.keys())):
    for j in xrange(len(primer_to_fragment.keys())):
        if i == j:
            continue
        if primer_to_fragment[primer_to_fragment.keys()[i]] == primer_to_fragment[primer_to_fragment.keys()[j]]:
            primer_to_fwdrev[primer_to_fragment[primer_to_fragment.keys()[i]]] = (primer_to_fragment.keys()[i], primer_to_fragment.keys()[j])


# print primer_to_fwdrev
# print primer_to_pos


name_to_pos = {}
name_to_len = {}
f = open(twist_order_file).readlines()
print "Input plate\tInput well\tFragment name\tDestination plate\tDestination well\tSize"
for line in f:
    line = line.strip()
    # print line
    plate, well, fragment, length = line.split("\t")
    # print plate, well, fragment
    name_to_pos[fragment] = (plate, well)
    name_to_len[fragment] = length
    # input plate, input well, fragment name, dest plate, dest well, fragment length
    print ("%s\t%s\t%s\t%s\t%s\t%s\t") % (name_to_pos[fragment][0], name_to_pos[fragment][1],
                                          fragment, dest_plate, name_to_pos[fragment][1],
                                          name_to_len[fragment])
    p1, p2 = primer_to_fwdrev[fragment]
    p1 = int(p1)
    p2 = int(p2)
    fwd = min(p1, p2)
    rev = max(p1, p2)
    print ("%s\t%s\t%s\t%s\t%s\t%s\t") % (primers_plate, primer_to_pos[fwd], primer_to_fullname[fwd],
                                          dest_plate, name_to_pos[fragment][1], name_to_len[fragment])
    print ("%s\t%s\t%s\t%s\t%s\t%s\t") % (primers_plate, primer_to_pos[rev], primer_to_fullname[rev],
                                          dest_plate, name_to_pos[fragment][1], name_to_len[fragment])



# print name_to_pos
