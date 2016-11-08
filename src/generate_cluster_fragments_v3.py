#!/usr/bin/python
import sys
import os
from Bio import SeqIO
from Bio import Restriction
from Bio.Restriction import *
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna, generic_protein
from Bio.SeqUtils import MeltingTemp as mt
from subprocess import call
from StringIO import StringIO

# v2working for cluster #5.

TRIM_SIZE = 1800  # maximum allowable length of a gene for ordering from Twist
OVERHANG_LEN = 60
OPERON_SIZE = 6000


# From Colin's script, modified by GV
# Use this function to chop up sequences and add name to each chopped fragment
def split_seq(this_seq, seq_name, trim_size, oh_length):
    seq = str(this_seq)
    split_seqs = {}
    seq_len = len(seq)
    # if the whole sequence length is smaller than the trim size, return the same sequence
    if seq_len < trim_size:
        split_seqs[seq_name] = seq
    else:
        split_num = int(seq_len/(trim_size - 2 * oh_length))
        if seq_len % (split_num * (trim_size - 2 * oh_length)) > 0:
            split_num += 1
        lfrag = int(seq_len / split_num)
        split_seqs[seq_name] = seq[0: (lfrag + oh_length / 2)]
        if split_num > 2:
            for i in range(1, split_num - 1):
                this_start = i * lfrag - (oh_length / 2)
                this_end = (i + 1) * lfrag + (oh_length / 2)
                seq_name_new = seq_name + "-" + str(i)
                split_seqs[seq_name_new] = seq[this_start: this_end]
                seq_name_new = seq_name
        split_seqs[seq_name] = seq[(seq_len - lfrag - (oh_length / 2)): seq_len]
    return split_seqs


# Convert fasta file into genbank
def fasta_to_gb(fasta_file_name):

    gbfile_list = []
    # in thefasta filename are the two constructs: pET and pRSF
    fasta_file = open(fasta_file_name, "rU")
    sequences = list(SeqIO.parse(fasta_file, "fasta"))
    # Asign generic_dna or generic_protein
    for seq in sequences:
        gbfile_name = fasta_file_name + "_" + seq.id + ".gb"
        gbfile = open(gbfile_name, "w")
        seq.seq.alphabet = generic_dna
        SeqIO.write(seq, gbfile, "genbank")
        gbfile.close()
        fasta_file.close()
        gbfile_list.append(gbfile_name)
    return gbfile_list


def gb_to_fullgb(gbfile_name, cluster_name, backbone_name):

    gbf = open(gbfile_name, 'r')
    gbfile = gbf.read()
    start_file, end_file = gbfile.split("FEATURES             Location/Qualifiers")

    # Find coordinates of all elelments in pET and pRSF constructs
    gbfile_full = cluster_name + "-" + backbone_name + ".gb"

    ff = open(os.path.join(cluster_dir, gbfile_full), "w")
    ff.write("%s\nFEATURES             Location/Qualifiers\n" % start_file)

    if "pETDuet" in gbfile_name:
        for element in pETDuet_construct_list:
            element_start = pETDuet_construct_seq_with_re.find(elements_to_seq[element]) + 1
            element_end = element_start + len(elements_to_seq[element]) - 1
            print element, element_start, element_end, len(elements_to_seq[element])
            f = open(os.path.join(cluster_dir, gbfile_full), "a")
            f.write("    gene            %s..%s\n" % (element_start, element_end))
            f.write("                    /locus_tag='%s'\n" % element)
        f.write("%s" % end_file)

    if "pRSFDuet" in gbfile_name:
        for element in pRSFDuet_construct_list:
            element_start = pRSFDuet_construct_seq_with_re.find(elements_to_seq[element]) + 1
            element_end = element_start + len(elements_to_seq[element]) - 1
            print element, element_start, element_end, len(elements_to_seq[element])
            f = open(os.path.join(cluster_dir, gbfile_full), "a")
            f.write("    gene            %s..%s\n" % (element_start, element_end))
            f.write("                    /locus_tag='%s'\n" % element)
        f.write("%s" % end_file)


elements_to_seq = {"T7_tetO": "taatacgactcactataggTCTATCATTGATAGGgtttccctctagata",
                   "UTRB": "tctagaaataattttgtttaactttaagaaggagatatacc",
                   "UTRA": "ataatttaaaaaacagacctcatatcgaaataaaagaaggagatatacc",
                   "UTRd1": "tgtagaaataattttgtttaactttaataaggagatatacc",
                   "UTRd2": "atcttagtatattagttaagtataagaaggagatatacata",
                   "insulator_wt": "agggctccggactgcgctgtatagt",
                   "insulator_m6": "cccgtacgcgagataaactgctagg",
                   "pheA1": "gacgaacaaTAAGGCCTCCCAAATCGGGGGGCCTTTTTTATTgaTaacaaaa",
                   "ECK120033737": "ggaaacacagAAAAAAGCCCGCACCTGACAGTGCGGGCTTTTTTTTTcgaccaaagg",
                   "ECK120029600": "TTCAGCCAAAAAACTTAAGACCGCCGGTCTTGTCCACTACCTTGCAGTAATGCGGTGGACAGGATCGGCGGTTTTCTTTTCTCTTCTCAA",
                   "NotI": "GCGGCCGC",
                   "SrfI": "GCCCGGGC",
                   "SwaI": "ATTTAAAT",
                   "PmlI": "CACGTG",
                   "PmeI": "GTTTAAAC",
                   "SalI": "GTCGAC",
                   "SexAI": "ACCAGGT",
                   "SgrAI": "CACCGGCG",
                   "NheI": "GCTAGC",
                   "HindIII": "AAGCTT",
                   "AscI": "GGCGCGCC",
                   "SpeI": "ACTAGT",
                   "ura3_extended_upoh": "tgtaagcggatgccgggagcagacaagcccgtcagggcgcgtcagcgggtgttggcgggt",
                   "ura3_extended_downoh": "TATTACCCTATGCGGTGTGAAATACCGCACAGATGCGTAAGGAGAAAATACCGCATCAGG"}
prom_list = ["T7_tetO"]
utr_list = ["UTRB", "UTRA", "UTRd1", "UTRd2"]
term_list = ["ECK120029600", "ECK120033737", "pheA1"]
ins_list = ["insulator_wt", "insulator_m6"]
re_list = ["NotI", "SrfI", "SwaI", "PmlI", "PmeI", "SalI", "SexAI", "SgrAI", "NheI", "HindIII", "AscI", "SpeI"]

# Input fasta file of backbones and ura3
backbones_list_file = sys.argv[1]
# backbone_name_to_seq = {}
backbone_list = []
for record in SeqIO.parse(backbones_list_file, "fasta"):
    seq_name = record.id
    seq = record.seq
    # add backbone sequences to dictionary with all genetic elements
    elements_to_seq[seq_name] = seq
    backbone_list.append(seq_name)

# Input fasta file of PKS orfs:
orf_list_file = sys.argv[2]
cluster_name = sys.argv[3]
total_cluster_len = 0
orf_list = []
name_to_seq = {}  # need this dict to calculate size of operon

for record in SeqIO.parse(orf_list_file, "fasta"):
    seq_name = record.id
    seq = record.seq
    seq_len = len(record.seq)
    total_cluster_len += seq_len
    print seq_name, seq_len, seq[:10]
    elements_to_seq[seq_name] = seq
    name_to_seq[seq_name] = seq  # need this dict to calculate size of operon
    orf_list.append(seq_name)


# list of all pks or nrps genes
pks_orfs_list = []
other_orfs_list = []
for orf in orf_list:

    # Manually specify with which genes each operon should start
    # (PKS genes should be first) and PKS list stored in a list
    if "PKS" in orf:
        pks_orfs_list.append(orf)
    else:
        other_orfs_list.append(orf)

num_genes = len(orf_list)
print "\nTotal clulster len", total_cluster_len
total_plasmid_size = total_cluster_len/2
print "Total plasmid size", total_plasmid_size

if total_plasmid_size > OPERON_SIZE:
    num_operons = 2
print "Number of operons is 2"
if total_plasmid_size > OPERON_SIZE * 2:
    print "increase size of operon"

# Not using it as of now, comment
# num_genes_in_operon = float(total_cluster_len) / num_genes

# print num_genes_in_operon

lista = []
start = 0  # counter for start of operon
k = 0  # counter for pks genes
l = 0  # counter for other genes
p = 0  # counter for promoters
u = 0  # counter for utrs
t = 0  # counter for terminators
i = 0  # counter for insulators
operon_len = 0
count = 0  # counter for initializing operon for the first time, every other time, an insulator sequence has to be added
# for i in xrange(num_operons):
o = 0  # counter for operon number


while k + l <= len(orf_list):
    print "\ngens added (k + l): ", k + l + 1
    print "start: ", start, "\npks gene index k: ", k, "\nother gene index l: ", \
          l, "\nindex utrs u:", u, "\n"

    # When initializing an Operon, start=0, add promoter first
    if start == 0:
        o += 1
        if o == 3:  # will work if I have 4 operons, 2 for each plasmid
            lista.append('pRSFDuet_backbone')
        if count != 0 and o != 3:  # add ura3 only when in the middle of 2 operon,s not when adding pRSF
            lista.append('ura3_extended')
            lista.append(ins_list[i])
            if i + 1 == len(ins_list):
                i = 0
            else:
                i += 1
        elif count == 0:
            lista.append('pETDuet_backbone')
        lista.append(prom_list[p])
        if p + 1 == len(prom_list):
            p = 0
        else:
            p += 1
        lista.append(utr_list[u])
        if k + 1 > len(pks_orfs_list):
            print "yes, there are no more PKS genes, so start adding the other genes"
            # print "k+1", k+1, "len pks orfs list: ", len(pks_orfs_list)
            lista.append(other_orfs_list[l])
            print "old operon len: ", operon_len
            operon_len += len(name_to_seq[other_orfs_list[l]])
            print "added operon len of : ", l, other_orfs_list[l], " gene ", len(name_to_seq[other_orfs_list[l]])
            print "new operon len: ", operon_len
            print "##### will have to increase l!!! #####"
            l += 1

        else:
            lista.append(pks_orfs_list[k])
            print "added gene", k, pks_orfs_list[k]

            print "old operon len: ", operon_len
            operon_len += len(name_to_seq[pks_orfs_list[k]])
            print "added operon len of pks : ", k, pks_orfs_list[k], " gene ", len(name_to_seq[pks_orfs_list[k]])
            print "new operon len: ", operon_len

        # if operon_len + len(name_to_seq[other_orfs_list[l]]) > OPERON_SIZE:
        #     print "operon size will exceed MAX SIZE"
        #     print "operon len", operon_len
        #     operon_len = 0
        #     lista.append('terminator')
        #     print "Reinitialized operon len: ", operon_len
        #     start = 0
        #     continue
        start = 1
        k += 1
        if u + 1 == len(utr_list):
            u = 0
        else:
            u += 1
        print lista

    # When not initializing an operon, start with UTR
    else:
        if operon_len + len(name_to_seq[other_orfs_list[l]]) > OPERON_SIZE:
            print "\n\n###############\n"
            print operon_len, len(name_to_seq[other_orfs_list[l]]), operon_len + len(name_to_seq[other_orfs_list[l]])
            print "operon size will exceed MAX SIZE if ", other_orfs_list[l], "is added with len ", len(name_to_seq[other_orfs_list[l]])
            print "Final operon len", operon_len, "operon len exceeding max: ", operon_len + len(name_to_seq[other_orfs_list[l]])
            operon_len = 0
            lista.append(term_list[t])
            if t + 1 == len(term_list):
                t = 0
            else:
                t += 1
            # print "operon len", operon_len
            print lista
            print "\n###############\n\n"
            start = 0
            continue

        lista.append(utr_list[u])
        if u + 1 == len(utr_list):
            print "range of utrs reached end, u = 0"
            u = 0
        else:
            u += 1
        lista.append(other_orfs_list[l])
        print "gene added: ", l, other_orfs_list[l]
        print "old operon len: ", operon_len
        operon_len += len(name_to_seq[other_orfs_list[l]])
        print "added operon len of : ", l, other_orfs_list[l], " gene ", len(name_to_seq[other_orfs_list[l]])
        print "new operon len: ", operon_len
        # print lista

        l += 1
        start = 1
        # print "u + 1: ", u + 1, "len utrlist:", len(utr_list)
        # print "operon len: ", operon_len, "utr index: ", u
    count += 1
lista.append(term_list[t])
print lista

# Concatenate all sequence elements into one big giant sequence with both backbones
construct_seq = ''
for element in lista:
    seq = elements_to_seq[element]
    # print seq
    construct_seq += seq
construct_seq.lower()
# print construct_seq


print "\n".join(lista)

# Split sequence into 2 plasmids
start_pRSF_construct = construct_seq.find(elements_to_seq['pRSFDuet_backbone'])
start_pRSF_construct_element = lista.index('pRSFDuet_backbone')
print "start_pRSF_construct_element", start_pRSF_construct_element

pETDuet_construct_seq = construct_seq[:start_pRSF_construct]
pETDuet_construct_list = lista[:start_pRSF_construct_element]
pRSFDuet_construct_seq = construct_seq[start_pRSF_construct:]
pRSFDuet_construct_list = lista[start_pRSF_construct_element:]

print "pET construct list: ", pETDuet_construct_list
print "pET construct ends with: ", pETDuet_construct_seq[-20:]
print "pRSF construct list: ", pRSFDuet_construct_list
print "pRSF construct starts with: ", pRSFDuet_construct_seq[:20], pRSFDuet_construct_seq[-20:]

# Search for restriction enzymes in sequences and add them
allowed_re = []
for re in re_list:
    # print re, NotI.site  # not working
    re_seq = elements_to_seq[re].lower()
    construct_seq = construct_seq.lower()
    print re, construct_seq.find(re_seq), "in the whole sequence"
    if construct_seq.find(re_seq) == -1:
        allowed_re.append(re)

if len(allowed_re) < 4:
    print "Add more enzymes to search for in the sequence"
allowed_re = sorted(allowed_re)
print sorted(allowed_re)

# Insert RE in pETDuet plasmid
pETDuet_construct_list.insert(1, allowed_re[0])
pETDuet_construct_list.insert(3, allowed_re[1])
ura3_index = pETDuet_construct_list.index('ura3_extended')
pETDuet_construct_list.insert(ura3_index - 1, allowed_re[1])
prom_indexes = [i for i, x in enumerate(pETDuet_construct_list) if x == 'T7_tetO']
last_prom_index = prom_indexes[-1]
print "prom index", last_prom_index
pETDuet_construct_list.insert(last_prom_index, allowed_re[2])
pETDuet_construct_list.insert(last_prom_index + 2, allowed_re[3])
last_term_index = len(pETDuet_construct_list)
pETDuet_construct_list.insert(last_term_index - 1, allowed_re[3])
print "\n###### pETDUET construct with RE:\n", pETDuet_construct_list

# Insert RE in pRSF construct list
pRSFDuet_construct_list.insert(1, allowed_re[0])
pRSFDuet_construct_list.insert(3, allowed_re[1])
ura3_index = pRSFDuet_construct_list.index('ura3_extended')
pRSFDuet_construct_list.insert(ura3_index - 1, allowed_re[1])
prom_indexes = [i for i, x in enumerate(pRSFDuet_construct_list) if x == 'T7_tetO']
last_prom_index = prom_indexes[-1]
print "prom index", last_prom_index
pRSFDuet_construct_list.insert(last_prom_index, allowed_re[2])
pRSFDuet_construct_list.insert(last_prom_index + 2, allowed_re[3])
last_term_index = len(pRSFDuet_construct_list)
pRSFDuet_construct_list.insert(last_term_index - 1, allowed_re[3])
print "\n###### pRSFDUET construct with RE:\n", pRSFDuet_construct_list


# Concatenate all seqeuce elements into two pET and pRSF constructs
pETDuet_construct_seq_with_re = ""
for element in pETDuet_construct_list:
    seq = elements_to_seq[element]
    pETDuet_construct_seq_with_re += seq

pRSFDuet_construct_seq_with_re = ""
for element in pRSFDuet_construct_list:
    seq = elements_to_seq[element]
    pRSFDuet_construct_seq_with_re += seq

print "pET construct len: ", len(pETDuet_construct_seq_with_re)
print "pRSF construct len: ", len(pRSFDuet_construct_seq_with_re)

# ########## Step # Write all sequences in files

# Write in files
cluster_dir = cluster_name
if not os.path.exists(cluster_dir):
    os.makedirs(cluster_dir)

# Write seq of two plasmids in a file
with open(os.path.join(cluster_dir, cluster_name + "_constructs.out.fasta"),"w") as f:
    # with open("constructs.out.fasta", "w") as f:
    f.write(">pETDuet_c\n{}\n>pRSFDuet_c\n{}\n".format(pETDuet_construct_seq_with_re, pRSFDuet_construct_seq_with_re))
f.close()

# Split sequence into Full fragments
start_index_ura = pETDuet_construct_seq_with_re.find(elements_to_seq['ura3_extended'])
fragment_1 = pETDuet_construct_seq_with_re[:start_index_ura + OVERHANG_LEN]
end_index_ura = start_index_ura + len(elements_to_seq['ura3_extended']) - OVERHANG_LEN
fragment_2 = pETDuet_construct_seq_with_re[end_index_ura:]
start_index_ura = pRSFDuet_construct_seq_with_re.find(elements_to_seq['ura3_extended'])
fragment_3 = pRSFDuet_construct_seq_with_re[:start_index_ura + OVERHANG_LEN]
end_index_ura = start_index_ura + len(elements_to_seq['ura3_extended']) - OVERHANG_LEN
fragment_4 = pRSFDuet_construct_seq_with_re[end_index_ura:]

# Put full fragments in a dictionary
fragment_name_1 = cluster_name + "_" + "Fr1"
fragment_name_2 = cluster_name + "_" + "Fr2"
fragment_name_3 = cluster_name + "_" + "Fr3"
fragment_name_4 = cluster_name + "_" + "Fr4"
full_fragments = {fragment_name_1: fragment_1, fragment_name_2: fragment_2,
                  fragment_name_3: fragment_3, fragment_name_4: fragment_4}

print "fragment 1", fragment_1[:20], fragment_1[-60:]
print "fragment 2", fragment_2[:120], fragment_2[-20:]
print "fragment 3", fragment_3[:20], fragment_3[-60:]
print "fragment 4", fragment_4[:120], fragment_4[-20:]

# Write full fragments in a file
with open(os.path.join(cluster_dir, cluster_name + "_full_fragments.out.fasta"),"w") as f:
    f.write(">{}\n{}".format(fragment_name_1, fragment_1))
    f.write("\n>{}\n{}".format(fragment_name_2, fragment_2))
    f.write("\n>{}\n{}".format(fragment_name_3, fragment_3))
    f.write("\n>{}\n{}".format(fragment_name_4, fragment_4))
    f.close()

w = open(os.path.join(cluster_dir, cluster_name + "_twist_fragments.out.fasta"), "w")
# Split Full fragments into twist fragments
for fragment_name in full_fragments.keys():
    twist_fragments = split_seq(full_fragments[fragment_name], fragment_name, TRIM_SIZE, OVERHANG_LEN)
    # print twist_fragments
    with open(os.path.join(cluster_dir, cluster_name + "_twist_fragments.out.fasta"), "a") as f:
        for fragment in twist_fragments.keys():
            f.write(">{}\n{}\n".format(fragment, twist_fragments[fragment]))





# Make gb file  NZ_CM003601_constructs.out.fasta_pRSFDuet_c.gb
gbfile_list = fasta_to_gb((os.path.join(cluster_dir, cluster_name + "_constructs.out.fasta")))

# ['NZ_CM003601/NZ_CM003601_constructs.out.fasta_pETDuet_c.gb', 'NZ_CM003601/NZ_CM003601_constructs.out.fasta_pRSFDuet_c.gb']
print gbfile_list

for gbfile_name in gbfile_list:
    if "pETDuet" in gbfile_name:
        backbone_name = "pETDuet"
    if "pRSFDuet" in gbfile_name:
        backbone_name = "pRSFDuet"
    gb_to_fullgb(gbfile_name, cluster_name, backbone_name)

#
#
# gbfile_name = "NZ_CM003601/NZ_CM003601_constructs.out.fasta_pETDuet_c.gb"
#
# gbf = open(gbfile_name, 'r')
# gbfile = gbf.read()
# start_file, end_file = gbfile.split("FEATURES             Location/Qualifiers")
#
# # Find coordinates of all elelments in pET and pRSF constructs
# gbfile_full = "NZ_CM003601-pETDuet.gb"
# ff = open(os.path.join(cluster_dir, gbfile_full), "w")
# ff.write("%s\nFEATURES             Location/Qualifiers\n" % start_file)
# for element in pETDuet_construct_list:
#     # print element, elements_to_seq[element]
#     element_start = pETDuet_construct_seq_with_re.find(elements_to_seq[element]) + 1
#     element_end = element_start + len(elements_to_seq[element]) - 1
#     print element, element_start, element_end, len(elements_to_seq[element])
#     f = open(os.path.join(cluster_dir, gbfile_full), "a")
#     f.write("\n    gene            %s..%s\n" % (element_start, element_end))
#     f.write("                    /locus_tag='%s'\n" % element)
#
# # f.write("%s" % end_file)
#
#

    #
    #
    #
    #
    #
    #
    #



































# # List of restriction enzymes to add to sequence
# allowed_re_pET = []
# for re in re_list:
#     # print re, NotI.site  # not working
#     re_seq = elements_to_seq[re].lower()
#     print re, pETDuet_construct_seq.find(re_seq), " in pETDuet"
#     if pETDuet_construct_seq.find(re_seq) == -1:
#         allowed_re_pET.append(re)
#
# allowed_re_pET = sorted(allowed_re_pET)
# print sorted(allowed_re_pET)
# if len(allowed_re_pET) < 4:
#     print "Add more enzymes to search for in the sequence"
# pETDuet_construct_list.insert(1, allowed_re_pET[0])
# pETDuet_construct_list.insert(3, allowed_re_pET[1])
# ura3_index = pETDuet_construct_list.index('ura3_extended')
# pETDuet_construct_list.insert(ura3_index - 1, allowed_re_pET[1])
# prom_indexes = [i for i, x in enumerate(pETDuet_construct_list) if x == 'T7_tetO']
# last_prom_index = prom_indexes[-1]
# print "prom index", last_prom_index
# pETDuet_construct_list.insert(last_prom_index, allowed_re_pET[2])
# pETDuet_construct_list.insert(last_prom_index + 2, allowed_re_pET[3])
# last_term_index = len(pETDuet_construct_list)
# pETDuet_construct_list.insert(last_term_index - 1, allowed_re_pET[3])
# print "\n###### pETDUET construct with RE:\n", pETDuet_construct_list
#
# allowed_re_pRSF = []
# for re in re_list:
#     # print re, NotI.site  # not working
#     re_seq = elements_to_seq[re].lower()
#     print re, pRSFDuet_construct_seq.find(re_seq), " in pRSFDuet"
#     if pRSFDuet_construct_seq.find(re_seq) == -1:
#         allowed_re_pRSF.append(re)
# allowed_re_pRSF = sorted(allowed_re_pRSF)
# print allowed_re_pRSF
# if len(allowed_re_pRSF) < 4:
#     print "Add more enzymes to search for in the sequence"
#
# # Add enzymes to construct list
# pRSFDuet_construct_list.insert(1, allowed_re_pRSF[0])
# pRSDFuet_construct_list.insert(3, allowed_re_pRSF[1])
# ura3_index = pDuet_construct_list.index('ura3_extended')
# pRSFDuet_construct_list.insert(ura3_index - 1, allowed_re_pRSF[1])
# prom_indexes = [i for i, x in enumerate(pRSFDuet_construct_list) if x == 'T7_tetO']
# last_prom_index = prom_indexes[-1]
# print "prom index", last_prom_index
# pRSFDuet_construct_list.insert(last_prom_index, allowed_re_pRSF[2])
# pRSFDuet_construct_list.insert(last_prom_index + 2, allowed_re_pRSF[3])
# last_term_index = len(pRSFDuet_construct_list)
# pRSFDuet_construct_list.insert(last_term_index - 1, allowed_re_pRSF[3])
# print "\n###### pRSFDUET construct with RE:\n", pRSFDuet_construct_list
