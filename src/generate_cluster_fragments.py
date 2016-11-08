#!/usr/bin/python
import sys
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt
from subprocess import call
from StringIO import StringIO

T7_tetO = "taatacgactcactataggTCTATCATTGATAGGgtttccctctagata"
UTRB = "tctagaaataattttgtttaactttaagaaggagatatacc"
UTRA = "ataatttaaaaaacagacctcatatcgaaataaaagaaggagatatacc"
UTRd1 = "tgtagaaataattttgtttaactttaataaggagatatacc"
UTRd2 = "atcttagtatattagttaagtataagaaggagatatacata"
insulator_wt = "agggctccggactgcgctgtatagt"
insulator_m6 = "cccgtacgcgagataaactgctagg"
pheA1 = "gacgaacaaTAAGGCCTCCCAAATCGGGGGGCCTTTTTTATTgaTaacaaaa"
ECK120033737 = "ggaaacacagAAAAAAGCCCGCACCTGACAGTGCGGGCTTTTTTTTTcgaccaaagg"
ECK120029600 = "TTCAGCCAAAAAACTTAAGACCGCCGGTCTTGTCCACTACCTTGCAGTAATGCGGTGGACAGGATCGGCGGTTTTCTTTTCTCTTCTCAA"
NotI = "GCGGCCGC"
SrfI = "GCCCGGGC"
SwaI = "ATTTAAAT"
PmlI = "CACGTG"
PmeI = "GTTTAAAC"
SalI = "GTCGAC"
SexAI = "ACCAGGT"
SgrAI = "CACCGGCG"
NheI = "GCTAGC"
HindIII = "AAGCTT"

utr_list = ["UTRB", "UTRA", "UTRd1", "UTRd2"]
term_list = ["pheA1", "ECK120033737", "ECK120029600"]
ins_list = ["insulator_wt", "insulator_m6"]
re_list = ["NotI", "SrfI", "SwaI", "PmlI", "PmeI", "SalI", "SexAI", "SgrAI", "NheI", "HindIII"]

orf_list_file = sys.argv[1]
primers = []
orf_list = []
name_to_seq = {}

OPERON_SIZE = 6100
MAX_GENES_PER_OPERON = 4
total_cluster_len = 0

for record in SeqIO.parse(orf_list_file, "fasta"):
    seq_name = record.id
    seq = record.seq
    seq_len = len(record.seq)
    total_cluster_len += seq_len
    print seq_name, seq_len, seq[:10]
    name_to_seq[seq_name] = seq
    orf_list.append(seq_name)

# list of all pks or nrps genes
pks_orfs_list = []
other_orfs_list = []
for orf in orf_list:
    if "PKS" in orf:
        pks_orfs_list.append(orf)
    else:
        other_orfs_list.append(orf)

num_genes = len(orf_list)
# print total_cluster_len
total_plasmid_size = total_cluster_len/2
# print total_plasmid_size

if total_plasmid_size > OPERON_SIZE:
    num_operons = 2
if total_plasmid_size > OPERON_SIZE * 2:
    print "increase size of operon"

num_genes_in_operon = float(total_cluster_len) / num_genes

# print num_genes_in_operon

lista = []
start = 0  # counter for start of operon
k = 0  # counter for pks genes
l = 0  # counter for other genes
u = 0  # counter for utrs
operon_len = 0

# for i in xrange(num_operons):


while k + l + 1 <= len(orf_list):
    print "\ngens added (k + l): ", k + l + 1
    print "start: ", start, "\npks gene index k: ", k, "\nother gene index l: ", \
          l, "\nindex utrs u:", u, "\n"
    if start == 0:
        lista.append('prmoter')
        lista.append(utr_list[u])
        if k + 1 > len(pks_orfs_list):
            print "yes, there are no more PKS genes, so start adding the other genes"
            # print "k+1", k+1, "len pks orfs list: ", len(pks_orfs_list)
            lista.append(pks_orfs_list[l])
            print "old operon len: ", operon_len
            operon_len += len(name_to_seq[other_orfs_list[l]])
            print "added operon len of : ", l, other_orfs_list[l], " gene ", len(name_to_seq[other_orfs_list[l]])
            print "new operon len: ", operon_len
            print  "##### will have to increase l!!! #####"
            l += 1
        else:
            lista.append(pks_orfs_list[k])
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
        print "operon len: ", operon_len, "utr index: ", u
    else:
        if operon_len + len(name_to_seq[other_orfs_list[l]]) > OPERON_SIZE:
            print "operon size will exceed MAX SIZE"
            print lista
            print "operon len", operon_len
            operon_len = 0
            lista.append('terminator')
            print "operon len", operon_len
            start = 0
            continue
        lista.append(utr_list[u])
        if u + 1 == len(utr_list):
            print "range of utrs reached end, u = 0"
            u = 0
        else:
            u += 1
        lista.append(other_orfs_list[l])
        l += 1

        print "old operon len: ", operon_len
        operon_len += len(name_to_seq[other_orfs_list[l]])
        print "added operon len of : ", l, other_orfs_list[l], " gene ", len(name_to_seq[other_orfs_list[l]])
        print "new operon len: ", operon_len
        start = 1
        print lista
        # print "u + 1: ", u + 1, "len utrlist:", len(utr_list)
        # print "operon len: ", operon_len, "utr index: ", u
        print "added operon len: ", len(name_to_seq[other_orfs_list[l]])
        print "new operon len: ", operon_len
