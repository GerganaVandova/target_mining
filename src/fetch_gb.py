#!/usr/bin/env python
import subprocess
import os
import os.path
import time
from urllib2 import HTTPError
from Bio import Entrez

# Fetches genbank ids
# ref.: http://wilke.openwetware.org/Parsing_Genbank_files_with_Biopython.html
Entrez.email = 'gergana.vandova@gmail.com'
# accession id works, returns genbank format, looks in the 'nucleotide' database:

gbids = []
gbidfile = "griselimycin.txt"
#gbidfile = "gbids.unique.84233.txt"
f = open(gbidfile, "r")
for line in f:
    gbid = line.strip()
    gbids.append(gbid)


count = 0
for gbid in gbids:
	#if count == 10000:
	#	break
        # store locally:
        gbfolder = "/mnt/gnpn/gnpn/projects/orphanpks/TargetMining/Genbank/gbdir"
        local_filename = os.path.join(gbfolder, gbid + ".gb")
        if os.path.exists(local_filename) == True:
		print local_filename, " already downloaded"
		continue
	count += 1
	try:
		handle = Entrez.efetch(db='nucleotide', id=gbid, rettype='gb')
		print count, gbid, local_filename
        	local_file = open(local_filename, 'w')
        	local_file.write(handle.read())
        	handle.close()
        	local_file.close()
	except:
		time.sleep(20)
		continue

