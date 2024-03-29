#Usage "python3 NAlbertaPre.py3 input.fas"

import sys, os
inFile = open(sys.argv[1], "r")
outfile = open("NAdated.fasta", "a+")

for line in inFile: 
	if line.startswith(">"):
		splitLn = line.strip('>').strip('\n').split()

		seqId = splitLn[0].strip('.1')
		seqDate = splitLn[3].split('_')[3]
		#seqYear = int(seqDate.split('-')[2]) + 2000 

		s = (">" + seqId + "_" + str(seqDate) + '\n')
		outfile.write(s)
	else: 
		outfile.write(line)
