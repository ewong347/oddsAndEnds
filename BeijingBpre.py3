#Usage "python3 BeijingBpre.py3 input.fas"

import sys, os
inFile = open(sys.argv[1], "r")
outfile = open("NAdated.fasta", "a+")

for line in inFile: 
	if line.startswith(">"):
		splitLn = line.strip('>').strip('\n').split('_')

		seqId = splitLn[0]
		seqDate = splitLn[6]


		s = (">" + seqId + "_" + str(seqDate) + '\n')
		outfile.write(s)
	else: 
		outfile.write(line)
