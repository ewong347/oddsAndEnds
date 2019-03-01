#Usage "python3 dateMerge.py3 input.fas"

import sys, os
inFile = open(sys.argv[1], "r")
idFile = open('ids.txt', "r")
dateFile = open('dates.txt', "r")
outfile = open("mergedDates.fasta",  "a+")

dates={}

for line in idFile:
	dates[line.strip('\n')] = dateFile.readline()


for line in inFile: 
	if line.startswith(">"):
		seqId = line.strip('>')
		seqId = seqId.strip('\n')
		seqDate = dates[seqId]

		s = (">" + seqId + "_" + seqDate)
		outfile.write(s)
	else: 
		outfile.write(line)
	
