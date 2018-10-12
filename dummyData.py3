#USAGE: python3 dummyData.py3 
#--n [Number of sequences in file] --m [Rate of mutation (a decimal number times the length)] 
#--l [length of each sequence] --f [Name of output file]

import argparse, os, random, textwrap

#Get arguments or initialize default arguments (300 sequences of 600 bases with up to 100 base substitutions each, fasta file)
parser=argparse.ArgumentParser()

parser.add_argument('--n', type=int)
parser.add_argument('--l', type=int)
parser.add_argument('--m', type=int)
parser.add_argument('--f', type=str)
args = parser.parse_args()

if (args.n==None): seqNum = 300
else: seqNum = args.n
if (args.l==None): seqLen = 600 
else: seqLen = args.l
if (args.m==None): mutRate = 100 
else: mutRate = seqLen*args.m
if (args.f==None): fileName = "dummyData.fas" 
else: fileName = args.f	

outfile = open(fileName, 'a+')
bases = ["G", "T", "A", "C"] #Altering base list can alter GC content	

#Create a base sequence from which all other sequences will be derived
baseSequence = list()
for i in range(0,(seqLen-1)):
	baseSequence.append(bases[random.randint(0,(len(bases)-1))])

#Initialize Text Output
outtext = ""

sequence = baseSequence

#Create sequence entry's for a defined number of id's
for x in range(11000,(11000 + seqNum)):
	outtext += (">" + str(x) + "_" + str(random.randint(2001,2018)) + "\n")	

	#Induce mutation (alter range to alter number of random base point mutations	
	for x in range(random.randint(0,mutRate)):
		index = random.randint(0,(len(sequence)-1))
		sequence[index] = bases[random.randint(0,(len(bases)-1))]
	sequenceText = textwrap.fill(''.join(sequence),60)	
	outtext += (sequenceText + "\n")

#write text output and 
outfile.write(outtext)
outfile.close()
