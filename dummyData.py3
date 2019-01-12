#USAGE: python3 dummyData.py3 
#--n [Number of sequences in file] --m [Rate of mutation (a decimal number times the length)] 
#--l [length of each sequence] --f [Name of output file]

import argparse, os, random, textwrap

#Get arguments or initialize default arguments (300 sequences of 600 bases with up to 100 base substitutions each, fasta file)
parser=argparse.ArgumentParser()

parser.add_argument('--n', type=int)
parser.add_argument('--l', type=int)
parser.add_argument('--m', type=float)
parser.add_argument('--f', type=str)
parser.add_argument('--y', type=int)
args = parser.parse_args()

if (args.n==None): seqNum = 720
else: seqNum = args.n
if (args.l==None): seqLen = 600 
else: seqLen = args.l
if (args.m==None): mutRate = 0.01*seqLen
else: mutRate = seqLen*args.m
if (args.f==None): fileName = "dummyData.fas" 
else: fileName = args.f	
if (args.y==None): years = 18 
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
year = 2019 - years
idNum = 11000

#Create sequence entry's for a defined number of id's
while year != 2019:
	for x in range(idNum,(idNum + seqNum//years)):
		outtext += (">" + str(idNum) + "_" + str(year) + "\n")	

		#Induce mutation (alter range to alter number of random base point mutations	
		for x in range(random.randint(0,mutRate)):
			index = random.randint(0,(len(sequence)-1))
			sequence[index] = bases[random.randint(0,(len(bases)-1))]

		sequenceText = textwrap.fill(''.join(sequence),60)	
		outtext += (sequenceText + "\n\n")
		idNum=idNum+1
	year = year+1


#write text output and 
outfile.write(outtext)
outfile.close()
