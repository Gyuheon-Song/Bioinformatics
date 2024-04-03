#!/usr/bin/python3 -u
#  Sungho Lee, 2017 Bioinformatics class
#  usage: ./code.py [sequence file]
#  output: [sequence file name]_aa_freq.txt
import sys

aa_list = {"A", "C", "D", "E", "F",
			"G", "H", "I", "K", "L",
			"M", "N", "P", "Q", "R",
			"S", "T", "V", "W", "Y"}

aa_count = dict()
sequence_length = int()

with open(sys.argv[1]) as LIST:
	for line in LIST:
		sequence_length += len(line.strip())
		for x in line.strip():
			if x not in aa_count:
				aa_count[x] = 0
			aa_count[x] += 1

output = open(sys.argv[1][:-4]+"_aa_freq.txt", "w")

for x in aa_list:
	print("Frequency of {AA} is {FREQ}".format(AA=x, FREQ=(aa_count[x]/sequence_length)), file=output)

output.close()
