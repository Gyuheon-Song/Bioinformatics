#!/usr/bin/python3 -u
#  Sungho Lee, 2017 Bioinformatics class
#  usage: ./code.py [sequence file]
#  output: [sequence file name]_tr_freq.txt

import sys

lines = list()
word_count = dict()
freq_base = {"T":["TS", "TT"], "S":["ST", "SS"]}

with open(sys.argv[1]) as LIST:
	for line in LIST:
		lines.append(line.strip())
		for pos in range(len(line.strip())-1):
			twin = line.strip()[pos:pos+2]
			if twin not in word_count:
				word_count[twin] = 0
			word_count[twin] += 1

whole_sequence = "".join(lines)

output = open(sys.argv[1][:-4]+"_tr_freq.txt", "w")

for x in ["S","T"]:
	print("Start probability of {ST} is {FREQ}".format(ST=x, FREQ=whole_sequence.count(x)/len(whole_sequence)), file=output)

print(file=output)

for x in word_count:
	print("Transition probability of {TR} is {FREQ}".format(TR="{} -> {}".format(x[0], x[1]), FREQ=word_count[x]/sum([word_count[k] for k in freq_base[x[0]]])), file=output)

output.close()
