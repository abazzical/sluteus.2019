#!/usr/local/pkg/MCBPython2010/epd-6.2-2-rh5-x86_64/bin/python
import sys
import egglib
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import itertools

path="/Users/abazzical/Documents/suillusFastSTRUCTURE/20190130_onlyDikaryons/calculateDxy/"

dir=sys.argv[1]
KB=int(sys.argv[2])

GP1=["1178750", "1178751", "1178744", "1178745", "1178746", "1178747", "1178748", "1178749", "1143059", "1143060", "1151569", "1143062", "1143054", "1143056", "1143075", "1143076", "1151566", "1151567"]
GP2=["1151576", "1151572", "1143066", "1143067", "1143068", "1151574", "1143052", "1143050", "1143051", "1143078", "1151580", "1151581", "1151582", "1143070", "1143073", "1151559", "1151560", "1151561", "1151562", "1151564"]
pops={'GP1':GP1,'GP2':GP2}
pop_ids=['GP1','GP2']
print "scaffold\twindow\tgroup\tS\tthetaW\tPi\tD\tHe\tK"

dataset=dir.replace("2019-01-28dikconcatscafs_fasta_","")#

for file in os.listdir(path+dir):
#  	print file
	data=egglib.Align(path+dir+"/"+file,groups=False)
	popdata={}
	for pop in pop_ids:
		popdata[pop]=[]
	scaf=file.replace(dataset+"_scaffold_","")
	scaf=scaf.replace(".fasta","")
	print scaf
	for seq in data:
		seq.name=(seq.name).replace("scaffold_"+scaf+"_","")
# 		print seq.name
		for pop in pop_ids:
			if seq.name in pops[pop]:
				popdata[pop].append(seq)

# 	for name,seq in popdata['GP1']:
# 		print name,seq
	for pop in popdata:
		align=egglib.Align.create(popdata[pop])
# 		print str(align[0])
# 		print align
		L=len(align[0][1])
		print L
		K=L/(KB*1000)
		for i in range(0,K):
# 			win=align.extract(range(i*KB*1000,(i+1)*KB*1000))
			stats=egglib.Align.polymorphism(align)
			print scaf,i,pop,stats["S"],stats["thetaW"],stats["Pi"],stats["D"],stats["He"],stats["K"]
