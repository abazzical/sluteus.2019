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

GP1=["502931_1178750", "502931_1178751", "502931_1178744", "502931_1178745", "502931_1178746", "502931_1178747", "502931_1178748", "502931_1178749", "502931_1143059", "502931_1143060", "502931_1151569", "502931_1143062", "502931_1143054", "502931_1143056", "502931_1143075", "502931_1143076", "502931_1151566", "502931_1151567"]
GP2=["502931_1151576", "502931_1151572", "502931_1143066", "502931_1143067", "502931_1143068", "502931_1151574", "502931_1143052", "502931_1143050", "502931_1143051", "502931_1143078", "502931_1151580", "502931_1151581", "502931_1151582", "502931_1143070", "502931_1143073", "502931_1151559", "502931_1151560", "502931_1151561", "502931_1151562", "502931_1151564"]
pops={'GP1':GP1,'GP2':GP2}
pop_ids=['GP1','GP2']
print "scaffold\twindow\tgroup\tS\tthetaW\tPi\tD\tHe\tK"

dataset=dir.replace("concatscafs_fasta","")#

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
			
			stats=egglib.Align.polymorphism(align)
			print scaf,i,pop,stats["S"],stats["thetaW"],stats["Pi"],stats["D"],stats["He"],stats["K"]
