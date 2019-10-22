#!/usr/local/pkg/MCBPython2010/epd-6.2-2-rh5-x86_64/bin/python
import sys
import subprocess
from Bio import SeqIO
#import numpy
#import scipy
#import pylab
#import matplotlib.pyplot as plt

fasta = 'Suilu4_AssemblyScaffolds.fasta'#use location for ref genome as argument.

handle = open(fasta, "rU")
Sbre = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
handle.close()

samples=["502931_1178750", "502931_1178751", "502931_1178744", "502931_1178745", "502931_1178746", "502931_1178747", "502931_1178748", "502931_1178749", "502931_1143059", "502931_1143060", "502931_1151569", "502931_1143062", "502931_1143054", "502931_1143056", "502931_1143075", "502931_1143076", "502931_1151566", "502931_1151567", "502931_1151576", "502931_1151572", "502931_1143066", "502931_1143067", "502931_1143068", "502931_1151574", "502931_1143052", "502931_1143050", "502931_1143051", "502931_1143078", "502931_1151580", "502931_1151581", "502931_1151582", "502931_1143070", "502931_1143073", "502931_1151559", "502931_1151560", "502931_1151561", "502931_1151562", "502931_1151564"]

ambig = {'A/G':'R','G/A':'R','C/T':'Y','T/C':'Y','G/C':'S','C/G':'S','A/T':'W','T/A':'W','G/T':'K','T/G':'K','A/C':'M','C/A':'M','A/A':'A','C/C':'C','G/G':'G','T/T':'T','./.':'N'}
         
##########
#read SNPs
tab=open("gvcf.filtered.2019-01-30dik.tab",'rU')
snps ={}#dictionary with scaffolds as keys and (SNP position + individual genotypes) as items
for line in tab:
    if not line.startswith('#'):
        parts=line.split('\t')
        if parts[0] not in snps:
            snps[parts[0]]={}
        snps[parts[0]][parts[1]]=parts[3:]

l=0
fastas={}
for scaf in snps:
    c=0
    for sample in samples:
        fastas[sample]=[]
        name=scaf+"_"+sample+"\n"
        fastas[sample].append('>'+name)
    for i in Sbre[scaf].seq:
        l+=1
        c+=1
        d=str(c)
        print l,scaf,c
        if d in snps[scaf]:
            ind=0
            for indiv in snps[scaf][d]:
                ind+=1
                fastas[samples[ind-1]].append(ambig[indiv.strip('\n')])
        else:
            for sample in samples:
                fastas[sample].append(i)
    for sample in fastas:
        out=open("2019-01-28dik_%s_%s.fasta" %(scaf,sample),'w')
        out.write("".join(fastas[sample]))
