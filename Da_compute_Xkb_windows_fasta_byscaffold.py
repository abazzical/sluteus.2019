#!/usr/local/pkg/MCBPython2010/epd-6.2-2-rh5-x86_64/bin/python
import sys
import egglib
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import itertools

path="/Users/abazzical/Documents/suillusFastSTRUCTURE/20190130_onlyDikaryons/calculateDxy"

dataset=sys.argv[1]
win=int(sys.argv[2])#window size in bp (1000, 5000, etc)

scaffolds=dataset+"concatscafs_fasta"#sys.argv[1]
fastas=dataset+"_scaffold_"

GP1=["502931_1178750", "502931_1178751", "502931_1178744", "502931_1178745", "502931_1178746", "502931_1178747", "502931_1178748", "502931_1178749", "502931_1143059", "502931_1143060", "502931_1151569", "502931_1143062", "502931_1143054", "502931_1143056", "502931_1143075", "502931_1143076", "502931_1151566", "502931_1151567"]
GP2=["502931_1151576", "502931_1151572", "502931_1143066", "502931_1143067", "502931_1143068", "502931_1151574", "502931_1143052", "502931_1143050", "502931_1143051", "502931_1143078", "502931_1151580", "502931_1151581", "502931_1151582", "502931_1143070", "502931_1143073", "502931_1151559", "502931_1151560", "502931_1151561", "502931_1151562", "502931_1151564"]
pops={'GP1':GP1,'GP2':GP2}
pop_ids=['GP1','GP2']
# GP1=NOmetal
# GP2=metal

print "SCAFFOLD\tPAIR\tSTART\tEND\tDxy\tDa\tPi1\tPi2"

scafs=["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", "34", "35", "36", "37", "38", "39", "40", "41", "42", "43", "44", "45", "46", "47", "48", "49", "50", "51", "52", "53", "54", "55", "56", "57", "58", "59", "60", "61", "62", "63", "64", "65", "66", "67"]
#listed by doing ls HiQ_Sb_woNsSin_woSb32Sb7_indivscafs_fasta_concat/ | while read f; do g=`echo $f|sed 's/HiQ_Sb_woNsSin_woSb32Sb7_scaffold_//'| sed 's/.fasta//'`; echo -n $g"\",\""; done

#for fasta in os.listdir(path+"/"+scaffolds):
for scaf in scafs: 
    fasta=fastas+scaf+".fasta"
    data=egglib.Align(path+"/"+scaffolds+"/"+fasta,groups=False)
    popdata={}
    for pop in pop_ids:
        popdata[pop]=[]
    for seq in data:
        seq.name=seq.name.replace('scaffold_'+scaf+'_','')
        for pop in pop_ids:
            if seq.name in pops[pop]:
                popdata[pop].append(seq)
    for pair in itertools.combinations(pop_ids, 2):
        if sorted(pair) not in [sorted(['WAca','CAWAca']),sorted(['WAnm','NMWAnm']),sorted(['AKEU','AK'])]:#*****CHECK THIS!!!!
            align1=egglib.Align.create(popdata[pair[0]])
            align2=egglib.Align.create(popdata[pair[1]])
            L=len(align1[0][1])
            if L<win:
                print scaf,pair,0,0,0,0,0,0
            else:
                K=L/win
#                print popdata[pair[0]]
                for k in range(0,K):
                    seq1=align1.extract(range(k*win,(k+1)*win))
                    seq2=align2.extract(range(k*win,(k+1)*win))
                    polym1=egglib.Align.polymorphism(seq1)
                    polym2=egglib.Align.polymorphism(seq2)
                    DIFFsum=0
                    DIFFtotal=0
                    for i in seq1:
                        for j in seq2:
                            seqpair=[]
                            seqpair.append(i)
                            seqpair.append(j)
                            alignpair = egglib.Align.create(seqpair)
#                            align=egglib.Align.create(i,j)
                            stats=egglib.Align.polymorphism(alignpair)
                            DIFFsum+=stats["S"]
                            DIFFtotal+=1
                    dxy=float(DIFFsum)/float(DIFFtotal)/float(win)
                    if polym1["Pi"] is None:
                        Pi1=0
                    else:
                        Pi1=polym1["Pi"]#/float(win)
                    if polym2["Pi"] is None:
                        Pi2=0
                    else:
                        Pi2=polym2["Pi"]#/float(win)
                    Da=dxy-(Pi1+Pi2)/2
                    print scaf,pair,k*win,(k+1)*win,dxy,Da,Pi1,Pi2
