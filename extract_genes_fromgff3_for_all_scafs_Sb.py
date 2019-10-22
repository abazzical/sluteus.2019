#!/usr/local/pkg/MCBPython2010/epd-6.2-2-rh5-x86_64/bin/python
import sys
import subprocess
import operator
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna

##############################
gfffilename = sys.argv[1]#gff file
# gfffilename = "Suilu4_all_genes_20170628.gff"
# type = "NsSing"
type = sys.argv[2]#type of dataset: i.e. with clones, without clones, with Ns...
gff = open(gfffilename, "r")#gff file name
###############################

#####################
#reads the gff and store everything in dictionaries
gene_info = {}#CDS names, start, stop for each CDS
gene_strand = {}#strand of each gene
parent_scaffold = {}#parent scaffold of each gene
child_CDS_ids={}#list of CDS of each gene; I need one to reverse order of CDS parsing  when on negative strand

for gffline in gff:
    gffline = gffline.strip()
    splitline = gffline.split("\t")
    n=len(splitline)
    if ((n==9) and ((splitline[2] == "CDS") or (splitline[2] == "exon") or (splitline[2] == "start_codon"))):
        if (splitline[2] == "CDS"):
            info= splitline[-1].split(";")
            ids=info[-2].split(" ")
            proteinId=ids[-1]
            if proteinId not in gene_info:
                gene_info[proteinId]={}
                child_CDS_ids[proteinId]=[]
                gene_strand[proteinId]=splitline[6]
                parent_scaffold[proteinId]=splitline[0]
            ids2=info[-1].split(" ")
            CDS=ids2[-1]
            child_CDS_ids[proteinId].append(CDS)
            gene_info[proteinId][CDS] = {}#one dictionary per CDS; ids2[-1] is CDS name 
            gene_info[proteinId][CDS]["start"]=int(splitline[3])
            gene_info[proteinId][CDS]["stop"]=int(splitline[4])
            print "reading gff for "+proteinId+", scaffold "+splitline[0]

c=0
for gene in gene_info:#I use 'gene' instead of 'proteinId'
    c+=1
    print c,gene
    OUT=open("/Users/abazzical/Documents/suillusFastSTRUCTURE/20190130_onlyDikaryons/calculateDxy/2019-01-28dikconcatscafs_fasta/%s.fasta" %(gene), 'w')
    
   
    #handle = open("/home/sbranco/Sb_SNPs/11.1.Prep_files_with_N_singletons/HiQ_Sb_%s_woSb32Sb7_indivscafs_fasta_concat/HiQ_Sb_%s_woSb32Sb7_%s.fasta" %(type,type,parent_scaffold[gene]), "rU")


    handle = open("/Users/abazzical/Documents/suillusFastSTRUCTURE/20190130_onlyDikaryons/calculateDxy/2019-01-28dikconcatscafs_fasta/%s.fasta" %(parent_scaffold[gene]), "rU")
    scafseq = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
    handle.close()#            mRNA_info[mRNA]["CDSs"] = sorted(mRNA_info[mRNA]["CDSs"],key=lambda CDS:CDS_info[CDS]["start"], reverse=True)
    d=0
    for ind in scafseq:
        print "writing data for gene %s, individual %s my dear" %(gene,ind)
        d+=1
        out=Seq("", generic_dna)
        if gene_strand[gene] == "-":
            for CDS in child_CDS_ids[gene]:
                print CDS+"-",
                cds=Seq(str(scafseq[ind][gene_info[gene][CDS]["start"]-1:gene_info[gene][CDS]["stop"]].seq),generic_dna)
                cds=cds.reverse_complement()
                out=cds+out
        else:
            for CDS in child_CDS_ids[gene]:
                print CDS+"+",
                cds=Seq(str(scafseq[ind][gene_info[gene][CDS]["start"]-1:gene_info[gene][CDS]["stop"]].seq),generic_dna)
                out=out+cds
        parts=scafseq[ind].id.split('_')
        name=parts[-1]
        a=SeqRecord(out,id=name,description="")
        OUT.write(a.format("fasta"))
