## python build_rf_matrix.py  hapA_lst.txt hapB_lst.txt  dataset/simulated_data/trio_data/long-reads/kidONT.R94.30X.p1.fasta.gz ttt 

from Bio.Seq import Seq
import gzip
from Bio import SeqIO
import sys
import statistics as st

def buildKdict(filename):
 files=open(filename,"r")
 klist=list()
 kstats=list()
 for fn in files:
  hash={}
  cov=list()
  f=open(fn.rstrip(),"r")
  for l in f:
   a=l.rstrip().split("\t")
   hash[a[0]]=int(a[1])
   cov.append(int(a[1]))
  
  stats=dict(mean=int(st.mean(cov)),std=int(st.stdev(cov)))
  kstats.append(stats)
  klist.append(hash)
 

 return klist,kstats

def printndict(n,d):
 tmp=dict(zip(list(d.keys())[:n], list(d.values())[:n]))
 ksize=len(list(d.keys())[0])
 print("Kmer:",ksize,"Size:",len(d.keys()),"items:",tmp)
 return ksize

if len(sys.argv) != 5:
 print("build_rf_matrix.py kmers_hap1 kmers_hap2 long_reads output\n")
 exit(1)

#Step1: we create the dictionaries for haplotypes A and B  and show 5 elements for each of them

lhka,lhkas=buildKdict(sys.argv[1])
ksizesA=list()
print("Haplotype A")
for i in lhka:
 ksizesA.append(printndict(5,i))

lhkb,lhkbs=buildKdict(sys.argv[2])
ksizesB=list()
print("Haplotype B")
for i in lhkb:
 ksizesB.append(printndict(5,i))

if ksizesA != ksizesB:
 print(ksizesA,ksizesB,"kmers for haplotype A and B differs")
 exit(1)



#Step2: we load the long-read fastq compressed file
out=open(sys.argv[4],"w")
print("id ha_k15 hb_k15 ha_k18 hb_k18 ha_k21 hb_k21 ha_k24 hb_k24 uha_k15 uhb_k15 uha_k18 uhb_k18 uha_k21 uhb_k21 uha_k24 uhb_k24",file=out)
with gzip.open(sys.argv[3], "rt") as handle:
 for record in SeqIO.parse(handle, "fasta"):
  #print(record.id,record.seq[0:10])
  index=0
  l=[0,0,0,0,0,0,0,0]
  #filtering repetitive k-mers one sd from average
  lw=[0,0,0,0,0,0,0,0]
  for k in ksizesA:
   for i in range(0, len(record.seq)-k+1):
    s=Seq(record.seq[i:i+k])
    if str(s) in lhka[index]:
     #print(str(s),"hapA",lhka[index][str(s)])
     l[index*2]+=1
     if lhka[index][str(s)] > lhkas[index]["mean"]-lhkas[index]["std"] and lhka[index][str(s)] < lhkas[index]["mean"]+lhkas[index]["std"]:
      lw[index*2]+=1
    elif str(s.reverse_complement()) in lhka[index]:
     #print(str(s.reverse_complement()),"hapA",lhka[index][str(s.reverse_complement())])
     l[index*2]+=1
     if lhka[index][str(s.reverse_complement())] > lhkas[index]["mean"]-lhkas[index]["std"] and lhka[index][str(s.reverse_complement())] < lhkas[index]["mean"]+lhkas[index]["std"]:
      lw[index*2]+=1
    elif str(s) in lhkb[index]:
     #print(str(s),"hapB",lhkb[index][str(s)])
     l[index*2+1]+=1
     if lhkb[index][str(s)] > lhkbs[index]["mean"]-lhkbs[index]["std"] and lhkb[index][str(s)] < lhkbs[index]["mean"]+lhkbs[index]["std"]:
      lw[index*2+1]+=1
    elif str(s.reverse_complement()) in lhkb[index]:
     #print(str(s.reverse_complement()),"hapB",lhkb[index][str(s.reverse_complement())])
     l[index*2+1]+=1
     if lhkb[index][str(s.reverse_complement())] > lhkbs[index]["mean"]-lhkbs[index]["std"] and lhkb[index][str(s.reverse_complement())] < lhkbs[index]["mean"]+lhkbs[index]["std"]:
      lw[index*2+1]+=1
   index+=1

#Step3: we output the variables per long read
  print(record.id,*l,*lw,file=out)

