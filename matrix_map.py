from Bio.Seq import Seq
import gzip
from Bio import SeqIO
import sys




def parser(line):
    return Seq(line[0:15])
    



print(len(sys.argv))

if len(sys.argv) != 5:
 print("matrix_map.py hap1 hap2 long_reads output\n")
 exit(1)

	
uk15_A=open(sys.argv[1],"r")
uk15_A=uk15_A.readlines()
uk15_A=(map(parser,uk15_A))

#print(uk15_A)

uk15_B=open(sys.argv[2],"r")
uk15_B=uk15_B.readlines()
uk15_B=list(map(parser,uk15_B))

#print(uk15_B[0:5])

long_reads=sys.argv[3]

pos=0
pos_lect=0
hapA_k15=[]
hapB_k15=[]
countA_k15=0
countB_k15=0

with gzip.open(long_reads, "rt") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        while pos < len(record.seq)-15:
            seq=record.seq[pos:pos+15]
            seq_rev=seq.reverse_complement()
            x=(map(lambda x: 1 if seq==x or seq_rev==x else 0, uk15_A))
            y=(map(lambda y: 1 if seq==y or seq_rev==y else 0, uk15_B))
            if sum(x)>=1:
                countA_k15+=sum(x)
                hapA_k15.append(pos+pos_lect)
                print(sum(x))
            if sum(y)>=1:
                countB_k15+=sum(y)
                hapB_k15.append(pos+pos_lect)
                print(sum(y))
        pos_lect+=pos+15
        pos=0

cov_A=(hapA_k15[0]+hapA_k15[-1])/(pos_lect+15)
cov_B=(hapB_k15[0]+hapB_k15[-1])/(pos_lect+15)

prom_A=[]
prom_B=[]

for i in range(len(hapA_k15)-1):
    prom_A.append((hapA_k15[i]+hapA_k15[i+1])/2)

for i in range(len(hapB_k15)-1):
    prom_B.append((hapB_k15[i]+hapB_k15[i+1])/2)

promA_K15=sum(prom_A)/len(prom_A)
promB_K15=sum(prom_B)/len(prom_B)

print(countA_k15, countB_k15, cov_A, cov_B)


res=open(sys.argv[4],"w")
res.write(str(countA_k15)) #count A
res.write(str(countB_k15)) #count B
res.write(str(cov_A)) #cov A
res.write(str(cov_B)) #cov B
res.write(str(promA_K15))
res.write(str(promB_K15))

res.close()
            
            
                

