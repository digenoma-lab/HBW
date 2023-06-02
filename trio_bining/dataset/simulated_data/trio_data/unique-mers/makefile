.DELETE_ON_ERROR:
#default k-mer size
KMER:=15

uk${KMER}/hapA_only_kmers.txt.gz:
	mkdir -p uk${KMER}
	find-unique-kmers  -k${KMER} -p 4 -o uk${KMER} paternal.fwd.fasta.gz,paternal.rev.fasta.gz maternal.fwd.fasta.gz,maternal.rev.fasta.gz	> uk${KMER}.log 2>uk${KMER}.err
	find uk${KMER} -type f | grep -v "only_kmers" | xargs rm -f

all: uk${KMER}/hapA_only_kmers.txt.gz

