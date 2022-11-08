.DELETE_ON_ERROR:
KMER=15

uk${KMER}.class.log:
	mkdir -p classified.${KMER}.p1
	classify-by-kmers kidONT.R94.30X.p1.fasta.gz  uk${KMER}/hapA_only_kmers.txt uk${KMER}/hapB_only_kmers.txt --haplotype-a-out-prefix classified.${KMER}.p1/paternal --haplotype-b-out-prefix classified.${KMER}.p1/maternal --unclassified-out-prefix classified.${KMER}.p1/unclassified > $@.p1 2> $@.p1.err
	mkdir -p classified.${KMER}.p2
	classify-by-kmers kidONT.R94.30X.p2.fasta.gz  uk${KMER}/hapA_only_kmers.txt uk${KMER}/hapB_only_kmers.txt --haplotype-a-out-prefix classified.${KMER}.p2/paternal --haplotype-b-out-prefix classified.${KMER}.p2/maternal --unclassified-out-prefix classified.${KMER}.p2/unclassified > $@.p2 2> $@.p2.err
	cat $@.p1 $@.p2 > $@
all : uk${KMER}.class.log
	
