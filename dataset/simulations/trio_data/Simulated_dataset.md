
# Trio simulated data



## Paternal, maternal and kid haplotypes

The strains selected to represent paternal, maternal, and kids haplotypes are:

1. Maternal Haplotypes (99.98% identity)
	- Escherichia coli BW2952 (CP001396.1) maternal A
	- Escherichia coli ER2796 (CP009644.1) maternal B

2. Paternal Haplotypes (99.9% identity)
	- Escherichia coli APEC O78 (CP004009) paternal A
	- Escherichia coli ACN001 (CP007442.1) paternal B
	
3. kids haplotypes (98.61% identity)
	- Escherichia coli BW2952 (CP001396.1) kid A (maternal)
	- Escherichia coli ACN001 (CP007442.1) Kid B (paternal)

## Short reads data

1. Paternal read set (50X coverage)
	- short-reads/paternal.fwd.fasta.gz
	- short-reads/paternal.rev.fasta.gz
1. Maternal read set (50X coverage)
	- short-reads/maternal.fwd.fasta.gz
	- short-reads/maternal.rev.fasta.gz
1. Kid read set (55X coverage)
	- short-reads/kid.fwd.fasta.gz
	- short-reads/kid.rev.fasta.gz

The short-read data was simulated using ***wgsim***.

## Long reads data

- Kid read set (60X coverage)

1. 15% average error rate.
	- long-reads/kidONT.R94.30X.p1.fasta.gz
	- long-reads/kidONT.R94.30X2.p2.fasta.gz
	
2. 5% average error rate.
	- long-reads/kidONT.R103.30X.p1.fasta.gz
	- long-reads/kidONT.R103.30X.p2.fasta.gz

the long-read data was simulated using ***pbsim2***
