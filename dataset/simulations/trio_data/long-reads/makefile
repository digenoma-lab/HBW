.DELETE_ON_ERROR:


#nanopore reads model R94 average identity 85%
kidONT.R94.30X.p1.fasta.gz:
	pbsim --id-prefix K --depth 30  --hmm_model ~/micromamba/envs/mummner/data/R94.model --prefix kidONT.R94.30X kid.fasta --length-min 1000  2> kidONT.R94.30X.log
	seqtk seq -A kidONT.R94.30X_0001.fastq | gzip > $@
	seqtk seq -A kidONT.R94.30X_0002.fastq | gzip > kidONT.R94.30X2.p2.fasta.gz

#nanopore reads model R103 average identity 95%
kidONT.R103.30X.p1.fasta.gz:
	pbsim --id-prefix K --depth 30 --accuracy-mean 0.95  --hmm_model ~/micromamba/envs/mummner/data/R103.model --prefix kidONT.R103.30X kid.fasta --length-min 1000  2> kidONT.R103.30X.log
	seqtk seq -A kidONT.R103.30X_0001.fastq | gzip > $@
	seqtk seq -A kidONT.R103.30X_0002.fastq | gzip > kidONT.R103.30X.p2.fasta.gz

all: kidONT.R94.30X.p1.fasta.gz kidONT.R103.30X.p1.fasta.gz
