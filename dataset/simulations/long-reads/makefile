.DELETE_ON_ERROR:


#nanopore reads model R94 average identity 85%
patONT.R94.30X_0001.fasta.gz:
	pbsim --id-prefix P --depth 30  --hmm_model ~/micromamba/envs/mummner/data/R94.model --prefix patONT.R94.30X paternal.fasta --length-min 1000  2> patONT.R94.30X.log
	seqtk seq -A patONT.R94.30X_0001.fastq | gzip > $@

matONT.R94.30X_0001.fasta.gz:
	pbsim --id-prefix M --depth 30  --hmm_model ~/micromamba/envs/mummner/data/R94.model --prefix matONT.R94.30X maternal.fasta --length-min 1000  2> matONT.R94.30X.log
	seqtk seq -A matONT.R94.30X_0001.fastq | gzip > $@

#nanopore reads model R103 average identity 95%
patONT.R103.30X_0001.fasta.gz:
	pbsim --id-prefix P --depth 30 --accuracy-mean 0.95  --hmm_model ~/micromamba/envs/mummner/data/R103.model --prefix patONT.R103.30X paternal.fasta --length-min 1000  2> patONT.R103.30X.log
	seqtk seq -A patONT.R103.30X_0001.fastq | gzip > $@

matONT.R103.30X_0001.fasta.gz:
	pbsim --id-prefix M --depth 30 --accuracy-mean 0.95 --hmm_model ~/micromamba/envs/mummner/data/R103.model --prefix matONT.R103.30X maternal.fasta --length-min 1000  2> matONT.R103.30X.log
	seqtk seq -A matONT.R103.30X_0001.fastq | gzip > $@


all: patONT.R94.30X_0001.fasta.gz matONT.R94.30X_0001.fasta.gz patONT.R103.30X_0001.fasta.gz matONT.R103.30X_0001.fasta.gz
