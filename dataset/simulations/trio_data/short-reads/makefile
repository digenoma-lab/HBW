.DELETE_ON_ERROR:


#simulating paternal reads 50X
paternal.fwd.fasta.gz:
	gzip -dc  paternalA.fasta.gz paternalB.fasta.gz > paternal.fasta
	wgsim -e 0.01 -1 150 -2 150 -N 800000 -h paternal.fasta paternal.fwd.fq paternal.rev.fq > paternal.wgsim.log
	seqtk seq -A paternal.rev.fq | gzip > paternal.rev.fasta.gz  
	seqtk seq -A paternal.fwd.fq | gzip > $@  
	-rm -f paternal.rev.fq paternal.fwd.fq

#simulating maternal reads 50X
maternal.fwd.fasta.gz:
	gzip -dc  maternalA.fasta.gz maternalB.fasta.gz > maternal.fasta
	wgsim -e 0.01 -1 150 -2 150 -N 800000 -h maternal.fasta maternal.fwd.fq maternal.rev.fq > maternal.wgsim.log
	seqtk seq -A maternal.rev.fq | gzip > maternal.rev.fasta.gz  
	seqtk seq -A maternal.fwd.fq | gzip > $@  
	-rm -f maternal.rev.fq maternal.fwd.fq

#simulating kid reads 55X
kid.fwd.fasta.gz:
	gzip -dc  maternalA.fasta.gz paternalB.fasta.gz > kid.fasta
	wgsim -e 0.01 -1 150 -2 150 -N 820000 -h kid.fasta kid.fwd.fq kid.rev.fq > kid.wgsim.log
	seqtk seq -A kid.rev.fq | gzip > kid.rev.fasta.gz  
	seqtk seq -A kid.fwd.fq | gzip > $@  
	-rm -f kid.rev.fq kid.fwd.fq


all : paternal.fwd.fasta.gz maternal.fwd.fasta.gz kid.fwd.fasta.gz

