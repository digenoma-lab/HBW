.DELETE_ON_ERROR:

CPU=44

MM=/home/adigenova/HSCAFF/CANCER/CELLLINE/analysis/SVs/software/minimap2/minimap2
#SM=/home/adigenova/binaries/samtools-0.1.19/samtools
SM=/home/adigenova/HSCAFF/CANCER/CELLLINE/analysis/SVs/software/samtools-1.10/SM/bin/samtools
SN=/home/adigenova/HSCAFF/CANCER/CELLLINE/analysis/SVs/software/Sniffles-1.0.11/bin/sniffles-core-1.0.11/sniffles
#REF=/home/adigenova/HSCAFF/CANCER/CELLLINE/analysis/SVs/reference/hs38DH.fa
REF=wengan.hg002.fa
#READS=A01_promethion_hac.fastq.gz A02_promethion_hac.fastq.gz A03_promethion_hac.fastq.gz B01_promethion_hac.fastq.gz B02_promethion_hac.fastq.gz B03_promethion_hac.fastq.gz

READS=HG002_promethion_hac.fastq.gz

VPATH=reads

TIME=/home/adigenova/HSCAFF/CANCER/CELLLINE/SKBR3/time

map_files=$(patsubst %.fastq.gz,%.srt.bam,$(READS))

%.srt.bam: %.fastq.gz
	${MM} -t${CPU} -x map-ont --MD -Y -L -a ${REF} $<  2>$(subst .srt.bam,,$@).mm.err | ${SM} view -Sbh -q 20 - | ${SM} sort -@${CPU} -o  $@ - 

map2ref: $(map_files)

index_files=$(addsuffix .bai,${map_files})
%.srt.bam.bai:%.srt.bam
	${SM} index -@${CPU} $< $@ 

index: $(index_files)


#/sniffles -m A01_promethion_hac.srt.bam -v  A01_promethion_hac_def

svs_def=$(patsubst %.srt.bam,%.sv_def.vcf,$(map_files))
#${SM} index -t ${CPU} --tmp_file ${PWD}/$@.tmp -m $< -v $@ 
%.sv_def.vcf:%.srt.bam 
	${SN} --report_seq -t ${CPU} -m $< -v $@ 2>$@.err >$@.log

calldef: $(svs_def) index

all: map2ref index calldef
