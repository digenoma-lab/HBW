.DELETE_ON_ERROR:

CPU=44

MM=/home/adigenova/HSCAFF/CANCER/CELLLINE/analysis/SVs/software/minimap2/minimap2
SM=/home/adigenova/HSCAFF/CANCER/CELLLINE/analysis/SVs/software/samtools-1.10/SM/bin/samtools
SN=/home/adigenova/HSCAFF/CANCER/CELLLINE/analysis/SVs/software/Sniffles-1.0.11/bin/sniffles-core-1.0.11/sniffles
REF=wengan.hg002.fa
READS=HG002_promethion_hac.fastq.gz
TIME=/home/adigenova/HSCAFF/CANCER/CELLLINE/SKBR3/time

#map and index
$(PREF).srt.bam:
	${MM} -t${CPU} -x map-ont --MD -Y -L -a ${REF} ${READS}  2>$(subst .srt.bam,,$@).mm.err | ${SM} view -Sbh -q 20 - | ${SM} sort -@${CPU} -o  $@ - 
	${SM} index -@${CPU} $@

$(PREF).sv_def.vcf:$(PREF).srt.bam
	${SN} --report_seq -t ${CPU} -m $< -v $@ 2>$@.err >$@.log

all:  $(PREF).srt.bam $(PREF).sv_def.vcf
