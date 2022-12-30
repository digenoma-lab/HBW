



.DELETE_ON_ERROR:

#HOME=/Users/adigenova
#BWA=${HOME}/Programs/bwa/bwa
#PILON=${HOME}/Programs/Pilon/pilon-1.24.jar

MM=/home/adigenova/HSCAFF/CANCER/CELLLINE/analysis/SVs/software/minimap2/minimap2
SM=/home/adigenova/HSCAFF/CANCER/CELLLINE/analysis/SVs/software/samtools-1.10/SM/bin/samtools
RACON=/home/adigenova/HSCAFF/CANCER/CELLLINE/analysis/SVs/software/racon/build/bin/racon

CPU=44
ASM=wengan.hg002.fa
BED=wengan.hg002.bed
#ASM=multi_w.fasta
#BED=multi_w.bed
READS=HG002_promethion_hac.fastq.gz
#PREF=multi_w
PREF=asm1
#RACON=${HOME}/Programs/racon/build/bin/racon
#RACON=${HOME}/Git/racon/build/bin/racon

#first round of long-read polishing
$(PREF).r1.paf:
	${MM} -x map-ont -t ${CPU} ${ASM} ${READS} > $@
$(PREF).r1.fa:$(PREF).r1.paf
	${RACON} -u -b -t ${CPU} ${READS} $< ${ASM} ${BED} > $@ 2> $(PREF).r1.racon.log
	egrep "^NWD" $(PREF).r1.racon.log | awk '{print $$2"\t"$$3"\t"$$4"\t"$$4}' > $(PREF).r1.racon.bed

#second round
$(PREF).r2.paf:$(PREF).r1.fa
	${MM} -x map-ont -t ${CPU} $< ${READS} > $@
$(PREF).r2.fa:$(PREF).r2.paf
	${RACON} -u -b -t ${CPU} ${READS} $< $(PREF).r1.fa $(PREF).r1.racon.bed > $@ 2> $(PREF).r2.racon.log
	egrep "^NWD" $(PREF).r2.racon.log | awk '{print $$2"\t"$$3"\t"$$4"\t"$$4}' > $(PREF).r2.racon.bed

#third round
$(PREF).r3.paf:$(PREF).r2.fa
	${MM} -x map-ont -t ${CPU} $< ${READS} > $@
$(PREF).r3.fa:$(PREF).r3.paf
	${RACON} -u -b -t ${CPU} ${READS} $< $(PREF).r2.fa $(PREF).r2.racon.bed > $@ 2> $(PREF).r3.racon.log
	egrep "^NWD" $(PREF).r3.racon.log | awk '{print $$2"\t"$$3"\t"$$4"\t"$$4}' > $(PREF).r3.racon.bed

#short read polishing	
#$(PREF).r3.fa.sa:  $(PREF).r3.fa
#	 ${BWA} index $<
#$(PREF).sort.bam:$(PREF).r3.fa.sa
#	 ${BWA} mem -t ${CPU} $(PREF).r3.fa  reads/EC.50X.R1.fastq.gz reads/EC.50X.R2.fastq.gz 2>$(PREF).bwa.log | \
#	 ${SM} view -b - | ${SM} sort -o $@ -
#	 ${SM} index $@ 
#all: ec_Wd_or1.SPolished.asm.wengan.fasta.sa w.sort.bam $(PREF).r1.fa $(PREF).r2.fa $(PREF).r3.fa
#all: $(PREF).r1.fa $(PREF).r2.fa $(PREF).r3.fa $(PREF).sort.bam
all: $(PREF).r1.fa $(PREF).r2.fa $(PREF).r3.fa
