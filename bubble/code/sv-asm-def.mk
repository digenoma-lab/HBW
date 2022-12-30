.DELETE_ON_ERROR:

CPU=8
REF=ref.fa
HAP1=asm1.hap1.fa
HAP2=asm1.hap2.fa

PREF=asm

BIN=/miniconda/envs/callers/bin

${PREF}.aln.hap1.sam:
	${BIN}/minimap2 --paf-no-hit -a -x asm5 --cs -r2k -t ${CPU} ${REF} ${HAP1} > ${PREF}.aln.hap1.sam
${PREF}.aln.hap2.sam:${PREF}.aln.hap1.sam
	${BIN}/minimap2 --paf-no-hit -a -x asm5 --cs -r2k -t ${CPU} ${REF} ${HAP2} > ${PREF}.aln.hap2.sam
${PREF}.aln.hap1.sorted.bam:${PREF}.aln.hap2.sam
	${BIN}/samtools sort -m1G -@${CPU} -o ${PREF}.aln.hap1.sorted.bam ${PREF}.aln.hap1.sam
	${BIN}/samtools index ${PREF}.aln.hap1.sorted.bam
${PREF}.aln.hap2.sorted.bam:${PREF}.aln.hap1.sorted.bam
	${BIN}/samtools sort -m1G -@${CPU} -o ${PREF}.aln.hap2.sorted.bam ${PREF}.aln.hap2.sam
	${BIN}/samtools index ${PREF}.aln.hap2.sorted.bam

#run sv-asm caller using values from the paper
${PREF}-sv-asm-dip/variants.vcf: ${PREF}.aln.hap2.sorted.bam
	mkdir ${PREF}-sv-asm-dip
	${BIN}/svim-asm diploid ${PREF}-sv-asm-dip ${PREF}.aln.hap1.sorted.bam ${PREF}.aln.hap2.sorted.bam ${REF} \
	 --min_sv_size 30 --tandem_duplications_as_insertions --interspersed_duplications_as_insertion
	
#run sv-asm caller using the haploid mode (primary contigs)
${PREF}-sv-asm-hap/variants.vcf: ${PREF}.aln.hap2.sorted.bam
	mkdir ${PREF}-sv-asm-hap
	${BIN}/svim-asm haploid ${PREF}-sv-asm-hap ${PREF}.aln.hap1.sorted.bam ${REF} \
	 --min_sv_size 30 --tandem_duplications_as_insertions --interspersed_duplications_as_insertions

all: ${PREF}-sv-asm-dip/variants.vcf ${PREF}-sv-asm-hap/variants.vcf
