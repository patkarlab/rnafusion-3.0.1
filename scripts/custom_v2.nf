#!/usr/bin/env nextflow
nextflow.enable.dsl=2

log.info """
STARTING PIPELINE
=*=*=*=*=*=*=*=*=

Sample list: ${params.input}
"""
process coverage {
	publishDir "$PWD/Final_Output/${sampleId}/", mode: 'copy', pattern: '*.counts.bed'
	input:
		tuple val(sampleId), path(read1)
	output:
		tuple val (sampleId), file ("*.counts.bed")
	script:
	"""
	#${params.bedtools} bamtobed -i $PWD/fusioninspector/${sampleId}.consolidated.bam > ${sampleId}.bed
	${params.bedtools} bamtobed -i $PWD/star_for_arriba/${sampleId}.Aligned.out.bam | awk 'BEGIN{OFS="\t"}{ \$1="chr"\$1; print }' > ${sampleId}.bed
	#${params.bedtools} bamtobed -i $PWD/picard/${sampleId}.bam | awk 'BEGIN{OFS="\t"}{ \$1="chr"\$1; print }' > ${sampleId}.bed
	#${params.bedtools} bamtobed -i $PWD/star_for_squid/${sampleId}.Aligned.sortedByCoord.out.bam | awk 'BEGIN{OFS="\t"}{ \$1="chr"\$1; print }' > ${sampleId}.bed
	#${params.bedtools} bamtobed -i $PWD/star_for_starfusion/${sampleId}.Aligned.sortedByCoord.out.bam | awk 'BEGIN{OFS="\t"}{ \$1="chr"\$1; print }' > ${sampleId}.bed
	#${params.bedtools} bamtobed -i $PWD/samtools/${sampleId}_chimeric.bam | awk 'BEGIN{OFS="\t"}{ \$1="chr"\$1; print }' > ${sampleId}.bed
	${params.bedtools} coverage -counts -a ${read1} -b ${sampleId}.bed > ${sampleId}.counts.bed

	"""
}

process bam {
	input:
		tuple val(sampleId), path(read1)
	output:
		tuple val (sampleId), file ("*.sorted.bam"), file ("*.sorted.bam.bai")
	script:
	"""
	${params.samtools} sort $PWD/star_for_arriba/${sampleId}.Aligned.out.bam -o ${sampleId}.sorted.bam
	${params.samtools} index ${sampleId}.sorted.bam
	cp ${sampleId}.sorted.bam* $PWD/Final_Output/${sampleId}/

	"""
}

process file_copy {
	input:
		tuple val(sampleId), file (coverage_bed), file (bamFile), file(bamBai)
	script:
	"""
	#if [ -f ${PWD}/arriba/${sampleId}.arriba.fusions.tsv ]; then
	#	cp ${PWD}/arriba/${sampleId}.arriba.fusions.tsv ${PWD}/Final_Output/${sampleId}/
	#fi

	if [ -f ${PWD}/arriba_visualisation/${sampleId}_combined_fusions_arriba_visualisation.pdf ]; then
		cp ${PWD}/arriba_visualisation/${sampleId}_combined_fusions_arriba_visualisation.pdf ${PWD}/Final_Output/${sampleId}/ 
	fi

	#if [ -f ${PWD}/fusioncatcher/${sampleId}.fusioncatcher.fusion-genes.txt ]; then
	#	cp ${PWD}/fusioncatcher/${sampleId}.fusioncatcher.fusion-genes.txt ${PWD}/Final_Output/${sampleId}/
	#fi

	if [ -f ${PWD}/fusioncatcher/${sampleId}.fusioncatcher.summary.txt ]; then
		${params.sed_sh} ${PWD}/fusioncatcher/${sampleId}.fusioncatcher.summary.txt
		#cp ${PWD}/fusioncatcher/${sampleId}.fusioncatcher.summary.txt ${PWD}/Final_Output/${sampleId}/
	fi

	#if [ -f ${PWD}/starfusion/${sampleId}.starfusion.fusion_predictions.tsv ]; then
	#	cp ${PWD}/starfusion/${sampleId}.starfusion.fusion_predictions.tsv ${PWD}/Final_Output/${sampleId}/
	#fi

	if [ -d ${PWD}/fusionreport/${sampleId} ]; then
		cp -r ${PWD}/fusionreport/${sampleId} ${PWD}/Final_Output/${sampleId}/${sampleId}_fusionreport
	fi

	if [ -f ${PWD}/fusioninspector/${sampleId}.fusion_inspector_web.html ]; then 
		cp -r ${PWD}/fusioninspector/${sampleId}.fusion_inspector_web.html ${PWD}/Final_Output/${sampleId}/
	fi

	if [ -f ${PWD}/stringtie/${sampleId}.transcripts.gtf ]; then
		echo -e "chr\tstringtie\tseqname source\tfeature\tstart\tend\tscore\tstrand\tframe attributes" > ${sampleId}_stringtie.tsv
		grep -v '^#' ${PWD}/stringtie/${sampleId}.transcripts.gtf >> ${sampleId}_stringtie.tsv		
	fi

	python3 ${params.merge_csvs_script} ${sampleId} ${PWD}/Final_Output/${sampleId}/${sampleId}.xlsx ${coverage_bed} ${PWD}/arriba/${sampleId}.arriba.fusions.tsv ${PWD}/fusioncatcher/${sampleId}.fusioncatcher.fusion-genes.txt ${PWD}/fusioncatcher/${sampleId}.fusioncatcher.summary.txt ${PWD}/starfusion/${sampleId}.starfusion.fusion_predictions.tsv ${sampleId}_stringtie.tsv
	"""
}

workflow COVERAGE {
	Channel
		.fromPath(params.input)
		.splitCsv(header:false)
		.view () { row -> "${row[0]},${row[1]}" }
		.set { samples_ch }

	main:
	coverage(samples_ch)
	bam(samples_ch)
	file_copy(coverage.out.join(bam.out))
}
