#!/usr/bin/env nextflow
nextflow.enable.dsl=2

log.info """
STARTING PIPELINE
=*=*=*=*=*=*=*=*=

Sample list: ${params.input}
"""

vardict = params.vardict

include { VAR_RNA } from '../workflows/var_rna.nf'
include { ANNOVAR as ANNOVAR_VARDICT } from '../modules/local/annovar/annotate/main'
include { FORMAT_VARDICT } from '../modules/local/python/format_vardict/main'

process COUNTS {
	tag "${sampleId}"
	publishDir "${PWD}/Final_Output/${sampleId}/", mode: 'copy'
	input:
		tuple val(sampleId), path(bedfile), path(arriba_bam)
	output:
		tuple val (sampleId), file("${sampleId}.counts.bed")
	script:
	"""
	${params.bedtools} coverage -counts -a ${bedfile} -b ${arriba_bam} > ${sampleId}.counts.bed
	"""
}

process BAM {
	tag "${sampleId}"
	label 'process_inter'
	publishDir "${PWD}/Final_Output/${sampleId}/", mode: 'copy'
	input:
		tuple val(sampleId), path(bedfile), path(arriba_bam)
	output:
		tuple val (sampleId), file ("${sampleId}.sorted.bam"), file ("${sampleId}.sorted.bam.bai")
	script:
	"""
	${params.samtools} sort -@ ${task.cpus} ${arriba_bam} -o ${sampleId}.sorted.bam
	${params.samtools} index -@ ${task.cpus} ${sampleId}.sorted.bam
	"""
}

process FILE_COPY {
	input:
		tuple val(sampleId), file (coverage_bed), file(vardict_csv), file(haplotypecaller_csv)
	output:
		val (sampleId)
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

	python3 ${params.merge_csvs_script} ${sampleId} ${PWD}/Final_Output/${sampleId}/${sampleId}.xlsx ${coverage_bed} ${PWD}/arriba/${sampleId}.arriba.fusions.tsv ${PWD}/fusioncatcher/${sampleId}.fusioncatcher.fusion-genes.txt ${PWD}/fusioncatcher/${sampleId}.fusioncatcher.summary.txt ${PWD}/starfusion/${sampleId}.starfusion.fusion_predictions.tsv ${vardict_csv} ${haplotypecaller_csv}

	"""
}

process VARDICT {
	tag "${sampleId}"
	label 'process_inter'
	input:
		tuple val (sampleId), file (sorted_bam), file (sorted_bam_bai), path(bedfile), path(squid_bam)
	output:
		tuple val(sampleId), file ("${sampleId}.vardict.vcf.gz"), file("${sampleId}.vardict.vcf.gz.tbi")
	script:
	"""
	VarDict -G ${params.genome} -f 0.05 -N ${sampleId} -b ${sorted_bam} -c 1 -S 2 -E 3 -g 4 -th ${task.cpus} -L 10000000 --verbose ${bedfile} | teststrandbias.R | var2vcf_valid.pl | gzip > ${sampleId}.vardict.vcf.gz
	touch ${sampleId}.vardict.vcf.gz.tbi
	"""
}

workflow COVERAGE {
	Channel
		.fromPath(params.input)
		.splitCsv(header:false)
		.view () { row -> "${row[0]},${row[1]}" }
		.set { samples_ch }

	main:
	COUNTS(samples_ch)
	BAM(samples_ch)
	haplotypecaller_csv = VAR_RNA(samples_ch)
	VARDICT(BAM.out.join(samples_ch))
	ANNOVAR_VARDICT(VARDICT.out, vardict)
	FORMAT_VARDICT(ANNOVAR_VARDICT.out)
	FILE_COPY(COUNTS.out.join(FORMAT_VARDICT.out.join(haplotypecaller_csv)))
}
