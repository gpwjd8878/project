#!usr/bin/python

import argparse
import coloredlogs
import logging
import os
import subprocess
import sys

def parse_argument():
	"""Parsing arguments"""
	parser = argparse.ArgumentParser()
	parser.add_argument(
			'-s',
			'--sample',
			dest='sample',
			help='sample name',
			required=True
			)
	parser.add_argument(
			'-c',
			'--configure',
			dest='configure_file',
			default='/data/CLINICAL/CLINICAL_conf.txt',
			help=(
					'configure file.'
					'default: /data/CLINICAL/CLINICAL_conf.txt'
					)
			)
	parser.add_argument(
			'-o',
			'--output_path',
			dest='output_path',
			default=os.getcwd(),
			help=(
					'output path'
					'default: current working directory'
					)
			)
	parser.add_argument(
			'-t',
			'--target_bed',
			dest='target_bed',
			default='',
			help='target bed'
			)
	parser.add_argument(
			'-vc',
			'--variant_calling',
			dest='variant_calling',
			default='sm',
			help=(
					'type of variant calling'
					'default: somatic calling'
					)
			)

	args = parser.parse_args()
	sample = args.sample
	output_path = args.output_path
	configure_file = args.configure_file
	target_bed = args.target_bed
	variant_calling = args.variant_calling
	return sample, configure_file, output_path, target_bed, variant_calling


def load_configure(configure_file):
	"""Loading a configure file"""
	conf_dict = {}
	with open(configure_file, 'r') as read_conf:
		for line in read_conf:
			if not line:
				break
			if line.startswith('#'):
				continue
			S = line.rstrip().split('=')
			if len(S) > 1:
				key = S[0].strip()
				value = S[1].strip()
				conf_dict[key] = value
	read_conf.close()
	return conf_dict


def make_logger(sample):
	"""Prepare logging system"""
	log_file = (
			'{sample}.analysis.log'
			.format(sample=sample)
			)
	logger = logging.getLogger(log_file)
	analysis_logger = logging.FileHandler(log_file)
	analysis_logger.setLevel(logging.DEBUG)
	formatter = coloredlogs.ColoredFormatter(
			'%(asctime)s - %(name)s - %(levelname)s - %(message)s'
			)
	analysis_logger.setFormatter(formatter)
	logger.addHandler(analysis_logger)
	coloredlogs.install(level='DEBUG')
#	logger = logging.getLogger(log_file)
#	coloredlogs.install(level='DEBUG')
#	logger.setLevel(logging.DEBUG)
#	file_handler = logging.FileHandler(log_file)
#	logger.addHandler(file_handler)
	return logger, log_file


def fastq_qc(sample, output_path, order):
	"""QC: running fastqc"""
	logger.info('---QC---')
	fastq_suffix = ''
	if order == 1:
		fastq_suffix = 'fastq.gz'
	elif order == 2:
		fastq_suffix = 'adapterTrimmed.fastq.gz'
	cmd = (
		'fastqc '
		'-t 2 '
		'--nogroup '
		'-o {output_path} '
		'{sample}_1.{fastq_suffix} {sample}_2.{fastq_suffix} '
		'>> {log_file} 2>&1 '
		.format(
				fastq_suffix=fastq_suffix,
				log_file=log_file,
				output_path=output_path,
				sample=sample
				)
		)
	logger.info(cmd)
	os.system(cmd)
	return True


def fastq_adtrim(sample, output_path):
	"""Trimming adapters """
	logger.info('---trimming---')
	cmd = (
		'fastp '
		'-i {output_path}/{sample}_1.fastq.gz '
		'-I {output_path}/{sample}_2.fastq.gz '
		'-o {output_path}/{sample}_1.adapterTrimmed.fastq.gz '
		'-O {output_path}/{sample}_2.adapterTrimmed.fastq.gz '
		'-x -y -c -3 -p -w 2 '
		'-j {output_path}/{sample}_fastp.json '
		'-h {output_path}/{sample}_fastp.html '
		'-R "{sample} fastp" '
		'--adapter_sequence {adapter1} '
		'--adapter_sequence_r2 {adapter2} '
		'>> {log_file} 2>&1 '
		.format(
				adapter1=conf_dict['adapter1'],
				adapter2=conf_dict['adapter2'],
				log_file=log_file,
				output_path=output_path,
				sample=sample
				)
		)
	logger.info(cmd)
	os.system(cmd)
	return True


def fastq_mapping(sample, output_path):
	"""Mapping fastq files to reference genome"""
	logger.info('---mapping---')
	cmd = (
		'{bwa} mem '
		'-R "@RG\\tID:RGID\\tPL:illumina\\tLB:Target_Panel\\tSM:{sample}" '
		'-M -t 1 '
		'{reference} '
		'{output_path}/{sample}_1.adapterTrimmed.fastq.gz '
		'{output_path}/{sample}_2.adapterTrimmed.fastq.gz '
		'> {output_path}/{sample}.sam '
		'\n'
		'{sambamba} view '
		'-t 2 -S -f bam '
		'{output_path}/{sample}.sam | '
		'{sambamba} sort '
		'-t 2 -o {output_path}/{sample}.sorted.bam /dev/stdin'
		'\n'
		'{sambamba} index '
		'-t 2 '
		'{output_path}/{sample}.sorted.bam '
		'>> {log_file} 2>&1 '
		.format(
				bwa=conf_dict['bwa'],
				log_file=log_file,
				reference=conf_dict['ref'],
				output_path=output_path,
				sambamba=conf_dict['sambamba'],
				sample=sample
				)
		)
	logger.info(cmd)
	os.system(cmd)
	return True


def mark_dup(sample, output_path, target_bed):
	"""Marking duplicates"""
	logger.info('---Markduplication---')
	cmd = (
		'{sambamba} markdup '
		'-t 2 '
		'{output_path}/{sample}.sorted.bam '
		'{output_path}/{sample}.markdup.bam '
		'\n'
		'{samtools} view '
		'-@ 2 -b -q 20 '
		.format(
				output_path=output_path,
				sambamba=conf_dict['sambamba'],
				samtools=conf_dict['samtools'],
				sample=sample
				)
		)
	if target_bed == '':
		pass
	else:
		cmd += (
		'-L {target_bed} '
		.format(target_bed=conf_dict[target_bed])
		)
	cmd += (
		'{output_path}/{sample}.markdup.bam '
		'> {output_path}/{sample}.HQ.target.bam '
		'\n'
		'{sambamba} index '
		'-t 2 '
		'{output_path}/{sample}.HQ.target.bam '
		'>> {log_file} 2>&1 '
		.format(
				log_file=log_file,
				output_path=output_path,
				sambamba=conf_dict['sambamba'],
				sample=sample
				)
		)
	logger.info(cmd)
	os.system(cmd)
	return True


def make_recal(sample):
	"""Making a recal bam"""
	logger.info('---BaseRecalibration---')
	cmd = (
		'{gatk4v1} '
		'--java-options "-XX:ParallelGCThreads=4 -Xmx6g" '
		'BaseRecalibrator '
		'-R {reference} '
		'-I {sample}.HQ.target.bam '
		'--known-sites {knwon_dbsnp138} '
		'--known-sites {known_1k_snps} '
		'--known-sites {knwon_Mills_1k} '
		'-O {sample}.recal.table '
		'\n'
		'{gatk4v1} '
		'--java-options "-XX:ParallelGCThreads=4 -Xmx6g" '
		'ApplyBQSR '
		'-R {reference} '
		'-I {sample}.HQ.target.bam '
		'-bqsr {sample}.recal.table '
		'-O {sample}.recal.bam '
		.format(
				gatk4v1=conf_dict['gatk4v1'],
				known_1k_snps=conf_dict['known_1k_snps'],
				knwon_dbsnp138=conf_dict['known_dbsnp138'],
				knwon_Mills_1k=conf_dict['known_Mills_1k'],
				reference=conf_dict['ref'],
				sample=sample
				)
		)
	logger.info(cmd)
	os.system(cmd)
	return True


def somatic_calling(sample, target_bed):
	"""Somatic variants calling: MuTect2 of GATK4v0, GATK4v1"""
	logger.info('---SomaticVariantsCalling---')
	cmd = (
		'{gatk4v0} '
		'--java-options "-XX:ParallelGCThreads=4 -Xmx6g" '
		'Mutect2 '
		'-R {reference} '
		'-I {sample}.recal.bam '
		'-tumor {sample} '
		.format(
				gatk4v0=conf_dict['gatk4v0'],
				reference=conf_dict['ref'],
				sample=sample
				)
		)
	if target_bed == '':
		pass
	else:
		cmd += (
		'-L {target_bed} '
		.format(target_bed=conf_dict[target_bed])
		)
	cmd += (
		'-O {sample}.m2.raw.v0.vcf '
		'-new-qual true '
		'-RF OverclippedReadFilter '
		'--dont-use-soft-clipped-bases '
		'\n'
		'{gatk4v0} '
		'--java-options "-XX:ParallelGCThreads=4 -Xmx6g" '
		'FilterMutectCalls '
		'-V {sample}.m2.raw.v0.vcf '
		'-O {sample}.m2.f1.v0.vcf '
		'--min-pcr-slippage-size 6 '
		'-new-qual '
		'\n'
		'{gatk4v0} '
		'--java-options "-XX:ParallelGCThreads=4 -Xmx6g" '
		'FilterAlignmentArtifacts '
		'-R {reference} '
		'-I {sample}.recal.bam '
		'--bwa-mem-index-image {reference}.img '
		'--dont-skip-filtered-variants '
		'-V {sample}.m2.f1.v0.vcf '
		'-O {sample}.m2.f2.v0.vcf '
		'\n'
		'{gatk4v0} '
		'--java-options "-XX:ParallelGCThreads=4 -Xmx6g" '
		'CollectSequencingArtifactMetrics '
		'-I {sample}.recal.bam '
		'-O {sample}_artifact_v0 '
		'--FILE_EXTENSION ".v0.txt" '
		'-R {reference} '
		'\n'
		'{gatk4v0} '
		'--java-options "-XX:ParallelGCThreads=4 -Xmx6g" '
		'FilterByOrientationBias '
		'-AM G/T -AM C/T '
		'-P {sample}_artifact_v0.pre_adapter_detail_metrics.v0.txt '
		'-V {sample}.m2.f2.v0.vcf '
		'-O {sample}.m2.f3.v0.vcf '
		'\n'
		'{vt} '
		'decompose '
		'-s -o {sample}.m2.f3.decompose.v0.vcf '
		'{sample}.m2.f3.v0.vcf '
		'\n'
		'{vt} '
		'normalize '
		'-r {reference} '
		'-o {sample}.m2.f3.decompose.norm.v0.vcf '
		'{sample}.m2.f3.decompose.v0.vcf '
		.format(
				sample=sample,
				gatk4v0=conf_dict['gatk4v0'],
				reference=conf_dict['ref'],
				vt=conf_dict['vt']
				)

		)
	logger.info(cmd)
	os.system(cmd)

	cmd2 = (
		'{gatk4v1} '
		'--java-options "-XX:ParallelGCThreads=4 -Xmx6g" '
		'Mutect2 '
		'-R {reference} '
		'-I {sample}.recal.bam '
		'-tumor {sample} '
		.format(
				gatk4v1=conf_dict['gatk4v1'],
				reference=conf_dict['ref'],
				sample=sample
				)
		)
	if target_bed == '':
		pass
	else:
		cmd2 += (
		'-L {target_bed} '
		.format(target_bed=conf_dict[target_bed])
		)
	cmd2 += (
		'-O {sample}.m2.raw.v1.vcf '
#		'-new-qual true '
		'-RF OverclippedReadFilter '
		'--dont-use-soft-clipped-bases '
		'\n'
		'{gatk4v1} '
		'--java-options "-XX:ParallelGCThreads=4 -Xmx6g" '
		'FilterMutectCalls '
		'-V {sample}.m2.raw.v1.vcf '
		'-O {sample}.m2.f1.v1.vcf '
#		'--min-pcr-slippage-size 6'
		'--min-slippage-length 6 '
#		'-new-qual '
		'-R {reference} '
		'\n'
		'{gatk4v1} '
		'--java-options "-XX:ParallelGCThreads=4 -Xmx6g" '
		'FilterAlignmentArtifacts '
		'-R {reference} '
		'-I {sample}.recal.bam '
		'--bwa-mem-index-image {reference}.img '
		'--dont-skip-filtered-variants '
		'-V {sample}.m2.f1.v1.vcf '
		'-O {sample}.m2.f2.v1.vcf '
		'\n'
		'{gatk4v1} '
		'--java-options "-XX:ParallelGCThreads=4 -Xmx6g" '
		'CollectSequencingArtifactMetrics '
		'-I {sample}.recal.bam '
		'-O {sample}_artifact_v1 '
		'--FILE_EXTENSION ".v1.txt" '
		'-R {reference} '
		'\n'
		'{gatk4v1} '
		'--java-options "-XX:ParallelGCThreads=4 -Xmx6g" '
		'FilterByOrientationBias '
		'-AM G/T -AM C/T '
		'-P {sample}_artifact_v1.pre_adapter_detail_metrics.v1.txt '
		'-V {sample}.m2.f2.v1.vcf '
		'-O {sample}.m2.f3.v1.vcf '
		'\n'
		'{vt} '
		'decompose '
		'-s -o {sample}.m2.f3.decompose.v1.vcf '
		'{sample}.m2.f3.v1.vcf '
		'\n'
		'{vt} '
		'normalize '
		'-r {reference} '
		'-o {sample}.m2.f3.decompose.norm.v1.vcf '
		'{sample}.m2.f3.decompose.v1.vcf '
		.format(
				sample=sample,
				gatk4v1=conf_dict['gatk4v1'],
				reference=conf_dict['ref'],
				vt=conf_dict['vt']
				)
		)
	logger.info(cmd2)
	os.system(cmd2)
	return True


def germline_calling(sample, target_bed):
	"""Germline variants calling: Haplotypecaller of GATK4v1, and deepvariant"""
	logger.info('---SomaticVariantsCalling---')
	cmd = (
		'{gatk4v1} '
		'--java-options "-XX:ParallelGCThreads=4 -Xmx6g" '
		'HaplotypeCaller '
		'-R {reference} '
		'-I {sample}.recal.bam '
		'-ERC GVCF '
		'--dbsnp {dbsnp} '
		'-o {sample}.g.vcf.gz '
		'-new-qual true '
		.format(
				gatk4v1=conf_dict['gatk4v1'],
				reference=conf_dict['ref'],
				dbsnp=conf_dict['dbsnp'],
				sample=sample
				)
		)
	if target_bed == '':
		pass
	else:
		cmd += (
			'-L {target_bed} '
			.format(target_bed=conf_dict[target_bed])
			)
	cmd += (
		'\n'
		'{gatk4v1} '
		'--java-options "-XX:ParallelGCThreads=4 -Xmx6g" '
		'HaplotypeCaller '
		'-R {reference} '
		'-I {sample}.recal.bam '
		'--dbsnp {dbsnp} '
		'-stand-call-conf 30 '
		'-O {sample}.raw.vcf '
		'-new-qual true '
		'--max-mnp-distance 1 '
		.format(
				gatk4v1=conf_dict['gatk4v1'],
				reference=conf_dict['ref'],
				dbsnp=conf_dict['dbsnp'],
				sample=sample
				)
		)
	if target_bed == '':
		pass
	else:
		cmd += (
			'-L {target_bed} '
			.format(target_bed=conf_dict[target_bed])
			)
	cmd += (
		'\n'
		'{vt} '
		'decompose '
		'-s -o {sample}.raw.decompose.vcf '
		'{sample}.raw.vcf '
		'\n'
		'{vt} '
		'normalize '
		'-r {reference} '
		'-o {sample}.raw.decompose.norm.vcf '
		'{sample}.raw.decompose.vcf '
		'\n'
		'{gatk4v1} '
		'--java-options "-XX:ParallelGCThreads=4 -Xmx6g" '
		'VariantFiltration '
		'-R {reference} '
		'-V {sample}.raw.decompose.vcf '
		'-O {sample}.hardfilter.vcf '
		'--filter-expression \'!vc.hasAttribute(\"DP\")\' '
		'--filter-name "noCoverage" '
		'--filter-expression "(vc.isSNP() && '
		'(vc.hasAttribute(\'ReadPosRankSum\') && ReadPosRankSum < -8.0)) || '
		'((vc.isIndel() || vc.isMixed()) && '
		'(vc.hasAttribute(\'ReadPosRankSum\') && '
		'ReadPosRankSum < -20.0)) || (vc.hasAttribute(\'QD\') && QD < 2.0) " '
		'--filter-name "badSeq" '
		'--filter-expression "(vc.isSNP() && '
		'((vc.hasAttribute(\'FS\') && FS > 60.0) || '
		'(vc.hasAttribute(\'SOR\') &&  SOR > 3.0))) || '
		'((vc.isIndel() || vc.isMixed()) && '
		'((vc.hasAttribute(\'FS\') && FS > 200.0) || '
		'(vc.hasAttribute(\'SOR\') &&  SOR > 10.0)))" '
		'--filter-name "badStrand" '
		'--filter-expression "vc.isSNP() && '
		'((vc.hasAttribute(\'MQ\') && MQ < 40.0) || '
		'(vc.hasAttribute(\'MQRankSum\') && MQRankSum < -12.5))" '
		'--filter-name "badMap" '
		'\n'
		'{gatk4v1} '
		'--java-options "-XX:ParallelGCThreads=4 -Xmx6g" '
		'SelectVariants '
		'-R {reference} '
		'-V {sample}.hardfilter.vcf '
		'-O {sample}.filtered.vcf '
		'--exclude-filtered '
		.format(
				vt=conf_dict['vt'],
				gatk4v1=conf_dict['gatk4v1'],
				reference=conf_dict['ref'],
				dbsnp=conf_dict['dbsnp'],
				sample=sample
				)
		)
	logger.info(cmd)
	os.system(cmd)
	return True


def annotation(sample):
	"""Annotate vcf files"""
	logger.info('---AnnotationVariants---')
	for vcf_type in ['m2.f3.decompose.norm.v1', 'm2.f3.decompose.norm.v0']:
		cmd = (
			'java -XX:ParallelGCThreads=4 -Xmx6g -jar {snpeff} '
			'hg19 '
			'-noNextProt -noMotif -noInteraction -csvStats '
			'{sample}_snpEff.csv '
			'-noLog '
			'{sample}.{vcf_type}.vcf > {sample}.{vcf_type}.ann.vcf '
			'\n'
			'java -XX:ParallelGCThreads=4 -Xmx6g -jar {snpsift} '
			'annotate '
			'-name EXAC_ '
			'-info AC,AC_AFR,AC_AMR,AC_Adj,AC_EAS,AC_FIN,AC_Hemi,AC_Het,AC_Hom,'
			'AC_NFE,AC_OTH,AC_SAS,AC_MALE,AC_FEMALE,AC_CONSANGUINEOUS,AC_POPMAX,'
			'AN,AN_AFR,AN_AMR,AN_Adj,AN_EAS,AN_FIN,AN_NFE,AN_OTH,AN_SAS,AN_MALE,'
			'AN_FEMALE,AN_CONSANGUINEOUS,AN_POPMAX '
			'{exac} '
			'{sample}.{vcf_type}.ann.vcf > {sample}.{vcf_type}.ann.exac.vcf '
			'\n'
			'java -XX:ParallelGCThreads=4 -Xmx6g -jar {snpsift} '
			'annotate '
			'-name GNOMAD_ '
			'-info AF,AF_AFR,AF_AMR,AF_ASJ,AF_EAS,AF_FIN,AF_NFE,AF_OTH,AF_SAS,'
			'AF_Male,AF_Female,AF_raw,AF_POPMAX,AF_AMR_Male,AF_FIN_Female,AF_FIN_Male,'
			'AF_AFR_Male,AF_SAS_Male,AF_OTH_Male,AF_NFE_Female,AF_EAS_Female,AF_EAS_Male,'
			'AF_SAS_Female,AF_AFR_Female,AF_AMR_Female,AF_ASJ_Male,AF_ASJ_Female,AF_OTH_Female,AF_NFE_Male '
			'{gnomad} '
			'{sample}.{vcf_type}.ann.exac.vcf > {sample}.{vcf_type}.ann.exac.gnomad.vcf '
			'\n'
			'ln -s {sample}.{vcf_type}.ann.exac.gnomad.vcf {sample}.{vcf_type}.A.vcf '
			.format(
					snpeff=conf_dict['snpeff'],
					snpsift=conf_dict['snpsift'],
					exac=conf_dict['exac'],
					gnomad=conf_dict['gnomad'],
					vcf_type=vcf_type,
					sample=sample
					)
			)
		logger.info(cmd)
		os.system(cmd)
	return True
		

def filteration(sample):
	"""Filtering variants"""
	logger.info('---FilteringVariants---')
	for vcf_type in ['m2.f3.decompose.norm.v1.A', 'm2.f3.decompose.norm.v0.A']:
		cmd = (
			'python {foxog_filter} '
			'-v {sample}.{vcf_type}.vcf --gatk_version gatk4 > {sample}.{vcf_type}.F.oxoG.vcf '
			'\n'
			'python {exacAdj_filter} '
			'-v {sample}.{vcf_type}.F.oxoG.vcf > {sample}.{vcf_type}.F.oxoG.exacAjd.vcf '
			'\n'
			'python {gnomad_filter} '
			'-v {sample}.{vcf_type}.F.oxoG.exacAjd.vcf > {sample}.{vcf_type}.F.oxoG.exacAjd.gnomad.vcf '
			'\n'
			'python {af_filter} '
			'-v {sample}.{vcf_type}.F.oxoG.exacAjd.gnomad.vcf > {sample}.{vcf_type}.F.oxoG.exacAjd.gnomad.af.vcf '
			'\n'
			'python {alt_filter} '
			'-v {sample}.{vcf_type}.F.oxoG.exacAjd.gnomad.af.vcf > {sample}.{vcf_type}.F.oxoG.exacAjd.gnomad.af.minAlt.vcf '
			'\n'
			'python {dp_filter} '
			'-v {sample}.{vcf_type}.F.oxoG.exacAjd.gnomad.af.minAlt.vcf > {sample}.{vcf_type}.F.oxoG.exacAjd.gnomad.af.minAlt.dp.vcf '
			'\n'
			'ln -s {sample}.{vcf_type}.F.oxoG.exacAjd.gnomad.af.minAlt.dp.vcf {sample}.{vcf_type}.F.vcf'
			'\n'
			'python {apply_filter} -v {sample}.{vcf_type}.F.vcf '
			'--filter t_lod,orientation_bias,read_position,strand_artifact,contamination,str_contraction,exac_adj0.01,gnomAD0.01,oxoG,vaf0.02_1.0,dp100 '
			'> {sample}.filtered.vcf '
			.format(
					sample=sample,
					foxog_filter=conf_dict['foxog_filter'],
					exacAdj_filter=conf_dict['exacAdj_filter'],
					gnomad_filter=conf_dict['gnomad_filter'],
					af_filter=conf_dict['af_filter'],
					alt_filter=conf_dict['alt_filter'],
					dp_filter=conf_dict['dp_filter'],
					apply_filter=conf_dict['apply_filter'],
					vcf_type=vcf_type
					)
			)
		logger.info(cmd)
		os.system(cmd)

	return True




def main():
	sample, configure_file, output_path, target_bed, variant_calling = parse_argument()
	global logger
	global log_file
	global conf_dict
	conf_dict = load_configure(configure_file)
	logger, log_file = make_logger(sample)
	logger.info('--START ANALYSIS--')
	fastq_qc(sample, output_path, 1)
	fastq_adtrim(sample, output_path)
	fastq_qc(sample, output_path, 2)
	fastq_mapping(sample, output_path)
	mark_dup(sample, output_path, target_bed)
	make_recal(sample)
	somatic_calling(sample, target_bed)
	germline_calling(sample, target_bed)
	annotation(sample)
	filteration(sample)


if __name__=='__main__':
	main()

