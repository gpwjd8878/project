#!usr/bin/python

import os
import sys

def help():
    print ''
    print '[USAGE] *.py [sample]'
    print ''


def main(s1):
        os.system('source CLINICAL_env.sh')

        print '---QC---'
#        os.system('fastqc -t 2 --nogroup -o ./ %s_1.fastq %s_2.fastq' % (s1,s1))
        os.system('fastqc -t 2 --nogroup -o ./ {sample}_1.fastq.gz {sample}_2.fastq.gz'.format(sample=s1))

        print '---trimming---'
        os.system('fastp -i {sample}_1.fastq.gz -I {sample}_2.fastq.gz -o {sample}_1.adapterTrimmed.fastq.gz -O {sample}_2.adapterTrimmed.fastq.gz -x -y -c -j {sample}_fastp.json -h {sample}_fastp.html -R "{sample} fastp" -3 -p --adapter_sequence AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC --adapter_sequence_r2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -w 12'.format(sample=s1))

        print '---RE QC---'
        os.system('fastqc -t 2 --nogroup -o ./ {sample}_1.adapterTrimmed.fastq.gz {sample}_2.adapterTrimmed.fastq.gz'.format(sample=s1))

        print '---mapping---'
        os.system('bwa mem -R "@RG\\tID:RGID\\tPL:illumina\\tLB:Target_Panel" -t 12 hg19.fa {sample}_1.adapterTrimmed.fastq.gz {sample}_2.adapterTrimmed.fastq.gz > {sample}.sam'.format(sample=s1))
        os.system('sambamba view -t 12 -S -f bam {sample}.sam | sambamba sort -t 12 -o {sample}.sorted.bam /dev/stdin'.format(sample=s1))
        os.system('sambamba index -t 12 {sample}.sorted.bam'.format(sample=s1))

        print '---Markduplication---'
        os.system('sambamba markdup -t 12 {sample}.sorted.bam {sample}.markdup.bam'.format(sample=s1))
        os.system('samtools view -@ 12 -b -q 20 -L truseq-exome-targeted-regions-manifest-v1-2.bed {sample}.markdup.bam > {sample}.HQ.target.bam'.format(sample=s1))
        os.system('sambamba index -t 12 {sample}.HQ.target.bam'.format(sample=s1))

        print '---BaseRecalibration---'
        os.system('gatk --java-options "-XX:ParallelGCThreads=4 -Xmx30g" BaseRecalibrator -R hg19.fa -I {sample}.HQ.target.bam --known-sites dbsnp_138.hg19.vcf --known-sites 1000G_phase1.snps.high_confidence.hg19.sites.vcf --known-sites Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -O {sample}.recal.table'.format(sample=s1))
        os.system('gatk --java-options "-XX:ParallelGCThreads=4 -Xmx30g" ApplyBQSR -R hg19.fa  -I {sample}.HQ.target.bam  -bqsr {sample}.recal.table  -O {sample}.recal.bam'.format(sample=s1))

        print '--Variant Calling---'

###### somatic version ######
        os.system('gatk --java-options "-XX:ParallelGCThreads=4 -Xmx30g" Mutect2 -R hg19.fa -I {sample}.recal.bam -tumor {sample} -L truseq-exome-targeted-regions-manifest-v1-2.bed -O {sample}.m2.raw.vcf -new-qual true -RF OverclippedReadFilter --dont-use-soft-clipped-bases'.format(sample=s1))
        os.system('gatk --java-options "-XX:ParallelGCThreads=4 -Xmx30g" FilterMutectCalls -V {sample}.m2.raw.vcf -O {sample}.m2.f1.vcf --min-pcr-slippage-size 6 -new-qual'.format(sample=s1))
        os.system('gatk --java-options "-XX:ParallelGCThreads=4 -Xmx30g" FilterAlignmentArtifacts -R hg19.fa -I {sample}.recal.bam --bwa-mem-index-image hg19.img --dont-skip-filtered-variants -V {sample}.m2.f1.vcf -O {sample}.m2.f2.vcf'.format(sample=s1))
        os.system('gatk --java-options "-XX:ParallelGCThreads=4 -Xmx30g" CollectSequencingArtifactMetrics -I {sample}.recal.bam -O {sample}_artifact --FILE_EXTENSION ".txt" -R hg19.fa'.format(sample=s1))
        os.system('gatk --java-options "-XX:ParallelGCThreads=4 -Xmx30g" FilterByOrientationBias -AM G/T -AM C/T -P {sample}_artifact.pre_adapter_detail_metrics.txt -V {sample}.m2.f2.vcf -O {sample}.m2.f3.vcf'.format(sample=s1))
        os.system('vt decompose -s -o {sample}.m2.f3.decompose.vcf {sample}.m2.f3.vcf'.format(sample=s1))
        os.system('vt normalize -r hg19.fa -o {sample}.m2.f3.decompose.norm.vcf {sample}.m2.f3.decompose.vcf'.format(sample=s1))

#      print '--VCF annotation, filtering--'
#        os.system('java -XX:ParallelGCThreads=4 -Xmx30g -jar snpEff.jar -c snpEff.config -v hg19 {sample}.vcf > {smaple}.snpEff.vcf')format.sample=s1))


###### germline ######
#        os.system('gatk --java-options "-XX:ParallelGCThreads=4 -Xmx30g" HaplotypeCaller -R hg19.fa -I {sample}.recal.bam -ERC GVCF -O {sample}.g.vcf.gz --dbsnp {dbsnp} -L {snv_bed} -new-qual true').format(sample=s1))




if __name__=='__main__':
    if len(sys.argv) < 2:
            sys.exit(help())
    main(sys.argv[1])
