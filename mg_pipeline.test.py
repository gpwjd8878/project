#!usr/bin/python

import os
import sys

def help():
    print ''
    print '[USAGE] *.py [sample]'
    print ''


def main(s1):
        print '---QC---'
        os.system('fastqc -t 2 --nogroup -o ./ {sample}_1.fastq {sample}_2.fastq'.format(sample=s1))

        print '---trimming---'
        os.system('fastp -i {sample}_1.fastq -I {sample}_2.fastq -o {sample}_1.adapterTrimmed.fastq.gz -O {sample}_2.adapterTrimmed.fastq.gz -x -y -c -j {sample}_fastp.json -h {sample}_fastp.html -R "{sample} fastp" -3 -p --adapter_sequence AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC --adapter_sequence_r2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -w 12'.format(sample=s1))

        print '---RE QC---'
        os.system('fastqc -t 2 --nogroup -o ./ {sample}_1.adapterTrimmed.fastq.gz {sample}_2.adapterTrimmed.fastq.gz'.format(sample=s1))

        print '---mapping---'
        os.system('bwa mem -R "@RG\\tID:RGID\\tPL:illumina\\tLB:Target_Panel" -t 12 hg19.fa {sample}_1.adapterTrimmed.fastq.gz {sample}_2.adapterTrimmed.fastq.gz > {sample}.sam'.format(sample=s1))
        os.system('sambamba view -t 12 -S -f bam {sample}.sam | sambamba sort -t 12 -o {sample}.sorted.bam /dev/stdin'.format(sample=s1))
        os.system('sambamba index -t 12 {sample}.sorted.bam'.format(sample=s1))

        print '---Markduplication---'
        os.system('sambamba markdup -t 12 {sample}.sorted.bam {sample}.markdup.bam'.format(sample=s1))
        os.system('samtools view -@ 12 -b -q 20 {sample}.markdup.bam > {sample}.HQ.target.bam'.format(sample=s1))
        os.system('sambamba index -t 12 {sample}.HQ.target.bam'.format(sample=s1))

        print '---BaseRecalibration---'
        os.system('gatk --java-options "-XX:ParallelGCThreads=4 -Xmx30g" BaseRecalibrator -R hg19.fa -I {sample}.HQ.target.bam --known-sites dbsnp_138.hg19.vcf --known-sites 1000G_phase1.snps.high_confidence.hg19.sites.vcf --known-sites Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -O {sample}.recal.table'.format(sample=s1))
        os.system('gatk --java-options "-XX:ParallelGCThreads=4 -Xmx30g" ApplyBQSR -R hg19.fa  -I {sample}.HQ.target.bam  -bqsr {sample}.recal.table  -O {sample}.recal.bam'.format(sample=s1))

#        print '--Variant Calling---'
#        os.system('gatk --java-options "-XX:ParallelGCThreads=4 -Xmx30g" Mutect2 -R hg19.fa -I {sample}.recal.bam -tumor {sample} -O {sample}.m2.raw.vcf -new-qual true -RF OverclippedReadFilter --dont-use-soft-clipped-bases').format(sample=s1))
        
if __name__=='__main__':
    if len(sys.argv) < 2:
            sys.exit(help())
    main(sys.argv[1])
