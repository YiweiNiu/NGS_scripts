#!/usr/bin/env python

'''
Usage:
python xxx.py rnaseq --species human --dataset GSE96800 --sample SRR5511204 --queue Fat --ppn 8
python xxx.py chipseq --species human --dataset GSE96800 --sample SRR5359608 --queue Fat --ppn 8
python xxx.py atacseq --species human --dataset GSE96800 --sample SRR5359592 --queue Fat --ppn 8
'''

import sys
import os
import argparse

VERSION = '20190116'

# https://deeptools.readthedocs.io/en/latest/content/feature/effectiveGenomeSize.html
effective_genome_size = {'human': 2913022398, 'mouse': 2652783500}

def run():
	"""The main function
	"""

	argparser = prepare_argparser()
	args = argparser.parse_args()

	subcommand  = args.subcommand_name

	if subcommand == "atacseq":
		atacseq_preprocess(args)

	elif subcommand == "rnaseq":
		rnaseq_preprocess(args)

	elif subcommand == "chipseq":
		chipseq_preprocess(args)


def prepare_argparser():
	"""Prepare optparser object. New options will be added in this function first.
	"""

	description = "%(prog)s -- submit jobs to the clusters in batch."
	epilog = "For command line options of each command, type: %(prog)s COMMAND -h"

	argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description=description, epilog=epilog) #, usage = usage )
	argparser.add_argument("--version", action="version", version="%(prog)s "+ VERSION)
	subparsers = argparser.add_subparsers(dest='subcommand_name') #help="sub-command help")

	# command for 'ATAC-seq'
	add_atacseq_parser(subparsers)

	# command for 'RNA-seq'
	add_rnaseq_parser(subparsers)

	# command for 'ChIP-seq'
	add_chipseq_parser(subparsers)

	return argparser


def add_common_parser(parser):
	"""Add common argument parsers.
	"""

	# dataset and sample
	group_data = parser.add_argument_group( "Input files arguments")
	group_data.add_argument('--species', dest="species", required=True, choices=['human', 'mouse'], help='Species of the sample (for now only human and mouse are supported).')
	group_data.add_argument('--dataset', dest="dataset", required=True, help='Dataset identifier.')
	group_data.add_argument('--sample', dest="sample", required=True, help='Sample identifier.')
	group_data.add_argument('--single', dest="single", action='store_true', help='The sequencing data is single-end.')

	# resources
	group_resources = parser.add_argument_group( "Cluster resoures arguments")
	group_resources.add_argument('--queue', dest="queue", default='Fat', choices=['Fat', 'Blade', 'Test'], help='The queue to use.')
	group_resources.add_argument('--node', dest="node", default=1, help='Node number or node name.')
	group_resources.add_argument('--ppn', dest="ppn", default=2, type=int, help='Required threads.')

	# DIRs
	group_directory = parser.add_argument_group( "Directory arguments")
	group_directory.add_argument('--workdir', dest="workdir", default='/home/niuyw/Project/Data_processing', help='Working directory.')
	group_directory.add_argument('--logdir', dest="logdir", default='/home/niuyw/Project/Data_processing/logs', help='The directory for logs.')
	group_directory.add_argument('--pbsdir', dest="pbsdir", default='/home/niuyw/Project/Data_processing/jobfiles', help='The directory for pbs files.')

	# whether submit the job?
	group_other = parser.add_argument_group( "Other arguments")
	group_other.add_argument('--submit', dest="submit", action='store_true', help='Submit the job when the script is created.')

	# overwrite the existing pbs file?
	group_other.add_argument('--overwrite', dest="overwrite", action='store_true', help='Overwrite existing pbs file.')


def add_atacseq_parser(subparsers):
	"""Add main function 'atacseq' argument parsers.
	"""
	argparser_atacseq = subparsers.add_parser("atacseq", help="Main Function: preprocess ATAC-seq data.")

	add_common_parser(argparser_atacseq)


def add_rnaseq_parser(subparsers):
	"""Add main function 'rnaseq' argument parsers.
	"""
	argparser_rnaseq = subparsers.add_parser("rnaseq", help="Main Function: preprocess RNA-seq data.")

	add_common_parser(argparser_rnaseq)

	return


def add_chipseq_parser(subparsers):
	"""Add main function 'rnaseq' argument parsers.
	"""
	argparser_chipseq = subparsers.add_parser("chipseq", help="Main Function: preprocess ChIP-seq data.")

	add_common_parser(argparser_chipseq)

	return


def pbs_header(pbs_file_h, sample, queue, node, ppn, logdir):
	"""pbs header
	"""

	str2print = """
	#!/bin/bash
	#PBS -V
	#PBS -j eo
	#PBS -N %s
	#PBS -q %s
	#PBS -l nodes=%s:ppn=%s
	#PBS -d %s

	""" %(sample, queue, node, ppn, logdir)

	pbs_file_h.write(str2print)


def sample_infor(pbs_file_h, workdir, ppn, dataset, sample):
	"""sample information
	"""

	str2print = """
	# WORK DIR
	WORKDIR=%s

	# PPN
	PPN=%s

	# dataset and sample
	dataset=%s
	sample=%s
	""" %(workdir, ppn, dataset, sample)

	pbs_file_h.write(str2print)


def anno_tools(pbs_file_h):
	"""annotations and tools
	"""

	str2print = r"""
	echo Start time is `date +%Y/%m/%d--%H:%M`

	# genome and annotations
	ANNODIR=/home/niuyw/RefData
	GRCh38_REFERENCE=$ANNODIR/Homo_sapiens/GRCh38_no_alt/genome.fa
	GRCh38_GTF=$ANNODIR/Homo_sapiens/GENCODE_v27/gencode.v27.annotation.gtf
	GRCh38_bowtie2=$ANNODIR/Homo_sapiens/GRCh38_no_alt/Bowtie2Index/GRCh38_no_alt
	GRCh38_STAR=$ANNODIR/Homo_sapiens/GRCh38_no_alt/STARgenomes
	mm10_REFERENCE=$ANNODIR/Mus_musculus/Bowtie2Index/mm10.fa
	mm10_GTF=$ANNODIR/Mus_musculus/GENCODE_vM18/gencode.vM18.annotation.gtf
	mm10_bowtie2=$ANNODIR/Mus_musculus/Bowtie2Index/mm10
	mm10_STAR=$ANNODIR/Mus_musculus/STARgenomes

	# tools dir
	TOOLDIR=/home/niuyw/software
	path2java=$TOOLDIR/jre1.8.0_111/bin/java
	path2samtools=$TOOLDIR/samtools.1.9/bin/samtools
	path2picard=/home/work01/tools/picard-2.9.2/picard.jar
	"""

	pbs_file_h.write(str2print)


def sra_to_fastq(pbs_file_h):
	"""sra to fastq
	"""

	str2print = r"""
	# file check
	if [ ! -s $WORKDIR/rawsra/$dataset/$sample ]; then echo "$WORKDIR/rawsra/$dataset/$sample not exist or empty!"; exit 1; fi

	# sra to fastq
	if [ ! -d $WORKDIR/rawfastq/$dataset/$sample ]; then mkdir -p $WORKDIR/rawfastq/$dataset/$sample; fi
	$TOOLDIR/sratoolkit.2.8.1-ubuntu64/bin/fastq-dump --split-3 --gzip $WORKDIR/rawsra/$dataset/$sample -O $WORKDIR/rawfastq/$dataset/$sample
	$TOOLDIR/FastQC/fastqc -t $PPN -q $WORKDIR/rawfastq/$dataset/$sample/${sample}*.fastq.gz -o $WORKDIR/rawfastq/$dataset/$sample
	"""

	pbs_file_h.write(str2print)


def TrimGalore(pbs_file_h, single):
	"""TrimGalore
	"""

	str2print = r"""
	# file check
	if [ ! -s $WORKDIR/rawfastq/$dataset/$sample/${sample}_1.fastq.gz ]; then echo "$WORKDIR/rawfastq/$dataset/$sample/${sample}_1.fastq.gz not exist or empty!"; exit 1; fi
	if [ ! -s $WORKDIR/rawfastq/$dataset/$sample/${sample}_2.fastq.gz ]; then echo "$WORKDIR/rawfastq/$dataset/$sample/${sample}_2.fastq.gz not exist or empty!"; exit 1; fi

	# TrimGalore
	if [ ! -d $WORKDIR/TrimGalore/$dataset}/$sample ]; then mkdir -p $WORKDIR/TrimGalore/$dataset/$sample; fi
	$TOOLDIR/TrimGalore-0.4.5/trim_galore --paired $WORKDIR/rawfastq/$dataset/$sample/${sample}_1.fastq.gz $WORKDIR/rawfastq/$dataset/$sample/${sample}_2.fastq.gz --length 36 --trim-n -o $WORKDIR/TrimGalore/$dataset/$sample
	$TOOLDIR/FastQC/fastqc -t $PPN -q $WORKDIR/TrimGalore/$dataset/$sample/${sample}*.fq.gz -o $WORKDIR/TrimGalore/$dataset/$sample
	"""

	pbs_file_h.write(str2print)


def peak_align(pbs_file_h, bowtie2_index, genome_size, single):
	"""alignment, filter
	"""

	str2print = """
	# file check
	if [ ! -s $WORKDIR/TrimGalore/$dataset/$sample/${{sample}}*_1.fq.gz ]; then echo "$WORKDIR/TrimGalore/$dataset/$sample/${{sample}}*_1.fq.gz not exist or empty!"; exit 1; fi
	if [ ! -s $WORKDIR/TrimGalore/$dataset/$sample/${{sample}}*_2.fq.gz ]; then echo "$WORKDIR/TrimGalore/$dataset/$sample/${{sample}}*_2.fq.gz not exist or empty!"; exit 1; fi

	# bowtie2
	if [ ! -d $WORKDIR/bowtie2/$dataset/$sample ]; then mkdir -p $WORKDIR/bowtie2/$dataset/$sample; fi

	/home/software/bowtie2-2.2.8/bowtie2 --very-sensitive -X 2000 -x ${bowtie2_index} -1 $WORKDIR/TrimGalore/$dataset/$sample/${{sample}}*_1.fq.gz -2 $WORKDIR/TrimGalore/$dataset/$sample/${{sample}}*_2.fq.gz -p $PPN 2> $WORKDIR/bowtie2/$dataset/$sample/${{sample}}.bowtie2.log | $path2samtools sort -@ $PPN -O bam -o $WORKDIR/bowtie2/$dataset/$sample/${{sample}}.sorted.bam
	$path2samtools index -@ $PPN $WORKDIR/bowtie2/$dataset/$sample/${{sample}}.sorted.bam

	# mark duplicates
	$path2java -XX:ParallelGCThreads=$PPN -Djava.io.tmpdir=/tmp -jar $path2picard MarkDuplicates QUIET=true INPUT=$WORKDIR/bowtie2/$dataset/$sample/${{sample}}.sorted.bam OUTPUT=$WORKDIR/bowtie2/$dataset/$sample/${{sample}}.sorted.marked.bam METRICS_FILE=$WORKDIR/bowtie2/$dataset/$sample/${{sample}}.sorted.marked.metrics REMOVE_DUPLICATES=false CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT TMP_DIR=$WORKDIR/tmp
	$path2samtools flagstat -@ $PPN $WORKDIR/bowtie2/$dataset/$sample/${{sample}}.sorted.marked.bam > $WORKDIR/bowtie2/$dataset/$sample/${{sample}}.sorted.marked.bam.flagstat
	rm $WORKDIR/bowtie2/$dataset/$sample/${{sample}}.sorted.bam $WORKDIR/bowtie2/$dataset/$sample/${{sample}}.sorted.bam.bai

	# insert size
	$path2java -XX:ParallelGCThreads=$PPN -Djava.io.tmpdir=/tmp -jar $path2picard CollectInsertSizeMetrics QUIET=true I=$WORKDIR/bowtie2/$dataset/$sample/${{sample}}.sorted.marked.bam O=$WORKDIR/bowtie2/$dataset/$sample/${{sample}}.insert_size_metrics.txt H=$WORKDIR/bowtie2/$dataset/$sample/${{sample}}.insert_size_histogram.pdf

	# fragment size
	$TOOLDIR/anaconda2/bin/bamPEFragmentSize --numberOfProcessors $PPN --maxFragmentLength 2000 --bamfiles $WORKDIR/bowtie2/$dataset/$sample/${{sample}}.sorted.marked.bam --histogram $WORKDIR/bowtie2/$dataset/$sample/${{sample}}.fragment_size_histogram.pdf --samplesLabel $sample --plotTitle "Fragment size distribution" 1> $WORKDIR/bowtie2/$dataset/$sample/${{sample}}.bamPEFragmentSize

	# mtDNA
	mtReads=$($path2samtools idxstats $WORKDIR/bowtie2/$dataset/$sample/${{sample}}.sorted.marked.bam | grep 'chrM' | cut -f 3)
	totalReads=$($path2samtools idxstats $WORKDIR/bowtie2/$dataset/$sample/${{sample}}.sorted.marked.bam | awk '{{SUM += $3}} END {{print SUM}}')
	echo '==> mtDNA Content:' $(bc <<< "scale=2;100*$mtReads/$totalReads")'%' > $WORKDIR/bowtie2/$dataset/$sample/${{sample}}.mtDNA

	# post-alignment filtering
	# -q 30 to remove non-unique alignment
	# -F 1804 to remove: reads unmapped, mate unmapped, not primary alignment, read fails platform/vendor quality checks, read is PCR or optical duplicate
	# grep -v remove chrM
	$path2samtools view -@ $PPN -h -q 30 -F 1804 $WORKDIR/bowtie2/$dataset/$sample/${{sample}}.sorted.marked.bam | grep -v chrM | $path2samtools sort -@ $PPN -O bam -o $WORKDIR/bowtie2/$dataset/$sample/${{sample}}.sorted.marked.filtered.bam
	$path2samtools index -@ $PPN $WORKDIR/bowtie2/$dataset/$sample/${{sample}}.sorted.marked.filtered.bam
	$path2samtools flagstat -@ $PPN $WORKDIR/bowtie2/$dataset/$sample/${{sample}}.sorted.marked.filtered.bam > $WORKDIR/bowtie2/$dataset/$sample/${{sample}}.sorted.marked.filtered.bam.flagstat

	# bam to bigwig, normalize using 1x effective genome size
	# effective genome size: https://deeptools.readthedocs.io/en/latest/content/feature/effectiveGenomeSize.html
	$TOOLDIR/anaconda2/bin/bamCoverage --numberOfProcessors $PPN --binSize 10 --normalizeUsing RPGC --effectiveGenomeSize {genome_size} --bam $WORKDIR/bowtie2/$dataset/$sample/${{sample}}.sorted.marked.filtered.bam -o $WORKDIR/bowtie2/$dataset/$sample/${{sample}}.sorted.marked.filtered.bw
	"""

	str2print = str2print.format(bowtie2_index=bowtie2_index, genome_size=genome_size)
	pbs_file_h.write(str2print)


def atacseq_shift(pbs_file_h):
	"""shift reads for atac-seq
	"""

	str2print = r"""
	# shift the reads
	$TOOLDIR/anaconda2/bin/alignmentSieve --numberOfProcessors $PPN --ATACshift --label $sample --bam $WORKDIR/bowtie2/$dataset/$sample/${sample}.sorted.marked.filtered.bam -o - | $path2samtools sort -@ $PPN -O bam -o $WORKDIR/bowtie2/$dataset/$sample/${sample}.tmp.bam
	$path2samtools sort -@ $PPN -O bam -o $WORKDIR/bowtie2/$dataset/$sample/${sample}.sorted.marked.filtered.shifted.bam $WORKDIR/bowtie2/$dataset/${sample}/${sample}.tmp.bam
	$path2samtools index -@ $PPN $WORKDIR/bowtie2/$dataset/$sample/${sample}.sorted.marked.filtered.shifted.bam
	rm $WORKDIR/bowtie2/$dataset/$sample/${sample}.tmp.bam
	"""

	pbs_file_h.write(str2print)


def rnaseq_align(pbs_file_h, STARindex, GTF, single):
	"""alignment for RNA-seq
	"""
	str2print = """
	# file check
	if [ ! -s $WORKDIR/TrimGalore/$dataset/$sample/${{sample}}*_1.fq.gz ]; then echo "$WORKDIR/TrimGalore/$dataset/$sample/${{sample}}*_1.fq.gz not exist or empty!"; exit 1; fi
	if [ ! -s $WORKDIR/TrimGalore/$dataset/$sample/${{sample}}*_2.fq.gz ]; then echo "$WORKDIR/TrimGalore/$dataset/$sample/${{sample}}*_2.fq.gz not exist or empty!"; exit 1; fi

	# STAR
	if [ ! -d $WORKDIR/STAR/$dataset/$sample ]; then mkdir -p $WORKDIR/STAR/$dataset/$sample; fi
	$TOOLDIR/STAR-2.5.3a/bin/Linux_x86_64_static/STAR --runMode alignReads --quantMode GeneCounts --genomeDir ${STARindex} --sjdbGTFfile ${GTF} --outSAMtype BAM SortedByCoordinate --outSAMattributes All --outSAMunmapped Within --twopassMode Basic --readFilesIn $WORKDIR/TrimGalore/$dataset/$sample/${{sample}}*_1.fq.gz $WORKDIR/TrimGalore/$dataset/$sample/${{sample}}*_2.fq.gz --readFilesCommand zcat --outFileNamePrefix "$WORKDIR/STAR/$dataset/$sample/" --runThreadN $PPN
	$path2samtools index -@ $PPN $WORKDIR/STAR/$dataset/$sample/Aligned.sortedByCoord.out.bam
	$path2samtools flagstat -@ $PPN $WORKDIR/STAR/$dataset/$sample/Aligned.sortedByCoord.out.bam > $WORKDIR/STAR/$dataset/$sample/Aligned.sortedByCoord.out.bam.flagstat

	# insert size
	$path2java -XX:ParallelGCThreads=$PPN -Djava.io.tmpdir=/tmp -jar $path2picard CollectInsertSizeMetrics QUIET=true I=$WORKDIR/STAR/$dataset/$sample/Aligned.sortedByCoord.out.bam O=$WORKDIR/STAR/$dataset/$sample/${{sample}}.insert_size_metrics.txt H=$WORKDIR/STAR/$dataset/$sample/${{sample}}.insert_size_histogram.pdf

	# fragment size
	$TOOLDIR/anaconda2/bin/bamPEFragmentSize --numberOfProcessors $PPN --maxFragmentLength 2000 --bamfiles $WORKDIR/STAR/$dataset/$sample/Aligned.sortedByCoord.out.bam --histogram $WORKDIR/STAR/$dataset/$sample/${{sample}}.fragment_size_histogram.pdf --samplesLabel $sample --plotTitle "Fragment size distribution" 1> $WORKDIR/STAR/$dataset/$sample/${{sample}}.bamPEFragmentSize
	"""

	str2print = str2print.format(STARindex=STARindex, GTF=GTF)
	pbs_file_h.write(str2print)


def pbs_footer(pbs_file_h):
	"""pbs header
	"""
	str2print = r"""
	echo Finish time is `date +%Y/%m/%d--%H:%M`
	"""

	pbs_file_h.write(str2print)
	pbs_file_h.close()


def atacseq_preprocess(args):
	"""preprocess the ATAC-seq data, from sra to bam
	"""
	species = args.species
	sample = args.sample
	dataset = args.dataset
	single = args.single

	queue = args.queue
	node = args.node
	ppn = args.ppn

	workdir = args.workdir
	logdir = args.logdir
	pbsdir = args.pbsdir

	submit = args.submit
	overwrite = args.overwrite

	pbs_file = pbsdir + '/' + sample + '.pbs'

	if os.path.exists(pbs_file):
		if overwrite:
			try:
				pbs_file_h = open(pbs_file, 'w')
			except IOError as err:
				print 'Open %s error: ' %(pbs_file) + str(err)
				sys.exit(1)
		else:
			print 'PBS file exists. Use --overwrite to overwrite it.'
			sys.exit(1)
	else:
		try:
			pbs_file_h = open(pbs_file, 'w')
		except IOError as err:
			print 'Open %s error: ' %(pbs_file) + str(err)
			sys.exit(1)

	if species == 'human':
		bowtie2_index = 'GRCh38_bowtie2'
		genome_size = effective_genome_size['human']
	elif species == 'mouse':
		bowtie2_index = 'mm10_bowtie2'
		genome_size = effective_genome_size['mouse']
	else:
		print 'Unsupported species.'
		sys.exit(1)

	pbs_header(pbs_file_h, sample, queue, node, ppn, logdir)
	sample_infor(pbs_file_h, workdir, ppn, dataset, sample)
	anno_tools(pbs_file_h)
	sra_to_fastq(pbs_file_h)
	TrimGalore(pbs_file_h, single)
	peak_align(pbs_file_h, bowtie2_index, genome_size, single)
	atacseq_shift(pbs_file_h)
	pbs_footer(pbs_file_h)

	if submit:
		os.system('qsub %s' %(pbs_file))


def chipseq_preprocess(args):
	"""preprocess the ChIP-seq data, from sra to bam
	"""
	species = args.species
	sample = args.sample
	dataset = args.dataset
	single = args.single

	queue = args.queue
	node = args.node
	ppn = args.ppn

	workdir = args.workdir
	logdir = args.logdir
	pbsdir = args.pbsdir

	submit = args.submit
	overwrite = args.overwrite

	pbs_file = pbsdir + '/' + sample + '.pbs'

	if os.path.exists(pbs_file):
		if overwrite:
			try:
				pbs_file_h = open(pbs_file, 'w')
			except IOError as err:
				print 'Open %s error: ' %(pbs_file) + str(err)
				sys.exit(1)
		else:
			print 'PBS file exists. Use --overwrite to overwrite it.'
			sys.exit(1)
	else:
		try:
			pbs_file_h = open(pbs_file, 'w')
		except IOError as err:
			print 'Open %s error: ' %(pbs_file) + str(err)
			sys.exit(1)

	if species == 'human':
		bowtie2_index = 'GRCh38_bowtie2'
		genome_size = effective_genome_size['human']
	elif species == 'mouse':
		bowtie2_index = 'mm10_bowtie2'
		genome_size = effective_genome_size['mouse']
	else:
		print 'Unsupported species.'
		sys.exit(1)

	pbs_header(pbs_file_h, sample, queue, node, ppn, logdir)
	sample_infor(pbs_file_h, workdir, ppn, dataset, sample)
	anno_tools(pbs_file_h)
	sra_to_fastq(pbs_file_h)
	TrimGalore(pbs_file_h, single)
	peak_align(pbs_file_h, bowtie2_index, genome_size, single)
	pbs_footer(pbs_file_h)

	if submit:
		os.system('qsub %s' %(pbs_file))


def rnaseq_preprocess(args):
	"""preprocess the RNA-seq data, from sra to bam
	"""
	species = args.species
	sample = args.sample
	dataset = args.dataset
	single = args.single

	queue = args.queue
	node = args.node
	ppn = args.ppn

	workdir = args.workdir
	logdir = args.logdir
	pbsdir = args.pbsdir

	submit = args.submit
	overwrite = args.overwrite

	pbs_file = pbsdir + '/' + sample + '.pbs'

	if os.path.exists(pbs_file):
		if overwrite:
			try:
				pbs_file_h = open(pbs_file, 'w')
			except IOError as err:
				print 'Open %s error: ' %(pbs_file) + str(err)
				sys.exit(1)
		else:
			print 'PBS file exists. Use --overwrite to overwrite it.'
			sys.exit(1)
	else:
		try:
			pbs_file_h = open(pbs_file, 'w')
		except IOError as err:
			print 'Open %s error: ' %(pbs_file) + str(err)
			sys.exit(1)

	if species == 'human':
		STARindex = 'GRCh38_STAR'
		REFERENCE = 'GRCh38_REFERENCE'
		GTF = 'GRCh38_GTF'
	elif species == 'mouse':
		STARindex = 'mm10_STAR'
		REFERENCE = 'mm10_REFERENCE'
		GTF = 'mm10_GTF'
	else:
		print 'Unsupported species.'
		sys.exit(1)

	pbs_header(pbs_file_h, sample, queue, node, ppn, logdir)
	sample_infor(pbs_file_h, workdir, ppn, dataset, sample)
	anno_tools(pbs_file_h)
	sra_to_fastq(pbs_file_h)
	TrimGalore(pbs_file_h, single)
	rnaseq_align(pbs_file_h, STARindex, GTF, single)
	pbs_footer(pbs_file_h)

	if submit:
		os.system('qsub %s' %(pbs_file))


if __name__ == '__main__':
	run()
