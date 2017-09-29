"""
Author: Y. Ahmed-Braimah
--- RNA-seq snakemake workflow: Use to map reads to genome without
___ pre-existing annotation (New PacBio genomes).

"""

import json
import os
from os.path import join, basename, dirname
from os import getcwd
from subprocess import check_output

##--------------------------------------------------------------------------------------##
## Functions
##--------------------------------------------------------------------------------------##

# To print process messages
def message(x):
  print()

# To remove suffix from a string
def rstrip(text, suffix):
    if not text.endswith(suffix):
        return text
    return text[:len(text)-len(suffix)]

## define environment variables

##--------------------------------------------------------------------------------------##
## Global config files: 
##--------------------------------------------------------------------------------------##

configfile: 'config.yml'

# Full path to an uncompressed FASTA file with all chromosome sequences.
DNA = config['DNA']
INDEX = config['INDEX']

# Full path to a folder where final output files will be deposited.
OUT_DIR = config['OUT_DIR']
WORK_DIR = config['WORK_DIR']
HOME_DIR = config['HOME_DIR']  # the "launch_snakemake.sh" and "config.yml" files should be here

SHARED_DATABASE = config['SHARED_DATABASE']
CUSTOM_DATABASE = config['CUSTOM_DATABASE']
TRINOTATE_HOME = config['TRINOTATE_HOME']
SQLITE = config['SQLITE']

## set the usr and job environments for each job (specific for CBSU qsub jobs)
USER = os.environ.get('USER')
JOB_ID = os.environ.get('JOB_ID')

# Samples and their corresponding filenames.
# single-end
seFILES = json.load(open(config['SE_SAMPLES_JSON'])) 
seSAMPLES = sorted(seFILES.keys())                  
# paired-end:
peFILES = json.load(open(config['PE_SAMPLES_JSON'])) 
peSAMPLES = sorted(peFILES.keys())           

# read both
# FILES = json.load(open(config['SAMPLES_JSON']))
combinedSam = [peSAMPLES, seSAMPLES]
SAMPLES = [y for x in combinedSam for y in x]  

## Create the final output directory if it doesn't already exist
if not os.path.exists(OUT_DIR):
            os.makedirs(OUT_DIR)

##--------------------------------------------------------------------------------------##
## RULES
##--------------------------------------------------------------------------------------##

## Final expected output(s)
rule all: 
    input: 
        join(OUT_DIR, 'MultiQC', 'multiqc_report.html'),
        join(OUT_DIR, 'ballgown', 'gene_counts.csv'), 
        join(OUT_DIR, 'ballgown', 'transcript_counts.csv'),
        join(OUT_DIR, 'Trinotate_report.xls'),
        join(OUT_DIR, 'Trinotate_report.xls.gene_ontology'),
        join(OUT_DIR, 'stringtie_merged.gtf' + '_' + rstrip(os.path.basename(DNA), '.fa') + '_exons.fa.transdecoder.genome.gff3')
        
##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

## Rule to check raw SE read quality
rule fastqcSE:
    input:
        r1 = lambda wildcards: seFILES[wildcards.sample]['R1']
    output:
        r1 = join(OUT_DIR, 'fastQC', '{sample}' + '.R1_fastqc.html')
    log:
        join(OUT_DIR, 'fastQC', '{sample}' + '.fastQC_se.log')
    benchmark:
        join(OUT_DIR, 'fastQC', '{sample}' + '.fastQC_se.benchmark.tsv')
    message: 
        """--- Checking read quality of SE sample "{wildcards.sample}" with FastQC """
    run:
        if not os.path.exists(join(OUT_DIR, 'fastQC')):
            os.makedirs(join(OUT_DIR, 'fastQC'))

        shell('mkdir -p ' + join(WORK_DIR, USER, JOB_ID) +
                ' && cp {input.r1} ' + join(WORK_DIR, USER, JOB_ID) +
                ' && cd ' + join(WORK_DIR, USER, JOB_ID) + 
                ' && fastqc {wildcards.sample}.R1.fq.gz' 
                ' > {log} 2>&1'
                ' && mv ' + join(WORK_DIR, USER, JOB_ID) + '/*fastqc* ' + join(OUT_DIR, 'fastQC'))
        shell('rm -r ' + join(WORK_DIR, USER, JOB_ID))

##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

## Rule to check raw PE read quality
rule fastqcPE:
    input:
        r1 = lambda wildcards: peFILES[wildcards.sample]['R1'],
        r2 = lambda wildcards: peFILES[wildcards.sample]['R2']
    output:
        r1 = join(OUT_DIR, 'fastQC', '{sample}' + '.R1_fastqc.html'),
        r2 = join(OUT_DIR, 'fastQC', '{sample}' + '.R2_fastqc.html')
    log:
        join(OUT_DIR, 'fastQC', '{sample}' + '.fastQC_init_pe.log')
    benchmark:
        join(OUT_DIR, 'fastQC', '{sample}' + '.fastQC_init_pe.benchmark.tsv')
    message: 
        """--- Checking read quality of PE sample "{wildcards.sample}" with FastQC """
    run:
        if not os.path.exists(join(OUT_DIR, 'fastQC')):
            os.makedirs(join(OUT_DIR, 'fastQC'))

        shell('mkdir -p ' + join(WORK_DIR, USER, JOB_ID) +
                ' && cp {input.r1} {input.r2} ' + join(WORK_DIR, USER, JOB_ID) +
                ' && cd ' + join(WORK_DIR, USER, JOB_ID) + 
                ' && fastqc {wildcards.sample}.R1.fq.gz {wildcards.sample}.R2.fq.gz' 
                ' > {log} 2>&1'
                ' && mv ' + join(WORK_DIR, USER, JOB_ID) + '/*fastqc* ' + join(OUT_DIR, 'fastQC'))
        shell('rm -r ' + join(WORK_DIR, USER, JOB_ID))

##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

## Rule to map PE reads with HISAT2
rule hisat2_se_mapping:
    input:
        r1 = lambda wildcards: seFILES[wildcards.sample]['R1']
    output:
        bam = join(OUT_DIR, 'HISAT-2', '{sample}', '{sample}' + '.csorted.bowtie2.bam')
    log:
        join(OUT_DIR, 'HISAT-2', '{sample}', 'hs2_map_se.log')
    benchmark:
        join(OUT_DIR, 'HISAT-2', '{sample}', 'hs2_map_se.benchmark.tsv')
    message: 
        """--- Mapping SE reads for sample {wildcards.sample} to genome with HISAT-2 """
    run:
        shell('mkdir -p ' + join(WORK_DIR, USER, JOB_ID) +
                ' && cp ' + join(DNA + '*') + ' ' + join(WORK_DIR, USER, JOB_ID) +
                ' && cp ' + join(INDEX, '*') + ' ' + join(WORK_DIR, USER, JOB_ID) +
                ' && cp {input.r1} ' + join(WORK_DIR, USER, JOB_ID) +
                ' && cd ' + join(WORK_DIR, USER, JOB_ID) + 
                ' && (hisat2'
                ' -p 16'
                ' --dta'
                ' -x ' + join(rstrip(os.path.basename(DNA), '.fa') + '_tran') +
                ' -U {wildcards.sample}.R1.fq.gz) 2>{log}'
                ' | samtools sort -@ 8 -o csorted.bowtie2.bam -')
        shell('mv ' + join(WORK_DIR, USER, JOB_ID, 'csorted.bowtie2.bam') + ' ' + join(OUT_DIR, 'HISAT-2', '{wildcards.sample}', '{wildcards.sample}' + '.csorted.bowtie2.bam'))
        shell('rm -r ' + join(WORK_DIR, USER, JOB_ID))

##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

## Rule to check raw PE read quality
rule hisat2_pe_mapping:
    input:
        r1 = lambda wildcards: peFILES[wildcards.sample]['R1'],
        r2 = lambda wildcards: peFILES[wildcards.sample]['R2']
    output:
        bam = join(OUT_DIR, 'HISAT-2', '{sample}', '{sample}' + '.csorted.bowtie2.bam')
    log:
        join(OUT_DIR, 'HISAT-2', '{sample}', 'hs2_map_pe.log')
    benchmark:
        join(OUT_DIR, 'HISAT-2', '{sample}', 'hs2_map_pe.benchmark.tsv')
    message: 
        """--- Mapping PE reads for sample {wildcards.sample} to genome with HISAT-2 """
    run:
        shell('mkdir -p ' + join(WORK_DIR, USER, JOB_ID) +
                ' && cp ' + join(DNA + '*') + ' ' + join(WORK_DIR, USER, JOB_ID) +
                ' && cp ' + join(INDEX, '*') + ' ' + join(WORK_DIR, USER, JOB_ID) +
                ' && cp {input.r1} {input.r2} ' + join(WORK_DIR, USER, JOB_ID) +
                ' && cd ' + join(WORK_DIR, USER, JOB_ID) + 
                ' && (hisat2'
                ' -p 16'
                ' --dta'
                ' -x ' + join(rstrip(os.path.basename(DNA), '.fa') + '_tran') +
                ' -1 {wildcards.sample}.R1.fq.gz'
                ' -2 {wildcards.sample}.R2.fq.gz)'
                ' 2>{log}'
                ' | samtools sort -@ 8 -o csorted.bowtie2.bam -')
        shell('mv ' + join(WORK_DIR, USER, JOB_ID, 'csorted.bowtie2.bam') + ' ' + join(OUT_DIR, 'HISAT-2', '{wildcards.sample}', '{wildcards.sample}' + '.csorted.bowtie2.bam'))
        shell('rm -r ' + join(WORK_DIR, USER, JOB_ID))

##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

## Rule to assemble transcripts with StringTie
rule stringtie_assembly:
    input:
        bam = join(OUT_DIR, 'HISAT-2', '{sample}', '{sample}' + '.csorted.bowtie2.bam')
    output:
        asmbly = join(OUT_DIR, 'StringTie', '{sample}', '{sample}' + '.gtf')
    log:
        join(OUT_DIR, 'StringTie', '{sample}', 'st_asmbly.log')
    benchmark:
        join(OUT_DIR, 'StringTie', '{sample}', 'st_asmbly.benchmark.tsv')
    message: 
        """--- Assembling transcripts for sample {wildcards.sample} with StringTie """
    run:
        shell('mkdir -p ' + join(WORK_DIR, USER, JOB_ID) +
                ' && cp {input.bam} ' + join(WORK_DIR, USER, JOB_ID) +
                ' && cd ' + join(WORK_DIR, USER, JOB_ID) + 
                ' && stringtie'
                ' {wildcards.sample}.csorted.bowtie2.bam'
                ' -p 16'
                ' -o {wildcards.sample}.gtf'
                ' -l {wildcards.sample} > {log}')
        shell('mv ' + join(WORK_DIR, USER, JOB_ID, '{wildcards.sample}.gtf') + ' ' + join(OUT_DIR, 'StringTie', '{wildcards.sample}'))
        shell('rm -r ' + join(WORK_DIR, USER, JOB_ID))

##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

## Rule to merge StringTie assemblies
rule merge_assemblies:
    input:
        assemblies = expand(join(OUT_DIR, 'StringTie', '{sample}', '{sample}' + '.gtf'), sample = SAMPLES)
    output:
        asmbly = join(OUT_DIR, 'StringTie', 'stringtie_merged.gtf')
    log:
        join(OUT_DIR, 'StringTie', 'st_mrg.index.log')
    benchmark:
        join(OUT_DIR, 'StringTie', 'st_mrg.index.benchmark.tsv')
    message: 
        """--- Merging StringTie transcripts """
    run:
        shell('mkdir -p ' + join(WORK_DIR, USER, JOB_ID) +
                ' && ls -1 ' + join(OUT_DIR) + '/StringTie/*/*.gtf > ' + join(OUT_DIR, 'StringTie', 'assemblies.txt') +
                ' && cd ' + join(WORK_DIR, USER, JOB_ID) + 
                ' && stringtie'
                ' --merge'
                ' -p 8'
                ' -o stringtie_merged.gtf ' + join(OUT_DIR, 'StringTie', 'assemblies.txt'))
        shell('mv ' + join(WORK_DIR, USER, JOB_ID, 'stringtie_merged.gtf') + ' ' + join(OUT_DIR, 'StringTie'))
        shell('rm -r ' + join(WORK_DIR, USER, JOB_ID))

##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

## Rule to measure transcript abundances with Stringtie
rule abundances:
    input:
        bam = join(OUT_DIR, 'HISAT-2', '{sample}', '{sample}' + '.csorted.bowtie2.bam'),
        mrgd = rules.merge_assemblies.output.asmbly
    output:
        abundance = join(OUT_DIR, 'ballgown', '{sample}', '{sample}' + '_abundance.gtf')
    log:
        join(OUT_DIR, 'ballgown', '{sample}', 'st_abnd.log')
    benchmark:
        join(OUT_DIR, 'ballgown', '{sample}', 'st_abnd.benchmark.tsv')
    message: 
        """--- Estimating transcript abundances for sample {wildcards.sample} with StringTie"""
    run:
        shell('mkdir -p ' + join(WORK_DIR, USER, JOB_ID) +
                ' && cp {input.bam} {input.mrgd} ' + join(WORK_DIR, USER, JOB_ID) +
                ' && cd ' + join(WORK_DIR, USER, JOB_ID) + 
                ' && stringtie'
                ' -e -B -p 8'
                ' -G stringtie_merged.gtf' 
                ' -o {wildcards.sample}_abundance.gtf'
                ' {wildcards.sample}.csorted.bowtie2.bam'
                ' > {log}')
        shell('mv ' + join(WORK_DIR, USER, JOB_ID, '{wildcards.sample}_abundance.gtf') + ' ' + join(OUT_DIR, 'ballgown', '{wildcards.sample}'))
        shell('rm -r ' + join(WORK_DIR, USER, JOB_ID))

##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

## Rule to combine abundance counts for downstream analysis
rule collate_counts:
    input:
        abundances = expand(join(OUT_DIR, 'ballgown', '{sample}', '{sample}' + '_abundance.gtf'), sample = SAMPLES)
    output:
        geneCounts = join(OUT_DIR, 'ballgown', 'gene_counts.csv'),
        transcriptCounts = join(OUT_DIR, 'ballgown', 'transcript_counts.csv')
    message: 
        """--- Outputting count matrices """
    run:
        shell('prepDE.py'
                ' -i ' + join(OUT_DIR, 'ballgown') + 
                ' -g ' + join(OUT_DIR, 'ballgown', 'gene_counts.csv') +
                ' -t ' + join(OUT_DIR, 'ballgown', 'transcript_counts.csv'))

##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

## Rule to collate fastQC and HISAT2 outputs with multiQC
rule multiQC:
    input:
        expand(join(OUT_DIR, 'HISAT-2', '{sample}', '{sample}' + '.csorted.bowtie2.bam'), sample = SAMPLES),
        expand(join(OUT_DIR, 'fastQC', '{sample}' + '.R1_fastqc.html'), sample = SAMPLES),
        expand(join(OUT_DIR, 'fastQC', '{sample}' + '.R2_fastqc.html'), sample = peSAMPLES)
    output:
        join(OUT_DIR, 'MultiQC', 'multiqc_report.html')
    log:
        join(OUT_DIR, 'MultiQC', 'multiQC.log')
    benchmark:
        join(OUT_DIR, 'MultiQC', 'multiQC.benchmark.tsv')
    message: 
        """--- Running MultiQC """
    run:
        shell('ls -1 ' + join(OUT_DIR) + '/HISAT-2/*/*log > ' + join(OUT_DIR, 'MultiQC', 'summary_files.txt'))
        shell('ls -1 ' + join(OUT_DIR) + '/fastQC/*fastqc.zip >> ' + join(OUT_DIR, 'MultiQC', 'summary_files.txt'))
        shell('multiqc'
                ' -f'
                ' -o ' + join(OUT_DIR, 'MultiQC') + ' -d -dd 2 -l ' + join(OUT_DIR, 'MultiQC', 'summary_files.txt') +
                ' > {log} 2>&1')

##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

rule gffread_exons:
    input:
        gtf = rules.merge_assemblies.output.asmbly,
        dna = DNA
    output:
        exons = join(OUT_DIR, 'sequences', 'stringtie_merged.gtf' + '_' + rstrip(os.path.basename(DNA), '.fa') + '_exons.fa'),
        gff3 = join(OUT_DIR, 'sequences', 'stringtie_merged.gtf' + '.gff3')
    message: 
        """--- Extracting exon sequences from the genome using the GTF file """
    run:
        shell('mkdir -p ' + join(WORK_DIR, USER, JOB_ID) +
                ' && cp {input.gtf} {input.dna} ' + join(WORK_DIR, USER, JOB_ID) +
                ' && cd ' + join(WORK_DIR, USER, JOB_ID) + 
                ' && gffread ' + os.path.basename(GTF) + ' -g ' + os.path.basename(DNA) + ' -w ' +
                join('stringtie_merged.gtf' + '_' + rstrip(os.path.basename(DNA), '.fa') + '_exons.fa') +
                ' && cufflinks_gtf_to_alignment_gff3.pl ' + os.path.basename(GTF) + ' > ' + 'stringtie_merged.gtf' + '.gff3'
                ' && mv ' + join(WORK_DIR, USER, JOB_ID) + '/* ' + join(OUT_DIR, 'sequences'))
        shell('rm -r ' + join(WORK_DIR, USER, JOB_ID))

##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

rule Transdecoder_LongOrfs:
    input:
        exons = rules.gffread_exons.output.exons
    output:
        longOrfs = join(OUT_DIR, 'stringtie_merged.gtf' + '_' + rstrip(os.path.basename(DNA), '.fa') + '_exons.fa.transdecoder_dir/longest_orfs.pep')
    log:
        join(OUT_DIR, 'LOGS', 'Transdecoder_LongOrfs.log')
    benchmark:
        join(OUT_DIR, 'benchmarks', 'Transdecoder_LongOrfs.benchmark.tsv')
    message: 
        """--- Extracting long ORFs with TransDecoder """
    run:
        shell('mkdir -p ' + join(WORK_DIR, USER, JOB_ID) +
                ' && cp {input.exons} ' + join(WORK_DIR, USER, JOB_ID) +
                ' && cd ' + join(WORK_DIR, USER, JOB_ID) + 
                ' && TransDecoder.LongOrfs -t ' + join('stringtie_merged.gtf' + '_' + rstrip(os.path.basename(DNA), '.fa') + '_exons.fa') + 
                ' -m 50' 
                ' > {log} 2>&1'
                ' && mv ' + join(WORK_DIR, USER, JOB_ID) + '/*transdecoder_dir/* ' + join(OUT_DIR, 'stringtie_merged.gtf' + '_' + rstrip(os.path.basename(DNA), '.fa') + '_exons.fa.transdecoder_dir/'))
        shell('rm -r ' + join(WORK_DIR, USER, JOB_ID))

##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

rule BLASTp_init:
    input:
        longOrfs = rules.Transdecoder_LongOrfs.output.longOrfs
    output:
        blastpI = join(OUT_DIR, 'BLAST_results', 'BLASTp_init.outfmt6')
    benchmark:
        join(OUT_DIR, 'benchmarks', 'BLASTp_init.benchmark.tsv')
    message: 
        """--- Initial BLASTp for TransDecoder """
    run:
        shell('mkdir -p ' + join(WORK_DIR, USER, JOB_ID) +
                ' && cp {input.longOrfs} ' + join(WORK_DIR, USER, JOB_ID) +
                ' && cp ' + join(SHARED_DATABASE, 'uniprot*') + ' ' + join(WORK_DIR, USER, JOB_ID) +
                ' && cd ' + join(WORK_DIR, USER, JOB_ID) + 
                ' && blastp -query longest_orfs.pep -db uniprot_sprot.pep -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads 10 > BLASTp_init.outfmt6'
                ' && mv ' + join(WORK_DIR, USER, JOB_ID, 'BLASTp_init.outfmt6') + ' ' + join(OUT_DIR, 'BLAST_results'))
        shell('rm -r ' + join(WORK_DIR, USER, JOB_ID))

##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

rule Pfam_init:
    input:
        longOrfs = rules.Transdecoder_LongOrfs.output.longOrfs
    output:
        pfamI = join(OUT_DIR, 'Pfam_results', 'pfam_i.domtblout')
    log:
        join(OUT_DIR, 'LOGS', 'Pfam_init.log')
    benchmark:
        join(OUT_DIR, 'benchmarks', 'Pfam_init.benchmark.tsv')
    message: 
        """--- Initial Pfam search for TransDecoder """
    run:
        shell('mkdir -p ' + join(WORK_DIR, USER, JOB_ID) +
                ' && cp {input.longOrfs} ' + join(WORK_DIR, USER, JOB_ID) +
                ' && cp ' + join(SHARED_DATABASE, 'Pfam*') + ' ' + join(WORK_DIR, USER, JOB_ID) +
                ' && cd ' + join(WORK_DIR, USER, JOB_ID) + 
                ' && hmmscan --cpu 8 --domtblout pfam_i.domtblout Pfam-A.hmm longest_orfs.pep > {log} 2>&1'
                ' && mv ' + join(WORK_DIR, USER, JOB_ID, 'pfam_i.domtblout') + ' ' + join(OUT_DIR, 'Pfam_results'))
        shell('rm -r ' + join(WORK_DIR, USER, JOB_ID))

##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

rule Transdecoder_Predict:
    input:
        exons = rules.gffread_exons.output.exons,
        blastpI = rules.BLASTp_init.output.blastpI,
        pfamI = rules.Pfam_init.output.pfamI
    output:
        TransGff3 = join(OUT_DIR, 'stringtie_merged.gtf' + '_' + rstrip(os.path.basename(DNA), '.fa') + '_exons.fa.transdecoder.gff3'),
        peptides = join(OUT_DIR, 'stringtie_merged.gtf' + '_' + rstrip(os.path.basename(DNA), '.fa') + '_exons.fa.transdecoder.pep')
    log:
        join(OUT_DIR, 'LOGS', 'Transdecoder_Predict.log')
    benchmark:
        join(OUT_DIR, 'benchmarks', 'Transdecoder_Predict.benchmark.tsv')
    message: 
        """--- Final ORF prediction """
    run:
        shell('mkdir -p ' + join(WORK_DIR, USER, JOB_ID) +
                ' && cp {input.exons} {input.blastpI} {input.pfamI} ' + join(WORK_DIR, USER, JOB_ID) +
                ' && cp -r ' + join(OUT_DIR, '*transdecoder_dir') + ' ' + join(WORK_DIR, USER, JOB_ID) +
                ' && cd ' + join(WORK_DIR, USER, JOB_ID) + 
                ' && TransDecoder.Predict'
                ' -t ' + join('stringtie_merged.gtf' + '_' + rstrip(os.path.basename(DNA), '.fa') + '_exons.fa') +
                ' --retain_pfam_hits pfam_i.domtblout --retain_blastp_hits BLASTp_init.outfmt6' 
                ' > {log} 2>&1'
                ' && mv ' + join(WORK_DIR, USER, JOB_ID, 'stringtie_merged.gtf' + '_' + rstrip(os.path.basename(DNA), '.fa') + '_exons.fa.transdecoder.*') + ' ' + OUT_DIR)
        shell('rm -r ' + join(WORK_DIR, USER, JOB_ID))

##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

rule ORF_mapping_to_genome:
    input:
        TransGff3 = rules.Transdecoder_Predict.output.TransGff3,
        exons = rules.gffread_exons.output.exons,
        gff3 = rules.gffread_exons.output.gff3
    output:
        GenomeGff3 = join(OUT_DIR, 'stringtie_merged.gtf' + '_' + rstrip(os.path.basename(DNA), '.fa') + '_exons.fa.transdecoder.genome.gff3')
    benchmark:
        join(OUT_DIR, 'benchmarks', 'ORF_mapping_to_genome.tsv')
    message: 
        """--- Generate genome ORF coordinate """
    run:
        shell('mkdir -p ' + join(WORK_DIR, USER, JOB_ID) +
                ' && cp {input.TransGff3} {input.exons} {input.gff3} ' + join(WORK_DIR, USER, JOB_ID) +
                ' && cd ' + join(WORK_DIR, USER, JOB_ID) + 
                ' && cdna_alignment_orf_to_genome_orf.pl ' + 
                join(OUT_DIR, 'stringtie_merged.gtf' + '_' + rstrip(os.path.basename(DNA), '.fa') + '_exons.fa.transdecoder.gff3') + ' ' +
                'stringtie_merged.gtf' + '.gff3 ' + 
                'stringtie_merged.gtf' + '_' + rstrip(os.path.basename(DNA), '.fa') + '_exons.fa'
                ' > ' + 'stringtie_merged.gtf' + '_' + rstrip(os.path.basename(DNA), '.fa') + '_exons.fa.transdecoder.genome.gff3'
                ' && mv ' + join(WORK_DIR, USER, JOB_ID) + '/*genome.gff3 ' + OUT_DIR)
        shell('rm -r ' + join(WORK_DIR, USER, JOB_ID))

##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

rule BLASTx:
    input:
        exons = rules.gffread_exons.output.exons
    output:
        blastX = join(OUT_DIR, 'BLAST_results', 'swissprot.blastx.outfmt6')
    benchmark:
        join(OUT_DIR, 'benchmarks', 'BLASTx.tsv')
    message: 
        """--- Transcript search against SwissProt (BLASTx)"""
    run:
        shell('mkdir -p ' + join(WORK_DIR, USER, JOB_ID) +
                ' && cp {input.exons} ' + join(WORK_DIR, USER, JOB_ID) +
                ' && cp ' + join(SHARED_DATABASE, 'uniprot*') + ' ' + join(WORK_DIR, USER, JOB_ID) +
                ' && cd ' + join(WORK_DIR, USER, JOB_ID) + 
                ' && blastx'
                ' -query ' + 'stringtie_merged.gtf' + '_' + rstrip(os.path.basename(DNA), '.fa') + '_exons.fa'
                ' -db uniprot_sprot.pep -num_threads 16 -max_target_seqs 1 -outfmt 6 > swissprot.blastx.outfmt6'
                ' && mv swissprot.blastx.outfmt6 ' + join(OUT_DIR, 'BLAST_results'))
        shell('rm -r ' + join(WORK_DIR, USER, JOB_ID))

##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

rule BLASTp:
    input:
        peptides = rules.Transdecoder_Predict.output.peptides
    output:
        blastP = join(OUT_DIR, 'BLAST_results', 'swissprot.blastp.outfmt6')
    benchmark:
        join(OUT_DIR, 'benchmarks', 'BLASTp.tsv')
    message: 
        """--- Peptide search against SwissProt (BLASTp)"""
    run:
        shell('mkdir -p ' + join(WORK_DIR, USER, JOB_ID) +
                ' && cp {input.peptides} ' + join(WORK_DIR, USER, JOB_ID) +
                ' && cp ' + join(SHARED_DATABASE, 'uniprot*') + ' ' + join(WORK_DIR, USER, JOB_ID) +
                ' && cd ' + join(WORK_DIR, USER, JOB_ID) + 
                ' && blastp'
                ' -query ' + 'stringtie_merged.gtf' + '_' + rstrip(os.path.basename(DNA), '.fa') + '_exons.fa'+ '.transdecoder.pep'
                ' -db uniprot_sprot.pep -num_threads 16 -max_target_seqs 1 -outfmt 6 > swissprot.blastp.outfmt6'
                ' && mv swissprot.blastp.outfmt6 ' + join(OUT_DIR, 'BLAST_results'))
        shell('rm -r ' + join(WORK_DIR, USER, JOB_ID))

##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

rule custom_BLASTx:
    input:
        exons = rules.gffread_exons.output.exons
    output:
        blastX = join(OUT_DIR, 'BLAST_results', 'custom.blastx.outfmt6')
    benchmark:
        join(OUT_DIR, 'benchmarks', 'custom_BLASTx.tsv')
    message: 
        """--- Transcript search against Custom database (BLASTx)"""
    run:
        shell('mkdir -p ' + join(WORK_DIR, USER, JOB_ID) +
                ' && cp {input.exons} ' + join(WORK_DIR, USER, JOB_ID) +
                ' && cp ' + CUSTOM_DATABASE + ' ' + join(WORK_DIR, USER, JOB_ID) +
                ' && cd ' + join(WORK_DIR, USER, JOB_ID) + 
                ' && makeblastdb -in ' + os.path.basename(CUSTOM_DATABASE) +
                ' -dbtype prot'
                ' && blastx'
                ' -query ' + 'stringtie_merged.gtf' + '_' + rstrip(os.path.basename(DNA), '.fa') + '_exons.fa'
                ' -db ' + os.path.basename(CUSTOM_DATABASE) + ' -num_threads 16 -max_target_seqs 1 -outfmt 6 > custom.blastx.outfmt6'
                ' && mv custom.blastx.outfmt6 ' + join(OUT_DIR, 'BLAST_results'))
        shell('rm -r ' + join(WORK_DIR, USER, JOB_ID))

##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

rule custom_BLASTp:
    input:
        peptides = rules.Transdecoder_Predict.output.peptides
    output:
        blastP = join(OUT_DIR, 'BLAST_results', 'custom.blastp.outfmt6')
    benchmark:
        join(OUT_DIR, 'benchmarks', 'custom_BLASTp.tsv')
    message: 
        """--- Peptide search against Custom database (BLASTx)"""
    run:
        shell('mkdir -p ' + join(WORK_DIR, USER, JOB_ID) +
                ' && cp {input.peptides} ' + join(WORK_DIR, USER, JOB_ID) +
                ' && cp ' + CUSTOM_DATABASE + ' ' + join(WORK_DIR, USER, JOB_ID) +
                ' && cd ' + join(WORK_DIR, USER, JOB_ID) + 
                ' && makeblastdb -in ' + os.path.basename(CUSTOM_DATABASE) +
                ' -dbtype prot'
                ' && blastp'
                ' -query ' + 'stringtie_merged.gtf' + '_' + rstrip(os.path.basename(DNA), '.fa') + '_exons.fa.transdecoder.pep'
                ' -db ' + os.path.basename(CUSTOM_DATABASE) + ' -num_threads 16 -max_target_seqs 1 -outfmt 6 > custom.blastp.outfmt6'
                ' && mv custom.blastp.outfmt6 ' + join(OUT_DIR, 'BLAST_results'))
        shell('rm -r ' + join(WORK_DIR, USER, JOB_ID))

##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

rule Pfam:
    input:
        peptides = rules.Transdecoder_Predict.output.peptides
    output:
        pfam_out = join(OUT_DIR, 'Pfam_results', 'TrinotatePFAM.out')
    log:
        join(OUT_DIR, 'LOGS', 'Pfam.log')
    benchmark:
        join(OUT_DIR, 'benchmarks', 'Pfam.benchmark.tsv')
    message: 
        """--- Pfam search with HMMSCAN """
    run:
        shell('mkdir -p ' + join(WORK_DIR, USER, JOB_ID) +
                ' && cp {input.peptides} ' + join(WORK_DIR, USER, JOB_ID) +
                ' && cp ' + join(SHARED_DATABASE, 'Pfam*') + ' ' + join(WORK_DIR, USER, JOB_ID) +
                ' && cd ' + join(WORK_DIR, USER, JOB_ID) + 
                ' && hmmscan --cpu 16 --domtblout TrinotatePFAM.out Pfam-A.hmm '
                + 'stringtie_merged.gtf' + '_' + rstrip(os.path.basename(DNA), '.fa') + '_exons.fa' + '.transdecoder.pep > {log} 2>&1'
                ' && mv ' + join(WORK_DIR, USER, JOB_ID, 'TrinotatePFAM.out') + ' ' + join(OUT_DIR, 'Pfam_results'))
        shell('rm -r ' + join(WORK_DIR, USER, JOB_ID))

##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

rule signalP:
    input:
        peptides = rules.Transdecoder_Predict.output.peptides
    output:
        signalp = join(OUT_DIR, 'signalP.out')
    log:
        join(OUT_DIR, 'LOGS', 'signalp.log')
    benchmark:
        join(OUT_DIR, 'benchmarks', 'signalp.benchmark.tsv')
    message: 
        """--- Signal peptide earch with signalP"""
    run:
        shell('mkdir -p ' + join(WORK_DIR, USER, JOB_ID) +
                ' && cp {input.peptides} ' + join(WORK_DIR, USER, JOB_ID) +
                ' && cd ' + join(WORK_DIR, USER, JOB_ID) + 
                ' && signalp -f short -n signalP.out '
                + 'stringtie_merged.gtf' + '_' + rstrip(os.path.basename(DNA), '.fa') + '_exons.fa.transdecoder.pep > {log} 2>&1'
                ' && mv ' + join(WORK_DIR, USER, JOB_ID, 'signalP.out ') + OUT_DIR)
        shell('rm -r ' + join(WORK_DIR, USER, JOB_ID))

##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

rule TMHMM:
    input:
        peptides = rules.Transdecoder_Predict.output.peptides
    output:
        tmhmm = join(OUT_DIR, 'tmhmm.out')
    benchmark:
        join(OUT_DIR, 'benchmarks', 'tmhmm.benchmark.tsv')
    message: 
        """--- Transmembrane domain prediction """
    run:
        shell('mkdir -p ' + join(WORK_DIR, USER, JOB_ID) +
                ' && cp {input.peptides} ' + join(WORK_DIR, USER, JOB_ID) +
                ' && cd ' + join(WORK_DIR, USER, JOB_ID) + 
                ' && tmhmm --short < '
                + 'stringtie_merged.gtf' + '_' + rstrip(os.path.basename(DNA), '.fa') + '_exons.fa.transdecoder.pep > tmhmm.out'
                ' && mv ' + join(WORK_DIR, USER, JOB_ID, 'tmhmm.out ') + OUT_DIR)
        shell('rm -r ' + join(WORK_DIR, USER, JOB_ID))

##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

rule RNAmmer:
    input:
        exons = rules.gffread_exons.output.exons
    output:
        rnammer = join(OUT_DIR, 'stringtie_merged.gtf' + '_' + rstrip(os.path.basename(DNA), '.fa') + '_exons.fa.rnammer.gff')
    benchmark:
        join(OUT_DIR, 'benchmarks', 'rnammer.tsv')
    message: 
        """--- Find ribosomal RNA loci"""
    run:
        shell('mkdir -p ' + join(WORK_DIR, USER, JOB_ID) +
                ' && cp {input.exons} ' + join(WORK_DIR, USER, JOB_ID) +
                ' && cd ' + join(WORK_DIR, USER, JOB_ID) + 
                ' && ' + join(TRINOTATE_HOME, 'util', 'rnammer_support', 'RnammerTranscriptome.pl') +
                ' --transcriptome ' + 'stringtie_merged.gtf' + '_' + rstrip(os.path.basename(DNA), '.fa') + '_exons.fa'
                ' --path_to_rnammer /programs/rnammer-1.2/rnammer'
                ' && mv *rnammer.gff ' + OUT_DIR)
        shell('rm -r ' + join(WORK_DIR, USER, JOB_ID))

##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

rule Trinotate:
    input:
        exons = rules.gffread_exons.output.exons,
        peptides = rules.Transdecoder_Predict.output.peptides,
        blastX = rules.BLASTx.output.blastX,
        blastP = rules.BLASTp.output.blastP,
        cBlastX = rules.custom_BLASTx.output.blastX,
        cBlastP = rules.custom_BLASTp.output.blastP,
        pfam = rules.Pfam.output.pfam_out,
        signalp = rules.signalP.output.signalp,
        tmhmm = rules.TMHMM.output.tmhmm,
        rnammer = rules.RNAmmer.output.rnammer,
        sqlite = SQLITE
    output:
        join(OUT_DIR, 'Trinotate_report.xls'),
        join(OUT_DIR, 'Trinotate_report.xls.gene_ontology')
    benchmark:
        join(OUT_DIR, 'benchmarks', 'trinotate.tsv')
    message: 
        """--- Combining annotation outputs into SQLite database """
    run:
        shell('mkdir -p ' + join(WORK_DIR, USER, JOB_ID) +
                ' && cp {input.exons} {input.peptides} ' + join(WORK_DIR, USER, JOB_ID) +
                ' && cp {input.sqlite} ' + join(WORK_DIR, USER, JOB_ID) +
                ' && cd ' + join(WORK_DIR, USER, JOB_ID) + 
                ' && grep ">" {input.exons} | sed "s/>//g" | sed "s/gene=//g" | awk \'{{print $2"\t"$1}}\' | sort -u > ' + 
                'stringtie_merged.gtf' + '_' + rstrip(os.path.basename(DNA), '.fa') + '_gene_trans_map')
        shell('cd ' + join(WORK_DIR, USER, JOB_ID) +
                ' && Trinotate Trinotate.sqlite init'
                ' --gene_trans_map ' + 'stringtie_merged.gtf' + '_' + rstrip(os.path.basename(DNA), '.fa') + '_gene_trans_map'
                ' --transcript_fasta ' + 'stringtie_merged.gtf' + '_' + rstrip(os.path.basename(DNA), '.fa') + '_exons.fa'
                ' --transdecoder_pep ' + 'stringtie_merged.gtf' + '_' + rstrip(os.path.basename(DNA), '.fa') + '_exons.fa' + '.transdecoder.pep'
                ' && Trinotate Trinotate.sqlite LOAD_swissprot_blastp {input.blastP}'
                ' && Trinotate Trinotate.sqlite LOAD_swissprot_blastx {input.cBlastX}'
                ' && Trinotate Trinotate.sqlite LOAD_custom_blast --outfmt6 {input.cBlastP} --prog blastp --dbtype ' + os.path.basename(CUSTOM_DATABASE) +
                ' && Trinotate Trinotate.sqlite LOAD_custom_blast --outfmt6 {input.cBlastX} --prog blastx --dbtype ' + os.path.basename(CUSTOM_DATABASE) +
                ' && Trinotate Trinotate.sqlite LOAD_pfam {input.pfam}'
                ' && Trinotate Trinotate.sqlite LOAD_tmhmm {input.tmhmm}'
                ' && Trinotate Trinotate.sqlite LOAD_signalp {input.signalp}'
                ' && Trinotate Trinotate.sqlite LOAD_rnammer {input.rnammer}'
                ' && Trinotate Trinotate.sqlite report > Trinotate_report.xls'
                ' && ' + join(TRINOTATE_HOME, 'util', 'extract_GO_assignments_from_Trinotate_xls.pl') + ' --Trinotate_xls Trinotate_report.xls -G -I > Trinotate_report.xls.gene_ontology'
                ' && mv Trinotate.sqlite Trinotate_report.xls Trinotate_report.xls.gene_ontology ' + OUT_DIR)
        shell('rm -r ' + join(WORK_DIR, USER, JOB_ID))

