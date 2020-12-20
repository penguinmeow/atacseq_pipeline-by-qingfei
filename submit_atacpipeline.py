#!/usr/bin/env python

# Script to submit ATAC-seq jobs in parallel
# Author Qingfei and Zhen Xie

import glob
import subprocess
import re
import os
import sys
import argparse
sys.path.append("/home/xzhen/pipelines/atac_pipeline/bin")
from encode_lib_common import (
    copy_f_to_dir, copy_f_to_f, log, ls_l, mkdir_p, read_tsv, rm_f,
    run_shell_cmd, strip_ext_fastq)

def parse_arguments():
    parser = argparse.ArgumentParser(prog='Wrapper for BSUB job submission for ATAC-seq data.', description='')
    parser.add_argument('--path-to-fastqs', default='', type=str, help= 'Path to FASTQ files.')
    parser.add_argument('--adapters',default='/home/xzhen/adapters/NexteraPE-PE.fa', type=str, help='Fasta file of adapter sequences')
    parser.add_argument('--trim-reads',default=1,type=int,choices=[1,0] , help='Option to trim or not to trim the input reads')
    parser.add_argument('--memory', default='10GB', type=str, help= 'Memory requested to run the analysis.')
    parser.add_argument('--queue', default='standard', type=str,help='Queue to submit the job in HPCF (use bqueues to choose).')
    parser.add_argument('--out-dir', default='ATAC_out', type=str,help='Output Directory.')
    parser.add_argument('--log-level', default='INFO',
                        choices=['NOTSET', 'DEBUG', 'INFO',
                                 'WARNING', 'CRITICAL', 'ERROR',
                                 'CRITICAL'],
                        help='Log level')
    args = parser.parse_args()
    log.setLevel(args.log_level)
    log.info(sys.argv)
    return args
    

args = parse_arguments()
log.info("Making output directory...")
mkdir_p(args.out_dir)

fastqR1 = glob.glob(args.path_to_fastqs+'/*R1_001.fastq.gz')
fastqR2 = glob.glob(args.path_to_fastqs+'/*R2_001.fastq.gz')

## Sort the list to include same sample at same index of two list
fastqR1.sort()
fastqR2.sort()

hpcfsubmit = 'bsub ' + '-R ' + '"rusage[mem=' + args.memory + ']" ' + '-q ' + args.queue + ' < '
#print 'Job submitted with command: {}'.format(hpcfsubmit)

def create_job_file_pe(samplefile1, samplefile2, adapters, out_dir, trim_reads):
    basename1 = os.path.basename(strip_ext_fastq(samplefile1))
    basename2 = os.path.basename(strip_ext_fastq(samplefile2))
    prefix = os.path.join(out_dir,basename1)
    
    job_header = '#!/bin/bash\n'
    job_header += '#BSUB -P QFATACseq\n'
    job_header += '#BSUB -J {}_QFATACseq\n'
    job_header += '#BSUB -oo {}'+'/ATACseqJ1log.out\n'
    job_header += '#BSUB -eo {}'+'/ATACseqJ1log.err\n'
    job_header += '#BSUB -n 1\n'
    job_header += '#BSUB -N xzhen@stjude.org\n'
    job_header = job_header.format(basename1,prefix, prefix)
    
    
    ### Load all the required module for analysis:
    module1 = 'module load trimmomatic/0.36\n'
    module1 += 'module load conda3/5.1.0\n'
    module1 += 'module load R/3.6.3\n'
    module1 += 'source activate /research/rgs01/project_space/yu3grp/software_JY/yu3grp/conda_env/yulab_env_2.7\n'
    
    apps1 = 'fastqc=/research/rgs01/applications/hpcf/apps/fastqc/install/0.11.5/fastqc\n'
    apps1 += 'Bowtie2=/research/rgs01/applications/hpcf/apps/bowtie/install/2.2.9/bin/bowtie2\n'
    apps1 += 'samtools=/research/rgs01/applications/hpcf/apps/samtools/install/1.2/bin/samtools\n'
    apps1 += 'java=/research/rgs01/applications/hpcf/apps/java/jdk1.8.0_66/bin/java\n'
    apps1 += 'picard=/research/rgs01/applications/hpcf/apps/picard/install/2.16.0/picard.jar\n'
    apps1 += 'bedtools=/research/rgs01/project_space/yu3grp/software_JY/yu3grp/yulab_apps/apps/bedtools2/bin/bedtools\n'
    apps1 += 'macs2=/research/rgs01/project_space/yu3grp/software_JY/yu3grp/conda_env/yulab_env_2.7/bin/macs2\n'
    apps1 += 'genomeDir=/research/projects/yu3grp/scRNASeq/yu3grp/qpan/Database/References/mm10/Gencode/Bowtie2/mm10\n'
    apps1 += 'blacklist=/research/rgs01/project_space/yu3grp/software_JY/yu3grp/yulab_databases/ENCODE_blacklist/mm10-blacklist.v2.bed\n'
    apps1 += 'mm10=/research/rgs01/project_space/yu3grp/scRNASeq/yu3grp/qpan/Database/References/mm10/Gencode/GRCm38.primary_assembly.genome.fa.size\n'
    
    job_body1 = ''
    job_body2 = ''
    job_body3 = ''
    job_body4 = ''
    job_body5 = ''
    job_body6 = ''
    job_body7 = ''

    
    ### Write job body to run each wrapper for the sample:
    
    log.info("Trimming with trimmomatic...")
    job_body1 = 'python /home/xzhen/pipelines/atac_pipeline/bin/trimmomatic_run.py --fastqs {} {} --nth 8 --adapters {} --out-dir {} --paired-end --trim-reads {}\n'
    job_body1 = job_body1.format(samplefile1, samplefile2, adapters, prefix, trim_reads)
    
    log.info("Doing Fastqc...")
    job_body2 = '$fastqc {} -o {}\n'
    job_body2 += '$fastqc {} -o {}\n'
    job_body2.format(prefix +'/'+ basename1 + '.trim.paired.fastq.gz',prefix, prefix +'/'+ basename2 + '.trim.paired.fastq.gz',prefix)
    
    log.info("Mapping paired-end files to reference genome...")
    job_body3 = '$Bowtie2 -x $genomeDir --very-sensitive --phred33 -X 2000 --no-mixed --no-discordant -p 8 -1 {} -2 {} | $samtools view -bS - > {}tmp.aligned.bam'
    job_body3 = job_body3.format(prefix +'/'+ basename1 + '.trim.paired.fastq.gz',prefix +'/'+ basename2 + '.trim.paired.fastq.gz',prefix+'/')
    
    log.info("PICARD is processing......")
    
    job_body4 = '$java -jar $picard AddOrReplaceReadGroups I= {}tmp.aligned.bam O= {}tmp.RGadded.bam SO=coordinate RGID={} RGLB={} RGPL=ATAC RGPU=Illumina RGSM={}\n'
    job_body4 += '$java -jar $picard MarkDuplicates I= {}tmp.RGadded.bam O= {}02_alignment.raw.bam M= {}tmp.dupMark.txt'
    job_body4 = job_body4.format(prefix+'/',prefix+'/',basename1, basename1, basename1, prefix+'/', prefix+'/', prefix+'/')
    
    log.info("Filtration is processing...")
    job_body5 = '$samtools view -b -F 780 -q 1 {}02_alignment.raw.bam > {}02_alignment.mapped.bam && $samtools flagstat {}02_alignment.mapped.bam > {}02_alignment.mapped.bam.txt\n'
    job_body5 += '$samtools view -b -F 1024 {}02_alignment.mapped.bam > {}02_alignment.mapped_rmdup.bam && $samtools flagstat {}02_alignment.mapped_rmdup.bam > {}02_alignment.mapped_rmdup.bam.txt\n'
    job_body5 += '$bedtools intersect -v -abam {}02_alignment.mapped_rmdup.bam -b $blacklist > {}02_alignment.mapped_rmdup_rmBLK.bam\n'
    job_body5 += '$samtools index {}02_alignment.mapped_rmdup_rmBLK.bam\n'
    job_body5 += '$samtools view -b {}02_alignment.mapped_rmdup_rmBLK.bam chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX chrY > {}02_alignment.mapped_rmdup_rmBLK_mitoFree.bam\n'
    job_body5 = job_body5.format(prefix+'/', prefix+'/', prefix+'/', prefix+'/', prefix+'/', prefix+'/', prefix+'/', prefix+'/', prefix+'/', prefix+'/', prefix+'/', prefix+'/', prefix+'/')
    
    log.info("Tn5 shift...")
    job_body6 = "$samtools view -h {}02_alignment.mapped_rmdup_rmBLK_mitoFree.bam | awk 'substr($0,1,1)=="@" || ($9<=146 && $9>=-146)' | $samtools view -bS - > {}02_alignment.mapped_rmdup_rmBLK_mitoFree.nucleoFree.bam\n"
    job_body6 += "$samtools view -h {}02_alignment.mapped_rmdup_rmBLK_mitoFree.bam | awk 'substr($0,1,1)=="@" || ($9>146 && $9<=307) || ($9<-146 && $9>=-307)' | $samtools view -bS - > {}02_alignment.mapped_rmdup_rmBLK_mitoFree.nucleoMono.bam\n"
    job_body6 += "$samtools view -h {}02_alignment.mapped_rmdup_rmBLK_mitoFree.bam | awk 'substr($0,1,1)=="@" || ($9>307) || ($9<-307)' | $samtools view -bS - > {}02_alignment.mapped_rmdup_rmBLK_mitoFree.nucleoTwoplus.bam\n"
    job_body6 += "$samtools sort -n {}02_alignment.mapped_rmdup_rmBLK_mitoFree.nucleoFree.bam {}02_alignment.mapped_rmdup_rmBLK_mitoFree.nucleoFree.sortedByName\n"
    job_body6 += "$bedtools bamtobed -i {}02_alignment.mapped_rmdup_rmBLK_mitoFree.nucleoFree.sortedByName.bam | awk -v OFS="\t" '{if($6=="+"){print $1,$2+4,$3+4,$4,$5,$6}else if($6=="-"){print $1,$2-5,$3-5,$4,$5,$6}}' > {}02_alignment.mapped_rmdup_rmBLK_mitoFree.nucleoFree.sortedByName.shifted.bed\n"
    job_body6 += "$bedtools bedtobam -ubam -i {}02_alignment.mapped_rmdup_rmBLK_mitoFree.nucleoFree.sortedByName.shifted.bed -g $mm10 > {}02_alignment.mapped_rmdup_rmBLK_mitoFree.nucleoFree.sortedByName.shifted.bam\n"
    job_body6 += "$samtools sort {}02_alignment.mapped_rmdup_rmBLK_mitoFree.nucleoFree.sortedByName.shifted.bam {}02_alignment.mapped_rmdup_rmBLK_mitoFree.nucleoFree.sortedByName.shifted.sortedByPos\n"
    job_body6 += "$samtools index {}02_alignment.mapped_rmdup_rmBLK_mitoFree.nucleoFree.sortedByName.shifted.sortedByPos.bam\n"
    job_body6 = job_body6.format(prefix+'/', prefix+'/', prefix+'/', prefix+'/', prefix+'/', prefix+'/', prefix+'/', prefix+'/', prefix+'/', prefix+'/', prefix+'/', prefix+'/', prefix+'/', prefix+'/', prefix+'/')
    
    log.info("Footprinting......")
    job_body7 = 'conda deactivate\n'
    job_body7 += 'source activate atac \n'
    job_body7 += 'rgt-hint footprinting --atac-seq --paired-end --organism=mm10 --output-location={} --output-prefix=04_footPrint 02_alignment.mapped_rmdup_rmBLK_mitoFree.nucleoFree.sortedByName.shifted.sortedByPos.bam 03_NucleosomeFree_peaks.narrowPeak\n'
    job_body7 += 'rgt-hint tracks --bc --bigWig --organism=mm10 02_alignment.mapped_rmdup_rmBLK_mitoFree.nucleoFree.sortedByName.shifted.sortedByPos.bam 03_NucleosomeFree_peaks.narrowPeak --output-prefix={}/footprintingTracks\n'
    job_body7 += 'rgt-motifanalysis matching --motif-dbs /research/rgs01/project_space/yu3grp/software_JY/yu3grp/git_repo/ATACseq_pipeline/data/rgtdata/motifs/transfac_mouse --input-files {}/04_footPrint.bed --output-location={} --organism=mm10\n'
    job_body7 += 'Rscript /research/rgs01/project_space/yu3grp/software_JY/yu3grp/git_repo/ATACseq_pipeline/scripts/peakAnnotation.R {}/04_footPrint_mpbs.bed mm10\n'
    job_body7 += 'perl /research/rgs01/project_space/yu3grp/software_JY/yu3grp/git_repo/ATACseq_pipeline/scripts/annotationToNetwork.pl {}/04_footPrint_mpbs.annotate.txt {}/04_networkATACseq.txt\n'
    job_body7 = job_body7.format(prefix, prefix, prefix, prefix, prefix, prefix, prefix)
    
    
    jobfile = prefix+".sh"
    with open(jobfile,"w") as new_file:
        new_file.write(job_header+module1+'\n'+job_body1+'\n'+job_body2+'\n'+job_body3+'\n'+job_body4+'\n'+job_body5+'\n')
    return jobfile


def submit_job(jobf):
  os.system('{}'.format(hpcfsubmit) + jobf)

log.info('Generating job file for each sample and submitting jobs....')

for fastq in range(0,len(fastqR1)):
    submit_job(create_job_file_pe(fastqR1[fastq], fastqR2[fastq], args.adapters, args.out_dir, args.trim_reads))

