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
    parser.add_argument('--memory', default='10GB', type=str, help= 'Memory requested to run the analysis.')
    parser.add_argument('--queue', default='standard', type=str,help='Queue to submit the job in HPCF (use bqueues to choose).')
    parser.add_argument('--out-dir', type=str,help='Output Directory.')
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

def create_job_file_pe(samplefile1, samplefile2, out_dir):
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

    ########## Section 3: Footprinting
    ### Load all the required module for analysis:
    module1 = 'module load R/3.6.3 \n'
    module1 += 'module load conda3/5.1.0 \n'
    module1 += 'source activate /research/rgs01/project_space/yu3grp/software_JY/yu3grp/conda_env/yulab_env_2.7 \n'
    module1 += 'export RGTDATA=/research/rgs01/project_space/yu3grp/software_JY/yu3grp/yulab_databases/rgtdata \n'
    
    apps1 = 'hint=/research/rgs01/project_space/yu3grp/software_JY/yu3grp/conda_env/yulab_env_2.7/bin/rgt-hint \n'
    apps1 += 'motifanalysis=/research/rgs01/project_space/yu3grp/software_JY/yu3grp/conda_env/yulab_env_2.7/bin/rgt-motifanalysis \n'
    apps1 += 'Rscript=/research/rgs01/project_space/yu3grp/scRNASeq/yu3grp/qpan/Software/R/v3.6.1/bin/Rscript \n'
    apps1 += 'genomeDir=/research/projects/yu3grp/scRNASeq/yu3grp/qpan/Database/References/mm10/Gencode/Bowtie2/mm10 \n'
    apps1 += 'blacklist=/research/rgs01/project_space/yu3grp/software_JY/yu3grp/yulab_databases/ENCODE_blacklist/mm10-blacklist.v2.bed \n'
    apps1 += 'mm10=/research/rgs01/project_space/yu3grp/scRNASeq/yu3grp/qpan/Database/References/mm10/Gencode/GRCm38.primary_assembly.genome.fa.size \n'

    job_body8 = ''

    log.info("Foot Printing......")
    job_body8 = '$hint footprinting --atac-seq --paired-end --organism=mm10 --output-location={} --output-prefix=04_footPrint {}/02_alignment.mapped_rmdup_rmBLK_mitoFree.nucleoFree.sortedByName.shifted.sortedByPos.bam {}/03_NucleosomeFree_peaks.narrowPeak\n'
    job_body8 += '$hint tracks --bc --bigWig --organism=mm10 {}/02_alignment.mapped_rmdup_rmBLK_mitoFree.nucleoFree.sortedByName.shifted.sortedByPos.bam {}/03_NucleosomeFree_peaks.narrowPeak --output-prefix={}/footprintingTracks\n'
    job_body8 += '$motifanalysis matching --motif-dbs $RGTDATA/motifs/transfac_mouse --input-files {}/04_footPrint.bed --output-location={} --organism=mm10\n'
    job_body8 += '$Rscript /research/rgs01/project_space/yu3grp/software_JY/yu3grp/git_repo/ATACseq_pipeline/scripts/peakAnnotation.R {}/04_footPrint_mpbs.bed mm10\n'
    job_body8 += 'perl /research/rgs01/project_space/yu3grp/software_JY/yu3grp/git_repo/ATACseq_pipeline/scripts/annotationToNetwork.pl {}/04_footPrint_mpbs.annotate.txt {}/04_networkATACseq.txt\n'
    job_body8 = job_body8.format(prefix, prefix, prefix, prefix, prefix, prefix, prefix, prefix, prefix, prefix, prefix)
    
    jobfile = prefix+".sh"
    with open(jobfile,"w") as new_file:
        new_file.write(job_header+module1+'\n'+apps1 +'\n' + job_body8+'\n')
    return jobfile
    
def submit_job(jobf):
  os.system('{}'.format(hpcfsubmit) + jobf)

log.info('Generating job file for each sample and submitting jobs....')

for fastq in range(0,len(fastqR1)):
    submit_job(create_job_file_pe(fastqR1[fastq], fastqR2[fastq], args.out_dir))
