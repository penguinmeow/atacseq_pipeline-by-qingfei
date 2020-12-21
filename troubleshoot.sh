#!/bin/bash
#BSUB -P QFATACseq
#BSUB -J 1905371_WT_BMDM_ATACseq_1_S1_L004_R1_001_QFATACseq
#BSUB -oo ATAC_out/1905371_WT_BMDM_ATACseq_1_S1_L004_R1_001/ATACseqJ1log.out
#BSUB -eo ATAC_out/1905371_WT_BMDM_ATACseq_1_S1_L004_R1_001/ATACseqJ1log.err
#BSUB -n 1
#BSUB -N xzhen@stjude.org
module load trimmomatic/0.36
module load conda3/5.1.0
module load R/3.6.3
source activate /research/rgs01/project_space/yu3grp/software_JY/yu3grp/conda_env/yulab_env_2.7

fastqc=/research/rgs01/applications/hpcf/apps/fastqc/install/0.11.5/fastqc
Bowtie2=/research/rgs01/applications/hpcf/apps/bowtie/install/2.2.9/bin/bowtie2
samtools=/research/rgs01/applications/hpcf/apps/samtools/install/1.2/bin/samtools
java=/research/rgs01/applications/hpcf/apps/java/jdk1.8.0_66/bin/java
picard=/research/rgs01/applications/hpcf/apps/picard/install/2.16.0/picard.jar
bedtools=/research/rgs01/project_space/yu3grp/software_JY/yu3grp/yulab_apps/apps/bedtools2/bin/bedtools
macs2=/research/rgs01/project_space/yu3grp/software_JY/yu3grp/conda_env/yulab_env_2.7/bin/macs2
genomeDir=/research/projects/yu3grp/scRNASeq/yu3grp/qpan/Database/References/mm10/Gencode/Bowtie2/mm10
blacklist=/research/rgs01/project_space/yu3grp/software_JY/yu3grp/yulab_databases/ENCODE_blacklist/mm10-blacklist.v2.bed
mm10=/research/rgs01/project_space/yu3grp/scRNASeq/yu3grp/qpan/Database/References/mm10/Gencode/GRCm38.primary_assembly.genome.fa.size

python /home/xzhen/pipelines/atac_pipeline/bin/trimmomatic_run.py --fastqs pytest/1905371_WT_BMDM_ATACseq_1_S1_L004_R1_001.fastq.gz pytest/1905371_WT_BMDM_ATACseq_1_S1_L004_R2_001.fastq.gz --nth 8 --adapters /home/xzhen/adapters/NexteraPE-PE.fa --out-dir ATAC_out/1905371_WT_BMDM_ATACseq_1_S1_L004_R1_001 --paired-end --trim-reads 1

$fastqc ATAC_out/1905371_WT_BMDM_ATACseq_1_S1_L004_R1_001/1905371_WT_BMDM_ATACseq_1_S1_L004_R1_001.trim.paired.fastq.gz -o ATAC_out/1905371_WT_BMDM_ATACseq_1_S1_L004_R1_001
$fastqc ATAC_out/1905371_WT_BMDM_ATACseq_1_S1_L004_R1_001/1905371_WT_BMDM_ATACseq_1_S1_L004_R2_001.trim.paired.fastq.gz -o ATAC_out/1905371_WT_BMDM_ATACseq_1_S1_L004_R1_001

#bow2.sh in /ATAC
$Bowtie2 -x $genomeDir --very-sensitive --phred33 -X 2000 --no-mixed --no-discordant -p 8 -1 ATAC_out/1905371_WT_BMDM_ATACseq_1_S1_L004_R1_001/1905371_WT_BMDM_ATACseq_1_S1_L004_R1_001.trim.paired.fastq.gz -2 ATAC_out/1905371_WT_BMDM_ATACseq_1_S1_L004_R1_001/1905371_WT_BMDM_ATACseq_1_S1_L004_R2_001.trim.paired.fastq.gz | $samtools view -bS - > ATAC_out/1905371_WT_BMDM_ATACseq_1_S1_L004_R1_001/tmp.aligned.bam
$java -jar $picard AddOrReplaceReadGroups I= ATAC_out/1905371_WT_BMDM_ATACseq_1_S1_L004_R1_001/tmp.aligned.bam O= ATAC_out/1905371_WT_BMDM_ATACseq_1_S1_L004_R1_001/tmp.RGadded.bam SO=coordinate RGID=1905371_WT_BMDM_ATACseq_1_S1_L004_R1_001 RGLB=1905371_WT_BMDM_ATACseq_1_S1_L004_R1_001 RGPL=ATAC RGPU=Illumina RGSM=1905371_WT_BMDM_ATACseq_1_S1_L004_R1_001
$java -jar $picard MarkDuplicates I= ATAC_out/1905371_WT_BMDM_ATACseq_1_S1_L004_R1_001/tmp.RGadded.bam O= ATAC_out/1905371_WT_BMDM_ATACseq_1_S1_L004_R1_001/02_alignment.raw.bam M= ATAC_out/1905371_WT_BMDM_ATACseq_1_S1_L004_R1_001/tmp.dupMark.txt
$samtools view -b -F 780 -q 1 ATAC_out/1905371_WT_BMDM_ATACseq_1_S1_L004_R1_001/02_alignment.raw.bam > ATAC_out/1905371_WT_BMDM_ATACseq_1_S1_L004_R1_001/02_alignment.mapped.bam && $samtools flagstat ATAC_out/1905371_WT_BMDM_ATACseq_1_S1_L004_R1_001/02_alignment.mapped.bam > ATAC_out/1905371_WT_BMDM_ATACseq_1_S1_L004_R1_001/02_alignment.mapped.bam.txt
$samtools view -b -F 1024 ATAC_out/1905371_WT_BMDM_ATACseq_1_S1_L004_R1_001/02_alignment.mapped.bam > ATAC_out/1905371_WT_BMDM_ATACseq_1_S1_L004_R1_001/02_alignment.mapped_rmdup.bam && $samtools flagstat ATAC_out/1905371_WT_BMDM_ATACseq_1_S1_L004_R1_001/02_alignment.mapped_rmdup.bam > ATAC_out/1905371_WT_BMDM_ATACseq_1_S1_L004_R1_001/02_alignment.mapped_rmdup.bam.txt
$bedtools intersect -v -abam ATAC_out/1905371_WT_BMDM_ATACseq_1_S1_L004_R1_001/02_alignment.mapped_rmdup.bam -b $blacklist > ATAC_out/1905371_WT_BMDM_ATACseq_1_S1_L004_R1_001/02_alignment.mapped_rmdup_rmBLK.bam
$samtools index ATAC_out/1905371_WT_BMDM_ATACseq_1_S1_L004_R1_001/02_alignment.mapped_rmdup_rmBLK.bam
$samtools view -b ATAC_out/1905371_WT_BMDM_ATACseq_1_S1_L004_R1_001/02_alignment.mapped_rmdup_rmBLK.bam chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX chrY > ATAC_out/1905371_WT_BMDM_ATACseq_1_S1_L004_R1_001/02_alignment.mapped_rmdup_rmBLK_mitoFree.bam

$bedtools bamtobed -i ATAC_out/1905371_WT_BMDM_ATACseq_1_S1_L004_R1_001/02_alignment.mapped_rmdup_rmBLK_mitoFree.nucleoFree.sortedByName.bam | awk -v OFS="\t" '{if($6=="+"){print $1,$2+4,$3+4,$4,$5,$6}else if($6=="-"){print $1,$2-5,$3-5,$4,$5,$6}}' > ATAC_out/1905371_WT_BMDM_ATACseq_1_S1_L004_R1_001/02_alignment.mapped_rmdup_rmBLK_mitoFree.nucleoFree.sortedByName.shifted.bed
$bedtools bedtobam -ubam -i ATAC_out/1905371_WT_BMDM_ATACseq_1_S1_L004_R1_001/02_alignment.mapped_rmdup_rmBLK_mitoFree.nucleoFree.sortedByName.shifted.bed -g $mm10 > ATAC_out/1905371_WT_BMDM_ATACseq_1_S1_L004_R1_001/02_alignment.mapped_rmdup_rmBLK_mitoFree.nucleoFree.sortedByName.shifted.bam
$samtools sort ATAC_out/1905371_WT_BMDM_ATACseq_1_S1_L004_R1_001/02_alignment.mapped_rmdup_rmBLK_mitoFree.nucleoFree.sortedByName.shifted.bam ATAC_out/1905371_WT_BMDM_ATACseq_1_S1_L004_R1_001/02_alignment.mapped_rmdup_rmBLK_mitoFree.nucleoFree.sortedByName.shifted.sortedByPos
$samtools index ATAC_out/1905371_WT_BMDM_ATACseq_1_S1_L004_R1_001/02_alignment.mapped_rmdup_rmBLK_mitoFree.nucleoFree.sortedByName.shifted.sortedByPos.bam

conda deactivate
source activate atac
rgt-hint footprinting --atac-seq --paired-end --organism=mm10 --output-location=ATAC_out/1905371_WT_BMDM_ATACseq_1_S1_L004_R1_001 --output-prefix=04_footPrint 02_alignment.mapped_rmdup_rmBLK_mitoFree.nucleoFree.sortedByName.shifted.sortedByPos.bam 03_NucleosomeFree_peaks.narrowPeak
rgt-hint tracks --bc --bigWig --organism=mm10 02_alignment.mapped_rmdup_rmBLK_mitoFree.nucleoFree.sortedByName.shifted.sortedByPos.bam 03_NucleosomeFree_peaks.narrowPeak --output-prefix=ATAC_out/1905371_WT_BMDM_ATACseq_1_S1_L004_R1_001/footprintingTracks
rgt-motifanalysis matching --motif-dbs /research/rgs01/project_space/yu3grp/software_JY/yu3grp/git_repo/ATACseq_pipeline/data/rgtdata/motifs/transfac_mouse --input-files ATAC_out/1905371_WT_BMDM_ATACseq_1_S1_L004_R1_001/04_footPrint.bed --output-location=ATAC_out/1905371_WT_BMDM_ATACseq_1_S1_L004_R1_001 --organism=mm10
Rscript /research/rgs01/project_space/yu3grp/software_JY/yu3grp/git_repo/ATACseq_pipeline/scripts/peakAnnotation.R ATAC_out/1905371_WT_BMDM_ATACseq_1_S1_L004_R1_001/04_footPrint_mpbs.bed mm10
perl /research/rgs01/project_space/yu3grp/software_JY/yu3grp/git_repo/ATACseq_pipeline/scripts/annotationToNetwork.pl ATAC_out/1905371_WT_BMDM_ATACseq_1_S1_L004_R1_001/04_footPrint_mpbs.annotate.txt ATAC_out/1905371_WT_BMDM_ATACseq_1_S1_L004_R1_001/04_networkATACseq.txt
