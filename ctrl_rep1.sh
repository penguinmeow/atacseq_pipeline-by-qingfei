#BSUB -P ATAC-Seq
#BSUB -n 8
#BSUB -M 3000
#BSUB -oo ctrl_rep1.out -eo ctrl_rep1.err
#BSUB -J BAM_Manipulation
#BSUB -q priority



module load conda3/5.1.0
source activate /research/rgs01/project_space/yu3grp/software_JY/yu3grp/conda_env/yulab_env_2.7

export RGTDATA=/research/rgs01/project_space/yu3grp/software_JY/yu3grp/yulab_databases/rgtdata

fastqc=/research/rgs01/applications/hpcf/apps/fastqc/install/0.11.5/fastqc
cutadapt=/hpcf/apps/python/install/3.6.1/bin/cutadapt
Bowtie2=/research/rgs01/applications/hpcf/apps/bowtie/install/2.2.9/bin/bowtie2
samtools=/research/rgs01/applications/hpcf/apps/samtools/install/1.2/bin/samtools
java=/research/rgs01/applications/hpcf/apps/java/jdk1.8.0_66/bin/java
picard=/research/rgs01/applications/hpcf/apps/picard/install/2.16.0/picard.jar
bedtools=/research/rgs01/project_space/yu3grp/software_JY/yu3grp/yulab_apps/apps/bedtools2/bin/bedtools
macs2=/research/rgs01/project_space/yu3grp/software_JY/yu3grp/conda_env/yulab_env_2.7/bin/macs2
hint=/research/rgs01/project_space/yu3grp/software_JY/yu3grp/conda_env/yulab_env_2.7/bin/rgt-hint
motifanalysis=/research/rgs01/project_space/yu3grp/software_JY/yu3grp/conda_env/yulab_env_2.7/bin/rgt-motifanalysis

genomeDir=/research/projects/yu3grp/scRNASeq/yu3grp/qpan/Database/References/mm10/Gencode/Bowtie2/mm10
blacklist=/research/rgs01/project_space/yu3grp/software_JY/yu3grp/yulab_databases/ENCODE_blacklist/mm10-blacklist.v2.bed
mm10=/research/rgs01/project_space/yu3grp/scRNASeq/yu3grp/qpan/Database/References/mm10/Gencode/GRCm38.primary_assembly.genome.fa.size


dir=/research/rgs01/project_space/yu3grp/software_JY/yu3grp/git_repo/ATACseq_pipeline/test
sample=ctrl_rep1
echo  -e "$sample is processing...\n"
mkdir -p $dir/$sample

########## Section 0: Prepare RAW FASTQ Files
cat sample1_L001_R1.fq.gz sample1_L002_R1.fq.gz sample1_L003_R1.fq.gz sample1_L004_R1.fq.gz > $dir/$sample/00_fqRaw.R1.fq.gz
cat sample1_L001_R2.fq.gz sample1_L002_R2.fq.gz sample1_L003_R2.fq.gz sample1_L004_R2.fq.gz > $dir/$sample/00_fqRaw.R1.fq.gz

########## Section 1: Process FASTQ Files
## 1.1 Quality Control of RAW FASTQ Files by FastQC
$fastqc $dir/$sample/00_fqRaw.R1.fq.gz -o $dir/$sample
$fastqc $dir/$sample/00_fqRaw.R2.fq.gz -o $dir/$sample

## 1.2 Adaptor Trimming
$cutadapt -a CTGTCTCTTATA -A CTGTCTCTTATA --trim-n --max-n=0.5 --quality-base=33 -q 30,20 -m 30 -o $dir/$sample/01_fqClean.R1.fq.gz -p $dir/$sample/01_fqClean.R2.fq.gz $dir/$sample/00_fqRaw.R1.fq.gz $dir/$sample/00_fqRaw.R2.fq.gz > 01_fqClean.trimStat.log 2>&1
$fastqc $dir/$sample/01_fqClean.R1.fq.gz -o $dir/$sample
$fastqc $dir/$sample/01_fqClean.R2.fq.gz -o $dir/$sample

########## Section 2: Alignment and Peak Calling
## 2.1 Alignment to reference genome
$Bowtie2 -x $genomeDir --very-sensitive --phred33 -X 2000 --no-mixed --no-discordant -p 8 -1 $dir/$sample/01_fqClean.R1.fq.gz -2 $dir/$sample/01_fqClean.R2.fq.gz | $samtools view -bS - > $dir/$sample/tmp.aligned.bam

## 2.2 BAM Manipulation: add Read Groups, Mark Duplicates
$java -jar $picard AddOrReplaceReadGroups I=$dir/$sample/tmp.aligned.bam O=$dir/$sample/tmp.RGadded.bam SO=coordinate RGID=$sample RGLB=$sample RGPL=ATAC RGPU=Illumina RGSM=$sample
$java -jar $picard MarkDuplicates I=$dir/$sample/tmp.RGadded.bam O=$dir/$sample/02_alignment.raw.bam M=$dir/$sample/tmp.dupMark.txt

## 2.3 BAM Filtration: remove unmapped reads, remove duplicates, mask blacklist regions, remove Mitochondrial reads
$samtools view -b -F 780 -q 1 $dir/$sample/02_alignment.raw.bam > $dir/$sample/02_alignment.mapped.bam && $samtools flagstat $dir/$sample/02_alignment.mapped.bam > $dir/$sample/02_alignment.mapped.bam.txt
$samtools view -b -F 1024 $dir/$sample/02_alignment.mapped.bam > $dir/$sample/02_alignment.mapped_rmdup.bam && $samtools flagstat $dir/$sample/02_alignment.mapped_rmdup.bam > $dir/$sample/02_alignment.mapped_rmdup.bam.txt
$bedtools intersect -v -abam $dir/$sample/02_alignment.mapped_rmdup.bam -b $blacklist > $dir/$sample/02_alignment.mapped_rmdup_rmBLK.bam
$samtools index $dir/$sample/02_alignment.mapped_rmdup_rmBLK.bam
$samtools view -b $dir/$sample/02_alignment.mapped_rmdup_rmBLK.bam chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX chrY > $dir/$sample/02_alignment.mapped_rmdup_rmBLK_mitoFree.bam

## 2.4 BAM Split and Alignment Shift
$samtools view -h $dir/$sample/02_alignment.mapped_rmdup_rmBLK_mitoFree.bam | awk 'substr($0,1,1)=="@" || ($9<=146 && $9>=-146)' | $samtools view -bS - > $dir/$sample/02_alignment.mapped_rmdup_rmBLK_mitoFree.nucleoFree.bam
$samtools view -h $dir/$sample/02_alignment.mapped_rmdup_rmBLK_mitoFree.bam | awk 'substr($0,1,1)=="@" || ($9>146 && $9<=307) || ($9<-146 && $9>=-307)' | $samtools view -bS - > $dir/$sample/02_alignment.mapped_rmdup_rmBLK_mitoFree.nucleoMono.bam
$samtools view -h $dir/$sample/02_alignment.mapped_rmdup_rmBLK_mitoFree.bam | awk 'substr($0,1,1)=="@" || ($9>307) || ($9<-307)' | $samtools view -bS - > $dir/$sample/02_alignment.mapped_rmdup_rmBLK_mitoFree.nucleoTwoplus.bam

$samtools sort -n $dir/$sample/02_alignment.mapped_rmdup_rmBLK_mitoFree.nucleoFree.bam $dir/$sample/02_alignment.mapped_rmdup_rmBLK_mitoFree.nucleoFree.sortedByName
$bedtools bamtobed -i $dir/$sample/02_alignment.mapped_rmdup_rmBLK_mitoFree.nucleoFree.sortedByName.bam | awk -v OFS="\t" '{if($6=="+"){print $1,$2+4,$3+4,$4,$5,$6}else if($6=="-"){print $1,$2-5,$3-5,$4,$5,$6}}' > $dir/$sample/02_alignment.mapped_rmdup_rmBLK_mitoFree.nucleoFree.sortedByName.shifted.bed
$bedtools bedtobam -ubam -i $dir/$sample/02_alignment.mapped_rmdup_rmBLK_mitoFree.nucleoFree.sortedByName.shifted.bed -g $mm10 > $dir/$sample/02_alignment.mapped_rmdup_rmBLK_mitoFree.nucleoFree.sortedByName.shifted.bam
$samtools sort $dir/$sample/02_alignment.mapped_rmdup_rmBLK_mitoFree.nucleoFree.sortedByName.shifted.bam $dir/$sample/02_alignment.mapped_rmdup_rmBLK_mitoFree.nucleoFree.sortedByName.shifted.sortedByPos
$samtools index $dir/$sample/02_alignment.mapped_rmdup_rmBLK_mitoFree.nucleoFree.sortedByName.shifted.sortedByPos.bam

## 2.5 Peak Calling
$macs2 callpeak -f BED -n 03_NucleosomeFree -g mm -q 0.01 --nomodel --extsize 200 --shift -100 -t $dir/$sample/02_alignment.mapped_rmdup_rmBLK_mitoFree.nucleoFree.sortedByName.shifted.bed --outdir $dir/$sample

########## Section 3: Footprinting
$hint footprinting --atac-seq --paired-end --organism=mm10 --output-location=$dir/$sample --output-prefix=04_footPrint $dir/$sample/02_alignment.mapped_rmdup_rmBLK_mitoFree.nucleoFree.sortedByName.shifted.sortedByPos.bam $dir/$sample/03_NucleosomeFree_peaks.narrowPeak
$hint tracks --bc --bigWig --organism=mm10 $dir/$sample/02_alignment.mapped_rmdup_rmBLK_mitoFree.nucleoFree.sortedByName.shifted.sortedByPos.bam $dir/$sample/03_NucleosomeFree_peaks.narrowPeak --output-prefix=$dir/$sample/footprintingTracks
$motifanalysis matching --motif-dbs ~/rgtdata/motifs/transfac_mouse --input-files $dir/$sample/04_footPrint.bed --output-location=$dir/$sample --organism=mm10
Rscript $dir/../scripts/peakAnnotation.R $dir/$sample/04_footPrint_mpbs.bed mm10
perl $dir/../scripts/annotationToNetwork.pl $dir/$sample/04_footPrint_mpbs.annotate.txt $dir/$sample/04_networkATACseq.txt
