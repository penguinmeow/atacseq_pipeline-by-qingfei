# atacseq_pipeline-by-qingfei
Following is the overview of the pipeline
![image](https://github.com/penguinmeow/atacseq_pipeline-by-qingfei/blob/main/ATAC-seq%20(by%20Qingfei).png)

**Modifications I made to the original version by Qingfei**
1. adding a python script to submit tasks in parallel;
2. using trimmomatics rather than cutadapt to trim adapters in my data.
3. construct a conda env called atac, in which rgt is installed.
4. Problem1 : Rscript: cannot detect GenomicFeatures package.
Reason: R version in the conda env by qingfei cannot cannot install GenomicFeatures.
Solved: module load R/3.6.3
install package ChIPseeker in R(3.6.3). 
Run smoothly in lsf.
5. waiting to see if perl script runs good.
