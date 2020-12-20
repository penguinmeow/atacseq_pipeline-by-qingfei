########## Prepare TxDb from GTF file ##########

library(GenomicFeatures)
#TxDb.Mmusculus.GENCODE.mm10.ensemblGene <- makeTxDbFromGFF('/Volumes/yu3grp/scRNASeq/yu3grp/qpan/Database/References/mm10/Gencode/gencode.vM17.primary_assembly.annotation.gtf', format = "gtf",dataSource = "GENCODEvM17", organism = "Mus musculus", taxonomyId = 10090)
library(AnnotationDbi)
#saveDb(TxDb.Mmusculus.GENCODE.mm10.ensemblGene, "/Volumes/yu3grp/scRNASeq/yu3grp/qpan/Database/References/mm10/Gencode/TxDb.Mmusculus.GENCODE.mm10.ensemblGene")

########## Load Packages ##########
library(ChIPseeker)
library(ggplot2)

########## Input Information ##########
args = commandArgs(trailingOnly=TRUE)
#dir = "/Volumes/yu3grp/IO_JY/yu3grp/Obesity/Hurwitz/bulkATACseq/04_motif"
#samples = c("CD45CD19Bcell_B6_1", "CD45CD19Bcell_B6_2", "CD45CD19Bcell_DIO_1", "CD45CD19Bcell_DIO_2")
mm10 <- loadDb(file = "/research/rgs01/project_space/yu3grp/scRNASeq/yu3grp/qpan/Database/References/mm10/Gencode/TxDb.Mmusculus.GENCODE.mm10.ensemblGene")

########## Peak Annotation ##########

## Set the working directory

input_file <- args[1]
txdb <- args[2]

input_table <- read.table(input_file, header = F)
names(input_table) <- c("Chrom", "ChromStart", "ChromEnd")
input_table <- input_table[grep("chr", input_table$Chrom),]
input_table$Strand <- "*"; input_table$Score <- 1
input_obj <- GRanges(
    seqnames = Rle(input_table$Chrom),
    ranges = IRanges(start = input_table$ChromStart, end = input_table$ChromEnd, names = input_table$Peak_ID),
    strand = Rle(input_table$Strand), score = input_table$Score
)

peakAnno <- annotatePeak(input_obj, TxDb = mm10, annoDb = NULL, level = "transcript", tssRegion = c(-3000, 3000), flankDistance = 5000, addFlankGeneInfo = FALSE,
                         assignGenomicAnnotation = TRUE, genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron", "Downstream", "Intergenic"),
                         sameStrand = FALSE, ignoreOverlap = FALSE, ignoreUpstream = FALSE, ignoreDownstream = FALSE, overlap = "TSS", verbose = TRUE)
anno_table <- data.frame(peakAnno@anno)
masterTable <- cbind(input_table, anno_table)
write.table(masterTable, file = paste0(basename(dirname(args[1])), "/", "/04_footPrint_mpbs.annotate.txt"), quote = F, row.names = F, col.names = T, sep = "\t")
