library(dplyr)
su.te <- read.csv('su_ccRCC_te.csv')
rownames(su.te) <- su.te$X
su.te$X <- NULL
head(su.te)

# filter for different log2fc thresholds
su.te.06 <- su.te %>% filter((avg_log2FC > 0.6) & (p_val_adj < 0.01))
su.te.08 <- su.te %>% filter((avg_log2FC > 0.8) & (p_val_adj < 0.01)) 
su.te.1  <- su.te %>% filter((avg_log2FC > 1) & (p_val_adj < 0.01))

# prepare gtf's with different cutoffs
te.gtf <- rtracklayer::import('GRCh38_rmsk_TE.gtf')

te.gtf.06 <- te.gtf[te.gtf$gene_id %in% su.te.06$gene]
export(te.gtf.06, 'GRCh38_rmsk_TE_06_avglog2FC.gtf')

te.gtf.08 <- te.gtf[te.gtf$gene_id %in% su.te.08$gene]
export(te.gtf.08, 'GRCh38_rmsk_TE_08_avglog2FC.gtf')

te.gtf.1 <- te.gtf[te.gtf$gene_id %in% su.te.1$gene]
export(te.gtf.1, 'GRCh38_rmsk_TE_1_avglog2FC.gtf')

# run featureCounts
# last number and letter = type of quantification according cell paper
library(Rsubread)

##########################################################
#################### avglog2FC = 0.06 #################### 
##########################################################
s51.counts.06.2a <- featureCounts('s51_cellranger_count_outs/possorted_genome_bam.bam', 
                               annot.ext = 'GRCh38_rmsk_TE_06_avglog2FC.gtf',
                               isGTFAnnotationFile = TRUE,
                               GTF.attrType = "transcript_id", 
                               ignoreDup = TRUE, 
                               nthreads = 12,
                               tmpDir = 'tmp/')

# # stats
# s51.counts.06.2a$stat

# convert to df and sort 
df.06.s1 <- as.data.frame(s51.counts.06.2a$counts)
df.06.s1$names = rownames(df.06.s1)
head(df.06.s1[order(df.06.s1$possorted_genome_bam.bam, decreasing=TRUE),], 20)



s52.counts.06.2a <- featureCounts('s52_cellranger_count_outs/possorted_genome_bam.bam', 
                               annot.ext = 'GRCh38_rmsk_TE_06_avglog2FC.gtf',
                               isGTFAnnotationFile = TRUE,
                               GTF.attrType = "transcript_id", 
                               ignoreDup = TRUE, 
                               nthreads = 12,
                               tmpDir = 'tmp/')

# convert counts to df and sort 
df.06.s52 <- as.data.frame(s52.counts.06.2a$counts)
df.06.s52$names = rownames(df.06.s52)
head(df.06.s52[order(df.06.s52$possorted_genome_bam.bam, decreasing=TRUE),], 20)

# check if rownames are the same then combine 
all(rownames(df.06.s1) ==  rownames(df.06.s52))
combined.06.df <- df.06.s1
combined.06.df$possorted_genome_bam.bam <- df.06.s1$possorted_genome_bam.bam + df.06.s52$possorted_genome_bam.bam


head(combined.06.df[order(combined.06.df$possorted_genome_bam.bam, decreasing=TRUE),], 20)
##########################################################
#################### avglog2FC = 0.08 #################### 
##########################################################
s51.counts.08.2a <- featureCounts('s51_cellranger_count_outs/possorted_genome_bam.bam', 
                                  annot.ext = 'GRCh38_rmsk_TE_08_avglog2FC.gtf',
                                  isGTFAnnotationFile = TRUE,
                                  GTF.attrType = "transcript_id", 
                                  ignoreDup = TRUE, 
                                  nthreads = 12,
                                  tmpDir = 'tmp/')

# # stats
# s51.counts.08.2a$stat

# convert to df and sort 
df.08.s1 <- as.data.frame(s51.counts.08.2a$counts)
df.08.s1$names = rownames(df.08.s1)
head(df.08.s1[order(df.08.s1$possorted_genome_bam.bam, decreasing=TRUE),], 20)



s52.counts.08.2a <- featureCounts('s52_cellranger_count_outs/possorted_genome_bam.bam', 
                                  annot.ext = 'GRCh38_rmsk_TE_08_avglog2FC.gtf',
                                  isGTFAnnotationFile = TRUE,
                                  GTF.attrType = "transcript_id", 
                                  ignoreDup = TRUE, 
                                  nthreads = 12,
                                  tmpDir = 'tmp/')

# convert counts to df and sort 
df.08.s52 <- as.data.frame(s52.counts.08.2a$counts)
df.08.s52$names = rownames(df.08.s52)
head(df.08.s52[order(df.08.s52$possorted_genome_bam.bam, decreasing=TRUE),], 20)

# check if rownames are the same then combine 
all(rownames(df.08.s1) ==  rownames(df.08.s52))
combined.08.df <- df.08.s1
combined.08.df$possorted_genome_bam.bam <- df.08.s1$possorted_genome_bam.bam + df.08.s52$possorted_genome_bam.bam


head(combined.08.df[order(combined.08.df$possorted_genome_bam.bam, decreasing=TRUE),], 20)


##########################################################
#################### avglog2FC = 1 #################### 
##########################################################

s51.counts.1.2a <- featureCounts('s51_cellranger_count_outs/possorted_genome_bam.bam', 
                                  annot.ext = 'GRCh38_rmsk_TE_1_avglog2FC.gtf',
                                  isGTFAnnotationFile = TRUE,
                                  GTF.attrType = "transcript_id", 
                                  ignoreDup = TRUE, 
                                  nthreads = 12,
                                  tmpDir = 'tmp/')

# # stats
# s51.counts.1.2a$stat

# convert to df and sort 
df.1.s1 <- as.data.frame(s51.counts.1.2a$counts)
df.1.s1$names = rownames(df.1.s1)
head(df.1.s1[order(df.1.s1$possorted_genome_bam.bam, decreasing=TRUE),], 20)



s52.counts.1.2a <- featureCounts('s52_cellranger_count_outs/possorted_genome_bam.bam', 
                                  annot.ext = 'GRCh38_rmsk_TE_1_avglog2FC.gtf',
                                  isGTFAnnotationFile = TRUE,
                                  GTF.attrType = "transcript_id", 
                                  ignoreDup = TRUE, 
                                  nthreads = 12,
                                  tmpDir = 'tmp/')

# convert counts to df and sort 
df.1.s52 <- as.data.frame(s52.counts.1.2a$counts)
df.1.s52$names = rownames(df.1.s52)
head(df.1.s52[order(df.1.s52$possorted_genome_bam.bam, decreasing=TRUE),], 20)

# check if rownames are the same then combine 
all(rownames(df.1.s1) ==  rownames(df.1.s52))
combined.1.df <- df.1.s1
combined.1.df$possorted_genome_bam.bam <- df.1.s1$possorted_genome_bam.bam + df.1.s52$possorted_genome_bam.bam


head(combined.1.df[order(combined.1.df$possorted_genome_bam.bam, decreasing=TRUE),], 20)

# write all to csv
write.csv(combined.06.df[order(combined.06.df$possorted_genome_bam.bam, decreasing=TRUE),], 
          'featureCounts.su.0.01.adjpval.06.log2fc.csv')
write.csv(combined.08.df[order(combined.08.df$possorted_genome_bam.bam, decreasing=TRUE),], 
          'featureCounts.su.0.01.adjpval.08.log2fc.csv')
write.csv(combined.1.df[order(combined.1.df$possorted_genome_bam.bam, decreasing=TRUE),], 
          'featureCounts.su.0.01.adjpval.1.log2fc.csv')
