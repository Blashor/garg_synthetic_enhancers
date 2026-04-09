#Sum counts of duplicate genes from counts.csv file produced from script rna_featurecounts.sh; this makes the count list compatible with subsequent edgeR analysis
#Read raw count table
counts_raw <- read.table(
  "counts_genesymbol.txt", #save original csv file generated from "rna_featurecounts.sh" as txt file and input here
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

#Define first column as containing gene symbols
gene_ids <- counts_raw[, 1]

#Remove columns
counts_raw <- counts_raw[, -(2:6)]

#Define remaining columns as containing counts
count_matrix <- counts_raw[, -1]

#Trim ENSEMBL Gene IDs such that they do not contain a decimal point and subsequent number (e.g. x.7)
gene_ids <- trimws(as.character(gene_ids))

#Sum counts across duplicate genes
counts_summed <- rowsum(count_matrix, group = gene_ids)

#Convert rownames into a gene symbol column
counts_summed_df <- data.frame(
  Geneid = rownames(counts_summed),
  counts_summed,
  row.names = NULL
)

write.table(
  counts_summed_df,
  file = "counts_genesymbol_summed.txt",
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

########
#The following edgeR pipeline is described as that utilized to analyze data from gastruloids derived from Shh PEP lines.
#For each sample, of which there are six, their corresponding paired-end reads were individually processed, aligned, and counted as described in rna_featurecounts.sh prior to this analysis.
########

library(edgeR) #(v4.8.2)
library(ggplot2) #(v4.0.2)

#Read raw count table
counts <- read.table("counts_genesymbol_summed.txt", header=TRUE, row.names=1)

#Define conditions (three groups with 2 replicates each)
group <- factor(c("0", "0", "20", "20", "50", "50"))

#Create DGE object
y <- DGEList(counts = counts, group = group)

#Filter lowly expressed genes
keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes = FALSE]

#EdgeR normalization
y <- calcNormFactors(y)

#Design matrix that has sample and condition information
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)
design

#Estimate dispersion of data
y <- estimateDisp(y, design)
plotBCV(y)

#Fit model
fit <- glmQLFit(y, design)l

#Define comparisons
contrast_20_vs_0 <- makeContrasts(20 - 0, levels = design)
contrast_50_vs_0 <- makeContrasts(50 - 0, levels = design)

# Store comparisons in a list
contrasts_list <- list(
  20_vs_0 = contrast_20_vs_0,
  50_vs_0 = contrast_50_vs_0,
)

# Run differential expression analysis for comparisons and store results in a list
results <- lapply(contrasts_list, function(ct) {
  qlf <- glmQLFTest(fit, contrast = ct)
  topTags(qlf, n = Inf)$table 
})

#Write to csv
for(name in names(results)) {
    write.csv(results[[name]], paste0(name, "_DE_results.csv"))
}