library(ggplot2)
library(tidyverse)
library(plyr)
# find all the count table files
files <- list.files(pattern = "*.htseq.counts.gz", recursive = TRUE)
# read count tables into a list of tables
datalist <- lapply(files, function(x){read.table(file=x,header=FALSE,col.names=c("gene", sub(".htseq.counts.gz", "", x)))})
# merge the individual count tables into a dataframe
m <- Reduce(function(...) merge(..., by=1, all = TRUE), datalist)
rownames(m) <- m[,1]
# get rid of the first few rows, they are summaries of the count tables
m <- m[6:nrow(m),-1]

# convert counts to z-score
m_scaled <- as.data.frame(t(scale(t(m))))

# read in the gene expression signature
signature9 <- read.csv('/Users/callamartyn/Systems/data/sc_signatures/bxpc3_leiden9_logfoldchangeGT50pct_genes.csv', stringsAsFactors = FALSE)
signature1 <- read.csv('/Users/callamartyn/Systems/data/sc_signatures/bxpc3_leiden1_logfoldchangeGT50pct_genes.csv', stringsAsFactors = FALSE)

# import expression table just to use for translating ensemble ids to gene names
express_table <- read.csv('/Users/callamartyn/Systems/data/sc_signatures/BXPC3_pvals_monacle.csv', stringsAsFactors = FALSE)
gene_dict <- express_table %>% select(gene_ids, gene_short_name)

# strip suffix (isoform) from ensembl IDs
m_scaled$gene_ids <- unlist(lapply(rownames(m_scaled),function(x) unlist(strsplit(x,'\\.'))[[1]]))
# join genedict table with expression table by ensemble IDs, table now has gene name column
m_gene_names <- join(m_scaled, gene_dict, by = "gene_ids", type = "inner")
# change rownames from ensembl IDs to gene names for easier selection
rownames(m_gene_names) <- m_gene_names$gene_short_name
# select only the rows matching the gene signature rows
m_sig1 <- m_gene_names[signature1$names,]
m_sig9 <- m_gene_names[signature9$names,]

write.csv(m_sig1, '/Users/callamartyn/Systems/Outputs/clin_zscores_sig1.csv')
write.csv(m_sig9, '/Users/callamartyn/Systems/Outputs/clin_zscores_sig9.csv')
