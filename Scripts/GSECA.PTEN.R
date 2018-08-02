# PRAD - PTEN loss analysis

# load configurations
source("Scripts/config.R")

# Gene expression matrix & sample type
M = read.delim("Examples/PRAD.ptenloss.M.tsv")
L = read.delim("Examples/PRAD.ptenloss.L.tsv")[,1]

# Pathway list
pl <- read.gmt.file("gene_sets/cereda.158.KEGG.gmt")

res = GSECA_executor(  M
                     , L
                     , symbol="ensembl_gene_id"
                     , pl
                     , outdir = "Results"
                     , analysis = "PTEN"
                     , p_adj_th = 0.1
                     , nClass = 7
                     , N.CORES = 2
                     , EMPIRICAL = T
                     , BOOTSTRP = T
                     , nsim = 2
                     , PSUMLOG = 0.25
                     , PADJ    = 0.1
                     , PEMP    = 1
                     , SRATE   = 0.7
                     , toprank = 20
)


