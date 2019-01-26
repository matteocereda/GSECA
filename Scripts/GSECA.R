# PRAD - PTEN loss example analysis ========

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
                     , N.CORES = 2
                     , EMPIRICAL = T
                     , BOOTSTRP = T
                     , nsim = 2
                     , AS = 0.25
                     , PEMP    = 1
                     , SR   = 0.7
                     , toprank = 20
)



