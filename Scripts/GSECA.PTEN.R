# PRAD - PTEN loss analysis
#
# 0. Resources ====
source("Scripts/config.R")

cat("Loading Data ...\n")
# Gene expression matrix & sample type
load("~/Dropbox (HuGeF)/GSECA_Paper/PRAD_PTEN.expr.M.L.Rdata")

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
                     , EMPIRICAL = F
                     , BOOTSTRP = F
                     , nsim = 2
                     , PSUMLOG = 0.01
                     , PADJ    = 0.1
                     , PEMP    = 1e-3
                     , SRATE   = 0.9
)


