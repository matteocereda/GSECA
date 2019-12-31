#  Toy analysis ========

# load configurations
source("Scripts/config.R")

# Gene expression matrix & sample type
M = read.delim("Examples/PRAD.ptenloss.M.tsv")
L = read.delim("Examples/PRAD.ptenloss.L.tsv")[,1]

# Gene set list
pl = read.gmt.file("gene_sets/cereda.158.KEGG.gmt")

# Run GSECA
res = GSECA_executor(  M # Gene Expression matrix
                     , L # Sample label list
                     , pl # gene set list
                     , outdir = "Results" #outdir folder
                     , analysis = "my_analysis"# analysis name
                     , N.CORES = 2 # number of cores
                     , EMPIRICAL = T # true if empirical p-value is requested
                     , BOOTSTRP = T
                     , nsim = 2 # number of bootstrapping
                     , AS = 0.25 # AS threshold
                     , PEMP = 1 # p.emp threshold
                     , SR   = 0.7 # success rate threshold
                     , toprank = 10 # top ranked pathways
                     , iphen = c("CASE", "CNTR") #  phenotype lables
)

devtools::session_info()

