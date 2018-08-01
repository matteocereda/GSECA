source("Scripts/config.R")

# load("Examples/toy_dataset.Rdata")
# load("Examples/toy_dataset.positional.Rdata")
# load("Examples/PRAD_ERG_fusion.Rdata")
# load("~/Dropbox (HuGeF)//GSECA/Examples/PRAD_ERG_fusion.Rdata")

load("/sto1/matteo/GSECA/PRAD_PTEN_loss/Rdata/TCGA/PRAD/PTEN_loss_sample.Rdata")
load("/sto1/matteo/GSECA/PRAD_PTEN_loss/Rdata/TCGA/PRAD/PRADexp_FPKM_complete.Rdata")
expr = subset(PRAD_exp, sample=='TP'  & gene_type=='protein_coding')

expr$type = factor("CNTR", levels=c("CNTR","CASE"))
expr$type[which(substr(expr$Tumor_Sample_Barcode, 1, 12)%in%pten_loss)] = "CASE"

save(expr, file="~/Dropbox (HuGeF)//GSECA/Examples/PRAD_PTEN.expr.Rdata")

M = dcast(expr, ensembl_gene_id~Tumor_Sample_Barcode, value.var = 'value')
L = unique(expr[,c('Tumor_Sample_Barcode','type')])$type
colnames(M)
save(M,L, file="Rdata/PRAD_PTEN.expr.M.L.Rdata")

x=lapply(pl, paste, collapse="\t")
y=mapply(function(a,b){paste0(a, '\t',b)}, names(x), x, SIMPLIFY = F)
z=do.call(rbind,y)

write(z, file="gene_sets/cereda.158.KEGG.gmt", ncolumns = 1)

# load("~/Dropbox (HuGeF)//GSECA/Examples/PRAD_PTEN_fusion.Rdata")

# 01. PREPARE DATASET ====

load("~/Dropbox (HuGeF)//GSECA/Examples/PRAD_PTEN.expr.M.L.Rdata")

expr = get_expression_dataset(M, L, 'ensembl_gene_id') # expr = get_expression_dataset(M, L, 'symbol')

expr = subset(expr, gene_type=='protein_coding')

# 02. MIXTURE MODEL ====

expr <- get_mixture_expr_class(expr
                               , normalization = "FPKM"
                               , ne_value=0.01
                               , nc = 3
)

# 03. LOAD GENE SETS ====

# pl = read.gmt.file("gene_sets/ncomm.cereda.186.KEGG.gmt")
# # pl = read.gmt.file("gene_sets/msigdb/c1.all.v6.0.symbols.gmt")
# sp = pl['KEGG_SPLICEOSOME']
# load("~/Lavoro/Prabs/Rdata/KEGG_gene_sets.157.Rdata")
# pl = c(pl, sp)

load("/sto1/matteo/GSECA/PRAD_PTEN_loss/Rdata/gene_sets/Andrea.Rdata")

# 04. SEA ====

gcr  = gene_class_representation(expr$expr_class
                                 , pl
                                 , method ='fisher'# 'fisher'
                                 , TAIL = 'two.sided'
                                 , PW_SIM = 1
                                 , correction='fdr'
)

# 05. GSECA ====

gseca = GSECA(  pl
                , gcr
                , method ='fisher' # fisher, chisq, wilcoxon
                , PW_SIM = 100
                , TAIL = "two.sided"
                , POWER = 0.9
                , correction='fdr')


save(expr, gcr, gseca, pl, file="Results/prad.PTEN.Rdata")

# 06. EMPIRICAL PVALUES ====

gseca.emp = GSECA.Bootstrap.empirical( expr$expr_class
                                       , gseca
                                       , pl
                                       , NSIM=1000
                                       , nc=3
                                       , method = 'fisher'
                                       , PW_SIM = 1
                                       , TAIL = "two.sided"
                                       , POWER = 0.9
                                       , correction = 'fdr'
                                       , cutoff = 0.01
)

gseca$p.emp = gseca.emp$p.emp[match(gseca$gene_set, gseca.emp$gene_set)]

save(expr, gcr, gseca, gseca.emp, pl, file="Results/prad.PTEN.Rdata")


# 06. CORRECTION FOR SAMPLE SIZE ====

# gseca.boot=GSECA.Bootstrap.sample_size(   expr$expr_class
#                              ,gseca
#                              , pl
#                              , phen = c("CASE","CNTR")
#                              , NSIM=1000
#                              , nc=3
#                              # , sig.threshold=0.01
#                              , method = 'fisher'
#                              , PW_SIM = 1
#                              , TAIL = "two.sided"
#                              , POWER = 0.9
#                              , correction = 'fdr'
#                              , cutoff = 0.01
# )

gseca.boot=GSECA.Bootstrap.sample_size(   expr$expr_class
                                          ,gseca
                                          , pl
                                          , phen = c("CASE","CNTR")
                                          , NSIM=1000
                                          , nc=3
                                          # , sig.threshold=0.01
                                          , method = 'fisher'
                                          , PW_SIM = 1
                                          , TAIL = "two.sided"
                                          , POWER = 0.9
                                          , correction = 'fdr'
                                          , cutoff = 0.05
)

gseca$success_rate = gseca.boot$success_rate[match(gseca$gene_set, gseca.boot$gene_set)]

save(expr, gcr, gseca, gseca.boot, gseca.emp, pl, file="Results/prad.PTEN.sr_0.05.Rdata")


View(unique(gseca[,c("gene_set","sumlog","p.emp","success_rate")]))

# 06. GSECA ====

# PSUMLOG = 5e-5
PSUMLOG = 0.01
PADJ    = 0.1
PW      = 0
PEMP    = 1e-3
# PEMP    = 5e-3
# SRATE   = 0.7
# SRATE   = 0.8
SRATE   = 0.9

ecmap = GSECA.ECmap( gseca
                     , filename=paste0("Figures/gseca_PTEN.sumlog.psumlog_",PSUMLOG,"_padj_",PADJ,"_pw_",PW,"_pemp_",PEMP,"_srate_",SRATE,".pdf")
                     , p_adj = PADJ
                     , pow = PW
                     , psumlog=PSUMLOG
                     , pemp = PEMP
                     , srate=SRATE)

kegg = read.csv("/sto1/matteo/KEGG.csv")
gseca = cbind.data.frame(kegg[match(gseca$gene_set, kegg$msigDB_158),c('category','superclass','kegg_id','gene_set')], fixed = fix_names(gseca$gene_set), gseca)

write.csv(gseca, file='/sto1/matteo/GSECA/Results/gseca_PTEN.csv')

prad = unique(gseca[,c('category','superclass','kegg_id','gene_set','sumlog','p.emp','success_rate')])
write.csv(prad, file=paste0('/sto1/matteo/GSECA/Results/gseca_PTEN.simple.csv'))

# 07. Assess gene expression distributions ====

plot_expression(expr, pl, gene_set="chr8q24",
                my_sample="TCGA-V1-A8WV-01A-11R-A37L-07")
