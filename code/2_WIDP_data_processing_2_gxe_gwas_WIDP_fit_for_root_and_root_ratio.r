#install necessary packages
#install.packages('devtools',repos = "http://cran.us.r-project.org", lib = '/storage/home/m/mul826/R/x86_64-redhat-linux-gnu-library/3.6')
library(devtools)
#install.packages('data.table', repos = "http://cran.us.r-project.org", lib = '/storage/home/m/mul826/R/x86_64-redhat-linux-gnu-library/3.6')
library(data.table)
#install_version("nloptr", version = "1.2.2", repos = "http://cran.us.r-project.org", lib = '/storage/home/m/mul826/R/x86_64-redhat-linux-gnu-library/3.6')
library(nloptr)
#install.packages('lme4', repos = "http://cran.us.r-project.org", lib = '/storage/home/m/mul826/R/x86_64-redhat-linux-gnu-library/3.6')
library(lme4)
#install.packages('glmnet',repos = "http://cran.us.r-project.org", lib = '/storage/home/m/mul826/R/x86_64-redhat-linux-gnu-library/3.6')
library(glmnet)
#devtools::install_github('deruncie/GridLMM',upgrade = FALSE)
library(GridLMM)
#install.packages('KRLS',repos = "http://cran.us.r-project.org", lib = '/storage/home/m/mul826/R/x86_64-redhat-linux-gnu-library/3.6')
#library(KRLS)
#install.packages('qqman',repos = "http://cran.us.r-project.org", lib = '/storage/home/m/mul826/R/x86_64-redhat-linux-gnu-library/3.6')
library(qqman)
library(dplyr)

#read phenotypic data
pheno <- read.csv('Genotype_add_ratio_blup_with_subpop_information.csv',header = TRUE, row.names = 'X')

pheno_trt_blup <- pheno %>%
    filter(source != 'root_blup') %>% # use ww_blup and ws_blup
    mutate(Treatment = case_when(source == 'ws_blup' ~ -1,
                                 source == 'ww_blup' ~ 1)) %>%
    droplevels()

#names(pheno_trt_blup)
#[1] "Taxa"          "MEDMETVA"      "MEDMETVD"      "TMETVA"        "NOMETV"        "RXSA"          "TCA"          
#[8] "TSA"           "AA"            "RCA"           "CCFN"          "CCS"           "ANG_TOP"       "BIOVEG"       
#[15] "BIOREPRO"      "YIELD"         "DIA_STM"       "TSA.RXSA"      "TSA.TCA"       "Subpopulation" "source"       
#[22] "Treatment"  

layers<-c("MEDMETVA", "MEDMETVD", "TMETVA", "NOMETV", "RXSA", "TCA", "TSA",  "AA", "RCA", "CCFN", "CCS", "ANG_TOP", "BIOVEG","BIOREPRO","YIELD","DIA_STM","TSA.RXSA", "TSA.TCA")

#genotype data was obtained from this paper: https://doi.org/10.1186/s12870-019-1653-x
#read genotypic data
myG_numeric <- fread('widiv_2015_175g_MAF_0.05_Tassel_numeric.txt')
mat<-data.frame(myG_numeric) # convert to data.frame
row.names(mat)<-mat[,1] #row.names = genotypes
mat<-mat[,-1] #remove the marker column
X<-mat
X<-X*2
X<-as.matrix(X)

X_1<-X[,-c(which(is.na(X),arr.ind=TRUE)[,2])] #dropping any column w/ NA for calculating kinship matrix, also the gxe_gwas can't work with NAs; this step dropped ~15% columns.

#kinship matrix
X_centered = sweep(X_1,2,colMeans(X_1),'-') # center marker genotypes; na.rm=TRUE allows NA to be omitted from calculation
K = tcrossprod(X_centered) / ncol(X_centered)
rownames(K) = colnames(K) = rownames(X)

K = K/mean(diag(K)) # this step is from Dan's code

#create the map file
map<-data.frame(snp=colnames(X),
                Chr=as.numeric(sapply(colnames(X),function(x) gsub(".*?([0-9]+).*", "\\1", x))), # extract first number of string
                pos=as.numeric(sapply(colnames(X),function(x) strsplit(x,'_')[[1]][2]))) #extract numbers after '_'

data<-pheno_trt_blup[pheno_trt_blup$Taxa %in% rownames(X),] #keep taxa the same as the genotype file

#make a loop to run through the manhattans:
sapply(layers,function(layer){
    trait<-data[,layer]
    #names(trait)<-trait$Taxa
    
    data_trait <- cbind.data.frame(trait, data[,c(1,22)])
    
    gxe_gwas = GridLMM_GWAS(
        formula = trait ~ 1 + Treatment + (1 + Treatment|Taxa), # the same error model is used for each marker. It is specified similarly to lmer
        test_formula = ~1 + Treatment, # this is the model for each marker. ~1 means an intercept for the marker. ~1 + cov means an intercept plus a slope on `cov` for each marker
        reduced_formula = ~1, # This is the null model for each test. ~0 means just the error model. ~1 means the null includes the intercept of the marker, but not anything additional
        data = data_trait, # The dataframe to look for terms from the 3 models
        weights = NULL, # optional observation-specific weights
        X = X_1, # The matrix of markers. Note: This can be of dimension n_g x p, where n_g is the number of unique genotypes.
        X_ID = 'Taxa', # The column of the data that identifies each genotype. Each level of data$Geno should match a rowname of X
        h2_start = NULL, # An optional vector of h2s to use as starting values for each random effect in the model. If NULL, will be calculated from the error model using GridLMM_ML
        h2_step = 0.01, # step size per random effect for testing alternate values of h2
        max_steps = 100, # maximum number of steps of size h2_step to take from h2_start
        X_map = map, # Optional. The marker positions.
        relmat = list(Taxa = K), # A list of Kernel matrices for the random effects. If X_ID (here Geno) is not included in this list, then it is calculated as tcrossprod(Xc)/ncol(Xc) where Xc is the centered (and optionally scaled) X. If any random effects are described in `error_model` but not provided here, the Kernel is assumed to be the identity matrix
        centerX = TRUE, # Should the markers be centered when calculating the GRM (only will be done if needed for calculating the GRM),
        scaleX = FALSE, # Should the markers be scaled to have constant variance when calculating the GRM?
        fillNAX = FALSE, # Should missing marker data be filled in with the mean allele frequency?
        method = 'REML', # Should the best model be selected by REML (if False, will be selected by ML)
        mc.cores = 1, #my_detectCores(), # How many cores should be used for parallel processing. Unless X is large, tends to actually be faster with mc.cores = 1
        verbose = FALSE # Should progress be printed to the screen?
    )
    
    #The GxE p-values are in the slot `p_value_ML`, or `p_value_REML.2`. The slot `p_value_REML.1` contains the p-values for the genotype main effect.
    
    colGxE<-data.frame(CHR=gxe_gwas$results$Chr,BP=gxe_gwas$results$pos,SNP=gxe_gwas$results$snp,
                       P_g=gxe_gwas$results$p_value_REML.1,P_gxe=gxe_gwas$results$p_value_REML.2,
                       fit_g = gxe_gwas$results$beta.3, fit_gxe = gxe_gwas$results$beta.4)
    #The GxE coefficient should be beta.4; the SNP main effect should be beta.3; the env covariate should be beta.2; and the intercept should be beta.1
    
    write.csv(colGxE, file = paste(layer,'G_GxE_GWAS_fitted_result.csv',sep = "_"))
    
    #colGxE_2 <- na.omit(colGxE)
    
    #manhattan plot G main effect only
    #png(filename = paste(layer,'G_GWAS_result.png',sep = "_"))
    #manhattan(colGxE_2, main = paste(layer,'G_GWAS',sep = '_'), chr="CHR", bp="BP", snp="SNP", p="P_g",suggestiveline = F, genomewideline = F)
    #dev.off()
    
})
