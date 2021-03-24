library(rcellminer)

# LOAD DATA ----
# RCELLMINER DATA ----
nci60MolDat <- rcellminer::getAllFeatureData(rcellminerData::molData)
nci60Exp <- nci60MolDat$exp
load("tcga_pancan_pathway_genes.rda")

# PARAMETERS ----
nrow <- 50 
ncol <- 1000

rows <- paste0("S", 1:nrow)
cols <- paste0("G", 1:ncol)

df <- matrix(rnorm(nrow*ncol), nrow=nrow, ncol=ncol, dimnames=list(row=rows, col=cols))

nsubset <- 7
t1 <- df[, 1:nsubset]

r1 <- cor(t1, df)
summary(as.vector(r1))

nsubset <- 1

df <- t(nci60Exp)

r2 <- NULL 
niter <- ncol(df) # 500
threshold <- 0.7

pb <- txtProgressBar(min=1, max=niter, style=3)
for(i in 1:niter) {
  setTxtProgressBar(pb, i)
  
  #selected_genes <- c("CDK4", "CDK6", "CCND1")
  #idx <- 1:nsubset
  
  #selected_genes <- sample(colnames(df), nsubset, replace=FALSE)
  #idx <- which(colnames(df) %in% selected_genes)
  
  idx <- i 

  t1 <- df[, idx]
  r1 <- cor(t1, df[, -idx])
  
  #t2 <- max(as.vector(r1), na.rm=TRUE)  
  t2 <- length(which(r1 >= threshold))
  #summary(as.vector(r1))
  r2 <- c(r2, t2)
}

r2 <- sort(r2, decreasing=TRUE)
hist(r2, breaks=287)
length(which(r2 != 0))
length(which(r2 == 1))

# X ----
genes_tcga <- sort(unique(unlist(tcga_pancan_pathway_genes)))
genes_tcga_available <- intersect(colnames(df), genes_tcga)

# GENERATE DATA ----
n_extra_genes <- 700
reltol <- 1e-2
extra_genes <- sample(colnames(df), n_extra_genes)
genes <- unique(c(genes_tcga_available, extra_genes))
x <- df[, genes] %>% as.matrix

tmp_beta0 <- rep(0, length(genes))
names(tmp_beta0) <- genes

## Selected features
genes_selected_features <- c("CDKN2A", "CDKN2C", "CDKN1A", "RB1", "CCND3", "CCND1")
tmp_beta0[genes_selected_features] <- 1*runif(length(genes_selected_features))

y_true <- x %*% tmp_beta0 %>% as.vector
tmp_beta0[tmp_beta0 > 0]
beta_true <- tmp_beta0

## Add noise 
beta_init <- 1*runif(length(tmp_beta0))
names(beta_init) <- names(beta_true)

# CORRELATION IN DATA.FRAME
r2 <- NULL 
niter <- ncol(x) # 500
threshold <- 0.8

pb <- txtProgressBar(min=1, max=niter, style=3)
for(i in 1:niter) {
  setTxtProgressBar(pb, i)
  
  idx <- i 
  
  t1 <- x[, idx]
  r1 <- cor(t1, x[, -idx])

  t2 <- length(which(r1 >= threshold))
  r2 <- c(r2, t2)
}

r2 <- sort(r2, decreasing=TRUE)
hist(r2, breaks=287)
length(which(r2 != 0))
length(which(r2 == 1))
sum(r2)
