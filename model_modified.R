library(numDeriv) 
library(fgsea)
library(magrittr)
#library(rcellminer)
library(impute)

source("tic_toc.R")

# LOAD DATA ---- 
## fgsea examples
#data(examplePathways)
#data(exampleRanks)

## Load rcellminer data 
#nci60MolDat <- rcellminer::getAllFeatureData(rcellminerData::molData)
#nci60Exp <- nci60MolDat$exp
#saveRDS(nci60Exp, "nci60Exp_bioc-2.12.0.rds")
nci60Exp <- readRDS("nci60Exp_bioc-2.12.0.rds")

# Load metadata
load("tcga_pancan_pathway_genes.rda")
hgnc <- read.table("hgnc_custom_20200106.txt", sep="\t", header=TRUE, comment.char="", quote="", stringsAsFactors=FALSE)

# PREPARE DATA ----
tcga_pancan_pathway_genes
tcga_genes <- sort(unique(unlist(tcga_pancan_pathway_genes)))

# # OPTIMIZE EXAMPLE ----
# # Gradient 
# #func <- function(x){ sin(x) }
# #x <- c(-pi, -0.5*pi, 0, 0.5*pi, pi)
# #numDeriv::grad(func, x)
# 
# loss_func <- function(x) {
#   (x[2] - x[1] * x[1])^2 + (1 - x[1])^2 + x[3]^2 + (1 - x[4])^3 + x[5]^2
# }
# beta_init <- c(-1.2, 1, 0.5, 0.2, -0.3)
# 
# # Find gradient of function and optimize function 
# results_grad <- numDeriv::grad(loss_func, beta_init)
# results_optim <- stats::optim(beta_init, loss_func, method="BFGS")
# 
# if(results_optim$convergence != 0) { stop("Convergence not reached") }
# 
# x <- results_optim$par
# loss_func(x)

# OPTIMIZE GSEA ----
loss_func_gsea <- function(beta_ranks, y_true, x, gsea_scale=1, log_file) {
  #beta_ranks <- beta_init; y_true <- y_true; x <- x; gsea_scale <- 1; log_file <- "del.txt"
  
  l1_scale <- 1e-1
  #gsea_scale <- 1 #1e2
  model_error_scale <- 1
  
  # REGULARIZATION ----
  # FIXME USE L0???
  l1_val <- l1_scale*sum(abs(beta_ranks))
  
  # FIXME DO NOT PUT THIS IN LOSS FUNCTION
  #beta_ranks[beta_ranks < 1e-6] <- 0
  
  # GSEA ----
  # NOTE: All genes must be given or sampling in GSEA will give incorrect p-values (p-values: not used here)
  # FIXME: Set pathways? 
  # NOTE: Sort the beta values and randomize 0 values so that similarly named genes of one pathway are not together (e.g., WNT)
  rand_beta_ranks <- c(sort(beta_ranks[beta_ranks > 0], decreasing=TRUE), sample(beta_ranks[beta_ranks <= 0]))
  results_fgsea <- fgsea(pathways=pathways, stats=rand_beta_ranks, minSize=3, maxSize=100, nperm=2)
  
  # Print statements do not output
  #min_padj <- min(results_fgsea$padj)
  #min_padj <- round(min_padj, 5)
  es_val <- max(results_fgsea$ES, na.rm=TRUE)
  es_size <- results_fgsea$size[results_fgsea$ES == es_val]
  es_pathway <- results_fgsea$pathway[results_fgsea$ES == es_val] 
  no_pathway <- setdiff(names(beta_ranks[beta_ranks > 0]), unique(unlist(pathways)))
  no_pathway_penalty <- sum(beta_ranks[no_pathway]^2)
  # GSEA penalty based on strength of un-normalized enrichment score (no permutations) and the size of gene set
  # Larger genesets are incentivized
  gsea_val <- gsea_scale*(((no_pathway_penalty/es_size)-es_val)+1) # Add +1 so lowest is 0
  #gsea_val <- ifelse(is.na(gsea_val), 1, gsea_val) # IGNORE
    
  # MODEL ERROR ----
  y_pred <- x %*% beta_ranks %>% as.vector
  model_error <- mean((y_true-y_pred)^2)
  model_error <- model_error_scale*model_error
    
  # RESULT ----  
  result_val <- l1_val + gsea_val + model_error
  
  # DEBUG ----
  sum_ypred <- sum(y_pred)
  sum_ytrue <- sum(y_true)
  sum_beta <- sum(beta_ranks)
  sum_x <- sum(x)
  rank_val <- paste(beta_ranks[abs(beta_ranks) > 0], collapse=";")
  
  out_str <- paste0(result_val, "|", 
                    gsea_val, "|", 
                    l1_val, "|", 
                    model_error, "|", 
                    sum_ypred, "|",
                    sum_ytrue, "|",
                    sum_beta, "|",
                    sum_x, "|", 
                    es_pathway)
  
  if(!is.null(log_file)) {
    cat(out_str, file=log_file, append=TRUE, sep="\n")    
  }
  
  return(result_val)
}

seed <- floor(runif(1, 0, 100000))
seed
#set.seed(234)
set.seed(seed)

# PARAMETERS ----
#n_genes <- 20
n_responses <- 50
n_extra_genes <- 700
reltol <- 1e-2
max_iter <- 1e4
optim_control_params <- list(trace=1, maxit=max_iter, reltol=reltol) # reltol=1e-8, 1e-2 faster
#optim_control_params <- list(trace=1, maxit=100, abstol=1e-2) # reltol=1e-8, 1e-2 faster
# 1e3 (default): From https://www.biostars.org/p/387492/ https://gsea-msigdb.github.io/gseapreranked-gpmodule/v6/index.html
nperm <- 1e4
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

# GENERATE DATA ----
## Get NCI60 data
### Impute data or the y_true will end up with NAs and there this error: initial value in 'vmmin' is not finite
nci60ExpImputed <- invisible(impute.knn(nci60Exp)$data)
df <- t(nci60ExpImputed)

## Set up gene lists 
pathways <- tcga_pancan_pathway_genes
genes_tcga <- sort(unique(unlist(tcga_pancan_pathway_genes)))
genes_tcga_available <- intersect(colnames(df), genes_tcga)
extra_genes <- sample(colnames(df), n_extra_genes)
genes <- unique(c(genes_tcga_available, extra_genes))

# Make starting matrix
x <- df[, genes] %>% as.matrix

## Make beta coefficient vectors
tmp_beta0 <- rep(0, length(genes))
names(tmp_beta0) <- genes

## Selected features
genes_selected_features <- c("CDKN2A", "CDKN2C", "CDKN1A", "RB1", "CCND3", "CCND1")
tmp_beta0[genes_selected_features] <- 1*runif(length(genes_selected_features))

tmp_beta0[tmp_beta0 > 0]
beta_true <- tmp_beta0
y_true <- x %*% tmp_beta0 %>% as.vector

## Add noise 
beta_init <- 1*runif(length(tmp_beta0))
names(beta_init) <- names(beta_true)
y_init <- x %*% beta_init %>% as.vector

## Correlation true to init
cor(y_true, y_init)

# RUN ANALYSIS ---
head(beta_init)
tmp_df <- data.frame(gene=names(beta_init), beta_true=beta_true, beta_init=beta_init, stringsAsFactors=FALSE)

# RUN OPTIMIZATION ----
method <- "BFGS" # "Nelder-Mead"
log_headers <- "result_val|gsea_val|l1_val|model_error|sum_ypred|sum_ytrue|sum_beta|sum_x|pathway"
## Run without GSEA
tic()
gsea_scale <- 0
verbose <- TRUE
log_file <- paste0("log_file_base_", timestamp, "_", seed, "seed.txt")
cat(log_headers, file=log_file, append=FALSE, sep = "\n")
results_optim_base <- stats::optim(par=beta_init, fn=loss_func_gsea, gr=NULL, y_true, x, gsea_scale, log_file, method=method, control=optim_control_params)
beta_final_base <- results_optim_base$par
toc()
## Run with GSEA
tic()
gsea_scale <- 1
verbose <- TRUE
log_file <- paste0("log_file_", timestamp, "_", seed, "seed.txt")
cat(log_headers, file=log_file, append=FALSE, sep = "\n")
results_optim <- stats::optim(par=beta_init, fn=loss_func_gsea, gr=NULL, y_true, x, gsea_scale, log_file, method=method, control=optim_control_params)
beta_final <- results_optim$par
toc()

# CHECK RESULTS ----
rand_beta_true <- c(sort(beta_true[beta_true > 0], decreasing=TRUE), sample(beta_true[beta_true <= 0]))
gsea_true <- fgsea(pathways=pathways, stats=rand_beta_true, minSize=3, maxSize=100, nperm=nperm)
gsea_true
min(gsea_true$pval)
min(gsea_true$ES)
rand_beta_init <- c(sort(beta_init[beta_init > 0], decreasing=TRUE), sample(beta_init[beta_init <= 0]))
gsea_init <- fgsea(pathways=pathways, stats=rand_beta_init, minSize=3, maxSize=100, nperm=nperm)
gsea_init
min(gsea_init$pval)
min(gsea_init$ES)
sum(gsea_init$ES)
rand_beta_final <- c(sort(beta_final[beta_final > 0], decreasing=TRUE), sample(beta_final[beta_final <= 0]))
gsea_final <- fgsea(pathways=pathways, stats=rand_beta_final, minSize=3, maxSize=100, nperm=nperm)
gsea_final
min(gsea_final$pval)
min(gsea_final$ES)
sum(gsea_final$ES)
rand_beta_final_base <- c(sort(beta_final_base[beta_final_base > 0], decreasing=TRUE), sample(beta_final_base[beta_final_base <= 0]))
gsea_final_base <- fgsea(pathways=pathways, stats=rand_beta_final_base, minSize=3, maxSize=100, nperm=nperm)
gsea_final_base
min(gsea_final_base$pval)
min(gsea_final_base$ES)
sum(gsea_final_base$ES)

# Final pathway matches initial? 
gsea_true$pval[gsea_true$pathway == "Cell Cycle"]
gsea_init$pval[gsea_init$pathway == "Cell Cycle"]
gsea_final$pval[gsea_final$pathway == "Cell Cycle"]
gsea_final_base$pval[gsea_final_base$pathway == "Cell Cycle"]

gsea_true$pathway[which.min(gsea_true$pval)]
gsea_init$pathway[which.min(gsea_init$pval)]
gsea_final$pathway[which.min(gsea_final$pval)]
gsea_final_base$pathway[which.min(gsea_final_base$pval)]

tmp_df$beta_final <- beta_final
tmp_df$abs_diff_init <- abs(tmp_df$beta_final - tmp_df$beta_init)
tmp_df$abs_diff_true <- abs(tmp_df$beta_final - tmp_df$beta_true)
tmp_df$in_pathway <- FALSE
tmp_df$in_pathway[tmp_df$gene %in% unique(unlist(pathways))] <- TRUE
#results_grad <- numDeriv::grad(loss_func, beta_ranks_final)

# PRED VS TRUE CORRELATION ----
y_pred <- x %*% beta_final %>% as.vector
y_pred_base <- x %*% beta_final_base %>% as.vector
plot(y_true, y_init)
plot(y_true, y_pred)
plot(y_true, y_pred_base)
cor(y_true, y_init)
cor(y_true, y_pred)
cor(y_true, y_pred_base)

sorted_tmp_df <- tmp_df[order(-tmp_df$abs_diff_true),]
sorted_tmp_df$diff_rank <- 1:nrow(sorted_tmp_df)
pathway <- "MYC"
sorted_tmp_df_pathway <- sorted_tmp_df[pathways[[pathway]], ]
hist(sorted_tmp_df_pathway$diff_rank, breaks=nrow(tmp_df), main=pathway)
pathway <- "Cell Cycle"
sorted_tmp_df_pathway <- sorted_tmp_df[pathways[[pathway]], ]
hist(sorted_tmp_df_pathway$diff_rank, breaks=nrow(tmp_df), main=pathway)

filename <- paste0("run_", 
                   length(genes), "selected_", 
                   n_extra_genes, "extra_", 
                   length(beta_true), "total_", 
                   n_responses, "samp_",
                   seed, "seed_",
                   timestamp, ".rda")
save.image(filename)
