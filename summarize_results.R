library(magrittr)

work_dir <- "run_o2"
# work_dir <- "./" 

files <- dir(work_dir, pattern="rda", full.names=TRUE)

results <- data.frame(pw_pval_true=numeric(0),
                      pw_pval_init=numeric(0),
                      pw_pval_final=numeric(0),
                      pw_pval_final_base=numeric(0),
                      
                      pw_min_true=character(0),
                      pw_min_init=character(0),
                      pw_min_final=character(0),
                      pw_min_final_base=character(0),
                      
                      cor_true_init=numeric(0),
                      cor_true_pred=numeric(0),
                      cor_true_pred_base=numeric(0),
                      
                      cor_pw_true_init=numeric(0),
                      cor_pw_true_pred=numeric(0),
                      cor_pw_true_pred_base=numeric(0),
                      cor_pw_is_final_better=logical(0), 
                      
                      seed=numeric(0),
                      
                      stringsAsFactors=FALSE)
for(i in 1:length(files)) {
  load(files[i])  
  
  # PRED VS TRUE CORRELATION ----
  y_init <- x %*% beta_init %>% as.vector
  y_pred <- x %*% beta_final %>% as.vector
  y_pred_base <- x %*% beta_final_base %>% as.vector
  
  pw_val_true <- x[, genes_selected_features] %*% beta_true[genes_selected_features]
  pw_val_init <- x[, genes_selected_features] %*% beta_init[genes_selected_features]
  pw_val_final <- x[, genes_selected_features] %*% beta_final[genes_selected_features]
  pw_val_final_base <- x[, genes_selected_features] %*% beta_final_base[genes_selected_features]
  
  tmp_results <- data.frame(
    pw_pval_true=gsea_true$pval[gsea_true$pathway == "Cell Cycle"],
    pw_pval_init=gsea_init$pval[gsea_init$pathway == "Cell Cycle"],
    pw_pval_final=gsea_final$pval[gsea_final$pathway == "Cell Cycle"],
    pw_pval_final_base=gsea_final_base$pval[gsea_final_base$pathway == "Cell Cycle"],
    
    pw_min_true=gsea_true$pathway[which.min(gsea_true$pval)],
    pw_min_init=gsea_init$pathway[which.min(gsea_init$pval)],
    pw_min_final=gsea_final$pathway[which.min(gsea_final$pval)],
    pw_min_final_base=gsea_final_base$pathway[which.min(gsea_final_base$pval)],
    
    cor_true_init=cor(y_true, y_init),
    cor_true_pred=cor(y_true, y_pred),
    cor_true_pred_base=cor(y_true, y_pred_base),
    
    cor_pw_true_init=cor(pw_val_true, pw_val_init),
    cor_pw_true_pred=cor(pw_val_true, pw_val_final),
    cor_pw_true_pred_base=cor(pw_val_true, pw_val_final_base),
    cor_pw_is_final_better=ifelse(cor(pw_val_true, pw_val_final) > cor(pw_val_true, pw_val_final_base), TRUE, FALSE), 
    
    seed=seed,
    
    stringsAsFactors=FALSE)

  results <- rbind(results, tmp_results)
}

#results <- unique(results)
saveRDS(results, "run_o2_results.rds")

results <- readRDS("run_o2_results.rds")
results <- unique(results[, -which(names(results) %in% c("seed"))])

# COR SUMMARY 
boxplot(results$cor_true_pred, results$cor_true_pred_base, xlim=c(0, 3), ylim=c(-0.1, 1.1), main="COR")
summary(results$cor_true_pred)
summary(results$cor_true_pred_base)

summary(results$cor_true_pred)
summary(results$cor_true_pred_base)

# PW SUMMARY
boxplot(results$pw_pval_final, results$pw_pval_final_base, xlim=c(0, 3), ylim=c(-0.1, 1.1), main="PVAL")
summary(results$pw_pval_final)
summary(results$pw_pval_final_base)

summary(results$pw_pval_final)
summary(results$pw_pval_final_base)

length(which(results$pw_min_final == "Cell Cycle"))
length(which(results$pw_min_final_base == "Cell Cycle"))

length(which(results$pw_min_final == "Cell Cycle" & results$cor_true_pred >= 0.7 & results$pw_pval_final <= 0.05))
length(which(results$pw_min_final_base == "Cell Cycle" & results$cor_true_pred_base >= 0.7 & results$pw_pval_final_base <= 0.05))


