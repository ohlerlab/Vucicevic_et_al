#' Running model in parallel
#'
#' This function will generate a cluster to run in parallel the linear model. Afterwards it will transform the output of the models into a data.frame containing the gene name
#' the p-value, the adjusted p-value according to the method selected ('fdr' by default), and the model estimate (or slope).
#' @param slim_df output of the slim function. Narrow, log-transformed, normalized data.frame
#' @param method p value adjustment method that wants to be used ('fdr' by default)
#' @param replicates indicates if the data contains biological replicates.
#' @return results with columns 'Gene','Estimate','adj_pvalue'
#' @export
model_testing <- function(slim_df, method ='fdr',replicates=FALSE){
  n_cores <- ifelse(parallel::detectCores()>8,8,parallel::detectCores()-2)
  cl <- parallel::makeCluster(n_cores, type="FORK")
  if (replicates==T) {
    # pvals <- parallel::parLapply(cl, unique(slim_df$Gene), function(gene) {
    pvals <- pbapply::pblapply(unique(slim_df$Gene), function(gene) {
      LM = lme4::lmer(Counts ~ Time + (1 | sgRNA) + (1 | Replicate), data = slim_df[slim_df$Gene == gene, ],REML = FALSE)
      LM_min = lme4::lmer(Counts ~ (1 | sgRNA) + (1 | Replicate), data = slim_df[slim_df$Gene ==  gene, ],REML = FALSE)
      A = anova(LM, LM_min)
      p = A$`Pr(>Chisq)`[2]
      estimate = (broom.mixed::tidy(LM)$estimate)[2]
      res <- list(p, estimate)
      return(res)
    }, cl=cl)
    parallel::stopCluster(cl)
    res_pvals <- data.frame(matrix(base::unlist(pvals), nrow = length(pvals), byrow = T))
    results <- data.frame(Gene = unique(slim_df$Gene), pvalue = res_pvals$X1, estimate = res_pvals$X2)
    results <- results %>% dplyr::mutate(adj_pvalue = p.adjust(pvalue, method = as.character(method)))
  }else{
    # pvals <- parallel::parLapply(cl,unique(slim_df$Gene),function(gene){
    pvals <- pbapply::pblapply(unique(slim_df$Gene), function(gene){
      LM = lme4::lmer(Counts ~ Time + (1|sgRNA)  ,data = slim_df[slim_df$Gene==gene,],REML = FALSE)
      LM_min = lme4::lmer(Counts ~ (1|sgRNA)  ,data = slim_df[slim_df$Gene==gene,],REML = FALSE)
      A=anova(LM,LM_min)
      p=A$`Pr(>Chisq)`[2]
      estimate=(broom.mixed::tidy(LM)$estimate)[2]
      res <- list(p,estimate)
      return(res)
    },cl=cl)
    parallel::stopCluster(cl)
    res_pvals <- data.frame(matrix(unlist(pvals), nrow=length(pvals), byrow=T))
    results <- data.frame('Gene'=unique(slim_df$Gene),'pvalue'=res_pvals$X1, 'estimate'=res_pvals$X2)
    results <- results %>%
      dplyr::mutate(adj_pvalue = p.adjust(pvalue, method=as.character(method)))
  }
}
