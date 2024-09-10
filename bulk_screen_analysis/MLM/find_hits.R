#' Running all functions
#'
#' Launcher that performs all data.frame transformations in the correct order.
#' @param df each guide as a row, and timepoints in the columns. Requires grouping (gene/region) column named Gene.
#' @param method pvalue adjustment method desired, default is FDR.
#' @param replicates indicates if the data.frame has biological replicates that need to be integrated.
#' @export
find_hits <- function(df, method = "fdr", replicates = FALSE) {
  results <- df %>%
    log_norm() %>%
    slim() %>%
    model_testing(method = as.character(method), replicates = replicates)
  return(results)
}
#' Filtering hits for shiny
#'
#' Function that filters hits according to FDR value.
#' @param fdr_value FDR threshold, default = 0.05
filter_hits <- function(results, threshold = 0.05) {
  dplyr::filter(results, adj_pvalue < threshold)
}