#' Normalize function
#'
#' This function normalized the count number by the library size (colSum) and multiplies by 10e6 to transform into
#' transcripts per million. The function is intended to be maped (tidyverse) or applied (base R).
#' @param x is the column that it wants to be applied to.
#' @return returns a normalized data.frame object.
#' @export
ls_norm <- function(x) {
  y <- numeric(length(x))
  for (i in 1:length(x)) {
    y[i] <- (x[i] / sum(x)) * 1e6
  }
  return(y)
}
#' Normalize and log transform
#'
#' This function will be applied to any double value in the dataframe that it is applied to. It will first ungroup
#' any grouped dataframe to prevent mistakes in the normalization. Then it will use the
#' function norm to normalize the column and then add a pseudocount and log2 transform the column.
#' @param input_df is the input data.frame following column name instructions (see Vignette)
#' @return norm_log_df is a Log2 transfromed, library size normalized dataframe
#' @export
log_norm <- function(input_df) {
  input_df %>%
    dplyr::ungroup() %>%
    dplyr::mutate_if(is.numeric, ls_norm) %>%
    dplyr::mutate_if(is.numeric, function(x) log2(x + 1))
}