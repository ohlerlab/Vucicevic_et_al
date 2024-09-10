#' Transforming into narrow data.frame
#'
#' This function will transform a a wide dataframe into a narrow one, in which columns containing timepoints will be
#' transformed into a variable named group and the count information into a variable named counts.
#' @param log_normed_df is a data.frame that has been normalized and log2 transformed by the log_norm function.
#' @return slim_df long form data.frame with single column for all Time-points.
#' @export
slim <- function(log_normed_df) {
  tidyr::pivot_longer(log_normed_df, cols = -c(sgRNA, Gene), names_to = "Group", values_to = "Counts") %>%
    dplyr::mutate("Replicate" = stringr::str_split_fixed(Group, "-", 2)[, 1]) %>%
    dplyr::mutate("Time" = as.numeric(stringr::str_split_fixed(Group, "_", 2)[, 2]))
}
