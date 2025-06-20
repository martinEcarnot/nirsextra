#' gridcv_mec
#'
#' Computes statistics to compare prediction and reference values
#'
#' @import dplyr
#' @export

mser = function(fm)  {


  result <- df %>%
    group_by(nlv) %>%
    summarise(
      npred = n(),
      sep = sep(y1),
      med = median(y1),
      ecart_type = sd(y1),
      .groups = "drop"
    ) %>%
    arrange(nlv)

}
