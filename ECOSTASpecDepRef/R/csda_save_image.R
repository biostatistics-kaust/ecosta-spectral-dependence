#' [internal]
#' @export
csda_save_image <- function(q, fname, height=15, width=30){
  ggsave(paste0("../FIGURES/", fname, ".pdf"), q, width=width, height=height, units="cm", device=cairo_pdf)
  ggsave(paste0("../FIGURES/", fname, ".svg"), q, width=width, height=height, units="cm")
}
