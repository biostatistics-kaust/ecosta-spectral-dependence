#' @export
csda_theme <- function(plain=FALSE) {
  #print("aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa")
  select_font <- function(name, alternative) if(name %in% extrafont::fonts()) name else alternative
  base_family = select_font("Segoe UI", select_font("Helvetica", select_font("Times New Roman", "")))
  if(plain){
    return(theme_void(base_family=base_family))
  }
  theme_minimal(base_family=base_family)
  #theme_minimal(base_family="Serif")
  #theme_minimal()
}
