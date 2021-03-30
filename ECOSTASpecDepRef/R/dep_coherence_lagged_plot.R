#' @param sign_lags    We have two interpretations of x(t-14)
#'                   lag is defined as x(t-h), lag is therefore 14 (sign_lags="-")
#'                   or lag is defined as x(t+h), and lag then is -14 (sign_lags="+")
#'                   By default, and following the paper's notation, sign_lags="-"
#' @export
dep_coherence_lagged_plot <- function(x, y, max_lags=10, label_id=1, xlabel="Lag", ylabel="Correlation", color="#e4a027", corr_var_name="R", x_var_name="X", y_var_name="Y", add_corr_definition_at_title=TRUE, use_dbl_line_title=TRUE, sign_lags="-") {
  #stopifnot(length(x) == length(y))
  #stopifnot(max_lags < 0.7 * length(x))
  #stopifnot(sign_lags %in% c("+", "-"))
  
  coh_info <- coherence_corr(x, y, max_lags=max_lags, sign_lags=sign_lags)
  
  if(as.numeric(label_id)) {
    label_id <- intToUtf8(utf8ToInt("A") + label_id - 1)
  }
  
  theme <- csda_theme()
  corr_lag <- coh_info$correlation_values
  corr_indices <- coh_info$correlation_indices
  corr_maximum <- coh_info$coherence
  corr_idx_maximum <- coh_info$max_corr_lag
    
  corr_def_title <- sprintf(
    "%s) $%s(h) = \\kappa \\sum %s(t)\\cdot %s(t) (t %s h)\\,$.  ",
    label_id,
    corr_var_name,
    x_var_name,
    y_var_name,
    sign_lags
  )
  corr_val_title <- sprintf(
    "$\\textbf{C}oh\\[%s, %s\\]%s = %.3f$ at $h=%d$",
    x_var_name,
    y_var_name,
    if(add_corr_definition_at_title) sprintf(" = \\max(%s(h))^2", corr_var_name) else "",
    corr_maximum,
    corr_idx_maximum
  )
  title <- paste0(if(add_corr_definition_at_title) corr_def_title else "", corr_val_title) 
  title <- sprintf("%s) %s", label_id, title)
  if(use_dbl_line_title){
    title <- sprintf(
      "$\\overset{%s)\\,\\textbf{C}oh\\[%s, %s\\]= \\max(%s(h)^2)}{ = %.4g$ at $h=%d$}",
      label_id,
      x_var_name,
      y_var_name,
      corr_var_name,
      corr_maximum,
      corr_idx_maximum
    )
  }
  corr_plot <- ggplot2::qplot(
    corr_indices,
    corr_lag,
    geom="line",
    xlab=TeX(xlabel),
    ylab=TeX(ylabel),
    color=I(color),
    main=TeX(title)
  ) + theme
}

