#######################
#' @export
run_tests <- function() {
  run_example_1()
  
  q <- run_example_2()
  csda_save_image(q, "coh-latent-delta-beta", height=11*3)
    
  q <- run_example_3()
  csda_save_image(q, "coh-lagged-latent-delta-beta", height=11*3)
  
  q <- run_example_4()
  csda_save_image(q, "time-sample-autospectrum-AR2")
  
  q <- run_example_4b()
  csda_save_image(q, "time-autospectrum-AR2")
  
  q <- run_example_4b()
  csda_save_image(q, "true-autospectrum-AR2", height=6)
  
  q <- run_example_5()
  csda_save_image(q, "coh-latent-3-comp", height=11*3, width=40)
  
  q <- run_example_5b()
  csda_save_image(q, "comparison-coh-pcoh", height=15)

  q <- run_example_5c()
  csda_save_image(q, "comparison-coh-pcoh-pdcoh", height=15 * 2.5)

  q <- run_example_5d()
  csda_save_image(q, "comparison-coh-pcoh-2", height=15)

  q <- run_example_6()
  csda_save_image(q, "phase-amplitude-coupling", height=15 * 2.5)
  
  q <- run_example_7()
  csda_save_image(q, "comparison-coh-pcoh-pdcoh-VAR", height=15 * 2.5)

  q <- run_example_8()
  csda_save_image(q, "filter-causal-noncausal", height=15)
  
  run_example_9()
  
}
