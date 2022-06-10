channels <- c("ACAA", "ACAC", "ACAG", "ACAT", "CCAA", "CCAC", "CCAG", "CCAT", "GCAA", "GCAC", "GCAG", "GCAT", "TCAA", "TCAC", "TCAG", "TCAT", "ACGA", "ACGC", "ACGG", "ACGT", "CCGA", "CCGC", "CCGG", "CCGT", "GCGA", "GCGC", "GCGG", "GCGT", "TCGA", "TCGC", "TCGG", "TCGT", "ACTA", "ACTC", "ACTG", "ACTT", "CCTA", "CCTC", "CCTG", "CCTT", "GCTA", "GCTC", "GCTG", "GCTT", "TCTA", "TCTC", "TCTG", "TCTT", "ATAA", "ATAC", "ATAG", "ATAT", "CTAA", "CTAC", "CTAG", "CTAT", "GTAA", "GTAC", "GTAG", "GTAT", "TTAA", "TTAC", "TTAG", "TTAT", "ATCA", "ATCC", "ATCG", "ATCT", "CTCA", "CTCC", "CTCG", "CTCT", "GTCA", "GTCC", "GTCG", "GTCT", "TTCA", "TTCC", "TTCG", "TTCT", "ATGA", "ATGC", "ATGG", "ATGT", "CTGA", "CTGC", "CTGG", "CTGT", "GTGA", "GTGC", "GTGG", "GTGT", "TTGA", "TTGC", "TTGG", "TTGT")

#' Plots point estimates and credible intervals of mutational
#'
#' @param Phi A 3D array of dimension V x K x S, where
#' V is the number of categories, K is the number of clusters (signatures) and
#' S is the number of posterior samples
#' @param dst The path of a directory where the signatures will be stored the allocation variables of N observations.
#'
#' @export
hdp_plot_sig_uncertainty <- function(Phi, dst = NULL) {

  # Phi is an array of dimension 96 x K x S
  K <- dim(Phi)[2] # Number of signatures

  ### contains_zero <- min(Z) == 0

  for(i in 1:K) {
    signatures <- coda::as.mcmc(t(Phi[,i,]))
    colnames(signatures) <- channels
    # png(paste0(dst, "/", "signature", i, ".png"), 1200, 300)
    sig_plot <- ggmcmc::ggs_caterpillar(
      ggmcmc::ggs(signatures, keep_original_order = TRUE),
      horizontal = FALSE,
      sort = FALSE) + ggplot2::aes(color = as.factor(rep(1:6, each = 16))) +
      ggplot2::scale_colour_manual(values = c("cyan", "black", "red", "grey", "green", "pink")) +
      ggplot2::ggtitle(paste("Signature", i)) +
      ggplot2::theme(plot.title = ggplot2::element_text(size = 24), axis.text.y = ggplot2::element_text(size = 16))
    print(sig_plot)
    # dev.off()
  }
}

