channels <- c("ACAA", "ACAC", "ACAG", "ACAT", "CCAA", "CCAC", "CCAG", "CCAT", "GCAA", "GCAC", "GCAG", "GCAT", "TCAA", "TCAC", "TCAG", "TCAT", "ACGA", "ACGC", "ACGG", "ACGT", "CCGA", "CCGC", "CCGG", "CCGT", "GCGA", "GCGC", "GCGG", "GCGT", "TCGA", "TCGC", "TCGG", "TCGT", "ACTA", "ACTC", "ACTG", "ACTT", "CCTA", "CCTC", "CCTG", "CCTT", "GCTA", "GCTC", "GCTG", "GCTT", "TCTA", "TCTC", "TCTG", "TCTT", "ATAA", "ATAC", "ATAG", "ATAT", "CTAA", "CTAC", "CTAG", "CTAT", "GTAA", "GTAC", "GTAG", "GTAT", "TTAA", "TTAC", "TTAG", "TTAT", "ATCA", "ATCC", "ATCG", "ATCT", "CTCA", "CTCC", "CTCG", "CTCT", "GTCA", "GTCC", "GTCG", "GTCT", "TTCA", "TTCC", "TTCG", "TTCT", "ATGA", "ATGC", "ATGG", "ATGT", "CTGA", "CTGC", "CTGG", "CTGT", "GTGA", "GTGC", "GTGG", "GTGT", "TTGA", "TTGC", "TTGG", "TTGT")

channels_dbs <- c("ACCA", "ACCG", "ACCT", "ACGA", "ACGG", "ACGT", "ACTA", "ACTG", "ACTT", "ATCA", "ATCC", "ATCG", "ATGA", "ATGC", "ATTA", "CCAA", "CCAG", "CCAT", "CCGA", "CCGG", "CCGT", "CCTA", "CCTG", "CCTT", "CGAT", "CGGC", "CGGT", "CGTA", "CGTC", "CGTT", "CTAA", "CTAC", "CTAG", "CTGA", "CTGC", "CTGG", "CTTA", "CTTC", "CTTG", "GCAA", "GCAG", "GCAT", "GCCA", "GCCG", "GCTA", "TAAT", "TACG", "TACT", "TAGC", "TAGG", "TAGT", "TCAA", "TCAG", "TCAT", "TCCA", "TCCG", "TCCT", "TCGA", "TCGG", "TCGT", "TGAA", "TGAC", "TGAT", "TGCA", "TGCC", "TGCT", "TGGA", "TGGC", "TGGT", "TTAA", "TTAC", "TTAG", "TTCA", "TTCC", "TTCG", "TTGA", "TTGC", "TTGG")

#' Plots point estimates and credible intervals of mutational
#'
#' @param Phi An array of size V x K x S, where V is the number of categories,
#' K is the number of clusters in the best partition and S is the number of
#' iterations of the MCMC sampler
#' @param Phi A 3D array of dimension V x K x S, where
#' V is the number of categories, K is the number of clusters (signatures) and
#' S is the number of posterior samples
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

#' Same as hdp_plot_sig_uncertainty but for DBS signatures
#'
#' @export
hdp_plot_sig_dbs_uncertainty <- function(Phi, dst = NULL) {

  # Phi is an array of dimension 96 x K x S
  K <- dim(Phi)[2] # Number of signatures

  ### contains_zero <- min(Z) == 0

  for(i in 1:K) {
    signatures <- coda::as.mcmc(t(Phi[,i,]))
    colnames(signatures) <- channels_dbs
    # png(paste0(dst, "/", "signature", i, ".png"), 1200, 300)
    sig_plot <- ggmcmc::ggs_caterpillar(
      ggmcmc::ggs(signatures, keep_original_order = TRUE),
      horizontal = FALSE,
      sort = FALSE) + ggplot2::aes(color = as.factor(rep(1:6, each = 13))) +
      ggplot2::scale_colour_manual(values = c("cyan", "black", "red", "grey", "green", "pink")) +
      ggplot2::ggtitle(paste("Signature", i)) +
      ggplot2::theme(plot.title = ggplot2::element_text(size = 24), axis.text.y = ggplot2::element_text(size = 16))
    print(sig_plot)
    # dev.off()
  }
}
