library(hdpExtra)
library(hdp)
library(coda)
library(salso)
library(lsa)
library(MCMCpack)

M <- read.csv("data/WGS_PCAWG.96.csv", header = TRUE, row.names = 1)
M <- M[, startsWith(colnames(M), "CNS.GBM")] # Keep glioblastoma samples only
M <- as.data.frame(t(M))
# initialise HDP
hdp_mut <- hdp_init(ppindex = c(0, rep(1, each = nrow(M))), # index of parental nodes
                    cpindex = c(1, rep(2, each = nrow(M))), # index of concentration param
                    hh=rep(1, ncol(M)), # prior is uniform over the 96 mutation categories
                    alphaa=rep(1, 2), # shape hyperparams for five different CPs
                    alphab=rep(1, 2)) # rate hyperparams for five different CPs

# add data to leaf nodes (one per cancer sample, in row order of mut_count)
hdp_mut <- hdp_setdata(hdp_mut,
                       dpindex = 2:numdp(hdp_mut), # index of nodes to add data to
                       M) # input data (mutation counts, sample rows match up with specified dpindex)

hdp_mut
hdp_extra_chains <- vector("list", 2)
burnin <- 100
numiter <- 100
space <- 1
for (i in 1:2){

  # activate DPs, 10 initial components
  hdp_activated <- dp_activate(hdp_mut, 1:numdp(hdp_mut), initcc=20, seed=i*200)

  hdp_extra_chains[[i]] <- hdpExtra_posterior(hdp_activated,
                               burnin=burnin,
                               n=numiter,
                               space=space,
                               cpiter=3,
                               seed=i*1e3)
}
hdp_extra_chains <- HdpExtraChainMulti(hdp_extra_chains)

Sys.time()
allocations_best <- salso(t(hdp_allocations(hdp_extra_chains)), loss = VI(a = 2/3))
Sys.time()
Phi_best <- hdp_postprocessing(hdp_extra_chains, allocations_best)
Sys.time()
hdp_plot_sig_uncertainty(Phi_best, "plots")

Sys.time()
allocations_best_binder <- salso(t(hdp_allocations(hdp_extra_chains)), loss = binder(a = 4/3))
Sys.time()
Phi_best_binder <- hdp_postprocessing(hdp_extra_chains, allocations_best_binder)
Sys.time()
hdp_plot_sig_uncertainty(Phi_best_binder, "plots")

sort(table(allocations_best_binder))

hdpExtra_diagnostic_plots(hdp_extra_chains, c("gamma", "alpha"))

which.max.cosine <- function(sig, sigs_matrix) {
  N <- ncol(sigs_matrix)
  similarities <- apply(sigs_matrix, 2, function(x) cosine(x, sig))
  if (max(similarities) > 0.85) names(which.max(similarities))
  else "other"
}

tolerance <- 0.001
cosmic_signatures <- read.table("data/cosmic_signatures.txt", header = TRUE, row.names = 1)
#S <- length(hdp_extra_chains@Beta)
S <- 100
weights_matrix <- matrix(0, nrow = ncol(cosmic_signatures) + 1, ncol = S)
rownames(weights_matrix) <- c(colnames(cosmic_signatures), "other")
dp_pi <- list()
for (s in 1:S) {
  Phi <- hdp_extra_chains@Phi[[s]]
  Beta <- hdp_extra_chains@Beta[[s]]
  K <- ncol(Phi)
  dp_alpha <- hdp_extra_chains@cp_values[[1]][burnin + (s-1) * space + 1, 2]
  dp_pi[[s]] <- rdirichlet(1, dp_alpha*Beta)
  for (k in 1:K) {
    which.signature <- which.max.cosine(Phi[, k], cosmic_signatures)
    weights_matrix[which.signature, s] <- weights_matrix[which.signature, s] + dp_pi[[s]][k]
  }
}
weights_matrix <- t(weights_matrix)
library(ggmcmc)
ggs_density(ggs(as.mcmc(weights_matrix[, c("SBS1", "SBS5", "SBS8", "SBS11", "SBS30", "SBS3", "other")])))


# 0.5 SBS1 + 0.3 SBS11 + 0.2 SBS30
simulated_weights <- c(0.5, 0.3, 0.2)
dim(as.matrix(simulated_weights))
simulated_patient <- as.matrix(cosmic_signatures[, c("SBS1", "SBS11", "SBS30")]) %*% as.matrix(simulated_weights)
simulated_patient <- round(1000 * simulated_patient)
simulated_patient
simulated_patient <- rep(1:96, simulated_patient)

setClass("DPPosteriorDraw",
  slots = c(
    Beta = "numeric",
    Phi = "matrix",
    dp_gamma = "numeric"
))

rg0 <- function(dpPosteriorDraw) {
  Beta <- dpPosteriorDraw@Beta
  k <- sample(1:length(Beta), size = 1, prob = Beta)
  if (k == length(Beta)) Phi[, k]
  else rdirichlet(1, rep(1, 96))
}

m <- 3
# for s = 1, ..., S do
s <- 1
alpha <- 1
g0 <- new(
  "DPPosteriorDraw",
  Beta = hdp_extra_chains@Beta[[s]],
  Phi = hdp_extra_chains@Phi[[s]],
  dp_gamma = hdp_extra_chains@cp_values[[1]][burnin + (s-1) * space + 1, 2]
)
n <- length(simulated_patient)
z <- sample(1:20, size = n, replace = TRUE)
Phi <- rdirichlet(20, rep(1, 96))
cluster_sizes <- sapply(1:20, function(k) sum(z == k))
# for i = 1, ..., n do
i <- 1
current_label <- z[i]
cluster_sizes[current_label] <- cluster_sizes[current_label] - 1
k.minus <- length(cluster_sizes)
h <- k.minus + m
if (cluster_sizes[current_label] > 0) {
  for (k in (k.minus+1):h) {
    Phi <- rbind(Phi, rg0(g0))
  }
} else {
  z[i] <- k.minus + 1
  cluster_sizes[k.minus + 1] <- 1
  for (k in (k.minus+2):h) {
    Phi <- rbind(Phi, rg0(g0))
  }
}
probs <- c(cluster_sizes, rep(alpha/m, m))
probs <- probs * Phi[, simulated_patient[i]]
z[i] <- sample(1:h, size = 1, prob = probs)
cluster_sizes[z[i]] <- cluster_sizes[z[i]] + 1
Phi <- Phi[which(cluster_sizes > 0), ]
for (k in 1:h) {
    if (cluster_sizes)
}
}
