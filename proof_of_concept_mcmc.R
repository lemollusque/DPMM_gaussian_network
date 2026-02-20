library(dirichletprocess)
library(gRbase)
library(BiDAG)
library(matrixStats)
library(cowplot)
library(pcalg)
library(rstan)
library(bridgesampling)

source("fns.R")
source("Fourier_fns.R")
source("BayesStanFns.R")
source("sampling_fns.R")

set.seed(101)
lambda <- 1
n <- 4  # number of nodes
N <- 100  # number of samples

myDAG <- pcalg::randomDAG(n, prob = 0.5, lB = 1, uB = 2) 
trueDAG <- as(myDAG, "matrix")
truegraph <- 1*(trueDAG != 0)
data <- Fou_nldata(truegraph, N, lambda = lambda, noise.sd = 1, standardize = T) 

vars  <- c("x1","x2","x3","x4")
colnames(data) = vars

#----------------------  overwrite functions ----------------------------------
# replace BIDAG functions
unlockBinding("usrscoreparameters", asNamespace("BiDAG"))
assign("usrscoreparameters", usrscoreparameters, envir = asNamespace("BiDAG"))
lockBinding("usrscoreparameters", asNamespace("BiDAG"))

unlockBinding("usrDAGcorescore", asNamespace("BiDAG"))
assign("usrDAGcorescore", usrDAGcorescore, envir = asNamespace("BiDAG"))
lockBinding("usrDAGcorescore", asNamespace("BiDAG"))
#----------------------  end overwrite functions ------------------------------
# Initiate params for DP and BGe
n <- ncol(data)
alpha_mu <- 1          
alpha_w  <- n + alpha_mu + 1      
t <- alpha_mu * (alpha_w - n - 1) / (alpha_mu + 1)

# DP
g0Priors <- list(
  mu0    = rep(0, n),
  Lambda = diag(n) / t,   # T = (1/t) I
  kappa0 = alpha_mu,
  nu     = alpha_w
)

scaled_data = scale(data) 
n_iter = 500
burnin = 300
L = 50 # sample to take

Gamma_list <- list()
vars  <- c("x1","x2","x3","x4")
for (child in vars){
  parents <- vars[vars != child]
  dp_data = scaled_data[,c(child, parents)]
  dp <-  DirichletProcessMvnormal(dp_data, g0Priors)
  dp <- Fit(dp, n_iter)
  
  Gamma_sample <- dp_membership_probs(dp, n_iter, burnin, L)
  Gamma_list <- add_membershipp(Gamma_list, Gamma_sample, child=child, parents=parents)
}

# scoring
usr_score_param <- BiDAG::scoreparameters(scoretype = "usr", 
                                          data = scaled_data, 
                                          usrpar = list(pctesttype = "bge",
                                                        membershipp_list = Gamma_list,
                                                        am = alpha_mu, 
                                                        aw = alpha_w, 
                                                        T0scale = t,
                                                        edgepf = 1
                                          )
)
bge_score_param <- scoreparameters("bge", 
                                   scaled_data, 
                                   bgepar = list(am = alpha_mu, 
                                                 aw = alpha_w, 
                                                 edgepf = 1)
)
#----------------------------------- posterior -------------------------------
all_dags <- list()
adj <- matrix(0, n, n)
dag.counter <- 0

all.comb <- rep(list(c(0,1)), n*(n-1))
all.comb <- expand.grid(all.comb)

for(i in 1:nrow(all.comb)) {
  adj[col(adj)!=row(adj)] <- as.numeric(all.comb[i, ])
  if(is.DAG(adj)) {
    dag.counter <- dag.counter + 1
    all_dags[[dag.counter]] <- adj
  }
}

true_score_usr <- sapply(all_dags, function(dag) BiDAG::DAGscore(usr_score_param, dag))
true_score_bge <- sapply(all_dags, function(dag) BiDAG::DAGscore(bge_score_param, dag))

true_p_usr <- exp(true_score_usr - matrixStats::logSumExp(true_score_usr))
true_p_bge <- exp(true_score_bge - matrixStats::logSumExp(true_score_bge))

# searchspace
searchspace_usr <- BiDAG::iterativeMCMC(scorepar = usr_score_param, scoreout = TRUE)
searchspace_bge <- BiDAG::iterativeMCMC(scorepar = bge_score_param, scoreout = TRUE)

# MCMC
toburn <- 250
iters <- c(30e1, 44e1, 65e1, 96e1, 14e2, 21e2, 31e2, 46e2, 67e2, 10e3)

results <- data.frame()
for(i in 1:length(iters)) {
  print(i)
  fit_usr <- BiDAG::partitionMCMC(usr_score_param, 
                                  alpha=0.2, 
                                  startDAG = searchspace_usr$DAG,
                                  scoretable = searchspace_usr$scoretable, 
                                  startspace = searchspace_usr$endspace,
                                  iterations = 2*(iters[i] + toburn), stepsave = 2)
  
  fit_bge <- BiDAG::partitionMCMC(bge_score_param, 
                                  alpha=0.2, 
                                  startDAG = searchspace_bge$DAG,
                                  scoretable = searchspace_bge$scoretable, 
                                  startspace = searchspace_bge$endspace,
                                  iterations = 2*(iters[i] + toburn), stepsave = 2)
  
  sampled_usr <- fit_usr$traceadd$incidence[-(1:toburn)]
  sampled_bge <- fit_bge$traceadd$incidence[-(1:toburn)]
  
  all_vecdags <- lapply(all_dags, c)
  
  post_indexes_usr <- sapply(sampled_usr, function(x)
    which(all_vecdags %in% list(as.numeric(x))))
  
  post_indexes_bge <- sapply(sampled_bge, function(x)
    which(all_vecdags %in% list(as.numeric(x))))
  
  tab_usr <- table(post_indexes_usr)
  p_hat_usr <- as.numeric(tab_usr) / sum(tab_usr)
  ind_usr <- as.numeric(names(tab_usr))
  
  tab_bge <- table(post_indexes_bge)
  p_hat_bge <- as.numeric(tab_bge) / sum(tab_bge)
  ind_bge <- as.numeric(names(tab_bge))
  
  
  KL_div <- function(est.post, est.ind, p_true) {
    q <- rep(NA, length(p_true))
    q[est.ind] <- est.post
    d <- q * log(q / p_true)
    sum(d, na.rm = TRUE)
  }
  
  kl_usr <- KL_div(p_hat_usr, ind_usr, true_p_usr)
  kl_bge <- KL_div(p_hat_bge, ind_bge, true_p_bge)
  
  
  # Save results
  results <- rbind(results, data.frame(kl = kl_usr, 
                                       iter = iters[i], 
                                       group = "usr"))
  results <- rbind(results, data.frame(kl = kl_bge, 
                                       iter = iters[i], 
                                       group = "BGe"))
}


# -----------------------------
# Plot KL divergences (usr vs BGe)
# -----------------------------
post_plots <- list()
color2 <- c("usr" = "#00acc7", "BGe" = "#db0000")
shape2 <- c("usr" = 4, "BGe" = 3)
linetype2 <- c("usr" = "dashed", "BGe" = "dotdash")

post_plots[[1]] <- ggplot(results, aes(x = iter, y = kl, group = group)) +
  geom_line(aes(linetype = group)) +
  geom_point(aes(shape = group, color = group)) +
  scale_shape_manual(values = shape2) +
  scale_color_manual(values = color2) +
  scale_linetype_manual(values = linetype2) +
  scale_x_continuous(trans = "log10") +
  xlab("Samples") + ylab("Reverse K-L divergence") +
  theme_light() +
  theme(legend.title = element_blank(), legend.position = "none")

# -----------------------------
# Posterior overlay (last setting)
# Choose what "true.p" means:
#   Option A: compare both estimates to BGe truth:
#     true.p <- true_p_bge
#   Option B: compare each to its own truth (not one bar plot):
#     (harder to show in one bar plot)
#
# Here I default to BGe truth in the bars (common baseline),
# and plot both usr and BGe estimates as points.
# -----------------------------
true.p <- true_p_bge

post_data <- rbind(
  data.frame(x = ind_usr, y = p_hat_usr, group = "usr"),
  data.frame(x = ind_bge, y = p_hat_bge, group = "BGe")
)

post_main <- ggplot() +
  geom_bar(
    data = data.frame(x = 1:dag.counter, y = true.p),
    aes(x, y),
    stat = "identity",
    width = 0.43,
    fill = "#5e5e5e"
  ) +
  geom_point(
    data = post_data,
    aes(x, y, shape = group, color = group),
    size = 1
  ) +
  scale_shape_manual(values = shape2) +
  scale_color_manual(values = color2) +
  xlab("DAG") + ylab("Posterior probability") +
  theme_light() +
  xlim(1, min(200, dag.counter)) +
  theme(legend.title = element_blank(), legend.position = "none") +
  scale_y_sqrt()

# Inset for the tail (201:dag.counter), only if there is a tail
if (dag.counter > 200) {
  post_inset <- ggplot() +
    geom_bar(
      data = data.frame(x = 1:dag.counter, y = true.p),
      aes(x, y),
      stat = "identity",
      width = 0.43,
      fill = "#5e5e5e"
    ) +
    geom_point(
      data = post_data,
      aes(x, y, shape = group, color = group),
      size = 1
    ) +
    scale_shape_manual(values = shape2) +
    scale_color_manual(values = color2) +
    xlab("DAG") + ylab("Posterior probability") +
    theme_light() +
    xlim(201, dag.counter) +
    theme(
      legend.title = element_blank(),
      legend.position = "none",
      panel.grid.major = element_blank(),
      axis.title = element_text(size = 8),
      axis.text  = element_text(size = 8)
    ) +
    scale_y_sqrt(limits = c(0, true.p[201]))
  
  post_plots[[2]] <- ggdraw() +
    draw_plot(post_main) +
    draw_plot(post_inset, x = 0.475, y = 0.45, width = 0.5, height = 0.5)
} else {
  # If <=200 DAGs, just show the main plot without inset
  post_plots[[2]] <- post_main
}

ggarrange(
  post_plots[[1]],
  post_plots[[2]],
  ncol = 2,
  widths = c(1, 2),
  common.legend = TRUE,
  legend = "bottom"
)
