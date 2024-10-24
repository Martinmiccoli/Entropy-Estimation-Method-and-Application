# Load libraries
library(gtools) 
library(ggplot2) 
library(entropy) 
library(gridExtra) 
library(ggraph)
library(igraph)
library(infotheo)
library(parmigene)
library(GeneNet)


# Shannon Entropy
shannon_entropy <- function(probabilities) {
  -sum(probabilities * log(probabilities + 1e-10))
}

# Maximum Likelihood Estimator
entropy_ml <- function(counts) {
  probabilities <- counts / sum(counts)
  shannon_entropy(probabilities)
}

# Miller-Madow Estimator
entropy_miller_madow <- function(counts) {
  n <- sum(counts)
  m_pos <- sum(counts > 0)
  entropy_ml(counts) + (m_pos - 1) / (2 * n)
}

# Bayesian Estimators with different priors
entropy_bayes <- function(counts, a) {
  n <- sum(counts)
  A <- sum(a)
  theta_bayes <- (counts + a) / (n + A)
  shannon_entropy(theta_bayes)
}

# Chao-Shen Estimator
entropy_chao_shen <- function(counts) {
  entropy.ChaoShen(counts)
}

# Shrinkage Estimator (James-Stein)
entropy_shrinkage <- function(counts) {
  entropy.shrink(counts)
}


generate_data <- function(p, sample_size, scenario) {
  if (scenario == 1) {
    alpha <- rep(0.0007, p)
    true_probs <- rdirichlet(1, alpha)
  } else if (scenario == 2) {
    alpha <- rep(1, p)
    true_probs <- rdirichlet(1, alpha)
  } else if (scenario == 3) {
    alpha <- rep(1, p/2)
    true_probs <- c(rdirichlet(1, alpha), rep(0, p/2))
    true_probs <- true_probs / sum(true_probs)
  } else if (scenario == 4) {
    zipf_probs <- 1 / (1:p)
    true_probs <- zipf_probs / sum(zipf_probs)
  } else {
    stop("Invalid scenario")
  }
  
  counts <- rmultinom(1, size = sample_size, prob = true_probs)
  list(true_probs = true_probs, counts = as.vector(counts))
}


calculate_mse_bias <- function(counts, true_probs) {
  true_entropy <- shannon_entropy(true_probs)
  
  # Entropy estimates
  entropy_estimates <- c(
    ML = entropy_ml(counts),
    MillerMadow = entropy_miller_madow(counts),
    Bayes_1 = entropy_bayes(counts, rep(1, length(counts))),
    Bayes_1_2 = entropy_bayes(counts, rep(0.5, length(counts))),
    Bayes_Minimax = entropy_bayes(counts, sqrt(counts)),
    Bayes_1_p = entropy_bayes(counts, rep(1/length(counts), length(counts))),
    Chao_Shen = entropy_chao_shen(counts),
    Shrinkage = entropy_shrinkage(counts)
  )
  
  
  mse_entropy <- (entropy_estimates - true_entropy)^2
  
  
  bias_entropy <- entropy_estimates - true_entropy
  
  return(list(mse_entropy = mse_entropy, bias_entropy = bias_entropy))
}

# Simulation Parameters
p <- 1000
sample_sizes <- c(10, 30, 100, 300, 1000, 3000, 10000)
scenarios <- 1:4
n_simulations <- 1000


results_list <- list()

for (scenario in scenarios) {
  for (sample_size in sample_sizes) {
    mse_list <- list()
    
    entropy_estimates <- matrix(NA, nrow = n_simulations, ncol = 8)
    bias_estimates <- matrix(NA, nrow = n_simulations, ncol = 8)
    colnames(entropy_estimates) <- c(
      "ML", "Miller-Madow", "Bayes_1", "Bayes_1/2", 
      "Bayes_Minimax", "Bayes_1/p", "Chao-Shen", "Shrinkage"
    )
    
    valid_simulations <- 0
    attempts <- 0
    max_attempts <- n_simulations * 2
    
    while (valid_simulations < n_simulations && attempts < max_attempts) {
      attempts <- attempts + 1
      data <- generate_data(p, sample_size, scenario)
      counts <- data$counts
      true_probs <- data$true_probs
      
      if (sum(counts > 0) >= 2) {
        valid_simulations <- valid_simulations + 1
        
        mse_bias <- calculate_mse_bias(counts, true_probs)
        
        entropy_estimates[valid_simulations, ] <- mse_bias$mse_entropy
        bias_estimates[valid_simulations, ] <- mse_bias$bias_entropy
      }
    }
    
    if (valid_simulations < n_simulations) {
      warning(paste(
        "Could not generate enough valid simulations for Scenario", scenario, 
        "and Sample Size", sample_size
      ))
      next
    }
    
    mse <- colMeans(entropy_estimates, na.rm = TRUE)
    bias <- colMeans(bias_estimates, na.rm = TRUE)
    
    results_list[[paste(
      "Scenario", scenario, "SampleSize", sample_size, sep = "_"
    )]] <- list(
      mse = mse,
      bias = bias,
      true_probs = data$true_probs,
      counts = counts
    )
  }
}

# Prepare Data for Plotting
plot_data <- data.frame()

for (name in names(results_list)) {
  parts <- strsplit(name, "_")[[1]]
  scenario <- parts[2]
  sample_size <- as.numeric(parts[4])
  mse <- results_list[[name]]$mse
  bias <- results_list[[name]]$bias
  
  temp_data <- data.frame(
    Scenario = scenario,
    SampleSize = sample_size,
    Estimator = names(mse),
    MSE = as.numeric(mse),
    Bias = as.numeric(bias)
  )
  
  plot_data <- rbind(plot_data, temp_data)
}

# Descriptive scenario labels
scenario_labels <- c(
  "Scenario 1: Dirichlet with small alpha",
  "Scenario 2: Dirichlet with alpha=1",
  "Scenario 3: Dirichlet with half zeros",
  "Scenario 4: Zipf distribution"
)

## PLOTTING FUNCTION
generate_plots <- function(scenario, plot_data) {
  scenario_label <- scenario_labels[as.numeric(scenario)]
  scenario_pattern <- paste0("Scenario_", scenario, "_")
  scenario_data <- results_list[grep(scenario_pattern, names(results_list))]
  
  if (length(scenario_data) == 0) {
    stop(paste("No data available for Scenario", scenario))
  }
  
  largest_sample_size_data <- scenario_data[[length(scenario_data)]]
  
 
  empirical_probs <- largest_sample_size_data$counts / sum(largest_sample_size_data$counts)
  

  prob_data <- data.frame(
    bin_number = 1:p,
    probability = empirical_probs
  )
  
  prob_density_plot <- ggplot(prob_data, aes(x = bin_number, y = probability)) +
    geom_bar(stat = "identity", color = "black", fill = "black") +
    labs(title = paste("Scenario", scenario, "- Empirical Probability Distribution (H =", 
                       round(shannon_entropy(empirical_probs), 2), ")"), 
         x = "Bin Number", y = "Probability") +
    scale_x_continuous(limits = c(0, p), breaks = seq(0, p, p/5)) +
    theme_minimal() +
    theme(
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      plot.title = element_text(size = 14, face = "bold")
    )
  

  line_types <- c(
    "ML" = "solid",
    "Miller-Madow" = "dashed",
    "Bayes_1" = "dashed",
    "Bayes_1/2" = "dashed",
    "Bayes_Minimax" = "dotted",
    "Bayes_1/p" = "dotted",
    "Chao-Shen" = "solid",
    "Shrinkage" = "solid"
  )
  
  # MSE Cell Frequencies Plot 
  filtered_plot_data_mse_cell_frequencies <- plot_data[!plot_data$Estimator %in% c("Miller-Madow", "Chao-Shen"), ]
  
  mse_cell_frequencies_plot <- ggplot(filtered_plot_data_mse_cell_frequencies[filtered_plot_data_mse_cell_frequencies$Scenario == scenario, ], 
                                      aes(x = SampleSize, y = MSE, color = Estimator, linetype = Estimator, group = Estimator)) +
    geom_line() +
    geom_point() +
    scale_x_log10(limits = c(10, 5000), breaks = c(10, 50, 500, 5000)) +
    scale_color_manual(values = c(
      "ML" = "red",
      "Miller-Madow" = "red",
      "Bayes_1" = "green",
      "Bayes_1/2" = "green",
      "Bayes_Minimax" = "blue",
      "Bayes_1/p" = "blue",
      "Chao-Shen" = "yellow",
      "Shrinkage" = "black"
    )) +
    scale_linetype_manual(values = line_types) +
    labs(title = paste("Scenario", scenario, "- MSE Cell Frequencies"), 
         x = "Sample Size (log scale)", y = "MSE") +
    theme_minimal()
  
  # MSE Entropy Plot 
  mse_entropy_plot <- ggplot(plot_data[plot_data$Scenario == scenario, ], 
                             aes(x = SampleSize, y = MSE, color = Estimator, linetype = Estimator, group = Estimator)) +
    geom_line() +
    geom_point() +
    scale_x_log10(limits = c(10, 5000), breaks = c(10, 50, 500, 5000)) +
    scale_color_manual(values = c(
      "ML" = "red",
      "Miller-Madow" = "red",
      "Bayes_1" = "green",
      "Bayes_1/2" = "green",
      "Bayes_Minimax" = "blue",
      "Bayes_1/p" = "blue",
      "Chao-Shen" = "yellow",
      "Shrinkage" = "black"
    )) +
    scale_linetype_manual(values = line_types) +
    labs(title = paste("Scenario", scenario, "- MSE Entropy"), 
         x = "Sample Size (log scale)", y = "MSE") +
    theme_minimal()
  
  # Bias Entropy Plot
  bias_entropy_plot <- ggplot(plot_data[plot_data$Scenario == scenario, ], 
                              aes(x = SampleSize, y = Bias, color = Estimator, linetype = Estimator, group = Estimator)) +
    geom_line() +
    geom_point() +
    scale_x_log10(limits = c(10, 5000), breaks = c(10, 50, 500, 5000)) +
    scale_color_manual(values = c(
      "ML" = "red",
      "Miller-Madow" = "red",
      "Bayes_1" = "green",
      "Bayes_1/2" = "green",
      "Bayes_Minimax" = "blue",
      "Bayes_1/p" = "blue",
      "Chao-Shen" = "yellow",
      "Shrinkage" = "black"
    )) +
    scale_linetype_manual(values = line_types) +
    labs(title = paste("Scenario", scenario, "- Bias Entropy"), 
         x = "Sample Size (log scale)", y = "Bias") +
    theme_minimal()
  
  return(list(
    prob_density_plot, 
    mse_cell_frequencies_plot, 
    mse_entropy_plot, 
    bias_entropy_plot
  ))
}


for (scenario in scenarios) {
  plots <- generate_plots(scenario, plot_data)
  grid.arrange(grobs = plots, ncol = 2)
}

####################
####################
####################


data(ecoli)
ecoli_data <- ecoli

# Function for mutual information 
mutual_information <- function(x, y, num_bins = 16) {
  if (!is.numeric(x) || !is.numeric(y)) {
    stop("Inputs x and y must be numeric")
  }
  
  x_disc <- cut(x, breaks = num_bins, labels = FALSE)
  y_disc <- cut(y, breaks = num_bins, labels = FALSE)
  
  mi <- mutinformation(x_disc, y_disc, method = "emp")
  
  return(mi)
}

discretize_data <- function(data, num_bins = 16) {
  apply(data, 2, function(x) cut(x, breaks = num_bins, labels = FALSE))
}

discretized_data <- discretize_data(ecoli_data)

num_genes <- ncol(discretized_data)

mi_matrix <- matrix(0, nrow = num_genes, ncol = num_genes)
colnames(mi_matrix) <- rownames(mi_matrix) <- colnames(discretized_data)

for (i in 1:(num_genes - 1)) {
  for (j in (i + 1):num_genes) {
    mi <- mutual_information(discretized_data[, i], discretized_data[, j])
    mi_matrix[i, j] <- mi
    mi_matrix[j, i] <- mi
  }
}

aracne_matrix <- aracne.a(mi_matrix)

mi_values <- mi_matrix[upper.tri(mi_matrix)]
aracne_values <- aracne_matrix[upper.tri(aracne_matrix)]

# Create histograms
p1 <- ggplot(data.frame(mi = mi_values), aes(x = mi)) +
  geom_histogram(bins = 30, fill = "white", color = "black") +
  labs(title = "MI Empirical Estimates", x = "Mutual Information", y = "Frequency") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

p2 <- ggplot(data.frame(mi = aracne_values), aes(x = mi)) +
  geom_histogram(bins = 30, fill = "white", color = "black") +
  labs(title = "ARACNE-Processed MIs", x = "Mutual Information", y = "Frequency") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

combined_plot <- grid.arrange(p1, p2, ncol = 2)

ggsave("mi_aracne_comparison.png", combined_plot, width = 10, height = 5)

##
##
##

shrinkage_entropy <- function(counts) {
  n <- sum(counts)
  p <- length(counts)
  
  ML <- counts / n
  target <- rep(1/p, p)
  
  var_ML <- ML * (1 - ML) / n
  lambda <- sum(var_ML) / sum((ML - target)^2)
  lambda <- min(1, lambda)
  
  p_shrink <- lambda * target + (1 - lambda) * ML
  
  H <- -sum(p_shrink * log2(p_shrink + 1e-10)) 
  
  return(list(entropy = H, lambda = lambda, probs = p_shrink))
}

mutual_information <- function(x, y, num_bins = 16) {
  if (!is.numeric(x) || !is.numeric(y)) {
    stop("Inputs x and y must be numeric")
  }
  
  x_disc <- cut(x, breaks = num_bins, labels = FALSE)
  y_disc <- cut(y, breaks = num_bins, labels = FALSE)
  
  joint_counts <- table(x_disc, y_disc)
  x_counts <- rowSums(joint_counts)
  y_counts <- colSums(joint_counts)
  
  H_x <- shrinkage_entropy(x_counts)$entropy
  H_y <- shrinkage_entropy(y_counts)$entropy
  H_xy <- shrinkage_entropy(as.vector(joint_counts))$entropy
  
  MI <- H_x + H_y - H_xy
  
  return(MI)
}

discretize_data <- function(data, num_bins = 16) {
  apply(data, 2, function(x) cut(x, breaks = num_bins, labels = FALSE))
}

discretized_data <- discretize_data(ecoli_data)

num_genes <- ncol(discretized_data)
mi_matrix <- matrix(0, nrow = num_genes, ncol = num_genes)
colnames(mi_matrix) <- rownames(mi_matrix) <- colnames(discretized_data)

for (i in 1:(num_genes-1)) {
  for (j in (i+1):num_genes) {
    mi <- mutual_information(discretized_data[,i], discretized_data[,j])
    mi_matrix[i,j] <- mi_matrix[j,i] <- mi
  }
}

# ARACNE Implementation
aracne <- function(mi_matrix, threshold) {
  n <- nrow(mi_matrix)
  for (i in 1:(n-2)) {
    for (j in (i+1):(n-1)) {
      for (k in (j+1):n) {
        if (mi_matrix[i,j] > threshold && mi_matrix[i,k] > threshold && mi_matrix[j,k] > threshold) {
          min_edge <- min(mi_matrix[i,j], mi_matrix[i,k], mi_matrix[j,k])
          if (min_edge == mi_matrix[i,j]) {
            mi_matrix[i,j] <- mi_matrix[j,i] <- 0
          } else if (min_edge == mi_matrix[i,k]) {
            mi_matrix[i,k] <- mi_matrix[k,i] <- 0
          } else {
            mi_matrix[j,k] <- mi_matrix[k,j] <- 0
          }
        }
      }
    }
  }
  return(mi_matrix)
}

# Apply ARACNE
threshold <- quantile(mi_matrix[upper.tri(mi_matrix)], 0.90)  
aracne_matrix <- aracne(mi_matrix, threshold)

# Create graph from adjacency matrix
g <- graph_from_adjacency_matrix(aracne_matrix, mode = "undirected", weighted = TRUE)

g <- simplify(g, remove.multiple = FALSE, remove.loops = TRUE)
weight_threshold <- 0.1  
g <- delete_edges(g, E(g)[which(E(g)$weight < weight_threshold)])


g <- delete_vertices(g, V(g)[degree(g) == 0])


largest_component <- induced_subgraph(g, which(components(g)$membership == which.max(components(g)$csize)))

# Create the plot object
plot_obj <- ggraph(largest_component, layout = "fr", niter = 1000) +  
  geom_edge_link(color = "black", alpha = 0.5, size = 0.3) +
  geom_node_point(color = "black", size = 12, stroke = 0.5, shape = 21, fill = "white") + 
  geom_node_text(aes(label = name), size = 3, color = "black") +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(hjust = 0.5, size = 14)
  ) +
  guides(color = "none", size = "none") +
  labs(title = "E. coli Gene Interaction Network") +
  coord_fixed(ratio = 1) +  
  scale_x_continuous(expand = expansion(mult = 0.2)) +  
  scale_y_continuous(expand = expansion(mult = 0.2))


pdf("ecoli_gene_network.pdf", width = 20, height = 20, family = "sans")
print(plot_obj)
dev.off()


cat("Number of nodes:", vcount(largest_component), "\n")
cat("Number of edges:", ecount(largest_component), "\n")
cat("Top 5 hub genes:", names(sort(degree(largest_component), decreasing=TRUE)[1:5]), "\n")
