library(tidyverse)
library(gridExtra)
library(plot3D)
library(ggrepel)

makeModelMat <- function(X, interaction) {
  N <- nrow(X)
  k <- ncol(X)
  if(k %in% c(1, 2, 5) & interaction == 0){
    F_mat <- cbind(rep(1, N), X)
  } else if (k == 1 & interaction == 2){
    F_mat <- cbind(rep(1, N), X, X^2)
  } else if (k == 2 & interaction == 1){
    F_mat <- cbind(rep(1, N), X, X[, 1] * X[, 2])
  } else if (k == 3) {
    F_mat <- cbind(rep(1, N), X, X[,1] * X[,2], X[,1] * X[,3], X[,2] * X[,3], X^2)
  } else if (k == 5 & interaction == 1){
    F_mat <- cbind(rep(1, N), X,
                   X[, 1] * X[, 2],
                   X[, 1] * X[, 3],
                   X[, 1] * X[, 4],
                   X[, 1] * X[, 5],
                   X[, 2] * X[, 3],
                   X[, 2] * X[, 4],
                   X[, 2] * X[, 5],
                   X[, 3] * X[, 4],
                   X[, 3] * X[, 5],
                   X[, 4] * X[, 5])
  } else if (k ==2 & interaction == 2){
    F_mat <- cbind(rep(1, N), X, X[, 1] * X[, 2], X^2)
  }
  return(F_mat)
}

K_1 <- matrix(c(1,    0,    1/3,
                0,    1/3,  0,
                1/3,  0,    1/5), 3, 3)

K_2 <- matrix(c(1, 0, 0, 0, (1/3), (1/3), 
                0, (1/3), 0, 0, 0, 0,
                0, 0, (1/3), 0, 0, 0,
                0, 0, 0,  (1/9), 0, 0,
                (1/3), 0, 0, 0, (1/5), (1/9),
                (1/3), 0, 0, 0, (1/9), (1/5)), 6, 6)

K_3 <- matrix(c(1,     0,     0,     0,     0,     0,     0,  (1/3),   (1/3),   (1/3),
                0, (1/3),     0,     0,     0,     0,     0,      0,      0,        0,
                0,     0, (1/3),     0,     0,     0,     0,      0,      0,        0,
                0,     0,     0, (1/3),     0,     0,     0,      0,      0,        0,
                0,     0,     0,     0, (1/9),     0,     0,      0,      0,        0,
                0,     0,     0,     0,     0, (1/9),     0,      0,      0,        0,
                0,     0,     0,     0,     0,     0, (1/9),      0,      0,        0,
                (1/3), 0,     0,     0,     0,     0,     0,  (1/5),  (1/9),    (1/9),
                (1/3), 0,     0,     0,     0,     0,     0,  (1/9),  (1/5),    (1/9),
                (1/3), 0,     0,     0,     0,     0,     0,  (1/9),  (1/9),    (1/5)),
              nrow = 10, byrow = TRUE)

###############
# I Criterion #
###############

I_crit <- function(X, interaction, k) {
  mat <- makeModelMat(X, interaction)
  det.inf <- det(t(mat)%*%mat)
  N <- nrow(X)
  M <- t(mat) %*% mat
  rank.inf <- qr(M)$rank
  if (rank.inf < ncol(X)){
    I.score <- 0
  } else{
    eff <- 1 / (sum(diag(solve(M, k))))
    I.score <- round(eff, 4)
  }
  return(I.score)
}

############
# Robust I #
############

rb_I <- function(X, interaction, k) {
  F_mat <- makeModelMat(X, interaction)
  N <- nrow(F_mat) - 1
  rank.inf <- qr(t(F_mat)%*%F_mat)$rank
  if (rank.inf < ncol(X)){
    I.score <- 0
  } else{
    eff <- vector("numeric", length = nrow(X))
    for (i in seq_len(nrow(X))) {
      mat <- F_mat[-i,]
      rank.inf2 <- qr(t(mat)%*%mat)$rank
      if(rank.inf2 < (ncol(F_mat))){
        eff[i] <- -Inf
      } else {
        M <- t(mat) %*% mat
        eff[i] <- 1 / (sum(diag(solve(M, k))))
      }
    }
  }
  EFF <- round(min(eff), 3)
  return(EFF)
}

# returns a list of drop-one scores
rb_I_l <- function(X, interaction, k) {
  F_mat <- makeModelMat(X, interaction)
  N <- nrow(F_mat) - 1
  rank.inf <- qr(t(F_mat)%*%F_mat)$rank
  if (rank.inf < ncol(X)){
    I.score <- 0
  } else{
    eff <- vector("numeric", length = nrow(X))
    for (i in seq_len(nrow(X))) {
      mat <- F_mat[-i,]
      rank.inf2 <- qr(t(mat)%*%mat)$rank
      if(rank.inf2 < (ncol(F_mat))){
        eff[i] <- 0
      } else {
        M <- t(mat) %*% mat
        eff[i] <- 1 / (sum(diag(solve(M, k))))
      }
    }
  }
  return(eff)
}

###############
# G Criterion #
###############

G_crit <- function(X, N, K, order){
  Xm <- makeModelMat(X, 2)
  p <- ncol(Xm)
  XpX <- t(Xm) %*% Xm
  det.inf <- qr(XpX)$rank
  if (K == 1) {
    Xpred <- matrix(c(1, 1, 1, 1, 1,
                      -1, -0.5, 0.0, 0.5, 1.0,
                      1, 0.25, 0.0, 0.25, 1.0), 5, 3)
  } else if (K == 2){
    part_1 <- matrix(c(rep(1, 25),
                       rep(c(-1, -.5, 0, .5, 1), 5),
                       c(rep(-1, 5), rep(-.5, 5), rep(0, 5), rep(.5, 5), rep(1, 5))
    ),25 ,3)
    part_2 <- matrix(c(part_1[,2] * part_1[,3], part_1[,2]^2, part_1[,3]^2), 25, 3)
    Xpred <- matrix(c(part_1, part_2), 25, 6)
  } else if (K == 3) {
    part_1 <- rep(1, 125)
    x1 <- c(-1, -0.5, 0.0, 0.5, 1.0)
    x2 <- expand.grid(x1, x1, x1)
    x3 <- cbind(x2[,1] * x2[,2], x2[,1] * x2[,3], x2[,2] * x2[,3], x2[,1]^2, x2[,2]^2, x2[,3]^2)
    Xpred <- as.matrix(cbind(part_1, x2, x3), 125, 10)
  }
  
  if (det.inf < K) {
    G_score <- -Inf
  } else {
    C <- chol(XpX, pivot = TRUE)
    C[lower.tri(C, diag = FALSE)] <- 0
    Z <- solve(t(C), t(Xpred))
    Tv <- Z^2
    D <- colSums(Tv)
    Mx <- max(D)
    result <- (p)/ Mx
  }
  return(round(result, 4))
}

############
# Robust G #
############

rb_G <- function(X, N, K, order){
  Xm <- makeModelMat(X, 2)
  p <- ncol(Xm)
  XpX <- t(Xm) %*% Xm
  det.inf <- qr(XpX)$rank
  if (det.inf < K) {
    G_score <- - Inf
  } else {
    sc <- vector(mode = 'numeric', length = N)
    if (K == 1) {
      Xpred <- matrix(c(1, 1, 1, 1, 1,
                        -1, -0.5, 0.0, 0.5, 1.0,
                        1, 0.25, 0.0, 0.25, 1.0), 5, 3)
    } else if (K == 2){
      part_1 <- matrix(c(rep(1, 25),
                         rep(c(-1, -.5, 0, .5, 1), 5),
                         c(rep(-1, 5), rep(-.5, 5), rep(0, 5), rep(.5, 5), rep(1, 5))),25 ,3)
      part_2 <- matrix(c(part_1[,2] * part_1[,3], part_1[,2]^2, part_1[,3]^2), 25, 3)
      Xpred <- matrix(c(part_1, part_2), 25, 6)
    } else if (K == 3) {
      part_1 <- rep(1, 1e+06)
      x1 <- c(seq(from = -1, to = 1, length.out = 100))
      x2 <- expand.grid(x1, x1, x1)
      x3 <- cbind(x2[,1] * x2[,2], x2[,1] * x2[,3], x2[,2] * x2[,3], x2[,1]^2, x2[,2]^2, x2[,3]^2)
      Xpred <- as.matrix(cbind(part_1, x2, x3), 1e+06, 10)
    }
    for (i in seq_len(N)){
      mat <- Xm[-i, ]
      det.inf_2 <- qr(t(mat)%*%mat)$rank
      if (det.inf_2 < ncol(mat)) {
        sc[i] <- 0
      } else {
        XpX2 <- t(mat) %*% mat
        C <- chol(XpX2, pivot = TRUE)
        C[lower.tri(C, diag = FALSE)] <- 0
        Z <- solve(t(C), t(Xpred))
        Tv <- Z^2
        D <- colSums(Tv)
        Mx <- max(D)
        sc[i] <- (p) / Mx
      }
    }
    G_score <- round(min(sc), 4)
  }
  return(G_score)
}

# returns a list of drop-one scores
rb_G_l <- function(X, N, K, order){
  Xm <- makeModelMat(X, 2)
  p <- ncol(Xm)
  XpX <- t(Xm) %*% Xm
  det.inf <- qr(XpX)$rank
  if (det.inf < K) {
    G_score <- - Inf
  } else {
    sc <- vector(mode = 'numeric', length = N)
    if (K == 1) {
      Xpred <- matrix(c(1, 1, 1, 1, 1,
                        -1, -0.5, 0.0, 0.5, 1.0,
                        1, 0.25, 0.0, 0.25, 1.0), 5, 3)
    } else if (K == 2){
      part_1 <- matrix(c(rep(1, 25),
                         rep(c(-1, -.5, 0, .5, 1), 5),
                         c(rep(-1, 5), rep(-.5, 5), rep(0, 5), rep(.5, 5), rep(1, 5))),25 ,3)
      part_2 <- matrix(c(part_1[,2] * part_1[,3], part_1[,2]^2, part_1[,3]^2), 25, 3)
      Xpred <- matrix(c(part_1, part_2), 25, 6)
    } else if (K == 3) {
      part_1 <- rep(1, 1e+06)
      x1 <- c(seq(from = -1, to = 1, length.out = 100))
      x2 <- expand.grid(x1, x1, x1)
      x3 <- cbind(x2[,1] * x2[,2], x2[,1] * x2[,3], x2[,2] * x2[,3], x2[,1]^2, x2[,2]^2, x2[,3]^2)
      Xpred <- as.matrix(cbind(part_1, x2, x3), 1e+06, 10)
    }
    for (i in seq_len(N)){
      mat <- Xm[-i, ]
      det.inf_2 <- qr(t(mat)%*%mat)$rank
      if (det.inf_2 < ncol(mat)) {
        sc[i] <- 0
      } else {
        XpX2 <- t(mat) %*% mat
        C <- chol(XpX2, pivot = TRUE)
        C[lower.tri(C, diag = FALSE)] <- 0
        Z <- solve(t(C), t(Xpred))
        Tv <- Z^2
        D <- colSums(Tv)
        Mx <- max(D)
        sc[i] <- (p) / Mx
      }
    }
  }
  return(sc)
}

#############
# RPG PLOTS #
#############

# The order the designs go in: I-Opt, RI-Opt, G-Opt, RG-Opt, CCD
# For the functions, func_1 is the robust function that returns
# a list while func_2 is the design I score without dropping values.

I_RPG <- function(..., func_1, func_2, order, k_v) {
  # Capture all design matrices as a list
  designs <- list(...)
  num_designs <- length(designs) # Number of designs provided
  
  N <- nrow(designs[[1]]) # Number of rows in each design matrix
  K <- ncol(designs[[1]]) # Number of columns in each design matrix
  
  # Compute lists for each design
  lists <- lapply(designs, function(X) {
    func1_result <- func_1(X, order, k_v)
    func2_result <- sort(func_2(X, order, k_v), decreasing = TRUE) # Sorted `func_2`
    c(func1_result, func2_result) # Combine `func_1` with sorted `func_2`
  })
  
  sorted_list <- lapply(lists, function(vec) {
    sort(vec, decreasing = TRUE)
  })
  
  # Find global min and max for normalization
  all_values <- unlist(sorted_list)
  v_min <- min(all_values)
  v_max <- max(all_values)
  
  # Create a sequence for plotting
  x_seq <- seq(0, N, 1)
  
  # Plot first design
  par(new = FALSE)
  plot(x_seq, (sorted_list[[1]] / v_max) * 100,
       main = sprintf("Robustness Profile Graph: I Design, K = %s, N = %s", K, N),
       col = "red", lwd = 3, ylab = "Efficiency",
       xlab = "Point Index", ylim = c(0, 100), type = "l")
  
  # Define colors and line types for designs
  colors <- c("red", "blue", "orange", "forestgreen", "black")
  lty <- c("solid", "dashed", "solid", "dashed", "twodash")
  
  # Add remaining designs to the plot
  for (i in 2:num_designs) {
    lines(x_seq, (sorted_list[[i]] / v_max) * 100, col = colors[i], lwd = 3, lty = lty[i])
  }
  
  # Add vertical line and legend
  abline(v = 1, col = "gray", lty = 2)
  legend_labels <- c("I-Opt", "RI-Opt", "G-Opt", "RG-Opt")
  if (num_designs == 5) legend_labels <- c(legend_labels, "CCD")
  legend("bottomleft", legend = legend_labels,
         col = colors[1:num_designs], lty = lty[1:num_designs])
}

# The order the designs go in: G-Opt, RG-Opt, I-Opt, RI-Opt, CCD
# For the functions, func_1 is the robust function that returns

G_RPG <- function(..., func_1, func_2, order) {
  # Capture all design matrices as a list
  designs <- list(...)
  num_designs <- length(designs) # Number of designs provided
  
  N <- nrow(designs[[1]]) # Number of rows in each design matrix
  K <- ncol(designs[[1]]) # Number of columns in each design matrix
  
  # Compute lists for each design
  lists <- lapply(designs, function(X) {
    func1_result <- func_1(X, N, K, order)
    func2_result <- sort(func_2(X, N, K, order), decreasing = TRUE) # Sorted `func_2`
    c(func1_result, func2_result) # Combine `func_1` with sorted `func_2`
  })
  
  sorted_list <- lapply(lists, function(vec) {
    sort(vec, decreasing = TRUE)
  })
  
  # Find global min and max for normalization
  all_values <- unlist(sorted_list)
  v_min <- min(all_values)
  v_max <- max(all_values)
  
  # Create a sequence for plotting
  x_seq <- seq(0, N, 1)
  
  # Plot first design
  par(new = FALSE)
  plot(x_seq, (sorted_list[[1]] / v_max) * 100,
       main = sprintf("Robustness Profile Graph: G Design, K = %s, N = %s", K, N),
       col = "red", lwd = 3, ylab = "Efficiency",
       xlab = "Point Index", ylim = c(0, 100), type = "l")
  
  # Define colors and line types for designs
  colors <- c("red", "blue", "orange", "forestgreen", "black")
  lty <- c("solid", "dashed", "solid", "dashed", "twodash")
  
  # Add remaining designs to the plot
  for (i in 2:num_designs) {
    lines(x_seq, (sorted_list[[i]] / v_max) * 100, col = colors[i], lwd = 3, lty = lty[i])
  }
  
  # Add vertical line and legend
  abline(v = 1, col = "gray", lty = 2)
  legend_labels <- c("G-Opt", "RG-Opt", "I-Opt", "RI-Opt")
  if (num_designs == 5) legend_labels <- c(legend_labels, "CCD")
  if (K > 1) {
  legend("topright", legend = legend_labels,
         col = colors[1:num_designs], lty = lty[1:num_designs])
  }else{
    legend("bottomleft", legend = legend_labels,
           col = colors[1:num_designs], lty = lty[1:num_designs])
    }
}

############
# 2D PLOTS #
############

I_2D <- function(X, l_crit, n_crit, interaction, k_type, title) {
  l_score <- round(l_crit(X, interaction = interaction, k = k_type), 4)
  n_score <- seq_len(nrow(X))
  min_v <- if (any(l_score > 0)) min(l_score[l_score > 0]) else NA
  max_v <- max(l_score, na.rm = TRUE)
  
  # Assign effect categories to each point
  effect_sc <- ifelse(
    l_score == 0, "Singular",
    ifelse(l_score == min_v, "Worst Case",
           ifelse(l_score == max_v, "Best Case", " "))
  )
  
  # Ensure all levels are included in the factor
  effect_sc <- factor(effect_sc, levels = c("Best Case", " ", "Worst Case", "Singular"))
  
  # Sort lists and rows of X by descending l_score
  order_indices <- order(l_score, decreasing = TRUE)
  sorted_l_score <- l_score[order_indices]
  sorted_effect_sc <- effect_sc[order_indices]
  new_X <- X[order_indices, , drop = FALSE]
  
  # Create data frame for ggplot
  plot_data <- data.frame(
    new_X,
    sorted_effect_sc,
    n_score,
    effect_sc
  )
  
  # Add dummy rows to ensure all legend levels appear
  dummy_rows <- data.frame(
    X1 = NA, X2 = NA, sorted_effect_sc = factor(c("Best Case", " ", "Worst Case", "Singular"), 
                                                levels = c("Best Case", " ", "Worst Case", "Singular")),
    n_score = NA, effect_sc = factor(c("Best Case", " ", "Worst Case", "Singular"), 
                                     levels = c("Best Case", " ", "Worst Case", "Singular"))
  )
  
  plot_data <- rbind(plot_data, dummy_rows)
  
  # Generate ggplot
  g <- ggplot(plot_data, aes(
    x = as.numeric(X1), 
    y = as.numeric(X2), 
    color = sorted_effect_sc
  )) +
    geom_point(size = 5, na.rm = TRUE) +
    geom_text_repel(aes(label = n_score), size = 6, color = "black", na.rm = TRUE) +
    labs(
      title = sprintf("%s: K = %s, N = %s", title, ncol(X), nrow(X)),
      subtitle = sprintf(
        "I-Efficiency: %s \n Robust Efficiency: %s",
        round(n_crit(X, interaction = interaction, k = k_type), 4),
        if (any(l_score == 0)) 0 else min_v
      ),
      x = "X1",
      y = "X2"
    ) +
    scale_color_manual(
      name = "Point Effect",
      values = c(
        "Best Case" = "green3", 
        " " = "black", 
        "Worst Case" = "orange", 
        "Singular" = "red"
      ),
      drop = FALSE # Ensure unused levels are displayed
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      plot.subtitle = element_text(hjust = 0.5),
      text = element_text(size = 12),
      panel.border = element_rect(color = "black", fill = NA, size = 2)
    )
  
  return(g)
}

G_2D <- function(X, l_crit, n_crit, title) {
  l_score <- round(l_crit(X, nrow(X), ncol(X), 2), 4)
  n_score <- seq_len(nrow(X))
  min_v <- if (any(l_score > 0)) min(l_score[l_score > 0]) else NA
  max_v <- max(l_score, na.rm = TRUE)
  
  # Assign effect categories to each point
  effect_sc <- ifelse(
    l_score == 0, "Singular",
    ifelse(l_score == min_v, "Worst Case",
           ifelse(l_score == max_v, "Best Case", " "))
  )
  
  # Ensure all levels are included in the factor
  effect_sc <- factor(effect_sc, levels = c("Best Case", " ", "Worst Case", "Singular"))
  
  # Sort lists and rows of X by descending l_score
  order_indices <- order(l_score, decreasing = TRUE)
  sorted_l_score <- l_score[order_indices]
  sorted_effect_sc <- effect_sc[order_indices]
  new_X <- X[order_indices, , drop = FALSE]
  
  # Create data frame for ggplot
  plot_data <- data.frame(
    new_X,
    sorted_effect_sc,
    n_score,
    effect_sc
  )
  
  first_row <- plot_data[1, , drop = FALSE]
  
  # Step 2: Remove the first row
  plot_data <- plot_data[-1, , drop = FALSE]
  
  # Step 3: Append the first row to the end
  plot_data <- rbind(plot_data, first_row)
  
  # Add dummy rows to ensure all legend levels appear
  dummy_rows <- data.frame(
    X1 = NA, X2 = NA, sorted_effect_sc = factor(c("Best Case", " ", "Worst Case", "Singular"), 
                                                levels = c("Best Case", " ", "Worst Case", "Singular")),
    n_score = NA, effect_sc = factor(c("Best Case", " ", "Worst Case", "Singular"), 
                                     levels = c("Best Case", " ", "Worst Case", "Singular"))
  )
  
  plot_data <- rbind(plot_data, dummy_rows)
  plot_data$effect_sc <- factor(plot_data$effect_sc, levels = c("Best Case", "Worst Case", "Singular", " "))
  
  g <- ggplot(plot_data, aes(
    x = as.numeric(X1), 
    y = as.numeric(X2), 
    color = sorted_effect_sc
  )) +
    geom_point(size = 5, na.rm = TRUE) +
    labs(
      title = sprintf("%s: K = %s, N = %s", title, ncol(X), nrow(X)),
      subtitle = sprintf(
        "G-Efficiency: %s \n Robust Efficiency: %s",
        round(n_crit(X, nrow(X), ncol(X), 2), 4),
        if (any(l_score == 0)) 0 else min_v
      ),
      x = "X1",
      y = "X2"
    ) +
    geom_text_repel(aes(label = n_score), size = 6, color = "black", na.rm = TRUE) +
    scale_color_manual(
      name = "Point Effect",
      values = c(
        "Best Case" = "green3", 
        " " = "black", 
        "Worst Case" = "orange", 
        "Singular" = "red"
      ),
      drop = FALSE # Ensure unused levels are displayed
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      plot.subtitle = element_text(hjust = 0.5),
      text = element_text(size = 12),
      panel.border = element_rect(color = "black", fill = NA, size = 2)
    )
  
  return(g)
}

############
# 3D PLOTS #
############

addLines <- function(){
  
  x0tt = c(0, -1, -1, 1, -1, -1)
  x1tt = c(0, 1, -1, 1, 1, 1)
  y0tt = c(-1, 0, -1, -1, -1, 1)
  y1tt = c(1, 0, 1, 1, -1, 1)
  z0tt = rep(0, 6)
  z1tt = rep(0, 6)
  
  x0tt = c(0, -1)
  x1tt = c(0, 1)
  y0tt = c(-1, 0)
  y1tt = c(1, 0)
  z0tt = rep(0, 2)
  z1tt = rep(0, 2)
  
  segments3D(x0 = x0tt,
             x1 = x1tt,
             y0 = y0tt,
             y1 = y1tt,
             z0 = z0tt,
             z1 = z1tt,
             add = T, lty = 2)
  
  x0tt = c(0, -1, -1, 1, -1, -1)
  x1tt = c(0, 1, -1, 1, 1, 1)
  y0tt = c(-1, 0, -1, -1, -1, 1)
  y1tt = c(1, 0, 1, 1, -1, 1)
  z0tt = rep(1, 6)
  z1tt = rep(1, 6)
  
  
  segments3D(x0 = x0tt,
             x1 = x1tt,
             y0 = y0tt,
             y1 = y1tt,
             z0 = z0tt,
             z1 = z1tt,
             add = T, lty = 1)
  
  x0tt = c(0, -1, -1, 1, -1, -1)
  x1tt = c(0, 1, -1, 1, 1, 1)
  y0tt = c(-1, 0, -1, -1, -1, 1)
  y1tt = c(1, 0, 1, 1, -1, 1)
  z0tt = rep(-1, 6)
  z1tt = rep(-1, 6)
  
  
  segments3D(x0 = x0tt,
             x1 = x1tt,
             y0 = y0tt,
             y1 = y1tt,
             z0 = z0tt,
             z1 = z1tt,
             add = T, lty = 1)
  
  x0tt = c(1, 0, 1, -1, -1)
  x1tt = c(1, 0,1, -1, -1)
  y0tt = c(-1, 0, 1, -1, 1)
  y1tt = c(-1, 0, 1, -1, 1)
  z0tt = c(-1, -1, -1, -1, -1)
  z1tt = c(1,1, 1, 1, 1)
  
  
  segments3D(x0 = x0tt,
             x1 = x1tt,
             y0 = y0tt,
             y1 = y1tt,
             z0 = z0tt,
             z1 = z1tt,
             add = T, lty = 2)
  
  
}

I_3D <- function(design, func_1, func_2, order, k_v, type, 
                                          xlim = c(-1, 1), ylim = c(-1, 1), zlim = c(-1, 1)) {
  # Calculate metrics
  func1_value <- func_1(design, order, k_v)
  func2_values <- round(func_2(design, order, k_v), 4) # Compute func_2 for each design point

  # Assign colors based on func_2 values
  colors <- rep("black", nrow(design)) # Default color
  # Identify all max and min values
  max_indices <- which(func2_values == max(func2_values))  # Indices of all max values
  min_indices <- which(func2_values == min(func2_values))  # Indices of all min values
  colors[max_indices] <- "green3"
  colors[min_indices] <- "orange"
  
  # Generate the scatter plot
  scatter3D(
    x      = design[, 1],
    y      = design[, 2],
    z      = design[, 3],
    xlim   = xlim,
    ylim   = ylim,
    zlim   = zlim,
    colvar = NULL,
    pch    = 16,
    xlab   = expression(x[1]),
    ylab   = expression(x[2]),
    zlab   = expression(x[3]),
    col    = colors,
    cex    = 2,
    ticktype = "detailed",
    bty      = "g",
    phi      = 22,
    theta    = 30,
    main     = sprintf("%s-Criterion: K = %s, N = %s", type, nrow(design), ncol(design)),
    cex.axis = 0.6,
    cex.lab  = 0.7
  )
  mtext(sprintf("I-Efficiency: %s", I_crit(design, 2, K_3)), side = 3, line = .6, cex = 0.8)
  mtext(sprintf("Robust Efficiency: %s", rb_I(design, 2, K_3)), side = 3, line = -0.15, cex = 0.8)
  addLines()
}

G_3D <- function(design, func_1, func_2, type, 
                 xlim = c(-1, 1), ylim = c(-1, 1), zlim = c(-1, 1)) {
  # Calculate metrics
  N <- nrow(design)
  K <- ncol(design)
  func1_value <- round(func_1(design, N, K, 2), 4)
  func2_values <- round(func_2(round(design, 3), N, K, 2), 3) # Compute func_2 for each design point
  
  # Assign colors based on func_2 values
  colors <- rep("black", nrow(design)) # Default color
  # Identify all max and min values
  max_indices <- which(func2_values == max(func2_values))  # Indices of all max values
  min_indices <- which(func2_values == min(func2_values))  # Indices of all min values
  colors[max_indices] <- "green3"
  colors[min_indices] <- "orange"
  
  # Generate the scatter plot
  scatter3D(
    x      = design[, 1],
    y      = design[, 2],
    z      = design[, 3],
    xlim   = xlim,
    ylim   = ylim,
    zlim   = zlim,
    colvar = NULL,
    pch    = 16,
    xlab   = expression(x[1]),
    ylab   = expression(x[2]),
    zlab   = expression(x[3]),
    col    = colors,
    cex    = 2,
    ticktype = "detailed",
    bty      = "g",
    phi      = 22,
    theta    = 30,
    main     = sprintf("%s-Criterion: K = %s, N = %s", type, K, N),
    cex.axis = 0.6,
    cex.lab  = 0.7
  )
  mtext(sprintf("G-Efficiency: %s", func1_value), side = 3, line = .6, cex = 0.8)
  mtext(sprintf("Robust Efficiency: %s", min(func2_values)), side = 3, line = -0.15, cex = 0.8)
  addLines()
}

################
# EXAMPLE WORK #
################

string_to_matrix <- function(input_string, N, K) {
  input_vector <- as.numeric(unlist(strsplit(input_string, " ")))
  matrix_result <- matrix(input_vector, nrow = N, ncol = K)
  return(matrix_result)
}

exp1 <- read.csv("Study_Collection/K=2/N=10/Optimal-N_10_K_2.csv")
E_I_2_10 <- string_to_matrix(exp1[1,3], 10, 2)
PSO_I_2_10 <- string_to_matrix(exp1[2,3], 10, 2)
E_G_2_10 <- string_to_matrix(exp1[3,3], 10, 2)
PSO_G_2_10 <- string_to_matrix(exp1[4,3], 10, 2)
ccf2_10 <- string_to_matrix(exp1[1,3], 10, 2)

exp2 <- read.csv("Study_Collection/K=2/N=9/Optimal-N_9_K_2.csv")
E_I_2_9 <- string_to_matrix(exp2[1,3], 9, 2)
PSO_I_2_9 <- string_to_matrix(exp2[2,3], 9, 2)
E_G_2_9 <- string_to_matrix(exp2[3,3], 9, 2)
PSO_G_2_9 <- string_to_matrix(exp2[4,3], 9, 2)
ccf2_9 <- string_to_matrix(exp2[5,3], 9, 2)

exp3 <- read.csv("Study_Collection/K=3/N=11/Optimal-N_11_K_3.csv")
E_I_3_11 <- string_to_matrix(exp3[1,3], 11, 3)

exp4 <- read.csv("Study_Collection/K=3/N=14/Optimal-N_14_K_3.csv")
PSO_G_3_14 <- string_to_matrix(exp4[4,3], 14, 3)

# example: RPG on the I-Criterion
I_RPG(E_I_2_10, PSO_I_2_10, E_G_2_10, PSO_G_2_10, ccf2_10,
      func_1 = rb_I_l, func_2 = I_crit, order = 2, k_v = K_2)

# example: RPG on the G-Criterion
G_RPG(E_G_2_9, PSO_G_2_9, E_I_2_9, PSO_I_2_9, ccf2_9, 
      func_1 = G_crit, func_2 = rb_G_l, order = 2)

# example: 2D plot of the K = 2 N = 10 Exact I-Optimal Design 
I_2D(E_I_2_10, rb_I_l, I_crit, 2, K_2, "I-Criterion")

# example: 2D plot of the K = 2 N = 9 Robust G-Optimal Design 
G_2D(PSO_G_2_9, rb_G_l, G_crit, "RG-Criterion")

# example: 3D plot of the K = 3 N = 11 Exact I-Optimal Design 
I_3D(E_I_3_11, I_crit, rb_I_l, order = 2, k_v = K_3, 
                              type = "I")
# example: 3D plot of the K = 3 N = 14 Robust G-Optimal Design 
G_3D(PSO_G_3_14, G_crit, rb_G_l, "RG")
