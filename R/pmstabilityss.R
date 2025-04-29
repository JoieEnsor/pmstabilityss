  #' Sample Size for Prediction Model Development Based on Prediction Stability
  #'
  #' @param data A data frame containing the variables for analysis
  #' @param varlist Character vector of variable names in the linear predictor
  #' @param prevalence Numeric value of prevalence for the outcome
  #' @param lp Character name of linear predictor variable
  #' @param subgroup Character name of subgroup variable (optional)
  #' @param nodraw Logical, if TRUE, suppresses drawing plots
  #' @param nocompare Logical, if TRUE, suppresses comparison plots
  #' @param pmss Logical, if TRUE, uses pmsampsize package to calculate sample size
  #' @param n Numeric vector of sample sizes to evaluate
  #' @param cstatistic Numeric C-statistic value (used for pmsampsize)
  #' @param simobs Number of observations for simulation (default 500000)
  #' @param seed Random seed value (default 1234)
  #' @param pciwidth Numeric vector of target prediction interval widths
  #' @param pcutpoints Numeric vector of probability cutpoints (must match pciwidth)
  #' @param threshold Numeric classification threshold (default 0)
  #' @param tolerance Numeric convergence tolerance (default 0.005)
  #' @param logitpvarincrement Numeric increment for logit p variance (default 0.001)
  #' @param color Plot color (default "darkgrey")
  #'
  #' @return A list containing tables and plot objects
  #' @export
  #'
  #' @importFrom stats logit invlogit lowess var median quantile sd
  #' @importFrom graphics plot lines points text
  #' @importFrom grDevices rgb
  #' @importFrom Matrix Matrix
  pmstabilityss <- function(data,
                          varlist,
                          prevalence,
                          lp,
                          subgroup = NULL,
                          nodraw = FALSE,
                          nocompare = FALSE,
                          color = "darkgrey",
                          pmss = FALSE,
                          n = NULL,
                          cstatistic = 0,
                          simobs = 500000,
                          seed = 1234,
                          pciwidth = NULL,
                          pcutpoints = NULL,
                          threshold = 0,
                          tolerance = 0.005,
                          logitpvarincrement = 0.001) {

    # Check for required packages
    required_packages <- c("dplyr", "ggplot2", "patchwork", "Matrix")

    # Install missing packages
    for (pkg in required_packages) {
      if (!requireNamespace(pkg, quietly = TRUE)) {
        install.packages(pkg)
      }
    }

    library(dplyr)
    library(ggplot2)
    library(gridExtra)
    library(Matrix)

    # Define functions needed
    logit <- function(p) log(p/(1-p))
    invlogit <- function(x) exp(x)/(1+exp(x))

    # Validate inputs
    if (is.null(lp) && is.null(betas)) {
      stop("Either lp or betas must be provided")
    }

    if (!is.null(pcutpoints) && !is.null(pciwidth)) {
      if (length(pcutpoints) != length(pciwidth)) {
        stop("Number of cutpoints must match number of target interval widths")
      }
      if (pcutpoints[length(pcutpoints)] != 1) {
        stop("pcutpoints numlist should end with 1")
      }
    }

    # Create variables for calculations
    data$pmsss_intercept <- 1  # Proxy variable for intercept

    # Create new lp_var
    data$lp_var <- data[[lp]]

    # Calculate Information matrix elements
    # Number of variables including intercept
    elements <- length(varlist) + 1

    # Store variable names for information matrix
    var_names <- c("pmsss_intercept", varlist)

    # Initialize information matrix
    I_matrix <- matrix(0, nrow = elements, ncol = elements)

    # Compute off-diagonals of information matrix
    for (col in 1:elements) {
      for (row in col:elements) {
        # Element-wise multiplication and taking mean
        var_product <- data[[var_names[row]]] * data[[var_names[col]]] *
          (exp(data$lp_var)/((1+exp(data$lp_var))^2))
        I_matrix[row, col] <- mean(var_product)
        # Fill in symmetric part
        if (row != col) {
          I_matrix[col, row] <- I_matrix[row, col]
        }
      }
    }

    # Inverse of information matrix
    invI <- solve(I_matrix)

    # Calculate individual unit variance
    unit_v <- numeric(nrow(data))
    for (i in 1:nrow(data)) {
      x_i <- c(data$pmsss_intercept[i], sapply(varlist, function(v) data[[v]][i]))
      unit_v[i] <- t(x_i) %*% invI %*% x_i
    }

    data$unit_v <- unit_v

    # Generate predicted probabilities
    data$p_true <- invlogit(data$lp_var)

    # Current dataset size
    cur_data_n <- nrow(data)
    cat("Fixed SS of input dataset =", cur_data_n, "\n\n")

    # Parse user cutpoints for p and desired widths
    if (!is.null(pcutpoints) && !is.null(pciwidth)) {
      widths <- data.frame(
        categories = pcutpoints,
        width = pciwidth,
        cat_code = 1:length(pcutpoints),
        avg_p = NA
      )

      # Assign categories and widths
      data$width <- NA
      data$cat_code <- NA

      for (i in 1:nrow(widths)) {
        idx <- which(is.na(data$width) & data$p_true <= widths$categories[i])
        data$width[idx] <- widths$width[i]
        data$cat_code[idx] <- widths$cat_code[i]
      }

      # Calculate average p by category
      for (i in 1:nrow(widths)) {
        widths$avg_p[i] <- median(data$p_true[data$cat_code == i], na.rm = TRUE)
      }

      # Matrix of values for p groupings
      p_values <- seq(0.01, 1.00, by = 0.01)
      pvars <- data.frame(
        pcode = 1:length(p_values),
        pcat = p_values,
        width = NA,
        logit_p_var = NA
      )

      pvars$pcat[length(pvars$pcat)] <- 0.99

      # Assign target widths
      widths_desc <- widths[order(widths$categories, decreasing = TRUE), ]

      for (j in 1:nrow(pvars)) {
        for (k in 1:nrow(widths_desc)) {
          if (round(pvars$pcat[j], 3) <= round(widths_desc$categories[k], 3)) {
            pvars$width[j] <- widths_desc$width[k]
          }
        }

        # Identify appropriate target var(logit(p)) given pciwidth
        width_p <- 0
        var_logit_p_loc <- 0


        while (width_p < pvars$width[j]) {
          var_logit_p_loc <- var_logit_p_loc + logitpvarincrement
          ub_p <- invlogit(logit(pvars$pcat[j]) + (1.96 * sqrt(var_logit_p_loc)))
          lb_p <- invlogit(logit(pvars$pcat[j]) - (1.96 * sqrt(var_logit_p_loc)))
          width_p <- ub_p - lb_p
        }

        pvars$logit_p_var[j] <- var_logit_p_loc
      }

      # Round p values for merging
      data$p_round <- round(data$p_true, 2)
      data$p_round[data$p_round == 0] <- 0.01
      data$p_round[data$p_round == 1] <- 0.99

      # Merge back to get var_logit_p
      data$var_logit_p <- NA
      for (i in 1:nrow(data)) {
        idx <- which(abs(pvars$pcat - data$p_round[i]) < 0.001)
        if (length(idx) > 0) {
          data$var_logit_p[i] <- pvars$logit_p_var[idx[1]]
        }
      }

      # SS required to meet desired precision
      data$ss_target_var <- (1/data$var_logit_p) * data$unit_v
      target_var_min_N <- ceiling(max(data$ss_target_var, na.rm = TRUE))

      cat("Minimum SS required to meet target UI widths =", target_var_min_N, "\n\n")
    }

    # Calculate minimum SS from pmsampsize if requested
    if (pmss) {
      # Check if pmsampsize is available
      if (!requireNamespace("pmsampsize", quietly = TRUE)) {
        install.packages("pmsampsize")
      }

      library(pmsampsize)

      # Use pmsampsize
      pmss_result <- pmsampsize::pmsampsize(
        type = "b",
        cstatistic = cstatistic,
        parameters = length(varlist),
        prevalence = prevalence
      )

      pmss_min_N <- pmss_result$sample_size
      cat("Minimum SS required by pmsampsize =", pmss_min_N, "\n\n")
    }

    # Build up the list of N's to investigate
    ss_tests <- cur_data_n

    if (!is.null(pcutpoints) && !is.null(pciwidth)) {
      ss_tests <- c(ss_tests, target_var_min_N)
    }

    if (pmss) {
      ss_tests <- c(ss_tests, pmss_min_N)
    }

    if (!is.null(n)) {
      ss_tests <- c(ss_tests, n)
    }

    # Remove duplicates and sort
    ss_tests <- sort(unique(ss_tests))


    # Create data frames for results
    overall_stats <- data.frame(
      N = numeric(),
      Mean = numeric(),
      Min = numeric(),
      Median = numeric(),
      Max = numeric()
    )

    if (threshold != 0) {
      threshold_stats <- data.frame(
        N = numeric(),
        Threshold = numeric(),
        Mean = numeric(),
        Min = numeric(),
        Median = numeric(),
        Max = numeric()
      )
    }

    if (!is.null(pcutpoints) && !is.null(pciwidth)) {
      cuts_stats <- data.frame(
        N = numeric(),
        P_category = numeric(),
        Target_width = numeric(),
        Mean = numeric(),
        Min = numeric(),
        Median = numeric(),
        Max = numeric(),
        Prop_target_width_met = numeric()
      )
    }

    # Lists for plots
    instability_plots <- list()
    classification_plots <- list()

    # Calculate statistics for each N
    for (num in ss_tests) {
      # Calculate variance for individual at sample size N
      data[[paste0("var_ind_", num)]] <- (1/num) * data$unit_v

      # Calculate confidence intervals on logit scale
      data[[paste0("lower_logitp_", num)]] <- data$lp_var - (1.96 * sqrt(data[[paste0("var_ind_", num)]]))
      data[[paste0("upper_logitp_", num)]] <- data$lp_var + (1.96 * sqrt(data[[paste0("var_ind_", num)]]))

      # Convert to probability scale
      data[[paste0("lower_p_", num)]] <- invlogit(data[[paste0("lower_logitp_", num)]])
      data[[paste0("upper_p_", num)]] <- invlogit(data[[paste0("upper_logitp_", num)]])

      # Calculate width of interval
      data[[paste0("width_", num)]] <- data[[paste0("upper_p_", num)]] - data[[paste0("lower_p_", num)]]

      # Summarize widths overall
      width_summary <- data[[paste0("width_", num)]]
      overall_stats <- rbind(overall_stats, data.frame(
        N = num,
        Mean = mean(width_summary, na.rm = TRUE),
        Min = min(width_summary, na.rm = TRUE),
        Median = median(width_summary, na.rm = TRUE),
        Max = max(width_summary, na.rm = TRUE)
      ))

      # Create instability plot
      plot_data <- data.frame(
        p_true = data$p_true,
        lower = data[[paste0("lower_p_", num)]],
        upper = data[[paste0("upper_p_", num)]],
        width = data[[paste0("width_", num)]]
      )

      insta_p <- suppressMessages(ggplot(plot_data, aes(x = p_true, y = p_true)) +
        geom_segment(aes(x = p_true, xend = p_true, y = lower, yend = upper),
                     color = color, alpha = 0.7) +
        geom_abline(slope = 1, intercept = 0, linetype = "solid") +
        geom_smooth(formula = y ~ x, aes(y = lower), method = "loess", se = FALSE,
                    linetype = "dashed", color = "black", span = 0.2) +
        geom_smooth(formula = y ~ x, aes(y = upper), method = "loess", se = FALSE,
                    linetype = "dashed", color = "black", span = 0.2) +
        annotate("text", x = 0.1, y = 0.9, label = paste("N =", num), hjust = 0) +
        labs(x = "True risk", y = "95% Uncertainty interval\nfor true risk") +
        theme_minimal() +
        theme(aspect.ratio = 1) +
        scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
        scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)))

      instability_plots[[as.character(num)]] <- insta_p

      # Handle cutpoints if specified
      if (!is.null(pcutpoints) && !is.null(pciwidth)) {
        data[[paste0("width_good_", num)]] <- ifelse(
          data[[paste0("width_", num)]] > data$width, 0, 1
        )

        # Calculate proportion meeting target width
        met_target_width <- mean(data[[paste0("width_good_", num)]], na.rm = TRUE)

        for (c in 1:nrow(widths)) {
          width_subset <- data[[paste0("width_", num)]][data$cat_code == c]

          cuts_stats <- rbind(cuts_stats, data.frame(
            N = num,
            P_category = widths$categories[c],
            Target_width = widths$width[c],
            Mean = mean(width_subset, na.rm = TRUE),
            Min = min(width_subset, na.rm = TRUE),
            Median = median(width_subset, na.rm = TRUE),
            Max = max(width_subset, na.rm = TRUE),
            Prop_target_width_met = met_target_width
          ))
        }
      }

      # Classification instability if threshold specified
      if (threshold != 0) {
        data[[paste0("z_", num)]] <- (logit(threshold) - data$lp_var) /
          sqrt(data[[paste0("var_ind_", num)]])
        data[[paste0("prob_above_", num)]] <- 1 - pnorm(data[[paste0("z_", num)]])
        data[[paste0("prob_below_", num)]] <- pnorm(data[[paste0("z_", num)]])

        data[[paste0("prob_different_", num)]] <- NA
        data[[paste0("prob_different_", num)]][data$p_true < threshold] <-
          data[[paste0("prob_above_", num)]][data$p_true < threshold]
        data[[paste0("prob_different_", num)]][data$p_true >= threshold] <-
          data[[paste0("prob_below_", num)]][data$p_true >= threshold]

        # Summarize misclassification
        prob_diff_summary <- data[[paste0("prob_different_", num)]]
        threshold_stats <- rbind(threshold_stats, data.frame(
          N = num,
          Threshold = threshold,
          Mean = mean(prob_diff_summary, na.rm = TRUE),
          Min = min(prob_diff_summary, na.rm = TRUE),
          Median = median(prob_diff_summary, na.rm = TRUE),
          Max = max(prob_diff_summary, na.rm = TRUE)
        ))

        # Create classification instability plot
        class_data <- data.frame(
          p_true = data$p_true,
          prob_different = data[[paste0("prob_different_", num)]]
        )

        class_p <- suppressMessages(ggplot(class_data, aes(x = p_true, y = prob_different)) +
          geom_point(color = color, alpha = 0.7, size = 1) +
          annotate("text", x = max(class_data$p_true, na.rm = TRUE),
                   y = max(class_data$prob_different, na.rm = TRUE),
                   label = paste("N =", num), hjust = 1) +
          labs(x = "True risk",
               y = "Probability of misclassification") +
          theme_minimal() +
          scale_x_continuous(breaks = seq(0, 1, 0.2)) +
          scale_y_continuous(breaks = seq(0, 1, 0.2)))

        classification_plots[[as.character(num)]] <- class_p
      }

      # Subgroup analysis if specified
      if (!is.null(subgroup)) {
        # To be implemented
      }
    }

    # Display tables

    
    # Generate combined plots if requested
    plot_output <- list()
    
    if (!nodraw && !nocompare) {
      if (length(instability_plots) > 0) {
        # Use patchwork to combine plots
        if (!requireNamespace("patchwork", quietly = TRUE)) {
          install.packages("patchwork")
        }
        library(patchwork)
        
        # Start with the first plot
        combined_instability <- instability_plots[[1]]
        
        # Add the rest of the plots
        if (length(instability_plots) > 1) {
          for (i in 2:length(instability_plots)) {
            combined_instability <- combined_instability + instability_plots[[i]]
          }
        }
        
        # Arrange in a grid
        combined_instability <- combined_instability + 
          plot_layout(ncol = min(3, length(instability_plots)))
        
        plot_output$instability_plots <- instability_plots
        plot_output$combined_instability <- combined_instability
        
        if (!is.null(dev.list())) {
          print(combined_instability)
        }
      }
      
      if (threshold != 0 && length(classification_plots) > 0) {
        # Start with the first plot
        combined_classification <- classification_plots[[1]]
        
        # Add the rest of the plots
        if (length(classification_plots) > 1) {
          for (i in 2:length(classification_plots)) {
            combined_classification <- combined_classification + classification_plots[[i]]
          }
        }
        
        # Arrange in a grid
        combined_classification <- combined_classification + 
          plot_layout(ncol = min(3, length(classification_plots)))
        
        plot_output$classification_plots <- classification_plots
        plot_output$combined_classification <- combined_classification
        
        if (!is.null(dev.list())) {
          print(combined_classification)
        }
      }
    }
    

    # Prepare return object
    result <- list(
      overall_stats = overall_stats,
      #unit_variance = data$unit_v,
      input_data = data,
      plots = plot_output
    )

    if (threshold != 0) {
      result$threshold_stats <- threshold_stats
    }

    if (!is.null(pcutpoints) && !is.null(pciwidth)) {
      result$cuts_stats <- cuts_stats
    }

    class(result) <- "pmstabilityss"
    return(result)
  }

#' Print method for pmstabilityss objects
#'
#' @param x A pmstabilityss object
#' @param ... Additional arguments
#'
#' @export
print.pmstabilityss <- function(x, ...) {
  cat("pmstabilityss object\n\n")
  cat("Overall summary UI widths:\n")
  print(x$overall_stats, row.names = FALSE)

  if (!is.null(x$threshold_stats)) {
    cat("\nSummary of proportion of misclassified:\n")
    print(x$threshold_stats, row.names = FALSE)
  }

  if (!is.null(x$cuts_stats)) {
    cat("\nSummary UI widths by probability categories:\n")
    print(x$cuts_stats, row.names = FALSE)
  }

  cat("\nUse $plots to access plot objects\n")
}

#' Plot method for pmstabilityss objects
#'
#' @param x A pmstabilityss object
#' @param type Type of plot to display: "instability" or "classification"
#' @param combined Whether to show combined plots (TRUE) or individual plots (FALSE)
#' @param ... Additional arguments
#'
#' @export

#' Plot method for pmstability objects
#'
#' @param x A pmstability object
#' @param type Type of plot to display: "instability" or "classification"
#' @param combined Whether to show combined plots (TRUE) or individual plots (FALSE)
#' @param ... Additional arguments
#'
#' @export
plot.pmstabilityss <- function(x, type = "instability", combined = TRUE, ...) {
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    install.packages("patchwork")
  }
  library(patchwork)
  
  if (type == "instability") {
    if (combined && !is.null(x$plots$combined_instability)) {
      return(x$plots$combined_instability)
    } else if (!is.null(x$plots$instability_plots) && combined) {
      # Use patchwork to combine plots
      plot_list <- x$plots$instability_plots
      # Start with the first plot
      if (length(plot_list) > 0) {
        combined_plot <- plot_list[[1]]
        # Add the rest of the plots
        if (length(plot_list) > 1) {
          for (i in 2:length(plot_list)) {
            combined_plot <- combined_plot + plot_list[[i]]
          }
        }
        # Arrange in a grid
        combined_plot <- combined_plot + 
          plot_layout(ncol = min(3, length(plot_list)))
        
        return(combined_plot)
      }
    } else if (!is.null(x$plots$instability_plots)) {
      return(x$plots$instability_plots)
    }
  } else if (type == "classification") {
    if (combined && !is.null(x$plots$combined_classification)) {
      return(x$plots$combined_classification)
    } else if (!is.null(x$plots$classification_plots) && combined) {
      # Use patchwork to combine plots
      plot_list <- x$plots$classification_plots
      # Start with the first plot
      if (length(plot_list) > 0) {
        combined_plot <- plot_list[[1]]
        # Add the rest of the plots
        if (length(plot_list) > 1) {
          for (i in 2:length(plot_list)) {
            combined_plot <- combined_plot + plot_list[[i]]
          }
        }
        # Arrange in a grid
        combined_plot <- combined_plot + 
          plot_layout(ncol = min(3, length(plot_list)))
        
        return(combined_plot)
      }
    } else if (!is.null(x$plots$classification_plots)) {
      return(x$plots$classification_plots)
    }
  }
  
  message("No plots of type '", type, "' available.")
  return(invisible(NULL))
}
                                               
