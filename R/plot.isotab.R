plot.isotab <- function(x, labels = FALSE, text.size = 15, title = NULL,
                        phi.min = "isotab", p.max = "isotab", ...) {

  if (!methods::is(x, 'isotab')) stop("x is not of class isotab")

  phi <- x$phi
  ft  <- x$fisher_p
  cluster_names <- colnames(x$phi)
  it_phi <- x$thresholds[1]

  # phi.min
  if(phi.min == "isotab") {
    phi.min <- x$thresholds["phi.min"]
  } else {
    if (phi.min < -1 | phi.min > 1) stop("phi.min must be between -1 and 1")
  }

  # p.max
  if(p.max == "isotab") {
    p.max <- x$thresholds["p.max"]
  } else {
    if (p.max < 0 | p.max > 1) stop("p.max must be between 0 and 1")
  }

  # Identify rows with max abs value greater than or equal to phi.min
  idx <- which(apply(phi, 1, max, na.rm = TRUE) >= abs(phi.min))
  nidx <- length(idx)
  if (nidx == 0) stop("phi.min not reached")

  # Assign values to clusters
  # Extract names
  nam <- names(idx)
  # Filter phi and p based on indices
  filtered_phi <- phi[idx, ]
  filtered_p <- ft[idx, ]

  if (nidx < 2) {
    stop("Not enough elements exceeding phi threshold")
  }

  grp <- max.col(filtered_phi)

  # Collect the relevant values
  pos <- data.frame(nam = nam,
                    grp = grp,
                    phi = filtered_phi[cbind(seq_along(grp), grp)],
                    p = filtered_p[cbind(seq_along(grp), grp)])

  # Filter by p
  pos_sig <- pos[pos$p <= p.max, ]
  if (nrow(pos_sig) < 2) {
    stop("Not enough elements exceeding p threshold")
  }
  
  # Sort the entire data frame by Group and Value
  df <- pos_sig[order(pos_sig$grp, -pos_sig$phi), ]

  # Add a width column
  df4plot <- df
  df4plot$order <- nrow(df4plot):1

  # Create a vector of colors
  # Determine the number of groups
  n_groups <- length(unique(df4plot$grp))
  colors <- rep(c("#606060", "#cecece"), length.out = n_groups)

  # Calculate label positions
  # Group the data by 'grp'
  grouped_df <- split(df4plot, df4plot$grp)

  # Map numeric group names to original cluster names
  name_map <- stats::setNames(cluster_names, seq_along(cluster_names))

  # Replace numeric group names with original cluster names in grouped_df
  grouped_df <- lapply(grouped_df, function(df) {
    df$grp <- name_map[as.character(df$grp)]
    df
  })

  # x position of group labels
  label_x <- max(df4plot$phi)

  # Initialize an empty list for the label positions
  label_positions_list <- vector("list", length(grouped_df))

  # Loop over groups
  for (i in seq_along(grouped_df)) {
    # Get the data for the current group
    group_data <- grouped_df[[i]]

    # Calculate the y label positions
    label_y <- min(group_data$order) + (max(group_data$order) -
                min(group_data$order)) / 2

    # Get the original group name
    original_grp_name <- name_map[names(grouped_df)[i]]

    # Add the label positions to the list
    label_positions_list[[i]] <- data.frame(grp = original_grp_name,
                                            label_y = label_y, 
                                            label_x = label_x)
  }
  # Combine the list into a data frame
  label_positions <- do.call(rbind, label_positions_list)

  # Significance column
  df4plot$sig <- ifelse(df4plot$p <= 0.001, "***",
                        ifelse(df4plot$p <= 0.01, "**",
                              ifelse(df4plot$p <= 0.05, "*", "")))
  # ggplot
  gg <- ggplot2::ggplot(df4plot, ggplot2::aes(x = phi, y = factor(nam,
                    levels = unique(df4plot$nam[order(df4plot$order)])))) +
                    ggplot2::geom_bar(ggplot2::aes(fill = factor(grp)),
                    stat = "identity") +
      ggplot2::scale_fill_manual(values = colors) +
      ggplot2::labs(title = title, x = "phi", y = NULL) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = text.size * 1.5,
                    face = "bold"),
        axis.title.x = ggplot2::element_text(size = text.size,
                    margin = ggplot2::margin(t = text.size / 2)),
        axis.text.x = ggplot2::element_text(size = text.size),
        axis.ticks.x = ggplot2::element_blank(),
        axis.text.y = if (labels) {
          ggplot2::element_text(size = text.size)
        } else { ggplot2::element_blank() },
        panel.grid.major.y = ggplot2::element_blank()#,
      ) +
      ggplot2::expand_limits(x = c(phi.min, 1)) +
      ggplot2::coord_cartesian(xlim = c(phi.min, 1)) +
      ggplot2::guides(fill = "none")

  # Get the maximum x limit of the plot
  plot_data <- ggplot2::ggplot_build(gg)
  max_x <- max(plot_data$layout$panel_params[[1]]$x.range)

  ## Add cluster labels
  gg <- gg + ggplot2::geom_text(data = label_positions,
      ggplot2::aes(x = max_x, y = label_y, label = grp),
      hjust = 1, vjust = 0.5, size = text.size * 0.352778)

  ## This is for later use when a proper p for phi is calculated
  #if (print.p) {
  #  gg <- gg + ggplot2::geom_text(ggplot2::aes(label = .data$sig,
  #                                             x = .data$phi - 0.02),
  #                                             hjust = 1,
  #                                             color = "white")
  #}

  plot(gg)
  return(invisible(gg))
}