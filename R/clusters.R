clusters <-  function(x, level = NULL, k = NULL, style = c("flat", "hierarchical")) {

  # Check if isopam object is provided
  if (!"isopam" %in% class(x)) {
    stop("x is not an object of class isopam")
  }

  # Handle style
  style <- match.arg(style)
  
  if (style == "flat") {
    flat <- TRUE
  } else if (style == "hierarchical") {
    flat <- FALSE
  } else {
    stop("style needs to be 'flat' or 'hierarchical'")
  }

  # Clustering solution hierarchical or not?
  hierarchical <- !is.null(x$dendro)

  # Adjust level if necessary
  if (hierarchical) {
    nlevels <- length(x$hier)

    if (!is.null(level)) {
      if (length(level) == 1) {
        if (level < 1) {
          level <- 1
          warning("Level adjusted to 1")
        } else if (level > nlevels) {
          level <- nlevels
          warning("Level adjusted to ", nlevels)
        }
      } else if (length(level) > 1) {
        if (min(level) > nlevels || max(level) < 1) {
          stop("No valid levels provided")
        }
        if (max(level) > nlevels) {
          level <- level[level <= nlevels]
          warning("Levels > ", nlevels, " are omitted")
        }
        if (min(level) < 1) {
          level <- level[level >= 1]
          warning("Levels < 1 are omitted")
        }
      }
    }
  }

  # Select and return by level
  if(!is.null(level)) {
    if (hierarchical) {
      if (length(level) > 1) {
        clustering <- as.data.frame(lapply(if (flat) x$flat[level]
                                           else x$hier[level], as.factor))
      } else {
        clustering <- as.factor(if (flat) x$flat[level][[1]]
                                else x$hier[level][[1]])
        names(clustering) <- rownames(x$hier)
      }
    } else { # if level is given but solution is not hierarchical
      clustering <- as.factor(x$flat)
      if (level > 1) {
        warning("Non-hierarchical cluster solution, level ignored")
      }
    }
  } else if(!is.null(k)) { # Select and return by k

    if(!is.null(level)) {
      stop("Provide either level or k, not both")
    }

    if(k < 2) {
      stop("k needs to be 2 or higher")
    }

    if (hierarchical) {

      # Numbers of clusters by level
      nc <- sapply(x$flat, function(x) length(unique(x)))

      # Which levels meet criterion k
      meeting_k <- nc %in% k

      if (any(meeting_k)) {
        # Return the levels that meet the criterion
        selected_levels <- which(meeting_k)
      } else {
        # Find the next lower value
        next_lower <- max(nc[nc < min(k)])
        selected_levels <- which(nc == next_lower)
      }

      # Return clustering matrix or vector
      clustering <- as.data.frame(lapply(if (flat) x$flat
                                         else x$hier, as.factor))
      clustering <- clustering[, selected_levels]
      names(clustering) <- rownames(x$hier)
    } else {  # if not hierarchical
      clustering <- as.factor(x$flat)
      warning("Non-hierarchical cluster solution, k ignored")
    }
  } else { # If neither level nor k are specified, output everything
    if (hierarchical) {
      clustering <- as.data.frame(lapply(if (flat) x$flat
                                         else x$hier, as.factor))
    } else {
      clustering <- as.factor(x$flat)
    }
  }

  return(clustering)
}
