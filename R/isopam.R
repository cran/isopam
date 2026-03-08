isopam <-  function(dat, c.fix = NULL, c.max = NULL, 
                    l.max = FALSE, stopat = c(1, 7), sieve = TRUE, 
                    Gs = 3.5, ind = NULL, centers = NULL, 
                    distance = "bray", k.max = 100, d.max = 7, 
                    juice = FALSE, polishing = c('strict', 'relaxed'),
                    ...) {

  # Check if the distance measure is available for vegdist or proxy
  vegdists <- c("manhattan", "euclidean", "canberra", "bray",
                "kulczynski", "gower", "morisita", "horn",
                "mountford", "jaccard", "raup", "binomial", "chao",
                "altGower", "cao", "mahalanobis", "clark", "chisq", 
                "chord", "hellinger", "aitchison", "robust.aitchison")
  method <- pmatch(distance, vegdists)
  distFunc <- "distFunc1"

  # If vegan::vegdist does not know the distance try with proxy::dist
  # and vegan::designdist
  if (is.na(method)) {
    if (requireNamespace("proxy", quietly = TRUE)) {
      proxydists <- rownames(as.data.frame(proxy::pr_DB))
      method <- distance %in% proxydists

      if (method) {
        distFunc <- "distFunc2"
        message("Using proxy::dist")
      } else {
        distFunc <- "distFunc3"
        message("Trying vegan::designdist")
      }    
    } else { 
        distFunc <- "distFunc3"
        message("proxy::dist not available, trying vegan::designdist")
    }
  }

  ## Make backwards compatible(this is mainly for use with Juice)
  if (!is.null(list(...)$fixed.number)) c.fix <- list(...)$fixed.number
  if (!is.null(list(...)$c.num)) c.fix <- list(...)$c.num ## old Juice vers.!
  if (!is.null(list(...)$opt.number)) c.opt <- list(...)$opt.number
  if (!is.null(list(...)$max.number)) c.max <- list(...)$max.number
  if (!is.null(list(...)$max.level)) l.max <- list(...)$max.level
  if (!is.null(list(...)$thresh)) Gs <- list(...)$thresh
  if (!is.null(list(...)$filtered)) sieve <- list(...)$filtered
  ## Note that Gs is set to a default of 30 in old Juice versions. The setting
  ## is ignored here

  # Check if the deprecated argument c.opt has been used
  args <- list(...)
  if ("c.opt" %in% names(args)) {
    warning("The 'c.opt' argument is deprecated. Use c.max = 2 instead.")
  }

  ## Prepare Juice session if applicable
  if (juice) {
    dir.create("isopam", showWarnings = FALSE)
  }

  ## Add fake "sample names" if necessary
  if (is.null(rownames(dat))) {
    rownames(dat) <- seq_len(nrow(dat))
  }

  ## Add fake "taxon names" if necessary
  if (is.null(colnames(dat))) {
    colnames(dat) <- seq_len(ncol(dat))
  }

  ## Convert to matrix if necessary while taking care of loosing rownames
  rownam <- rownames(dat)
  colnam <- colnames(dat)
  dat <- as.matrix(dat)
  rownames(dat) <- rownam
  dimdat_orig <- dim(dat)
  dat_orig <- dat # used later to extract unallocated sites and species

  # Match the polishing argument
  polishing <- match.arg(polishing)
  
  # Polishing
  if (polishing == 'strict') {
    repeat {
      IO <- ifelse(dat > 0, 1, 0)
      initial_rows <- nrow(dat)
      initial_cols <- ncol(dat)
      dat <- dat[, colSums(IO) > 1] # Species with > 1 occurrence
      dat <- dat[rowSums(IO) > 1, ] # ... and plots with > 1 species
      dat <- dat[, apply(dat, 2, var) > 0]
      if (nrow(dat) == initial_rows && ncol(dat) == initial_cols) {
        break
      }
    }
  } else if (polishing == 'relaxed') {
    dat <- dat[, colSums(dat) > 0] # Omit species without occurrence
    dat <- dat[, apply(dat, 2, var) > 0] # Omit species with no variance
    dat <- dat[rowSums(dat) > 0, ] # ... and rows without species
  }

  dimdat <- dim(dat)

  if (is.null(dimdat) | dimdat[1] < 3 | dimdat[2] < 2) {
      stop("Not enough valid rows or columns")
  }

  ## In case of predefined indicator species check their validity
  usr_ind <- FALSE # Initialize usr_ind
  if (!is.null(ind)) {
    if (all(ind %in% colnames(dat))) {
      usr_ind <- TRUE
    } else {
      stop("Predefined indicators not valid")
    }
  }

  ## In case of predefined centers check their validity
  semi_supervised <- FALSE
  semi_supervised_search <- FALSE
  n_fixed_centers <- 0
  
  ## Store whether c.fix/c.max were explicitly provided by user
  c.fix_user <- !is.null(c.fix)
  c.max_user <- !is.null(c.max)

  if (!is.null(centers)) {

    if (is.character(centers)) {
      if (!all(centers %in% rownames(dat))) {
        stop("Centers not valid")
      } else {
        # Transform into indices (preserving order)
        centers <- match(centers, rownames(dat))
      }
    }
    if (is.numeric(centers) && max(centers) > nrow(dat)) {
      stop("Centers not in range")
    }

    n_fixed_centers <- length(centers)
    
    # Validate: need at least 2 centers for clustering
    if (n_fixed_centers < 2 && !c.fix_user && !c.max_user) {
      stop("At least 2 centers required for supervised mode")
    }
  }

  ## ----------- Mode detection and validation ----------- ##
  
  ## Case 1: centers NULL, c.fix set - standard partitioning
  if (is.null(centers) && c.fix_user) {
    if (c.fix < 2) stop("c.fix must be >= 2")
    if (c.fix > nrow(dat) - 1) c.fix <- nrow(dat) - 1
    l.max <- 1
    c.min <- c.fix
    c.max <- c.fix
  
  ## Case 2: centers set, c.fix FALSE, c.max NULL - fully supervised
  } else if (!is.null(centers) && !c.fix_user && !c.max_user) {
    if (n_fixed_centers < 2) stop("At least 2 centers required")
    l.max <- 1
    c.min <- n_fixed_centers
    c.max <- n_fixed_centers
  
  ## Case 3: centers set, c.fix set - semi-supervised fixed (c.max ignored)
  } else if (!is.null(centers) && c.fix_user) {
    if (c.max_user) warning("Both c.fix and c.max provided; c.max ignored")
    if (c.fix < n_fixed_centers) stop("c.fix must be >= number of centers")
    if (c.fix < 2) stop("c.fix must be >= 2")
    if (c.fix > nrow(dat) - 1) c.fix <- nrow(dat) - 1
    l.max <- 1
    c.min <- c.fix
    c.max <- c.fix
    if (c.fix > n_fixed_centers) {
      semi_supervised <- TRUE
      message("Semi-supervised mode: ", c.fix - n_fixed_centers,
              " additional cluster(s)")
    }
  
  ## Case 4: centers set, c.max set - semi-supervised search
  } else if (!is.null(centers) && c.max_user) {
    if (c.max < n_fixed_centers) stop("c.max must be >= number of centers")
    if (c.max < 2) stop("c.max must be >= 2")
    if (c.max > nrow(dat) - 1) c.max <- nrow(dat) - 1
    l.max <- 1
    c.min <- max(2, n_fixed_centers)  # At least 2 clusters
    if (c.max > n_fixed_centers) {
      semi_supervised <- TRUE
      semi_supervised_search <- TRUE
      message("Semi-supervised search mode: up to ", c.max - n_fixed_centers,
              " additional cluster(s)")
    }
  
  ## Default case: no centers, no c.fix - hierarchical or default partitioning
  } else {
    c.min <- 2
    if (!c.max_user) c.max <- 6
  }

  ## Initiate count
  count <- 1

  ## Initialize parallel workers once (outside core to avoid repeated startup)
  parallel_initialized <- FALSE

  ## ----------- core function ---------------------------------------------- ##  
  
  core <- function(xdat) {

    IO.xdat <- ifelse(xdat > 0, 1L, 0L)

    ## Some useful descriptors
    N.xdat <- nrow(xdat)                         ## Total number of plots
    SP.xdat <- ncol(xdat)                        ## Total number of species
    frq.xdat <- t(as.matrix(colSums(IO.xdat)))   ## Species frequencies

    ## For Williams` correction
    w3 <- N.xdat *((1 / frq.xdat) +(1 /(N.xdat - frq.xdat))) - 1

    ## In case of predefined indicators: which columns?
    if (usr_ind) {
      xind <- which(colnames(xdat) %in% ind)
    }

    ## Distance matrix
    if (distFunc == "distFunc1") {
      dst.xdat <- vegan::vegdist(xdat, method = distance)
    }  
    if (distFunc == "distFunc2") {
      dst.xdat <- proxy::dist(xdat, method = distance)
    }
    if (distFunc == "distFunc3") {
      dst.xdat <- vegan::designdist(xdat, method = distance)
    }

    ## Exit with dignity if N.xdat < 3
    if (N.xdat < 3) {
      mess1 <- "Not enough data."
      if (juice) {
        write.table(mess1, "isopam/alert.txt", row.names = FALSE,
                     col.names = FALSE, quote = FALSE)
      }
      stop(mess1)
    }

    ## Determine maximum Isomap k
    if (k.max > N.xdat - 1) {
      k.max <- N.xdat - 1
    }

    ## Determine minimum Isomap k(from isomapdist.r, vegan)
    dmtr <- as.matrix(dst.xdat)
    diag(dmtr) <- NA

    k.min <- 2

    for (a in c(2:k.max)) {

      ## Check min. k
      dm <- dmtr

      is.na(dm) <- apply(dm, 2, function(xx) xx >
                         xx[order(xx, na.last = TRUE)[a]])
      dm <- pmax(as.dist(dm), as.dist(t(dm)), na.rm = TRUE)
      fragm <- vegan::distconnected(dm, toolong = 0,  trace = FALSE)
      if (length(unique(fragm)) > 1) {
        k.min <- k.min + 1
      } else break
    }

    ## Adjust d.max if necessary
    if (d.max > N.xdat - 1) {
      d.max <- N.xdat - 1
    }

    ## Adjust c.max if necessary
    if (c.max > N.xdat - 1) c.max <- N.xdat - 1
    if (c.max < 2) stop("c.max < 2")

    # Criteria for parallel processing (subject to further experiments)
    fut <- TRUE # Stays TRUE if parallel processing seems to be possible
    rg.k <- k.max - k.min
    # Check current memory limit
    future_memlim <- unlist(options("future.globals.maxSize"))
    # In case there is no explicit limit, the default applies:
    if (is.null(future_memlim)) {
      future_memlim <- 524288000
    }

    # Condition: Data set exceeds certain minimum size
    if (rg.k <= 50) {
      fut <- FALSE
    }

    # Reporting
    if (!already_shown) {
      already_shown <<- TRUE
    }

    if (count == 1) {
      message("Level 1: Partitioning .", appendLF = FALSE)
    }

    ## ----------- Note ------------------------------ ##

    ## A check for slow or fast mode in the inner loops
    ## would slow down the process by 10% -20%.

    ## ----------- Parallel processing --------------- ##

    if (fut) {

      ## Initialize workers only once (reuse across core() calls)
      if (!parallel_initialized) {
        future::plan(future::multisession)
        parallel_initialized <<- TRUE
      }

      ## b-loop: Isomap k
      ## Wrap entire future_sapply to suppress package loading warnings in workers
      ## Use future.packages to explicitly load dependencies in workers
      out.array <- suppressWarnings(array(future.apply::future_sapply(k.min:k.max, 
                         future.seed = TRUE, 
                         future.packages = c("vegan", "cluster", "fastkmedoids"),
                         function(b) {

          ## Isomap
          suppressMessages(isom <- isomap(dst.xdat, ndim = d.max, k = b))

          ## Fixing the maximum of dimensions considered when calculating
          ## the distance matrix for the isomap space
          d.max.new <- min(sum(isom$eig > 0), ncol(isom$points), d.max, na.rm = T)
          out.mat <- matrix(NA, nrow = d.max - 1, ncol = c.max - c.min + 1 )

          ## d-Loops
          if (d.max.new > 1) {
            
            # Optimization: Pre-standardize points to avoid repeated work in daisy
            # daisy(stand=TRUE) uses Mean Absolute Deviation (MAD)
            points_std <- apply(isom$points[, 1:d.max.new, drop = FALSE], 2, function(x) {
              mad_val <- mean(abs(x - mean(x)))
              if (mad_val < 1e-12) return(rep(0, length(x)))
              (x - mean(x)) / mad_val
            })
            
            # Initialize squared distance matrix with 1st dimension
            dist_sq <- outer(points_std[, 1], points_std[, 1], "-")^2

            for (d in 2:d.max.new) {

              # Incrementally add squared differences of the new dimension
              dist_sq <- dist_sq + outer(points_std[, d], points_std[, d], "-")^2
              isodiss <- as.dist(sqrt(dist_sq))

              # For fastpam
              isodiss_vector <- as.vector(isodiss)
              # For semi-supervised mode (computed once per d-level)
              dmat_iso <- if (semi_supervised) as.matrix(isodiss) else NULL
              # Cache for incremental medoid selection across e values
              prev_optimized_medoids <- NULL

              for (e in c.min:c.max) { ## e-loop: Cluster no.

                ## --------- Partitioning (PAM) ------------------------------- #

                if (!is.null(centers)) {
                  if (semi_supervised && e > n_fixed_centers) {
                    # Semi-supervised: fixed centers + additional free medoids
                    n_additional <- e - n_fixed_centers
                    if (is.null(prev_optimized_medoids)) {
                      # First semi-supervised iteration: full initialization
                      additional_medoids <- select_additional_medoids(dmat_iso, 
                                              centers, n_additional)
                      initial_medoids <- c(centers, additional_medoids)
                    } else {
                      # Incremental: add one more medoid to previous result
                      new_medoid <- select_additional_medoids(dmat_iso, 
                                      prev_optimized_medoids, 1L)
                      initial_medoids <- c(prev_optimized_medoids, new_medoid)
                    }
                    # Apply restricted swap (only free medoids can move)
                    optimized_medoids <- restricted_pam_swap(dmat_iso, 
                                          initial_medoids, n_fixed_centers)
                    prev_optimized_medoids <- optimized_medoids
                    # Fast assignment without full PAM call
                    cl <- assign_to_medoids(dmat_iso, optimized_medoids)
                    ci <- tabulate(cl, nbins = e)
                  } else {
                    # Standard supervised: only fixed centers
                    cl.iso <- cluster::pam(isodiss, k = e, medoids = centers,
                                  diss = TRUE, do.swap = FALSE)
                    cl <- cl.iso$clustering
                    ci <- cl.iso$clusinfo[, 1]
                  }
                } else {
                  cl.iso <- fastkmedoids::fastpam(isodiss_vector, k = e, n = N.xdat)
                  cl <- cl.iso@assignment
                  ci <- tabulate(cl, nbins = e)
                }

                ## ----------- Start vectorized calculation --------------- ##

                # Calculate G-values using the helper function
                gt <- calc_G_vectorized(IO.xdat, cl, ci, frq.xdat, N.xdat, w3, e)

                ## Standardization (Botta-Dukat et al. 2005)
                gt.ex <- e - 1                 ## Expected G
                gt.sd <- sqrt(2 * gt.ex)       ## Expected sd
                G <- (gt - gt.ex) / gt.sd

                ## Using predefined indicators
                if (usr_ind) {
                  glgth <- length(G[xind])
                  if (glgth == 0) {
                    out.mat[d - 1, e + 1 - c.min] <- NA
                  } else {
                    out.mat[d - 1, e + 1 - c.min] <- mean(G[xind])
                  }
                }

                if (!sieve && !usr_ind) {
                  out.mat[d - 1, e + 1 - c.min] <- mean(G)
                }

                ## Standard: Filtering and averaging
                if (sieve && !usr_ind) {
                  ## Filtering by G
                  glgth <- length(G[G >= Gs])
                  if (glgth == 0) {
                    out.mat[d - 1, e + 1 - c.min] <- NA
                  } else {
                    ## Averaging
                    out.mat[d - 1, e + 1 - c.min] <- mean(G[G >= Gs]) * glgth
                  }
                }
                
                ## ----------- End vectorized calculation --------------- ##

              }
            }
          }

        return(out.mat)

      }), dim = c(d.max - 1, c.max - c.min + 1, k.max - k.min + 1))) # End suppressWarnings

    ## ----------- End parallel mode --------------------- ##

    } else {

    ## ----------- Non-parallel processing --------------- ##

      ## b-loop: Isomap k
      ## Wrap entire sapply to suppress package loading warnings
      out.array <- suppressWarnings(array(sapply(k.min:k.max, function(b) {

          ## Isomap
          suppressMessages(isom <- isomap(dst.xdat, ndim = d.max, k = b))

          ## Fixing the maximum of dimensions considered when calculating
          ## the distance matrix for the isomap space
          d.max.new <- min(sum(isom$eig > 0), ncol(isom$points), d.max, na.rm = T)
          out.mat <- matrix(NA, nrow = d.max - 1, ncol = c.max - c.min + 1)

          ## d-Loops
          if (d.max.new > 1) {
            
            # Optimization: Pre-standardize points to avoid repeated work in daisy
            # daisy(stand=TRUE) uses Mean Absolute Deviation (MAD)
            points_std <- apply(isom$points[, 1:d.max.new, drop = FALSE], 2, function(x) {
              mad_val <- mean(abs(x - mean(x)))
              if (mad_val < 1e-12) return(rep(0, length(x)))
              (x - mean(x)) / mad_val
            })
            
            # Initialize squared distance matrix with 1st dimension
            dist_sq <- outer(points_std[, 1], points_std[, 1], "-")^2

            for (d in 2:d.max.new) {

              # Incrementally add squared differences of the new dimension
              dist_sq <- dist_sq + outer(points_std[, d], points_std[, d], "-")^2
              isodiss <- as.dist(sqrt(dist_sq))

              # For fastpam
              isodiss_vector <- as.vector(isodiss)
              # For semi-supervised mode (computed once per d-level)
              dmat_iso <- if (semi_supervised) as.matrix(isodiss) else NULL
              # Cache for incremental medoid selection across e values
              prev_optimized_medoids <- NULL

              for (e in c.min:c.max) { ## e-loop: Cluster no.

                ## --------- Partitioning(PAM) ------------------------------ #

                if (!is.null(centers)) {
                  if (semi_supervised && e > n_fixed_centers) {
                    # Semi-supervised: fixed centers + additional free medoids
                    n_additional <- e - n_fixed_centers
                    if (is.null(prev_optimized_medoids)) {
                      # First semi-supervised iteration: full initialization
                      additional_medoids <- select_additional_medoids(dmat_iso, 
                                              centers, n_additional)
                      initial_medoids <- c(centers, additional_medoids)
                    } else {
                      # Incremental: add one more medoid to previous result
                      new_medoid <- select_additional_medoids(dmat_iso, 
                                      prev_optimized_medoids, 1L)
                      initial_medoids <- c(prev_optimized_medoids, new_medoid)
                    }
                    # Apply restricted swap (only free medoids can move)
                    optimized_medoids <- restricted_pam_swap(dmat_iso, 
                                          initial_medoids, n_fixed_centers)
                    prev_optimized_medoids <- optimized_medoids
                    # Fast assignment without full PAM call
                    cl <- assign_to_medoids(dmat_iso, optimized_medoids)
                    ci <- tabulate(cl, nbins = e)
                  } else {
                    # Standard supervised: only fixed centers
                    cl.iso <- cluster::pam(isodiss, k = e, medoids = centers,
                                  diss = TRUE, do.swap = FALSE)
                    cl <- cl.iso$clustering
                    ci <- cl.iso$clusinfo[, 1]
                  }
                } else {
                  cl.iso <- fastkmedoids::fastpam(isodiss_vector, k = e, n = N.xdat)
                  cl <- cl.iso@assignment
                  ci <- tabulate(cl, nbins = e)
                }
                
                ## ----------- Start vectorized calculation --------------- ##

                # Calculate G-values using the helper function
                gt <- calc_G_vectorized(IO.xdat, cl, ci, frq.xdat, N.xdat, w3, e)

                ## Standardization (Botta-Dukat et al. 2005)
                gt.ex <- e - 1                 ## Expected G
                gt.sd <- sqrt(2 * gt.ex)       ## Expected sd
                G <- (gt - gt.ex) / gt.sd

                ## Using predefined indicators
                if (usr_ind) {
                  glgth <- length(G[xind])
                  if (glgth == 0) {
                    out.mat[d - 1, e + 1 - c.min] <- NA
                  } else {
                    out.mat[d - 1, e + 1 - c.min] <- mean(G[xind])
                  }
                }

                if (!sieve && !usr_ind) {
                  out.mat[d - 1, e + 1 - c.min] <- mean(G)
                }

                ## Standard: Filtering and averaging
                if (sieve && !usr_ind) {
                  ## Filtering by G
                  glgth <- length(G[G >= Gs])
                  if (glgth == 0) {
                    out.mat[d - 1, e + 1 - c.min] <- NA
                  } else {
                    ## Averaging
                    out.mat[d - 1, e + 1 - c.min] <- mean(G[G >= Gs]) * glgth
                  }
                }
                
                ## ----------- End vectorized calculation --------------- ##

              }
            }
          }

        return(out.mat)
      }),dim = c(d.max - 1, c.max - c.min + 1, k.max - k.min + 1))) # End suppressWarnings

    ## ----------- End non-parallel mode --------------- ##
    }

    ## ----------- End parameter search ------------------------------------- ##

    solution <- TRUE
    out.array[is.na(out.array)] <- 0

    if (length(out.array[out.array > 0]) == 0) {
      mess2 <- "No solution found with current settings"

      if (count == 1) {

        if (juice) {
          write.table(mess2, "isopam/alert.txt",
                       row.names = FALSE, col.names = FALSE, quote = FALSE)
        }
        stop("No solution found with current settings")

      } else {
        solution <- FALSE
      }
    }

    if (solution) {

      mn.iso <- max(out.array)

      ## ----------- Parameters for final run ------------------------------- ##

      wmx.iso <- which(out.array == mn.iso, arr.ind = TRUE) ## Cases with max
      colnames(wmx.iso) <- c("iso.dim", "clusters", "iso.k")

      ## In case of multiple best solutions select the one with max clusters,
      ## max dims and max k (in this rank order)

      ## Select cases with maximum number of clusters
      try(wmx.iso <- wmx.iso[which(wmx.iso[, 2] == max
                                     (wmx.iso[, 2])), ], silent = TRUE)
      ## Select cases with maximum number of Isomap dimensions
      try(wmx.iso <- wmx.iso[which(wmx.iso[, 1] == max
                                     (wmx.iso[, 1])), ], silent = TRUE)
      ## Select cases with maximum k
      try(wmx.iso <- wmx.iso[which(wmx.iso[, 3] == max
                                     (wmx.iso[, 3])), ], silent = TRUE)

      mc <- wmx.iso[2] - 1 + c.min
      md <- wmx.iso[1] + 1
      mk <- wmx.iso[3] - 1 + k.min

      ## ----------- Final run ---------------------------------------------- ##

      suppressMessages(isom <- isomap(dst.xdat, ndim = d.max, k = mk))
      d.iso <- daisy(isom$points[, 1:md], metric = "euclidean", stand = TRUE)

      if (!is.null(centers)) {
        if (semi_supervised && mc > n_fixed_centers) {
          # Semi-supervised: fixed centers + additional free medoids
          dmat_final <- as.matrix(d.iso)
          n_additional <- mc - n_fixed_centers
          additional_medoids <- select_additional_medoids(dmat_final, 
                                  centers, n_additional)
          initial_medoids <- c(centers, additional_medoids)
          # Apply restricted swap (only free medoids can move)
          optimized_medoids <- restricted_pam_swap(dmat_final, 
                                initial_medoids, n_fixed_centers)
          # Final assignment
          CLS <- assign_to_medoids(dmat_final, optimized_medoids)
          MDS <- optimized_medoids
          CLI <- tabulate(CLS, nbins = mc)
        } else {
          # Standard supervised: only fixed centers
          cl.iso <- cluster::pam(d.iso, k = mc, medoids = centers,
                        diss = TRUE, do.swap = FALSE)
          CLS <- cl.iso$clustering
          MDS <- cl.iso$medoids
          CLI <- t(cl.iso$clusinfo[, 1])
        }
      } else {
        cl.iso <- cluster::pam(d.iso, k = mc, diss = TRUE)
        CLS <- cl.iso$clustering
        MDS <- cl.iso$medoids
        CLI <- t(cl.iso$clusinfo[, 1])
      }

      ## ------------ method == G ------------------------ #

      # Calculate G-values using the vectorized helper function
      # CLS = cluster assignments, CLI = cluster sizes, mc = number of clusters
      g.1 <- calc_G_vectorized(IO.xdat, CLS, CLI, frq.xdat, N.xdat, w3, mc)

      ## Standardization (Botta-Dukat et al. 2005)
      gt.ex <- mc - 1                       ## Expected G
      gt.sd <- sqrt(2 * (mc - 1))           ## Expected sd
      sG <- (g.1 - gt.ex) / gt.sd

      ## Some analytical output
      ## Averaged G
      ivx <- round(mean(sG), 1)
      if (usr_ind) {
        ## Averaged G (only indicators)
        ivi <- round(mean(sG[xind]), 1)
        ## Number of indicators >= Gs
        noi <- length(sG[xind])
        ## Indicator identities
        INDN <- colnames(IO.xdat)[xind]
      }
      else if (sieve && !usr_ind) {
        ## Averaged G (only indicators)
        ivi <- round(mean(sG[sG >= Gs]), 1)
        ## Number of indicators >= Gs
        noi <- length(sG[sG >= Gs])
        ## Indicator identities
        INDN <- colnames(IO.xdat)[sG >= Gs]
      }
      else if (!sieve && !usr_ind) {
        ## Averaged G (only indicators)
        ivi <- "NA"
        ## Number of indicators >= Gs
        noi <- "NA"
        ## Indicator identities
        INDN <- "NA"
      }

      ## Was this a good partition?
      ## At least stopcrit[1] descriptors with g >= stopcrit[2]
      if (length(sG[sG >= stopat[2]]) >= stopat[1] * mc) {
        fine <- TRUE
      } else {
        fine <- FALSE
      }

      ## ---------------- preparing output ------------------#

      out <- list(
        medoids = MDS,
        clusters = CLS,
        sizes = CLI,
        is.ok = fine,
        k.min = k.min,
        k.max = k.max,
        k = mk,
        d = md,
        noi = noi,
        ivx = ivx,
        ivi = ivi,
        indnames = INDN)

    } else {

      out <- list(
        medoids = NULL,
        clusters = NULL,
        sizes = NULL,
        is.ok = FALSE,
        k.min = NULL,
        k.max = NULL,
        k = NULL,
        d = NULL,
        noi = NULL,
        ivx = NULL,
        ivi = NULL,
        indnames = NULL )
    }
    return(out)

  }

  ## ----------- End core function ----------------------------------------- ##
  ## ----------- Dendrogram function(code: J. Collison, 2009) -------------- ##

  create_dendro <- function(clust) {
    ## Expects a list of vectors containing cluster affiliations
    ## without info about hierarchy(running number).
    ## Returns an object of class "hclust"

    num_obs <- length(clust[[1]])

    dendro <- list(
      merge  = array(NA, dim = c(num_obs - 1, 2)),
      height = rep(0, times = num_obs - 1),
      order  = c(),
      labels = names(clust[[1]]),
      method = NULL,
      call   = NULL)

    class(dendro) <- "hclust"


    ## hclust requires a set of merge operations, one observation at a time.
    ## Subsequent operations [may] refer to the index of a previous merge,
    ## so we store the operation index after everything we do.

    opnum <- 0

    group_opnums <- c()

    for (level in length(clust):1) {

      groups <- clust[[level]]

      groupnum <- 0

      curlevel_group_opnums <- c()

      repeat {

        groupnum <- groupnum + 1

        log_in_group <- (groups == groupnum)

        num_in_clust <- sum(log_in_group);

        if (num_in_clust < 1) break


        if (level == length(clust)) {
          ## Bottom  level
          ## Join all these and add them to the order vector

          prev_opnum <- 0
          prev_index <- 0

          for (j in 1:num_obs) {

            if (!log_in_group[j]) next

            dendro$order <- c(dendro$order, j)

            if (prev_opnum) {
              opnum <- opnum + 1
              dendro$merge[opnum, ] <- c(-j, prev_opnum)
              dendro$height[opnum] <- 1
              prev_opnum <- opnum
            } else if (prev_index) {
              opnum <- opnum + 1
              dendro$merge[opnum,] <- c(-j, -prev_index)
              dendro$height[opnum] <- 1
              prev_opnum <- opnum
            } else {
              prev_index <- j
            }
          }

          ## Done merging this cluster
          ## Save the opnum for the last merge in this cluster
          curlevel_group_opnums <- c(curlevel_group_opnums,0)
          curlevel_group_opnums[groupnum] <- prev_opnum

          ## Special case for singletons
          if (!prev_opnum) {
            curlevel_group_opnums[groupnum] <- -prev_index
          }

        } else {

          ## higher levels

          ## For all members of this group, see what their groupid was
          ## one level deeper and store those(uniquely)

          groups_to_join <- c()

          for (j in 1:num_obs) {
            if (!log_in_group[j]) next

            subgroup <- clust[[level + 1]][j]

            if (sum(groups_to_join == subgroup) == 0) {
              ## Not there, add it
              groups_to_join <- c(groups_to_join,subgroup)
            }
          }

          if (length(groups_to_join) >= 2) {
            ## Merge them.
            prev_opnum <- group_opnums[groups_to_join[1]];
            for (r in 2:length(groups_to_join)) {
              opnum <- opnum + 1
              dendro$merge[opnum,] <- c(prev_opnum,
                                        group_opnums[groups_to_join[r]])
              dendro$height[opnum] <- length(clust) - level + 1
              prev_opnum <- opnum
            }

            curlevel_group_opnums <- c(curlevel_group_opnums, 0)
            curlevel_group_opnums[groupnum] <- prev_opnum
          } else {
            ## Nothing to merge this one with at this time;
            ## bubble up for next level
            prev_opnum <- group_opnums[groups_to_join[1]]
            curlevel_group_opnums <- c(curlevel_group_opnums, 0)
            curlevel_group_opnums[groupnum] <- prev_opnum
          }
        }
      } ## end groupnum loop

      ## Done processing this level
      group_opnums <- curlevel_group_opnums
    }

    ## All levels have been processed.
    ## Need to merge whatever remains

    if (length(group_opnums) > 1) {

      prev_opnum <- group_opnums[1]

      for (r in 2:length(group_opnums)) {
        opnum <- opnum + 1
        dendro$merge[opnum,] <- c(prev_opnum,group_opnums[r])
        dendro$height[opnum] <- length(clust) + 1
        prev_opnum <- opnum
      }
    }

    return(dendro)
  }
  
  ## ------------ End dendro function -------------------------------------- ##
  ##------------------ Helpers --------------------------------------------- ##

  ## Helper function to select additional medoids (deterministic, distance-based)
  select_additional_medoids <- function(dmat, fixed_centers, n_additional) {
    if (n_additional <= 0) return(integer(0))
    n_samples <- nrow(dmat)
    available <- setdiff(1:n_samples, fixed_centers)
    if (length(available) <= n_additional) return(available)
    
    # Select medoids that are far from fixed centers (greedy maximin)
    selected <- integer(n_additional)
    all_selected <- fixed_centers
    
    for (i in 1:n_additional) {
      # For each available point, find min distance to already selected
      min_dists <- rowMins(dmat[available, all_selected, drop = FALSE])
      # Select the one with maximum min-distance
      best_idx <- which.max(min_dists)
      new_medoid <- available[best_idx]
      selected[i] <- new_medoid
      all_selected <- c(all_selected, new_medoid)
      available <- available[-best_idx]
    }
    return(selected)
  }

  ## Fast row-wise minimum (avoid apply overhead)
  rowMins <- function(m) {
    if (ncol(m) == 1L) return(m[, 1L])
    result <- m[, 1L]
    for (j in 2:ncol(m)) result <- pmin(result, m[, j])
    result
  }

  ## Restricted swap for semi-supervised mode - optimized version
  ## Only swaps free medoids while keeping fixed centers locked
  restricted_pam_swap <- function(dmat, initial_medoids, n_fixed, max_iter = 100) {
    n <- nrow(dmat)
    k <- length(initial_medoids)
    n_free <- k - n_fixed
    
    if (n_free == 0L) return(initial_medoids)
    
    current_medoids <- initial_medoids
    non_medoid_set <- setdiff(1:n, current_medoids)
    
    # Pre-compute distances from all points to current medoids
    D <- dmat[, current_medoids, drop = FALSE]
    
    # For each point: index of closest medoid and the distance
    closest_idx <- max.col(-D, ties.method = "first")
    closest_dist <- D[cbind(1:n, closest_idx)]
    current_cost <- sum(closest_dist)
    
    for (iter in 1:max_iter) {
      improved <- FALSE
      
      # Try swapping each free medoid
      for (fi in 1:n_free) {
        pos <- n_fixed + fi  # position in medoid vector
        old_med <- current_medoids[pos]
        
        best_cand <- old_med
        best_cost <- current_cost
        
        # Pre-compute: for points whose closest is this medoid, 
        # what's their second-best distance?
        affected <- which(closest_idx == pos)
        second_best <- if (length(affected) > 0 && k > 1) {
          D_others <- D[affected, -pos, drop = FALSE]
          rowMins(D_others)
        } else numeric(0)
        
        # Vectorized: compute cost for all candidates at once
        n_cands <- length(non_medoid_set)
        if (n_cands == 0) next
        
        d_to_cands <- dmat[, non_medoid_set, drop = FALSE]  # n x n_cands
        
        # For unaffected points: pmin(closest_dist, d_to_cand) summed per candidate
        unaffected <- which(closest_idx != pos)
        n_unaff <- length(unaffected)
        if (n_unaff > 0) {
          d_unaff <- matrix(d_to_cands[unaffected, ], nrow = n_unaff, ncol = n_cands)
          cd_unaff <- closest_dist[unaffected]
          # pmin broadcasts vector across columns correctly
          contrib_unaffected <- colSums(pmin(d_unaff, cd_unaff))
        } else {
          contrib_unaffected <- rep(0, n_cands)
        }
        
        # For affected points: pmin(second_best, d_to_cand)
        n_aff <- length(affected)
        if (n_aff > 0) {
          d_aff <- matrix(d_to_cands[affected, ], nrow = n_aff, ncol = n_cands)
          contrib_affected <- colSums(pmin(d_aff, second_best))
        } else {
          contrib_affected <- rep(0, n_cands)
        }
        
        costs <- contrib_unaffected + contrib_affected
        
        best_idx <- which.min(costs)
        if (costs[best_idx] < best_cost) {
          best_cost <- costs[best_idx]
          best_cand <- non_medoid_set[best_idx]
        }
        
        # Apply best swap for this position if improvement found
        if (best_cand != old_med) {
          # Update medoids
          current_medoids[pos] <- best_cand
          non_medoid_set <- setdiff(1:n, current_medoids)
          
          # Update distance matrix and closest info
          D[, pos] <- dmat[, best_cand]
          closest_idx <- max.col(-D, ties.method = "first")
          closest_dist <- D[cbind(1:n, closest_idx)]
          current_cost <- best_cost
          improved <- TRUE
        }
      }
      
      if (!improved) break
    }
    
    return(current_medoids)
  }

  ## Fast assignment of points to medoids (avoids full PAM call)
  assign_to_medoids <- function(dmat, medoids) {
    D <- dmat[, medoids, drop = FALSE]
    max.col(-D, ties.method = "first")
  }

  ## Vectorized G-value calculation
  calc_G_vectorized <- function(IO.xdat, cl, ci, frq.xdat, N.xdat, w3, e) {
    SP.xdat <- ncol(IO.xdat)
    
    # Calculate fra1 matrix (rows=species, cols=clusters)
    # Optimization: Use crossprod for frequency calculation instead of rowsum
    # Construct indicator matrix for clusters (N.xdat x e)
    C_mat <- matrix(0, nrow = N.xdat, ncol = e)
    C_mat[cbind(seq_len(N.xdat), cl)] <- 1
    # crossprod(A, B) computes t(A) %*% B -> (S x N) * (N x e) -> S x e
    fra1_mat <- crossprod(IO.xdat, C_mat)
    
    # Expand vectors to matrices for element-wise operations
    Nj_mat <- matrix(ci, nrow = SP.xdat, ncol = e, byrow = TRUE)
    bom_vec <- as.vector(frq.xdat) / N.xdat
    bom_mat <- matrix(bom_vec, nrow = SP.xdat, ncol = e, byrow = FALSE)
    bim_mat <- 1 - bom_mat
    
    fra0_mat <- Nj_mat - fra1_mat
    
    # Calculate expected frequencies
    bum_mat <- fra1_mat / (Nj_mat * bom_mat)
    bam_mat <- fra0_mat / (Nj_mat * bim_mat)
    
    # Calculate Log-Likelihood terms (handling 0 * log(0) = 0)
    term1 <- matrix(0, nrow = SP.xdat, ncol = e)
    idx1 <- fra1_mat > 0 & bum_mat > 0 & !is.na(bum_mat)
    term1[idx1] <- fra1_mat[idx1] * log(bum_mat[idx1])
    
    term2 <- matrix(0, nrow = SP.xdat, ncol = e)
    idx2 <- fra0_mat > 0 & bam_mat > 0 & !is.na(bam_mat)
    term2[idx2] <- fra0_mat[idx2] * log(bam_mat[idx2])
    
    DDD_vec <- rowSums(term1 + term2) * 2
    
    # Williams correction
    w1 <- N.xdat * sum(1 / ci) - 1
    w2 <- 6 * N.xdat * (e - 1)
    willi_vec <- 1 + ((w1 * as.vector(w3)) / w2)
    
    # Return G-values as a column matrix
    as.matrix(DDD_vec / willi_vec)
  }

  ## ------------ End Helpers ----------------------------------------------- ##
  ## ------------ Isopam call ----------------------------------------------- ##

  ## Initiate container for results
  matr <- matrix(NA, nrow = nrow(dat), ncol = 1) ## Initiate cluster output
  rownames(matr) <- rownames(dat)
  colnames(matr) <- "lev.1"

  ## Initiate container for indicators
  indlist <- list()

  ## Initiate matrix for summary("analytics")
  summ <- matrix(NA, nrow = 9, ncol = 1)
  rownames(summ) <- c("Name",
                       "Subgroups",
                       "Isomap.dim",
                       "Isomap.k.min",
                       "Isomap.k",
                       "Isomap.k.max",
                       "Ind.N",
                       "Ind.Gs",
                       "Global.Gs")
  colnames(summ) <- "Part.1"

  # This initializes a variable that prevents a message from showing
  # up multiple times in the core function
  already_shown <- FALSE

  ## Run core function
  output <- core(dat)

  ## Fill cluster container
  matr[, 1] <- output$clusters

  ## Fill indicator container
  indlist$"Part.1" <- as.character(output$indnames)

  ## Fill summary matrix
  summ[1,1] <- 0                          ## Name
  summ[2,1] <- length(output$medoids)     ## No. of subgroups
  summ[3,1] <- output$d                   ## Isomap dimensions
  summ[4,1] <- output$k.min               ## Minimum k
  summ[5,1] <- output$k                   ## Isomap k
  summ[6,1] <- output$k.max               ## Maximum k
  summ[7,1] <- output$noi                 ## No. of indicators used
  summ[8,1] <- format(output$ivi, digits = 3)  ## Mean sG of indicators
  summ[9,1] <- format(output$ivx, digits = 3)  ## Mean sG of all descriptors

  ## Medoids
  med <- list()
  med[[1]] <- output$medoids
  names(med[[1]]) <- seq_along(output$medoids)

  ## stop if max.level is 1
  stepdown <- ifelse(l.max == 1, FALSE, TRUE)
  ## Which group is large enough?
  spl <- as.numeric(output$sizes > 2)
  ## stop if there is nothing to split
  if (sum(spl) == 0) {
    stepdown <- FALSE
  }

  ## Preparative stuff
  ctb <- matr ## Cluster table
  colnames(ctb) <- "lev.1"
  mtb <- matrix(NA, nrow = nrow(dat), ncol = 0) ## Medoid table
  rownames(mtb) <- rownames(dat)

  ## Create ctb.flat for the one level case:
  ## cluster affiliations without info about hierarchy(running number)
  ctb.flat <- as.numeric(as.factor(ctb))
  names(ctb.flat) <- rownames(ctb)

  count <- 2 ## Counter for cluster levels

  ## -------------- Follow-up runs ------------------------------------------ ##

  while(stepdown) {
    if (l.max & count > l.max) {
      stepdown <- FALSE
    }
    if (stepdown) {

      ifelse(sum(spl) == 1, cas <- "group ", cas <- "groups ")
      if (cas != 0) {
        message("\nLevel ", count, ": Partitioning ", sum(spl), " ", cas, appendLF = FALSE)
      }
      output.sub <- list()              ## Empty list for results
      count.2 <- 1
      for (j in 1:length(spl)) {          ## Loop through partitions

        ## Report progress
        message(".", appendLF = FALSE)

        if (spl[j] == 1) { ## splittable?
          x.sub <- dat[matr[, ncol(matr)] == j, ] ## Create data subset

          # Polishing for subsequent runs
          if (polishing == 'strict') {
            repeat {
              IO.x.sub <- x.sub
              IO.x.sub[IO.x.sub > 0] <- 1
              initial_rows <- nrow(x.sub)
              initial_cols <- ncol(x.sub)
              x.sub <- x.sub[, colSums(IO.x.sub) > 1] # Species with > 1 occurrence
              x.sub <- x.sub[rowSums(IO.x.sub) > 1, ] # ... and plots with > 1 species
              x.sub <- x.sub[, apply(x.sub, 2, var) > 0]
              if (nrow(x.sub) == initial_rows && ncol(x.sub) == initial_cols) {
                break
              }
            }
            if (is.null(dim(x.sub)) | nrow(x.sub) < 3 | ncol(x.sub) < 2) {
              output.sub[[j]] <- NA  # Nothing to split
            } else {                 
                output.sub[[j]] <- core(data.matrix(x.sub, rownames.force = NA))
                count.2 <- count.2 + 1
            }
          } else if (polishing == 'relaxed') {
            x.sub <- x.sub[, colSums(x.sub) > 0]
            x.sub <- x.sub[, apply(x.sub, 2, var) > 0] # + species with variance
            x.sub <- x.sub[rowSums(x.sub) > 0, ]       # + plots with species

            if (is.null(dim(x.sub)) | nrow(x.sub) < 3 | ncol(x.sub) < 2) {
              output.sub[[j]] <- NA  # Nothing to split
            } else {
                output.sub[[j]] <- core(data.matrix(x.sub, rownames.force = NA))
                count.2 <- count.2 + 1
            }
          }
        } else {
          output.sub[[j]] <- NA
        }
      }

      ## Look if some of the partitions should be rejected
      ok.vec <- vector()
      for (l in 1:length(spl)) {
        if (is.na(output.sub[[l]][1])) {
          ok.vec[[l]] <- FALSE
        } else {
          ok.vec[[l]] <- output.sub[[l]]$is.ok
        }
      }

      ## Report
      ok.n <- sum(as.numeric(ok.vec))

      ## Stop if all partitions are rejected
      if (ok.n == 0) stepdown <- FALSE

      if (stepdown) {
        ## Write new clusters to container
        matr <- cbind(matr, matrix(NA, nrow = nrow(matr), ncol = 1))
        colnames(matr)[ncol(matr)] <- paste("lev.", ncol(matr), sep = "")

        for (m in 1:length(spl)) { ## Fill in values

          if (ok.vec[m]) {
            sc <- output.sub[[m]]$clusters
            idx.sc <- rownames(matr) %in% names(sc)
            matr[idx.sc, ncol(matr)] <- sc
          }
        }

        ## Make matrix with names expressing hierarchy
        ctb <- cbind(ctb, matrix(0, nrow = nrow(dat), ncol = 1))
        colnames(ctb)[ncol(ctb)] <- paste("lev.", ncol(ctb), sep = "")
        ctb.red <- na.omit(ctb)
        blub <- matr[,ncol(matr)]
        blub[is.na(blub)] <- 0
        ctb.red[,ncol(ctb.red)] <- paste(ctb.red[,ncol(ctb.red)-1],
                                            blub, sep = ".")
        ctb[rownames(ctb) %in% rownames(ctb.red), ncol(ctb)] <-
          ctb.red[,ncol(ctb.red)]
        ctb <- as.data.frame(ctb)

        ## Overwrite ctb.flat in the multiple level case
        ## Cluster affiliations without info about hierarchy(running number)
        ctb.flat <- list()
        for (o in 1:ncol(ctb)) {
          ctb.flat[[o]] <- as.numeric(as.factor(ctb[,o]))
          names(ctb.flat[[o]]) <- rownames(ctb)
          names(ctb.flat)[o] <- paste("lev.", o, sep="")
        }

        ## Make working matrix
        matr.red <- na.omit(matr)
        matr.red [,ncol(matr.red)] <- as.numeric(as.factor(paste
            (matr.red[,ncol(matr)-1], matr.red[,ncol(matr.red)], sep = ".")))
        matr[rownames(matr) %in% rownames(matr.red),] <- matr.red
        matr[is.na(matr)] <- 0

        ## Medoids
        mtb <- cbind(mtb, matrix(0, nrow = nrow(dat), ncol = 1))
        rownames(mtb) <- rownames(dat)
        for (n in 1:length(spl)) {
          if (ok.vec[n]) {
            mtb[output.sub[[n]]$medoids, ncol(mtb)] <- 1
          }
        }
        med[[count]] <- names(mtb[mtb[,ncol(mtb)] == 1, ncol(mtb)])
        names(med[[count]]) <- ctb[mtb[,ncol(mtb)] == 1, ncol(ctb)]
        nam <- sort(names(med[[count]]))
        med[[count]] <- med[[count]][nam]

        ## Add new partitions to the summary matrix
        for (x in 1:length(spl)) {
          if (ok.vec[x]) {
            summ <- cbind(summ, matrix(NA, nrow = 9, ncol = 1))
            colnames(summ)[ncol(summ)] <- paste("Part.", ncol(summ), sep = "")
            summ[2, ncol(summ)] <- length(output.sub[[x]]$medoids) ## Subgroups
            summ[3, ncol(summ)] <- output.sub[[x]]$d             ## Isomap dims
            summ[4, ncol(summ)] <- output.sub[[x]]$k.min         ## Minimum k
            summ[5, ncol(summ)] <- output.sub[[x]]$k             ## Selected k
            summ[6, ncol(summ)] <- output.sub[[x]]$k.max         ## Maximum k
            summ[7, ncol(summ)] <- output.sub[[x]]$noi           ## Indicators
            summ[8, ncol(summ)] <- output.sub[[x]]$ivi           ## IV(Indic.)
            summ[9, ncol(summ)] <- output.sub[[x]]$ivx           ## IV(all)
          }
        }

        ## Add new indicators to indicator container
        for (m2 in 1:length(spl)) {
          if (ok.vec[m2]) {
            idx <- length(indlist) + 1
            indlist[[idx]] <- output.sub[[m2]]$indnames
          }
        }
        names(indlist) <- colnames(summ)

        ## OK, now we have a matrix with cluster affiliations of this level
        ## Repeat the search for splittable units:

        subtab <- table(matr[,ncol(matr)])
        subtab <- subtab[names(subtab) != "0"]
        spl <- as.numeric(subtab > 2)
        names(spl) <- names(subtab)

        count <- count + 1

        ## Create object of class "hclust"
        dendro <- create_dendro(ctb.flat)
      }
    }
  } ## End while-loop

  ## Cluster names
  if (ncol(ctb) > 1) {
    cnam <- vector()
    for (y in 1:length(med)) {
      cnam <- c(cnam, names(med[[y]]))
    }

    childs <- vector()
    for (z in 1:length(cnam)) {
      nch <- nchar(cnam[z])
      childs <- c(childs, length(grep(cnam[z], substr(cnam, 1, nch),
                                      value = TRUE)))
    }
    summ[1,2:ncol(summ)] <- cnam[childs > 1]
  }
  summ <- as.data.frame(summ)

  ## Build supervised description string
  supervised <- NULL
  if (!is.null(centers)) {
    if (semi_supervised_search) {
      supervised <- "cluster centroids partly provided, flexible number of additional clusters"
    } else if (semi_supervised) {
      supervised <- "cluster centroids partly provided, fixed number of additional clusters"
    } else {
      supervised <- "cluster centroids all provided"
    }
  }
  if (!is.null(ind)) {
    if (is.null(supervised)) {
      supervised <- "indicator species provided"
    } else {
      supervised <- paste0(supervised, ", indicator species provided")
    }
  }

  ## ------------ Output -------------------------------------------------- ##

  if (ncol(ctb) == 1) {
    OUT <- list(
      call        = sys.call(),
      distance    = distance,
      flat        = ctb.flat,
      hier        = NULL,
      medoids     = med,
      analytics   = summ,
      dendro      = NULL,
      centers_usr = centers,
      ind_usr     = ind,
      supervised  = supervised,
      indicators  = indlist,
      dat         = dat
    )
  }
  else if (ncol(ctb) > 1) {
    OUT <- list(
      call        = sys.call(),
      distance    = distance,
      flat        = ctb.flat,
      hier        = ctb,
      medoids     = med,
      analytics   = summ,
      dendro      = dendro,
      centers_usr = centers,
      ind_usr     = ind,
      supervised  = supervised,
      indicators  = indlist,
      dat         = dat
    )
  }

  class(OUT) <- "isopam"

  ## save output
  if (juice) {
    write.table(ctb, file = "isopam/juicein.txt",
                 col.names = FALSE, quote = FALSE)
    if (!is.null(OUT$dendro[1])) {
      wth <-(nrow(ctb) * 11) + 100
      bmp(filename = "isopam/juicetree.bmp", width = wth)
      plot(OUT$dendro)
      dev.off()
    }
  }

  if (ncol(ctb) == 1) {
    message("\nNon-hierarchical partition created")
  } else {
    message("\nCluster tree with ", ncol(ctb), " levels created")
  }

  # Report missing rows
  unused_sites <- length(setdiff(rownames(dat_orig), rownames(dat)))
  if (unused_sites > 0) {
    site_label <- ifelse (unused_sites == 1, "site", "sites")
    message(paste(unused_sites, "unallocated", site_label))
  }

  invisible(OUT)
}


## --------------------- S3 methods --------------------------------------- ##

plot.isopam <- function(x, ...) {

  if (!is.null(x$dendro[1])) {
    tree <- as.dendrogram(x$dendro)
    plot(tree, main = format(x$call), ...)
    return(invisible(tree))
  } else {
    message("No cluster hierarchy - nothing to plot")
  }
}

identify.isopam <- function(x, ...) {
  identify(x$dendro, MAXCLUSTER = nrow(x$dat), ...)
}

summary.isopam <- function(object, ...) {

  # Retrieving info and building output
  nobj <- nrow(object$dat)
  ana <- object$analytics[-c(4, 6), , drop = FALSE]
  supervised <- object$supervised
  r <- list(
    call      = object$call,
    clusters  = object$flat,
    hierarchy = object$hier,
    medoids   = object$medoids,
    analytics = object$analytics)
  class(r) <- "summary.isopam"

  # Reporting
  cat("Call:", format(object$call), "\n")
  cat("Distance measure:", object$distance, "\n")
  if (!is.null(supervised)) {
    cat("Supervised:", supervised, "\n")
  }

  if (is.null(object$hier)) {
    partition <- TRUE
    ana <- object$analytics[-1, , drop = FALSE]   
    cat(nobj, "items grouped in a non-hierarchical partition\n")
  } else {
    partition <- FALSE
    lvl <- length(object$hier)
    cat(nobj, "items arranged in a cluster tree with", lvl, "levels\n")
  }

  print(ana)
  invisible(r)
}

print.isopam <- function(x, ...) {

  # Retrieving info and building output
  nobj <- nrow(x$dat)
  supervised <- x$supervised

  # Reporting
  cat("Call:", format(x$call), "\n")
  cat("Distance measure:", x$distance, "\n")
  if (!is.null(supervised)) {
    cat("Supervised:", supervised, "\n")
  }

  if (is.null(x$hier)) {

    cat(nobj, "items grouped in a non-hierarchical partition\n")
    cat("\nClustering vector\n")
    print(unlist(x[["flat"]]))

    # Medoids
    cat("\nMedoids\n")
    medoids <- unlist(x[["medoids"]])
    medoids <- as.data.frame(medoids)
    names(medoids) <- "Medoid"
    print(medoids)

  } else {

    lvl <- length(x$hier)

    # Cluster assignments
    cat(nobj, "items arranged in a cluster tree with", lvl, "levels\n")
    hier <- x[["hier"]]
    hier <- tibble::rownames_to_column(hier, "Site")
    print(tibble::as_tibble(hier, ...))

    # Medoids
    cat("\nMedoids\n")
    medoids <- unlist(x[["medoids"]])
    medoids <- as.data.frame(medoids)
    names(medoids) <- "Medoid"
    print(medoids)
  }

  invisible(x)
}
