isopam <-  function(dat, c.fix = FALSE, c.max = 6, 
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

  ## In case of predefined centers check their validity,
  if (!is.null(centers)) {

    if (is.character(centers)) {
      if (!all(centers %in% rownames(dat))) {
        stop("Centers not valid")
      } else {
        # Transform into indices
        centers <- which(rownames(dat) %in% centers)
      }
    }
    if (is.numeric(centers) && max(centers) > nrow(dat)) {
      stop("Centers not in range")
    }

    c.fix <- length(centers)
  }

  ## Initiate count
  count <- 1

  ## Default minimum cluster number
  c.min <- 2

  if (is.numeric(c.fix)) {
    l.max <- 1
    if (c.fix < 2) stop("c.fix < 2")
    if (c.fix > nrow(dat) - 1) c.fix <- nrow(dat) - 1
    c.min <- c.fix
    c.max <- c.fix
  }

  ## ----------- core function ---------------------------------------------- ##

  core <- function(xdat) {

    IO.xdat <- ifelse(xdat > 0, TRUE, FALSE)

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

    slowmode <- FALSE # Stays FALSE if fast loops are possible

    ## -------------- Dealing with memory limitations ------------------ ##

    # Available RAM
    sys_info <- ps::ps_system_memory()
    available_mem <- sys_info$avail

    ## Try to prepare output array
    ## Need a place to store the calculations done in the g2 loop
    ## so they may be reused rather than recalculated.
    ## The +1 is due to zeros needing to be indexed as well.

    object_size <- 0
    slowmode <- tryCatch({
      DDD_lookup_table <- array(NA_real_, c(N.xdat, N.xdat + 1,
                                N.xdat + 1))
      object_size <- object.size(DDD_lookup_table)
      FALSE
    }, error = function(e) {
      TRUE
    })

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

    # Condition: In case of fastmode (with lookup array):
    # lookup fits in mem limit
    if (!slowmode & object_size > future_memlim) {
      fut <- FALSE
    }

    # Reporting
    if (!already_shown) {
      # Slowdown alert not only in slowmode but also if
      # lookup exceeds memory (causing memory swapping)
      if (slowmode | object_size > available_mem) {
        message("Hitting memory limit: Expect substantial slowdown")
      }
      # Conditional suggestion to increase future mem limit
      if (!slowmode & object_size > future_memlim & 
          object_size <= available_mem) {
        if (object_size >= 10 ^ 9) {
          obj_siz_unit <- "GB"
          obj_siz <- round(object_size / 2 ^ 30, 1)
        } else {
          obj_siz_unit <- "MB"
          obj_siz <- round(object_size / 2 ^ 20, 1)
        }
        if (available_mem >= 10 ^ 9) {
          mem_siz_unit <- "GB"
          mem_siz <- round(available_mem / 2 ^ 30, 1)
        } else {
          mem_siz_unit <- "MB"
          mem_siz <- round(available_mem / 2 ^ 20, 1)
        }
        message("At least ", object_size, " bytes (",
        obj_siz, " ", obj_siz_unit, ") required for the use of a",
        "  lookup array in parallel operation.")
        message("Consider increasing 'future.globals.maxSize' if RAM allows.")
        message("Available memory: ", available_mem, " bytes (", mem_siz,
                " ", mem_siz_unit, ")")
      }
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

      future::plan(future::multisession)

      ## b-loop: Isomap k
      out.array <- array(future.apply::future_sapply(k.min:k.max, 
                         future.seed = TRUE, function(b) {

        ## Isomap
        suppressMessages(isom <- isomap(dst.xdat, ndim = d.max, k = b))

        ## Fixing the maximum of dimensions considered when calculating
        ## the distance matrix for the isomap space
        d.max.new <- min(sum(isom$eig > 0), ncol(isom$points), d.max, na.rm = T)
        out.mat <- matrix(NA, nrow = d.max - 1, ncol = c.max - c.min + 1 )

        ## d-Loops
        if (d.max.new > 1) {
          for (d in 2:d.max.new) {

            isodiss <- suppressWarnings(daisy(isom$points[, 1:d], metric =
                                              "euclidean", stand = TRUE))
            # For fastpam
            isodiss_vector <- as.vector(isodiss)

            for (e in c.min:c.max) { ## e-loop: Cluster no.

              ## --------- Partitioning (PAM) ------------------------------- #

              if (!is.null(centers)) {
                cl.iso <- cluster::pam(isodiss, k = e, medoids = centers,
                              diss = TRUE, do.swap = FALSE)
                cl <- cl.iso$clustering                     ## Group affiliation
                ci <- cl.iso$clusinfo[, 1]                  ## Cluster size
              } else {
                cl.iso <- fastkmedoids::fastpam(isodiss_vector, k = e, n = N.xdat)
                cl <- cl.iso@assignment
                ci <- as.numeric(table(cl))
              }

              ## ----------- Start parallel fast mode --------------- ##

              if (!slowmode) {
                
                ## Catching a situation where memory limit is hit in fast mode
                tryCatch({
                  ## For Williams` correction
                  w1 <- N.xdat * sum(1 / ci) - 1
                  w2 <- 6 * N.xdat * (e - 1)

                  ## Compute G-values for species(code adapted from Lubomir Tichy)

                  gt <- matrix(NA, SP.xdat, 1) ## Matrix for G-test results

                  for (g1 in 1:SP.xdat) {        ## g1-loop through species

                    DDD <- 0
                    spec_frq <- frq.xdat[g1]

                    ## Precompute bom and bim for the current species
                    bom <- spec_frq / N.xdat
                    bim <- 1 - bom
                    
                    ## For Williams` correction
                    willi <- 1 + ((w1 * w3[g1]) / w2)

                    ## Mask out entries in cl using appropriate col in IO.xdat
                    groupids <- cl[IO.xdat[, g1]]

                    for (g.fast in 1:e) {      ## g.fast-loop through clusters

                      fra1 <- sum(groupids == g.fast)   ## Species occ. in cluster
                      Nj <- ci[g.fast]                  ## Cluster size

                      ## Have we calculated this before?
                      DDDadd <- DDD_lookup_table[spec_frq, Nj + 1, fra1 + 1]

                      if (!is.na(DDDadd)) {
                        ## Already existed in the lookup table; use it.
                        DDD <- DDD + DDDadd
                      } else {
                        ## so need to calculate it ...
                        fra0 <- Nj - fra1
                        bum <- fra1 / (Nj * bom)
                        bam <- fra0 / (Nj * bim)
                        DDDadd <- 0
                        if (!is.na(bum) && bum > 0) {
                          DDDadd <- DDDadd + (fra1 * log(bum))
                        }
                        if (!is.na(bam) && bam > 0) {
                          DDDadd <- DDDadd + (fra0 * log(bam))
                        }
                        DDD <- DDD + DDDadd
                        ## ... and store for next time
                        DDD_lookup_table[spec_frq, Nj + 1, fra1 + 1] <- DDDadd
                      }
                    }
                    DDD <- DDD * 2
                    gt[g1, ] <- DDD / willi ## Williams` correction
                  }

                  ## Standardization(Botta-Dukat et al. 2005)
                  gt.ex <- e - 1                 ## Expected G
                  gt.sd <- sqrt(2 * gt.ex)      ## Expected sd
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
                }, error = function(e) {
                  slowmode <- TRUE
                })
              }

              ## ----------- Start parallel slow mode -------------- ##

              if (slowmode) {

                ## For Williams` correction
                w1 <- N.xdat * sum(1 / ci) - 1
                w2 <- 6 * N.xdat * (e - 1)

                ## Compute G-values for species
                ## (code adapted from Lubomir Tichy)

                gt <- matrix(NA, SP.xdat, 1) ## Matrix for G-test results

                for (g1 in 1:SP.xdat) {         ## g1-loop through species

                  DDD <- 0
                  spec_frq <- frq.xdat[g1]

                  ## Precompute bom and bim for the current species
                  bom <- spec_frq / N.xdat
                  bim <- 1 - bom

                  ## For Williams` correction
                  willi <- 1 + ((w1 * w3[g1]) / w2)

                  ## Mask out entries in cl using appropriate col in IO.xdat
                  groupids <- cl[IO.xdat[, g1]]

                  for (g.slow in 1:e) {   ## g.slow-loop through clusters

                    fra1 <- sum(groupids == g.slow)   ## Species in cluster
                    Nj <- ci[g.slow]                  ## Cluster size
                    fra0 <- Nj - fra1
                    bum <- fra1 / (Nj * bom)
                    bam <- fra0 / (Nj * bim)
                    DDDadd <- 0
                    if (!is.na(bum) && bum > 0) {
                      DDDadd <- DDDadd + (fra1 * log(bum))
                    }
                    if (!is.na(bam) && bam > 0) {
                      DDDadd <- DDDadd + (fra0 * log(bam))
                    }
                    DDD <- DDD + DDDadd
                  }
                  DDD <- DDD * 2
                  gt[g1,] <- DDD / willi ## Williams` correction
                }

                ## Standardization(Botta-Dukat et al. 2005)
                gt.ex <- e - 1                 ## Expected G
                gt.sd <- sqrt(2 * gt.ex)      ## Expected sd
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

                ## Standard: Filtering and aver aging
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
              }

            ## ----------- End parallel slow mode --------------- ##

            }
          }
        }

        return(out.mat)

      }), dim = c(d.max - 1, c.max - c.min + 1, k.max - k.min + 1))

    ## ----------- End parallel mode --------------------- ##

    } else {

    ## ----------- Non-parallel processing --------------- ##

      # this line is used to reset to default state
      future::plan(future::sequential)

      ## b-loop: Isomap k
      out.array <- array(sapply(k.min:k.max, function(b) {

        ## Isomap
        suppressMessages(isom <- isomap(dst.xdat, ndim = d.max, k = b))

        ## Fixing the maximum of dimensions considered when calculating
        ## the distance matrix for the isomap space
        d.max.new <- min(sum(isom$eig > 0), ncol(isom$points), d.max, na.rm = T)
        out.mat <- matrix(NA, nrow = d.max - 1, ncol = c.max - c.min + 1)

        if (d.max.new > 1) {
          for (d in 2:d.max.new) {

            isodiss <- suppressWarnings(daisy(isom$points[, 1:d], metric =
                                        "euclidean", stand = TRUE))
            # For fastpam
            isodiss_vector <- as.vector(isodiss)

            for (e in c.min:c.max) { ## e-loop: Cluster no.

              ## --------- Partitioning(PAM) ------------------------------ #

              if (!is.null(centers)) {
                cl.iso <- cluster::pam(isodiss, k = e, medoids = centers,
                              diss = TRUE, do.swap = FALSE)
                cl <- cl.iso$clustering                     ## Group affiliation
                ci <- cl.iso$clusinfo[, 1]                  ## Cluster size
              } else {
                cl.iso <- fastkmedoids::fastpam(isodiss_vector, k = e, n = N.xdat)
                cl <- cl.iso@assignment
                ci <- as.numeric(table(cl))
              }
              
              ## ----------- Start non-parallel fastmode --------------- ##

              if (!slowmode) {

                ## Catching a situation where memory limit is hit in fast mode
                tryCatch({
                  ## For Williams` correction
                  w1 <- N.xdat * sum(1 / ci) - 1
                  w2 <- 6 * N.xdat * (e - 1)

                  ## Compute G-values for species(code adapted from Lubomir Tichy)

                  gt <- matrix(NA, SP.xdat, 1) ## Matrix for G-test results

                  for (g1 in 1:SP.xdat) {      ## g1-loop through species

                    DDD <- 0
                    spec_frq <- frq.xdat[g1]

                    ## Precompute bom and bim for the current species
                    bom <- spec_frq / N.xdat
                    bim <- 1 - bom

                    ## For Williams` correction
                    willi <- 1 + ((w1 * w3[g1]) / w2)

                    ## Mask out entries in cl using appropriate col in IO.xdat
                    groupids <- cl[IO.xdat[, g1]]

                    for (g.fast in 1:e) {   ## g.fast-loop through clusters

                      fra1 <- sum(groupids == g.fast)   ## Species in cluster
                      Nj <- ci[g.fast]                  ## Cluster size

                      ## Have we calculated this before?
                      DDDadd <- DDD_lookup_table[spec_frq, Nj + 1, fra1 + 1]

                      if (!is.na(DDDadd)) {
                        ## Already existed in the lookup table; use it.
                        DDD <- DDD + DDDadd

                      } else {
                        ## so need to calculate it ...
                        fra0 <- Nj - fra1
                        bum <- fra1 /(Nj * bom)
                        bam <- fra0 /(Nj * bim)
                        DDDadd <- 0
                        if (!is.na(bum) && bum > 0) {
                          DDDadd <- DDDadd + (fra1 * log(bum))
                        }
                        if (!is.na(bam) && bam > 0) {
                          DDDadd <- DDDadd + (fra0 * log(bam))
                        }
                        DDD <- DDD + DDDadd

                        ## ... and store for next time
                        ## This is the memory bottleneck
                        DDD_lookup_table[spec_frq, Nj + 1, fra1 + 1] <- DDDadd
                      }
                    }

                    DDD <- DDD * 2
                    gt[g1,] <- DDD / willi ## Williams" correction
                  }

                  ## Standardization(Botta-Dukat et al. 2005)
                  gt.ex <- e - 1                 ## Expected G
                  gt.sd <- sqrt(2 * gt.ex)      ## Expected sd
                  G <- (gt - gt.ex) / gt.sd

                  ## Using predefined indicators
                  if (usr_ind) {
                    glgth <- length(G[xind])
                    if (glgth == 0) out.mat[d - 1, e + 1 - c.min] <- NA
                    else out.mat[d - 1, e + 1 - c.min] <- mean(G[xind])
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
                }, error = function(e) {
                  slowmode <- TRUE
                })
              }

              ## ----------- Start non-parallel slowmode --------------- ##

              if (slowmode) {

                ## For Williams` correction
                w1 <- N.xdat * sum(1 / ci) - 1
                w2 <- 6 * N.xdat * (e - 1)

                ## Compute G-values for species(code adapted from Lubomir Tichy)

                gt <- matrix(NA, SP.xdat, 1) ## Matrix for G-test results

                for (g1 in 1:SP.xdat) {       ## g1-loop through species

                  DDD <- 0
                  spec_frq <- frq.xdat[g1]

                  ## Precompute bom and bim for the current species
                  bom <- spec_frq / N.xdat
                  bim <- 1 - bom

                  ## For Williams` correction
                  willi <- 1 + ((w1 * w3[g1]) / w2)

                  ## Mask out entries in cl using appropriate col in IO.xdat
                  groupids <- cl[IO.xdat[, g1]]

                  for (g.slow in 1:e) {    ## g.slow-loop through clusters

                    fra1 <- sum(groupids == g.slow)   ## Species in cluster
                    Nj <- ci[g.slow]                  ## Cluster size
                    fra0 <- Nj - fra1
                    bum <- fra1 / (Nj * bom)
                    bam <- fra0 / (Nj * bim)
                    DDDadd <- 0
                    if (!is.na(bum) && bum > 0) {
                      DDDadd <- DDDadd + (fra1 * log(bum))
                    }
                    if (!is.na(bam) && bam > 0) {
                      DDDadd <- DDDadd + (fra0 * log(bam))
                    }
                    DDD <- DDD + DDDadd
                  }
                  DDD <- DDD * 2
                  gt[g1,] <- DDD / willi ## Williams" correction
                }

                ## Standardization(Botta-Dukat et al. 2005)
                gt.ex <- e - 1                 ## Expected G
                gt.sd <- sqrt(2 * gt.ex)      ## Expected sd
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

                ## Averaging
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
              }

              ## ----------- End non-parallel slowmode --------------- ##

            }
          }
        }

        return(out.mat)
      }),dim = c(d.max - 1, c.max - c.min + 1, k.max - k.min + 1))

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
        cl.iso <- cluster::pam(d.iso, k = mc, medoids = centers,
                      diss = TRUE, do.swap = FALSE)
      } else {
        cl.iso <- cluster::pam(d.iso, k = mc, diss = TRUE)
      }

      CLS <- cl.iso$clustering                 ## Group affiliation
      MDS <- cl.iso$medoids                    ## Medoids
      CLI <- t(cl.iso$clusinfo[, 1])           ## Cluster size

      ## ------------ method == G ------------------------ #

      ## Contingency table
      tab <- t(aggregate(IO.xdat, by = list(CLS), FUN = sum)[, -1])
      ## Frequency table
      #inf <- matrix(rep(CLI, SP.xdat), nrow = SP.xdat, byrow = TRUE)
      #FRQ <- tab / inf

      ## Prepare Williams` correction
      w1 <- N.xdat * sum(1 / CLI) - 1
      w2 <- 6 * N.xdat * (mc - 1)

      ## Matrix for results
      g.1 <- matrix(NA, SP.xdat, 1)

      for (g3 in 1:SP.xdat) {

        willi <- 1 + ((w1 * w3[g3]) / w2)
        DDD <- 0
        bom <- frq.xdat[g3] / N.xdat
        bim <- 1 - bom
        rip <- tab[g3, ]

        for (g4 in 1:mc) {

          fra1 <- rip[g4]
          Nj <- CLI[g4]
          fra0 <- Nj - fra1
          bum <- fra1 / (Nj * bom)
          bam <- fra0 / (Nj * bim)

          ## adding up values
          if (!is.na(bum) && bum > 0) {
            DDD <- DDD + (fra1 * log(bum))
          }
          if (!is.na(bam) && bam > 0) {
            DDD <- DDD + (fra0 * log(bam))
          }  
        }
        DDD <- DDD * 2
        g.1[g3, ] <- DDD / willi
      }

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
  centers <- object$centers_usr
  ind <- object$ind_usr
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
  if (!is.null(centers) && is.null(ind)) {
    cat("Supervised mode with medians suggested by user\n")
  }
  else if (!is.null(ind) && is.null(centers)) {
    cat("Supervised mode with indicators suggested by user\n")
  }
  else if (!is.null(centers) & !is.null(ind)) {
    cat("Supervised mode with medians and indicators suggested by user\n")
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
  centers <- x$centers_usr
  ind <- x$ind_usr

  # Reporting
  cat("Call:", format(x$call), "\n")
  cat("Distance measure:", x$distance, "\n")
  if (!is.null(centers) && is.null(ind)) {
    cat("Supervised mode with medians suggested by user\n")
  }  
  else if (!is.null(ind) && is.null(centers)) {
    cat("Supervised mode with indicators suggested by user\n")
  }  
  else if (!is.null(centers) & !is.null(ind)) {
    cat("Supervised mode with medians and indicators suggested by user\n")
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
