isotab <- function(x, level = NULL, clusters = NULL, phi.min = "isotab",
                   p.max = .05) {

  ## ---------------- Functions --------------------- #

  # phi_equalized inspired by indicspecies (Miquel De Cáceres
  # and Florian Jansen)
  phi_equalized <- function(x, cluster, SP, N) {

    cluster_names  <- levels(as.factor(cluster))
    cluster_n <- length(cluster_names)
    cluster_indices <- apply(sapply(cluster_names, "==", cluster), 1, which)
    cluster_members <- diag(1, cluster_n, cluster_n)[cluster_indices, ]
    ni       <- colSums(cluster_members * cluster_members)
    tx       <- t(x)
    aisp     <- tx %*% cluster_members
    lisp     <- tx^2 %*% cluster_members
    aispni   <- sweep(aisp, 2, ni, "/")
    lispni   <- sweep(lisp, 2, ni, "/")
    #Corrected sum and length of species vectors
    aspK     <- (N / cluster_n) * rowSums(aispni)
    lspK     <- (N / cluster_n) * rowSums(lispni)
    #Corrected sum of species values in combinations
    aspC     <- matrix(0, nrow = SP, ncol = cluster_n)
    nC       <- vector(mode = "numeric", length = cluster_n)
    cnt      <- 1
    for(level in 1:cluster_n) {
      co <- utils::combn(1:cluster_n, level)
      for(j in 1:min(ncol(co), cluster_n - cnt + 1)) {
        if(nrow(co) > 1) {
          aspC[, cnt] <- rowSums(aispni[, co[, j]])
        } else {
          aspC[, cnt] <- aispni[, co[, j]]
        }
        nC[cnt] <- length(co[, j])
        cnt <- cnt + 1
        if(cnt > cluster_n) break
      }
      if(cnt > cluster_n) break
    }
    aspC <- (N / cluster_n) * aspC
    nC <- (N / cluster_n) * nC
    #Compute phi
    num <- N * aspC - aspK %o% nC
    den <- sqrt(((N * lspK) - aspK^2) %o% (N * nC - (nC^2)))
    str <- num / den
    colnames(str) <- cluster_names
    return(str)
  }

  ## Two-tailed Fisher's exact test for 2x2 tables
  fshtest <- function(x) {
    PVAL <- NULL
    m <- sum(x[, 1])
    n <- sum(x[, 2])
    k <- sum(x[1, ])
    x <- x[1, 1]
    lo <- max(0, k - n)
    hi <- min(k, m)
    support <- lo:hi
    logdc <- dhyper(support, m, n, k, log = TRUE)

    dnhyper <- function(ncp) {
      d <- logdc + log(ncp) * support
      d <- exp(d - max(d))
      d / sum(d)
    }

    pnhyper <- function(q, ncp = 1, upper.tail = FALSE) {
      if (ncp == 1) {
        if (upper.tail)
          return(phyper(x - 1, m, n, k, lower.tail = FALSE))
        else return(phyper(x, m, n, k))
      }
      else if (ncp == 0) {
        if (upper.tail)
          return(as.numeric(q <= lo))
        else return(as.numeric(q >= lo))
      }
      else if (ncp == Inf) {
        if (upper.tail)
          return(as.numeric(q <= hi))
        else return(as.numeric(q >= hi))
      }
      d <- dnhyper(ncp)
      if (upper.tail) {
          sum(d[support >= q])
      } else {
        sum(d[support <= q])
      }
    }
    PVAL <- switch("two.sided", less = pnhyper(x, 1),
                   greater = pnhyper(x, 1, upper.tail = TRUE),
                   two.sided = {
                     relErr <- 1 + 10^(-7)
                     d <- dnhyper(1)
                     sum(d[d <= d[x - lo + 1] * relErr])
                   })
    return(PVAL)
  }

  ## ---------- Function isotab starts here -------------- #

  # If isopam object is provided
  if ("isopam" %in% class(x)) {
    snam <- colnames(x$dat)
    bigtab <- x$dat
    IO <- ifelse(x$dat > 0, 1, 0)    # presence-absence x
    N <- nrow(IO)                    # site number
    SP <- ncol(IO)                   # species number
    n_occurrences <- colSums(IO)     # total species frequency
    if (is.null(x$hier)) {           # depth of hierarchy
      depth <- 1
    } else {
      depth <- ncol(x$hier)
    }

    # If level is NULL:
    # Set level to second level if possible
    if (is.null(level)) {
      if (depth == 1) {
        level <- 1
      } else {
        level <- 2
      }
    }

    # If level is set:
    if (is.null(x$hier)) {
      if (level > 1) message("No hierarchy levels available")
      clusters <- x$flat # extract clusters from data object
      tab <- t(aggregate(IO, by = list(clusters), FUN = sum))
      colnames(tab) <- tab[1, ]
    } else {
      if (level > depth) {
        level <- depth
        message("Switching to lowest level ", depth)
      }
      clusters <- x$hier[, level] # extract clusters from data object
      names(clusters) <- rownames(x$hier)
      tab <- t(aggregate(IO, by = list(clusters), FUN = sum))
      colnames(tab) <- tab[1, ]
    }

  } else {

  # In case non-isopam object is provided

    # Is this a tibble?
    if ("tbl_df" %in% class(x)) {
      # If the first column is <chr> and all other columns are
      # <dbl> or <int>, transform into dataframe
      if (is.character(x[[1]]) && all(sapply(x[, -1], is.numeric))) {
          rn <- x[[1]]
          x <- as.data.frame(x[, -1])
          rownames(x) <- rn
          message("First column of tibble used for naming observations")
      } else {
        stop("Name column (chr) or matrix or dataframe with rownames required")
      }
    }

    # Are there clusters?
    if (is.null(clusters)) {
      stop("No clusters provided")
    }
    # Are there plot names in clusters?
    if (is.null(names(clusters))) {
      stop("No site names in clusters")
    } else {
      plotnames <- names(clusters)
    }
    # Make sure that clusters are of class char
    clusters <- as.character(clusters)
    names(clusters) <- plotnames

    # Are there site names (row names) in data?
    if (is.null(rownames(x))) {
      stop("No rownames in data")
    }
    # Are there species names (column names) in data?
    if (is.null(colnames(x))) {
      stop("No colnames in data")
    }
    clusters <- clusters[rownames(x)]
    # Do sitenames match?
    if (!identical(rownames(x), names(clusters))) {
      stop("Names do not match")
    }
    # Are there species in all sites?
    if (!all(rowSums(x) > 0)) {
      stop("Remove empty rows")
    }  
    # Do all species have occurrences?
    if (!all(colSums(x) > 0)) {
      stop("Remove empty columns")
    }  


    snam <- colnames(x)
    bigtab <- x
    IO <- ifelse(x > 0, 1, 0)        # presence-absence data
    N <- nrow(IO)                    # site number
    SP <- ncol(IO)                   # species number
    n_occurrences <- colSums(IO)     # total species frequency
    names(n_occurrences) <- snam
    tab <- t(aggregate(IO, by = list(clusters), FUN = sum))
    colnames(tab) <- tab[1, ]
    depth <- 1
    level <- 1
  }

  ## 1) Contingency table
  cnam <- tab[1, ] # Cluster names are in the first row
  tab <- tab[-1, ] # Remove that row
  tab <- matrix(as.numeric(tab), nrow = nrow(tab))
  colnames(tab) <- cnam
  rownames(tab) <- snam
  nc <- ncol(tab)

  ## 2) Frequency table
  siz <- table(clusters) ## Cluster sizes
  spc <- t(as.matrix(siz))[rep(1, nrow(tab)),]
  frequency <- tab / spc ## Frequency
  frequency <- round(frequency * 100, 0) ## as percentage

  ## 3) Fisher's significance table
  ft <- tab 
  for (fsp in 1:SP) {       # Species loop
    spec_frq <- n_occurrences[fsp]
    spec_io <- IO[, fsp]

    Nj <- siz
    insd <- tab[fsp, ]
    absci <- Nj - insd
    outs <- spec_frq - insd
    absco <- N - Nj - outs
    for (fcl in 1:nc) {     # Cluster loop
      # Create the 2x2 matrix
      fshm <- matrix(c(insd[fcl], absci[fcl], outs[fcl], absco[fcl]), nrow = 2,
                     byrow = TRUE)
      # Perform the test and store the result
      ft[fsp, fcl] <- fshtest(fshm)
    }
  }
  
  ## 4) Significance symbols
  ft.symb <- ft
  ft.symb[ft > 0.05] <- ""
  ft.symb[ft <= 0.05] <- "*"
  ft.symb[ft <= 0.01] <- "**"
  ft.symb[ft <= 0.001] <- "***"

  ## 5) Combined frequency table with significance symbols
  frq_ft <- matrix(paste(frequency, ft.symb, sep = ""),
                   nrow(frequency), ncol(frequency))
  frq_ft <- data.frame(frq_ft)
  colnames(frq_ft) <- colnames(ft.symb)
  rownames(frq_ft) <- rownames(ft.symb)

  ## 6) Equalized phi (Tichy & Chitry 2006)
  phi <- phi_equalized(IO, clusters, SP, N)
  phi[is.na(phi)] <- 0  ## Replace NaN-values by 0

  ## 7) Calculate phi.min threshold using a sigmoidal function
  if (phi.min == "isotab") {
    phi.min <- round(0.7 / (1 + exp(0.008 * (N - 5))) + 0.14, 2)
  }

  ## Table sorting
  ## 8) .... by phi and frequency
  # Note: The group affiliation is defined by phi. Equal values are frequent
  # but this leads to lack of significance and species don't pass default 
  # p.max. Accordingly, these cases do not receive special treatment.
  spc_own_clust <- apply(phi, 1, which.max)
  # Vector with the maximum freqencies
  frequency_max <- apply(frequency, 1, max)
  # Retrieve the phi value from 'own' cluster
  phi_at_max <- phi[cbind(1:nrow(phi), spc_own_clust)]
  names(phi_at_max) <- rownames(phi)
  # Retrieve the significance from 'own' cluster
  ft_at_max <- ft[cbind(1:nrow(ft), spc_own_clust)]
  names(ft_at_max) <- rownames(ft)

  # Order complete frequency table first by cluster assignment and than, 
  # descending, by max. frequency and descending phi.
  order_idx <- order(spc_own_clust, -frequency_max, -phi_at_max)
  # Total species frequencies, sorted accordingly
  frq_ordered <- n_occurrences[order_idx]
  names_ordered <- names(frq_ordered)

  # Transfer this order:
  # ... to the phi table
  phi_ordered <- phi[names_ordered, ]
  # ... and to phi_at_max
  phi_at_max_ordered <- phi_at_max[names_ordered]
  # ... and to the vector with cluster assignments
  spc_own_clust_ordered <- spc_own_clust[names_ordered]
  # ... and to the significance table
  ft_ordered <- ft[names_ordered, ]
  # ... and to ft_at_max
  ft_at_max_ordered <- ft_at_max[names_ordered]
  # ... and to the combined frequency and significance
  frq_ft_ordered <- frq_ft[names_ordered, ]

  ## 9) Set of diagnostic species to be shifted up
  ## Only look at "own" cluster
  filter1 <- ft_at_max_ordered <= p.max
  filter2 <- phi_at_max_ordered >= phi.min
  spc_filtered <- filter1 & filter2
  spc_filtered_n <- sum(spc_filtered) ## how many diagnostic species
  if (spc_filtered_n == 0) diag <- "No diagnostic species with given thresholds"
  if (spc_filtered_n > 0) diag <- frq_ft_ordered[spc_filtered, ]

  ## 10) For later use in the bottom part of the tables
  names_ordered_by_frq <- names(n_occurrences[order(-n_occurrences)])
  frq_ft_ordered_by_frq <- frq_ft[names_ordered_by_frq, ]

  ## 11) Move diagnostic species to top
  if (spc_filtered_n > 0) {
    FRQ <- rbind(diag, frq_ft_ordered_by_frq[!rownames(frq_ft_ordered_by_frq) %in%
      rownames(diag), ])
  } else {FRQ <- frq_ft_ordered_by_frq}

  phi_final_order <- phi_ordered[rownames(FRQ), ]
  ft_final_order <- ft_ordered[rownames(FRQ), ]

  ## 12) Report cluster sizes
  siz <- t(as.matrix(siz))
  rownames(siz) <- "n"

  ## 13) Parameters
  param <- c(phi.min, p.max)
  names(param) <- c("phi.min", "p.max")

  ## 14) Info about diagnostic species
  dig2 <- spc_own_clust_ordered[names(spc_own_clust_ordered) %in% rownames(diag)]
  typ <- list()
  ispec <- vector()
  for (i in 1:nc) {
    if (length(names(dig2)[dig2 == i]) > 0) {
      typ[i] <- paste(names(dig2)[dig2 == i], collapse = ", ")
      ispec <- c(ispec, names(dig2)[dig2 == i])
    } else { typ[i] <- "Nothing particularly typical" }
  }
  names(typ) <- cnam

  ## 15) Ordered, comprehensive table
  # Cluster assignments
  clusters <- sort(clusters)

  # Sort original table
  bigtab <- bigtab[, rownames(FRQ)]   # Adapt species order
  bigtab <- bigtab[names(clusters), ] # 

  ## 17) Output
  isotab_out <- list(
    call           = match.call(),
    depth          = depth,
    level          = level,
    tab            = FRQ,
    phi            = phi_final_order,
    fisher_p       = ft_final_order,
    n              = siz,
    thresholds     = param,
    typical        = typ,
    typical_vector = ispec,
    sorted_table   = bigtab)

  class(isotab_out) <- "isotab"

  invisible(isotab_out)
}

print.isotab <- function(x, n = NA, ...) {

  tab_for_tibble <- tibble::rownames_to_column(x$tab, var = "Species")
  tabtibble <- tibble::as_tibble(tab_for_tibble)

  cat("\nThresholds used for arranging the table\n")
  cat("(adjust as needed while calling isotab)\n")
  print(x$thresholds)

  if (x$depth == 1) {
    cat("\nFrequencies,", ncol(x$tab), "clusters, dataframe in $tab\n")
  } else {
    cat("\nFrequencies at cluster level", x$level, "\n")
    cat(ncol(x$tab), "clusters, dataframe in $tab\n")
  }

  if (is.na(n)) {
    # n is set to the number of diagnostic species + 10
    print(tabtibble, n = length(x$typical_vector) + 10, ...)
  } else {
    print(tabtibble, n = n, ...)
  }

  cat("\nCluster sizes\n")
  print(x$n)

  cat("\nTypical with given settings\n")
  print(x$typical)

  invisible(x)
}
