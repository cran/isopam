isopam.2 <-
function (dat, c.num = FALSE, c.max = 10, filtered = TRUE, 
          distance = 'bray', g.min = 3.5, k.max = 100, 
          stopcrit = c(2,7), maxlev = FALSE,
          juice = FALSE)
{
  require (vegan) || stop ("needs package vegan")
  require (cluster) || stop ("needs package cluster")
  if (distance != "bray" & distance != "jaccard") 
    require (proxy) || stop ("needs package proxy")  

  ## Create directory for juice output
  if (juice == TRUE) dir.create ("isopam", showWarnings = FALSE)
  
  ## Add Sample names if necessary
  if (is.null (rownames (dat)) == TRUE) 
  {
    rownames (dat) <- c(1:nrow(dat))
  }
   
  ## Add taxon names if necessary
  if (is.null (colnames (dat)) == TRUE) 
  {
    colnames (dat) <- c(1:ncol(dat))
  }

  ## Convert to matrix if necessary
  nam <- rownames (dat)
  dat <- as.matrix (dat)
  rownames (dat) <- nam

  ## Remove empty columns
  dat <- dat [,colSums (dat) > 0] ## Remove empty columns
  dat <- dat [rowSums (dat) > 0,] ## ... and rows

  ## Initiate count
  count <- 1
  
  ## Set maxlev == 1 if c.num is a number
  if (is.numeric (c.num)) maxlev <- 1
  
  ## ----------- core function ---------------------------------------------- ##

  core <- function (xdat) 
  {
    IO.xdat <- xdat
    IO.xdat [IO.xdat > 0] <- 1

    ## Some useful descriptors
    N.xdat <- nrow (xdat)                         ## Total number of plots
    SP.xdat <- ncol (xdat)                        ## Total number of species
    frq.xdat <- t (as.matrix (colSums (IO.xdat))) ## Species frequencies
    w3 <- N.xdat * ((1 / frq.xdat) + (1 / (N.xdat - frq.xdat))) - 1

    ## Distance matrix
    ## ... using vegan
    if (distance == "bray") dst.xdat <- vegdist (xdat, "bray")     
    if (distance == "jaccard") dst.xdat <- vegdist (xdat, "jaccard")
    ## ... using proxy
    if (distance != "bray" & distance != "jaccard")    
      dst.xdat <- dist (xdat, method = distance) 
    
    ## Exit with dignity if N.xdat < 3 (applies to the first hierarchy level)
    if (N.xdat < 3)                                                      ## <-----
    {
      mess1 <- "Not enough data."
      if (juice == TRUE) 
      {
      write.table (mess1, "isopam/alert.txt", row.names = FALSE, 
        col.names = FALSE, quote = FALSE)
      }
      stop (mess1)
    }

    ## Determine maximum Isomap k
    if (k.max > N.xdat - 1) k.max <- N.xdat - 1

    ## Determine minimum Isomap k (from isomapdist.r, vegan)
    dmtr <- as.matrix (dst.xdat)
    diag (dmtr) <- NA

    k.min <- 2

    for (a in c(2:k.max))
    {
      ## Check min. k
      dm <- dmtr

      is.na (dm) <- apply (dm, 2, function (xx) xx > 
        xx [order (xx, na.last = TRUE) [a]])

      dm <- pmax (as.dist(dm), as.dist(t(dm)), na.rm = TRUE)
      fragm <- distconnected (dm, toolong=0,  trace=FALSE)
      if (length (unique (fragm)) > 1)
        k.min <- k.min + 1
      else
        break
    }

    ## Adjust d.max
    if (N.xdat < 11) d.max <- N.xdat - 1
    else d.max <- 10
    
    ## Set c.min
    c.min <- 2
    
    ## Adjust c.max if necessary
    if (c.max > N.xdat - 1) c.max <- N.xdat - 1
    if (c.max < 2) stop ("c.max < 2")
    
    ## Adjust c.num if necessary
    if (is.numeric (c.num))
    {
      if (c.num < 2) stop ("c.num< 2")
      if (c.num > N.xdat - 1) c.num <- N.xdat - 1
      c.min <- c.num
      c.max <- c.num
    }  
    
    ## Prepare progress bar
    b.loops <- (k.max - k.min) + 1
    d.loops <- d.max - 1
    pb.mx <- (b.loops * d.loops) - 1
    pb <- txtProgressBar (min = 0, max = pb.mx, char = ".",
      width = 45, style = 3)

    ## Prepare output array
    out.array <- array (NA, dim = c(10, k.max, c.max))

    ## Using a matrix of logical values is a _hair_ faster inside the loop
    IO.xdat.logical <- IO.xdat != 0

    ## ----------- Begin with the big loops .... ---------------------------- ##

    for (b in c(k.min:k.max))  ## b-loop: Isomap k
    {   

      suppressWarnings (isom <- isomap (dst.xdat, ndim = d.max, k = b)) ## Isomap   

      for (d in 2:d.max) 
      {

        isodiss <- suppressWarnings (daisy (isom$points[,1:d], metric =
          "euclidean", stand = TRUE))

        for (e in c.min:c.max) ## e-loop: Cluster no.
        {                      
          ## --------- Partitioning (PAM) ------------------------------------ #
          
          cl.iso <- pam (isodiss, k = e, diss = TRUE) ## PAM
          cl <- cl.iso$clustering                     ## Group affiliation
          ci <- cl.iso$clusinfo[,1]                   ## Cluster size

          ## ---------- Standardized G-test ---------------------------------- #

          ## For Williams' correction
          w1 <- N.xdat * sum (1 / ci) - 1 
          w2 <- 6 * N.xdat * (e - 1)
          
          ## Compute G-values for all species (code adapted from Lubomir Tichy)
          ## This is a G-test of association
          
          gt <- matrix (NA, SP.xdat, 1) ## Matrix for G-test results

          for (g1 in 1:SP.xdat)         ## g1-loop through species
          {                    
            DDD <- 0
            spec_frq <- frq.xdat [g1]

            ## For Williams' correction
            willi <- 1 + ((w1 * w3 [g1]) / w2)
            
            ## Mask out entries in cl using appropriate col in IO.xdat
            groupids <- cl[IO.xdat.logical[,g1]]; 

            for (g.fast in 1:e)             ## g.fast-loop through clusters
            {                  
              fra1 <- sum (groupids == g.fast)   ## Species occ. in cluster
              Nj <- ci [g.fast]                  ## Cluster size

              ## so need to calculate it and store for later.
              ## This does not get executed very often at all.

              bom <- spec_frq / N.xdat 
              bim <- 1 - bom           
              fra0 <- Nj - fra1        
              bum <- fra1 / (Nj * bom)  
              bam <- fra0 / (Nj * bim) 
              DDDadd <- 0
              
              if (!is.na (bum)) 
                if (bum > 0) DDDadd <- DDDadd + (fra1 * log (bum)) 
              
              if (!is.na (bam)) 
                if (bam > 0) DDDadd <- DDDadd + (fra0 * log (bam))             
              
              DDD <- DDD + DDDadd
              
            }
            DDD <- DDD * 2
            gt [g1,] <- DDD / willi ## Williams' correction
          }

          ## Standardization (Botta-Dukat et al. 2005)
          gt.ex <- e - 1                 ## Expected G
          gt.sd <- sqrt (2 * gt.ex)      ## Expected sd
          G <- (gt - gt.ex) / gt.sd
          
          ## ----------------------------------------------------------------- #
          
          ## Averaging
          if (filtered == FALSE) out.array [d,b,e] <- mean (G)
          
          ## Standard: Filtering and averaging
          if (filtered == TRUE) 
          {
            ## Filtering by G
            glgth <- length (G [G >= g.min]) 
            if (glgth == 0) out.array [d,b,e] <- NA

            ## Averaging
            else out.array [d,b,e] <- mean (G [G >= g.min]) * glgth
          }
        }
        setTxtProgressBar (pb, ((b - k.min) * (d.max - 1)) + d - 1)
      }
    }

    close (pb) ## Close progress bar

    ## ----------- End parameter search ------------------------------------- ##
    
    solution <- TRUE
    out.array [is.na (out.array)] <- 0
    if (sum (out.array) <= 0) 
    { 
      mess2 <- "No solution found with current settings"
      
      if (count == 1) 
      {
        if (juice == TRUE) 
        {
          write.table (mess2, "isopam/alert.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE)
        }
        stop ("No solution found with current settings")
      }
      else solution <- FALSE
    }  
    
    if (solution == TRUE) 
    {
      mn.iso <- max (out.array)

      ## ----------- Parameters for final run ------------------------------- ##

      wmx.iso <- which (out.array == mn.iso, arr.ind = TRUE) ## Cases with max
      colnames (wmx.iso) <- c("iso.dim", "iso.k", "clusters")
      
      ## In case of multiple best solutions select the one with max clusters, 
      ## max dims and max k (in this rank order)

      ## Select cases with maximum number of clusters
      try (wmx.iso <- wmx.iso [which (wmx.iso [,3] == max
        (wmx.iso [,3])),], silent = TRUE)
      ## Select cases with maximum number of Isomap dimensions
      try (wmx.iso <- wmx.iso [which (wmx.iso [,1] == max
        (wmx.iso [,1])),], silent = TRUE)
      ## Select cases with maximum k
      try (wmx.iso <- wmx.iso [which (wmx.iso [,2] == max
        (wmx.iso [,2])),], silent = TRUE)

      mc <- wmx.iso [3]; md <- wmx.iso [1]; mk <- wmx.iso [2]
  
      ## ----------- Final run ---------------------------------------------- ##

      suppressWarnings (isom <- isomap (dst.xdat, ndim = d.max, k = mk))          
      d.iso <- daisy (isom$points[,1:md], metric = "euclidean", stand = TRUE)
      cl.iso <- pam (d.iso, k = mc, diss = TRUE)

      CLS <- cl.iso$clustering                 ## Group affiliation
      MDS <- cl.iso$medoids                    ## Medoids
      CLI <- t (cl.iso$clusinfo [,1])          ## Cluster size
      ISO <- isom$points [,1:3]                ## Isomap scores    

      ## Contingency table
      tab <- t (aggregate (IO.xdat, by = list (CLS), FUN = sum)[,-1]) 
  
      ## Frequency table
      inf <- matrix (rep (CLI, SP.xdat), nrow = SP.xdat, byrow = TRUE)
      FRQ <- tab / inf 
  
      ## G-statistics
      w1 <- N.xdat * sum (1 / CLI) - 1 
      w2 <- 6 * N.xdat * (mc - 1)

      g.1 <- matrix (NA, SP.xdat, 1)  
      for (g3 in 1:SP.xdat)           
      {                             
        willi <- 1 + ((w1 * w3 [g3]) / w2)
        DDD <- 0
        bom <- frq.xdat [g3] / N.xdat 
        bim <- 1 - bom                
        rip <- tab [g3,]              
        
        for (g4 in 1:mc)              
        {                           
          fra1 <- rip [g4]            
          Nj <- CLI [g4]              
          fra0 <- Nj - fra1           
          bum <- fra1 / (Nj * bom)              
          bam <- fra0 / (Nj * bim)            
  
          ## adding up values
          if (!is.na (bum))   
            if (bum > 0) DDD <- DDD + (fra1 * log (bum))
          if (!is.na (bam))   
            if (bam > 0) DDD <- DDD + (fra0 * log (bam))
        }
        DDD <- DDD * 2
        g.1 [g3,] <- DDD / willi
      }

      ## Standardization (Botta-Dukat et al. 2005)
      gt.ex <- mc - 1                        ## Expected G
      gt.sd <- sqrt (2 * (mc - 1))           ## Expected sd
      sG <- (g.1 - gt.ex) / gt.sd
           
      ## Number of indicators >= g.min
      noi <- length (sG [sG >= g.min])      
     
      ## Was this a good partition? 
      ## At least stopcrit [1] descriptors with g >= stopcrit [2]
      if (length (sG [sG >= stopcrit [2]]) >= stopcrit [1]) fine <- TRUE                        
      else fine <- FALSE
      
      out <- list (
        medoids = MDS, 
        clusters = CLS,
        sizes = CLI,
        is.ok = fine,
        scores = ISO, 
        k.min = k.min,
        k.max = k.max,
        k = mk,
        d = md,
        noi = noi
        )
    }
    else 
    {    

      out <- list (
        medoids = NULL,
        clusters = NULL,
        sizes = NULL,
        is.ok = FALSE,
        scores = NULL, 
        k.min = NULL,
        k.max = NULL,
        k = NULL,
        d = NULL,
        noi = NULL
     )
    }

    return (out)
  }  
  
  ## ----------- dendrogram function (code: J. Collison, 2009) -------------- ##

  create_dendro <- function (clust)
  {
    ## Expects a list of vectors containing cluster affiliations
    ## without info about hierarchy (running number).
    ## Returns an object of class "hclust"
  
    num_obs <- length(clust[[1]])
  
    dendro <- list(
      merge=array(NA, dim=c(num_obs-1,2)),
      height=rep(0, times=num_obs-1),
      order=c(),
      labels=names(clust[[1]]),
      method=NULL,
      call=NULL);
  
    class(dendro) <- "hclust"
  
  
    ## hclust requires a set of merge operations, one observation at a time.
    ## Subsequent operations [may] refer to the index of a previous merge,
    ## so we store the operation index after everything we do.
  
    opnum <- 0
  
    group_opnums <- c();
  
    for (level in length(clust):1)
    {
      groups <- clust[[level]];
  
      groupnum <- 0
  
      curlevel_group_opnums <- c();
  
      repeat
      {
        groupnum <- groupnum+1
  
        log_in_group <- (groups == groupnum)
  
		num_in_clust <- sum(log_in_group);

		if (num_in_clust < 1)
		{
			break
		}

		if (level == length(clust))
		{
		  ## Bottom level
		  ## Join all these and add them to the order vector

		  prev_opnum <- 0
		  prev_index <- 0

		  for (j in 1:num_obs)
          {
            if (!log_in_group[j])
              next
  
            dendro$order <- c(dendro$order, j)
  
            if (prev_opnum)
            {
              opnum <- opnum+1
              dendro$merge[opnum,] <- c(-j,prev_opnum)
              dendro$height[opnum] <- 1
			  prev_opnum <- opnum
            }
            else if (prev_index)
            {
              opnum <- opnum+1
              dendro$merge[opnum,] <- c(-j,-prev_index)
              dendro$height[opnum] <- 1
              prev_opnum <- opnum
            }
            else
            {
              prev_index <- j
            }
          }
  
          ## Done merging this cluster
          ## Save the opnum for the last merge in this cluster
          curlevel_group_opnums <- c(curlevel_group_opnums,0)
		  curlevel_group_opnums[groupnum] <- prev_opnum

		  ## Special case for singletons
		  if (!prev_opnum)
		  {
			  curlevel_group_opnums[groupnum] <- -prev_index
		  }
        }
        else
        {
          ## higher levels
  
          ## For all members of this group, see what their groupid was
          ## one level deeper and store those (uniquely)
  
          groups_to_join <- c()
  
          for (j in 1:num_obs)
          {
            if (!log_in_group[j])
              next
  
            subgroup = clust[[level+1]][j]
  
            if (sum(groups_to_join == subgroup) == 0)
            {
              ## Not there, add it
              groups_to_join <- c(groups_to_join,subgroup)
            }
          }
  
          if (length(groups_to_join)>=2)
          {
            ## Merge them.
            prev_opnum = group_opnums[groups_to_join[1]];
            for (r in 2:length(groups_to_join))
            {
              opnum <- opnum+1
              dendro$merge[opnum,] <- c(prev_opnum, 
                group_opnums[groups_to_join[r]])
              dendro$height[opnum] <- length(clust) - level + 1
              prev_opnum <- opnum;
            }
  
            curlevel_group_opnums <- c(curlevel_group_opnums,0)
            curlevel_group_opnums[groupnum] <- prev_opnum
          }
          else
          {
            ## Nothing to merge this one with at this time;
            ## bubble up for next level
            prev_opnum = group_opnums[groups_to_join[1]];
            curlevel_group_opnums <- c(curlevel_group_opnums,0)
            curlevel_group_opnums[groupnum] <- prev_opnum
          }
        }
      } ## end groupnum loop
  
      ## Done processing this level
      group_opnums <- curlevel_group_opnums
    }
  
    ## All levels have been processed.
    ## Need to merge whatever remains
  
    if (length(group_opnums)>1)
    {
      prev_opnum = group_opnums[1];
  
      for (r in 2:length(group_opnums))
      {
        opnum <- opnum+1
        dendro$merge[opnum,] <- c(prev_opnum,group_opnums[r])
        dendro$height[opnum] <- length(clust) + 1
        prev_opnum <- opnum;
      }
    }
  
    return(dendro)
  }

  ## ------------ Isopam call ----------------------------------------------- ##

  print ("Processing level 1", quote = FALSE)

  ## Prepare first partition  
  IO <- dat
  IO [IO > 0] <- 1
  dat1 <- dat [,colSums (IO) > 1]  ## Omit species with < 2 (!) occurrences    
  if (is.null (dim (dat1))) stop ("Not enough occurrences")      
  dat1 <- dat1 [,apply (dat1, 2, var) > 0] ## Omit species without variance
  if (is.null (dim (dat1))) stop ("Not enough variance in species")      
  dat1 <- dat1 [apply (dat1, 1, var) > 0,] ## Omit plots without variance
  if (is.null (dim (dat1))) stop ("Not enough variance in plots")      

  ## Initiate container for results
  matr <- matrix (NA, nrow = nrow (dat1), ncol = 1) ## Initiate cluster output
  rownames (matr) <- rownames (dat1)
  colnames (matr) <- "lev.1"

  ## Initiate matrix for summary ("analytics")
  summ <- matrix (NA, nrow = 7, ncol = 1) 
  rownames (summ) <- c("Name", 
    "Subgroups", 
    "Isomap.dim", 
    "Isomap.k.min",
    "Isomap.k", 
    "Isomap.k.max",
    "Indicators") 
  colnames (summ) <- "Part.1"
                           
  ## Run core function
  output <- core (dat1)
 
  ## Fill cluster container
  matr [,1] <- output$clusters

  ## Fill summary matrix
  summ [1,1] <- 0                          ## Name
  summ [2,1] <- length (output$medoids)    ## No. of subgroups
  summ [3,1] <- output$d                   ## Isomap dimensions
  summ [4,1] <- output$k.min               ## Minimum k
  summ [5,1] <- output$k                   ## Isomap k
  summ [6,1] <- output$k.max               ## Maximum k
  summ [7,1] <- output$noi                 ## No. of indicators used

  ## Medoids
  med <- list ()
  med [[1]] <- output$medoids
  names (med [[1]]) <- c (1:length (output$medoids))

  ## stop if maxlev is 1
  if (maxlev == 1) stepdown <- FALSE
  else stepdown <- TRUE
  ## Which group is large enough?
  spl <- as.numeric (output$sizes > 2) 
  ## stop if there is nothing to split
  if (sum (spl) == 0) stepdown <- FALSE
                 
  ## Preparative stuff   
  ctb <- matr ## Cluster table 
  colnames (ctb) <- "lev.1"
  mtb <- matrix (NA, nrow = nrow (dat), ncol = 0) ## Medoid table
  rownames (mtb) <- rownames (dat)

  ## Create ctb.flat for the one level case:
  ## cluster affiliations without info about hierarchy (running number)
  ctb.flat <- as.numeric (as.factor (ctb))
  names (ctb.flat) <- rownames (ctb)

  count <- 2 ## Counter for cluster levels
  
  ## -------------- Follow-up runs ------------------------------------------ ##
  
  while (stepdown == TRUE) 
  {
    if (maxlev != FALSE & count > maxlev) stepdown <- FALSE
    if (stepdown == TRUE)
    {
    
      ifelse (sum (spl) == 1, cas <- "group", cas <- "groups")
      if (cas != 0) print (paste ("Level ", count, ": Processing ", sum (spl), 
        " ", cas, sep = ""), quote = FALSE)
      output.sub <- list ()              ## Empty list for results
      count.2 <- 1
      for (j in 1:length (spl))          ## Loop through partitions
      {                       
        if (spl [j] == 1)                ## splittable? 
        {                                            
          x.sub <- dat [matr [,ncol (matr)] == j,]    ## Create data subset

          ## Prepare second run
          IO.x.sub <- x.sub
          IO.x.sub [IO.x.sub > 0] <- 1
          x.sub <- x.sub [,colSums (IO.x.sub) > 1]
          if (is.null (dim (x.sub))) output.sub [[j]] <- NA  ## one species left 
          else
          {
            ## Omit species and plots without variance
            x.sub <- x.sub [,apply (x.sub, 2, var) > 0]
            x.sub <- x.sub [apply (x.sub, 1, var) > 0,]                        
            ## and go ...
            print (paste ("Group", count.2), quote = FALSE)
            output.sub [[j]] <- core (x.sub)            ## Clustering
            count.2 <- count.2 + 1
          }
        }
        else output.sub [[j]] <- NA
      }
  
      ## Look if some of the partitions should be rejected
      ok.vec <- vector ()
      for (l in 1:length (spl)) 
      {
        if (is.na (output.sub [[l]][1]) == TRUE) ok.vec [[l]] <- FALSE 
        else ok.vec [[l]] <- output.sub [[l]]$is.ok
      }
  
      ## Report
      ok.n <- sum (as.numeric (ok.vec))
      
      ## Stop if all partitions are rejected
      if (ok.n == 0) stepdown <- FALSE
          
      if (stepdown == TRUE)
      {
        ## Write new clusters to container 
        matr <- cbind (matr, matrix (NA, nrow = nrow (matr), ncol = 1))
        colnames (matr) [ncol (matr)] <- paste ("lev.", ncol (matr), sep = "")
        
        for (m in 1:length (spl))  ## Fill in values 
        {                     
          if (ok.vec [m] == TRUE) 
          {
            sc <- output.sub [[m]]$clusters
            idx.sc <- rownames (matr) %in% names (sc) 
            matr [idx.sc, ncol (matr)] <- sc
          }
        }            
    
        ## Make matrix with names expressing hierarchy
        ctb <- cbind (ctb, matrix (0, nrow = nrow (dat), ncol = 1))       
        colnames (ctb) [ncol (ctb)] <- paste ("lev.", ncol (ctb), sep = "")
        ctb.red <- na.omit (ctb)
        blub <- matr [,ncol (matr)]
        blub [is.na (blub)] <- 0
        ctb.red [,ncol (ctb.red)] <- paste (ctb.red [,ncol (ctb.red)-1],
          blub, sep = ".")  
        ctb [rownames (ctb) %in% rownames (ctb.red), ncol (ctb)] <-
          ctb.red [,ncol (ctb.red)]    
        ctb <- as.data.frame (ctb)
        
        ## Overwrite ctb.flat in the multiple level case
        ## Cluster affiliations without info about hierarchy (running number)
        ctb.flat <- list ()
        for (o in 1:ncol (ctb)) 
        {
          ctb.flat [[o]] <- as.numeric (as.factor (ctb [,o]))
          names (ctb.flat [[o]]) <- rownames (ctb)
          names (ctb.flat) [o] <- paste ("level.", o, sep="")
        }
    
        ## Make working matrix
        matr.red <- na.omit (matr)
        matr.red  [,ncol (matr.red)] <- as.numeric (as.factor ( paste
          (matr.red [,ncol (matr)-1], matr.red [,ncol (matr.red)],
          sep = ".")))
        matr [rownames (matr) %in% rownames (matr.red),] <- matr.red
        matr [is.na (matr)] <- 0
    
        ## Medoids
        mtb <- cbind (mtb, matrix (0, nrow = nrow (dat), ncol = 1))
        rownames (mtb) <- rownames (dat)
        for (n in 1:length (spl)) 
        {
            if (ok.vec [n] == TRUE) 
            {
              mtb [output.sub [[n]]$medoids, ncol (mtb)] <- 1
            }
        }
        med [[count]] <- names (mtb [mtb [,ncol (mtb)] == 1, ncol (mtb)])
        names (med [[count]]) <- ctb [mtb [,ncol (mtb)] == 1, ncol (ctb)]
        nam <- sort (names (med [[count]]))
        med [[count]] <- med [[count]][nam]
        
        ## Add new partitions to the summary matrix            
        for (x in 1:length (spl))
        {
          if (ok.vec [x] == TRUE)  
          {
            summ <- cbind (summ, matrix (NA, nrow = 7, ncol = 1))
            colnames (summ) [ncol(summ)] <- paste ("Part.", ncol(summ), sep = "")                                      
            summ [2,ncol(summ)] <- length (output.sub[[x]]$medoids) ## Subgroups
            summ [3,ncol(summ)] <- output.sub[[x]]$d             ## Isomap dims
            summ [4,ncol(summ)] <- output.sub[[x]]$k.min         ## Minimum k            
            summ [5,ncol(summ)] <- output.sub[[x]]$k             ## Selected k             
            summ [6,ncol(summ)] <- output.sub[[x]]$k.max         ## Maximum k            
            summ [7,ncol(summ)] <- output.sub[[x]]$noi           ## Indicators    
          }
        }
  
        ## OK, now we have a matrix with cluster affiliations of this level
        ## Repeat the search for splittable units:
    
        subtab <- table (matr [,ncol (matr)])                    
        subtab <- subtab [names (subtab) != '0'] 
        spl <- as.numeric (subtab > 2)                                
        names (spl) <- names (subtab)
               
        count <- count + 1
        
        ## Create object of class 'hclust'
        dendro <- create_dendro(ctb.flat)
      }
    }     
  } ## End while-loop    

  ## Cluster names for the summary matrix
  if (ncol (ctb) > 1)
  {
    cnam <- vector ()
    for (y in 1:length(med))
    {
      cnam <- c (cnam, names (med [[y]]))
    }

    childs <- vector ()
    for (z in 1:length (cnam))
    {
      nch <- nchar (cnam [z])
      childs <- c (childs, length (grep (cnam [z], substr (cnam, 1, nch), 
        value = TRUE)))
    }
    summ [1,2:ncol(summ)] <- cnam [childs > 1]  
  }
  summ <- as.data.frame (summ)

  ## ------------ Output ---------------------------------------------------- ##
  
  if (ncol (ctb) == 1)
    OUT <- list (
      call = sys.call (),
      distance = distance,
      flat = ctb.flat,
      hier = NULL,
      medoids = med,
      analytics = summ,
      dendro = NULL,
      dat = dat
    )
  if (ncol (ctb) > 1)
    OUT <- list (
      call = sys.call (),
      distance = distance,
      flat = ctb.flat,
      hier = ctb,
      medoids = med,
      analytics = summ,
      dendro = dendro,
      dat = dat
    )

    class (OUT) <- "isopam"
    
  ## save output
  if (juice == TRUE)
  {
    write.table (ctb, file = "isopam/juicein.txt", 
      col.names = FALSE, quote = FALSE)  
    if (!is.null (OUT$dendro [1]))
    {
      wth <- (nrow (ctb) * 11) + 100
      bmp (file = "isopam/juicetree.bmp", width = wth)
      plot (OUT$dendro)
      dev.off()
    }  
  }
  
  invisible (OUT)
}

plot.isopam <-
function (x, ...)
{
  if (!is.null (x$dendro [1])) 
  {  
    tree <- x$dendro
    plot (tree, main = format(x$call), ...)  
  }     
  else print ("No cluster hierarchy", quote = FALSE)
}

identify.isopam <-
function (x, ...)
{
  identify (x$dendro, MAXCLUSTER = nrow (x$dat), ...)  
}

   