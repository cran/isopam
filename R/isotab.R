isotab <-
function (ip, level = 1, phi.min = .5, p.max = .05)
{ 
  IO <- ip$dat
  IO [IO > 0] <- 1 
  N <- nrow (IO)
  SP <- ncol (IO)
  frq <- t (as.matrix (colSums (IO)))
  
  if (is.null (ip$hier)) 
  {
    if (level > 1) print ("No hierarchy levels available", quote = FALSE)
    level <- 1
    tab <- t (aggregate (IO, by = list (ip$flat), FUN = sum))
  }
  else
  { 
    depth <- ncol (ip$hier)
    if (level > depth) 
    {
      level <- depth 
      print (paste ("Switching to lowest level", depth), 
        quote = FALSE)
    }
    tab <- t (aggregate (IO, by = list (ip$hier [,level]), FUN = sum))
  }
  
  ## 1) Contingency table  
  cnam <- tab [1,]
  tab <- tab [-1,]
  rnam <- rownames (tab)
  tab <- matrix (as.numeric (tab), nrow = nrow (tab))
  colnames (tab) <- cnam
  rownames (tab) <- rnam
  nc <- ncol(tab)
  
  ## 2) Frequency table
  if (is.null (ip$hier)) siz <- table (ip$flat) ## Cluster sizes
  else siz <- table (ip$hier [,level]) ## Cluster sizes
  spc <- t (as.matrix (siz))[rep (1, nrow (tab)),]
  frq.2 <- tab / spc ## Frequency
  frq.2 <- round (frq.2 * 100, 0) ## as percentage

  ## 3) Fisher's significance table
  ft <- tab
  for (u in 1:SP) 
  {
    for (v in 1:nc) 
    {
      insd <- tab [u, v]
      outs <- tab [u,-v]
      x2 <- matrix (c(
              insd,                        ## a occ. of species in cluster v
              siz [v] - insd,              ## c abs. in cluster v
              sum (outs),                  ## b occ. outside cluster v
              sum (siz [-v]) - sum (outs)  ## d abs. outside v
              ), 2, 2)    
      ft [u,v] <- fisher.test (x2)$p.value
    }
  }

  ## 4) Significance symbols
  ft.symb <- ft
  ft.symb [ft > 0.05] <- ""
  ft.symb [ft <= 0.05] <- "*"
  ft.symb [ft <= 0.01] <- "**"
  ft.symb [ft <= 0.001] <- "***"

  ## 5) Combined frequency table with significance symbols 
  frq.ft <- matrix (paste (frq.2, ft.symb, sep = ""), 
      nrow(frq.2), ncol(frq.2))
  frq.ft <- data.frame (frq.ft)
  colnames (frq.ft) <- colnames (ft.symb)
  rownames (frq.ft) <- rownames (ft.symb)

  ## 6) Standardized phi table 
  S <- 1 / nc                                   ## Constant s (Tichy '06)
  cs <- S * N                                   ## new cluster sizes
  phi <- tab

  for (i in 1:SP) 
  {
    for (j in 1:nc) 
    {
      insd <- tab [i, j]                        ## original n in cluster j
      outs <- sum (tab [i,-j])                  ## original n outside cluster j
      oc <- cs * (insd / siz [j])                ## new n in cluster j
      on <- (N - cs) * (outs / (N - siz [j]))    ## new n outside cluster j
      total <- oc + on                          ## new total value
      phi.1 <- nv <- (N * oc - total * cs) 
      phi.2 <- sqrt (total * cs * (N - total) * (N - cs))            
      nv <- phi.1 / phi.2
      phi [i,j] <- nv
    }
  }
  phi [is.na (phi)] <- 0                        ## Replace NaN-values by 0

  ## Table sorting     
  ## 7) .... by phi and frequency
  phi.idx <- apply (phi, 1, which.max) ## Group affiliation by phi 
  frq.ord <- phi.idx  
  for (i in 1:length(frq.ord)) frq.ord [i] <- frq.2 [i, phi.idx [i]] 
  frq.top <- t (frq) [order (phi.idx, -frq.ord),] ## Sorting
  ord.top <- names (frq.top)
  frq.ft.top <- frq.ft [ord.top,]
  ft <- ft [ord.top,]
  phi <- phi [ord.top,]

  ## 8) Filter diagnostic species
  filter1 <- apply (ft, 1, min) <= p.max                                          
  filter2 <- apply (phi, 1, max) >= phi.min
  dia <- which (filter1 [filter2 == TRUE] == TRUE) ## diagnostic species 
  n.dia <- length (dia) ## how many diagnostic species
  if (n.dia == 0) diag <- "No diagnostic species with given thresholds." 
  if (n.dia > 0) diag <- frq.ft.top [names (dia),]

  ## 9) For later use in the bottom part of the tables
  ord.bot <- names (t (frq) [order (-frq),])
  frq.ft.b <- frq.ft [ord.bot,]

  ## 10) Move diagnostic species to top
  if (n.dia > 0) 
  {
    FRQ <- rbind (diag, frq.ft.b [rownames (frq.ft.b) %in% 
      rownames (diag) == FALSE,])
  }
  else FRQ <- frq.ft.b
  
  ## 11) Report cluster sizes
  siz <- t (as.matrix (siz))
  rownames (siz) <- "n"
  
  ## 12) Output
  isotab.out <- list (
   n = siz,
   tab = FRQ)
  
  return (isotab.out)  
}

