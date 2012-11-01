pkgname <- "isopam"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
library('isopam')

assign(".oldSearch", search(), pos = 'CheckExEnv')
cleanEx()
nameEx("isopam")
### * isopam

flush(stderr()); flush(stdout())

### Name: isopam
### Title: Isopam (Clustering)
### Aliases: isopam plot.isopam identify.isopam
### Keywords: cluster

### ** Examples

     ## load data to the current environment
     data(andechs)
     
     ## call isopam with the standard options
     ip<-isopam(andechs)
     
     ## examine cluster hierarchy
     plot(ip)
     
     ## examine grouping
     ip$flat 

     ## examine frequency table (second hierarchy level)
     isotab(ip, 2)
     
     ## non-hierarchical partitioning
     ip<-isopam(andechs,c.fix=3)
     ip$flat
     
     


cleanEx()
nameEx("isotab")
### * isotab

flush(stderr()); flush(stdout())

### Name: isotab
### Title: Ordered frequency table for Isopam clusters
### Aliases: isotab

### ** Examples

   ## load data to the current environment
   data(andechs)
     
   ## call isopam with the standard options
   ip<-isopam(andechs)
    
   ## build table (uppermost hierarchy level)
   isotab(ip)

   ## build table (lower hierarchy level)
   isotab(ip,2)



### * <FOOTER>
###
cat("Time elapsed: ", proc.time() - get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
