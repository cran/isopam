     \encoding{UTF-8}
     \name{isopam}
     \Rdversion{1.1}
     \alias{isopam}
     \alias{plot.isopam}
     \alias{identify.isopam}
     \title{Isopam (Clustering)}
     \description{
       Isopam classification is performed either as a 
       hierarchical, divisive method, or as non-hierarchical 
       partitioning. Optimizes clusters and optionally cluster 
       numbers for maximum performance of group indicators. 
       Developed for matrices representing species abundances 
       in plots.        
     }
     \usage{
     isopam (dat, c.fix = FALSE, c.opt = TRUE, c.max = 6, 
            l.max = FALSE, stopat = c(1,7), sieve = TRUE, 
            Gs = 3.5, ind = NULL, centers = NULL, distance = 'bray', 
            k.max = 100, d.max = 7, \dots, juice = FALSE) 
     
     \method{identify}{isopam}(x, \dots)
     \method{plot}{isopam}(x, \dots)
     }
     \arguments{
       \item{dat}{data matrix: each row corresponds to an 
         object (typically a plot), each column corresponds 
         to a descriptor (typically a species). All variables must 
         be numeric. Missing values (NAs) are not allowed. At least
         3 rows (plots) are required.}
       \item{c.fix}{number of clusters (defaults to \code{FALSE}). 
         If a number is given, non-hierarchical partitioning is 
         performed, \code{c.opt} and \code{c.max} are ignored and 
         \code{l.max} is set to one.}
       \item{c.opt}{if \code{TRUE} (the default) cluster 
         numbers are optimized in the range between 2 and 
         \code{c.max} (slow and thorough). If \code{FALSE} 
         groups are divided into two subgroups (quick and dirty).}
       \item{c.max}{maximum number of clusters per partition.
         Applies to all partitioning steps if \code{c.opt = TRUE}.}
       \item{l.max}{maximum number of hierarchy levels. Defaults 
         to \code{FALSE} (no maximum number). Note that divisions
         may stop well before this number is reached (see 
         \code{stopat}). Use \code{l.max = 1} for non-hierarchical 
         partitioning (or use \code{c.fix}).}
       \item{stopat}{vector with stopping rules for hierarchical
         clustering. Two values define if a partition should be 
         retained in hierarchical clustering: the first determines 
         how many indicators must be present per cluster, the second 
         defines the standardized G-value that must be reached by 
         these indicators. \code{stopat} is not effective at the first 
         hierarchy level or in non-hierarchical partitioning.}                                      
       \item{sieve}{logical. If \code{TRUE} (the deafult), only 
         descriptors (species) exceeding a threshold defined by 
         \code{Gs} are used in the search for a good clustering 
         solution. Their number is multiplied with their mean 
         standardized G-value. The product is used as optimality 
         criterion. If \code{FALSE} all descriptors are used for 
         optimization.}
       \item{Gs}{threshold (standardized G value) for descriptors 
         (species) to be considered in the search for a good 
         clustering solution. Effective with \code{sieve = TRUE}.}
       \item{ind}{optional vector of column names from \code{dat}
         defining descriptors (species) used as indicators. This 
         turns Isopam in an expert system. Replaces the automated
         selection of indicators with \code{sieve = TRUE}
         (\code{ind} overrules \code{sieve}).}
       \item{centers}{optional vector with observations used as 
	     cluster cores (supervised classification).}
	   \item{distance}{distance measure for the distance matrix
         used as a starting point for Isomap. Any distance measure
         implemented in packages \pkg{vegan} or \pkg{proxy} can be 
         used (see details).} 
       \item{k.max}{maximum Isomap \emph{k}.}
       \item{d.max}{maximum number of Isomap dimensions.}
       \item{...}{other arguments to S3 functions \code{plot} and 
         \code{identify} corresponding to \link[stats]{hclust}.}         
       \item{juice}{logical. If \code{TRUE} input files for Juice are
         generated.}
       \item{x}{an \code{isopam} result object.}
     }
     \value{
       \item{call}{generating call}
       \item{distance}{distance measure used by Isomap}
       \item{flat}{observations (plots) with group affiliation. 
         Running group numbers for each level of the hierarchy.}
       \item{hier}{observations (plots) with group affiliation. 
         Group identifiers reflect the cluster hierarchy. Not 
         present with only one level of partitioning.}
       \item{medoids}{observations (plots) representing the 
         medoids of the resulting groups.}
       \item{analytics}{table summarizing parameter settings for
         the final partitioning steps. \code{Name}: name of the 
         respective parent cluster (0 in case of the first partition); 
         \code{Subgroups}: number of subgroups; \code{Isomap.dim}:
         Isomap dimensions used; \code{Isomap.k.min}: minimum
         possible Isomap \emph{k}; \code{Isomap.k}: Isomap 
         \emph{k} used; \code{Isomap.k.max}: maximum possible 
         Isomap \emph{k}; \code{Ind.N}: number of indicators 
         reaching or exceeding \code{Gs}; \code{Ind.Gs}: the average 
         standardized G value of these indicators; and 
         \code{Global.Gs}: the average standardized G value of all 
         descriptors.} 
       \item{dendro}{an object of class \code{hclust} representing
         the clustering. Not present with only one level of 
         partitioning.}  
       \item{dat}{data used}         
     }        
     \details{
       Isopam is described in Schmidtlein et al. (2010). 
       It consists of dimensionality reduction (Isomap: Tenenbaum 
       et al. 2000; \link[vegan:isomap]{isomap} in \pkg{vegan}) and 
       partitioning of the resulting ordination space (PAM: Kaufman 
       & Rousseeuw 1990; \link[cluster:pam]{pam} in \pkg{cluster}). 
       The classification is performed either 
       as a hierarchical, divisive method, or as non-hierarchical 
       partitioning. Compared to other clustering methods, it has 
       the following features: (a) it optimizes partitions for the 
       performance of group indicators (typically species) or for 
       maximum average 'fidelity' of descriptors to groups; (b) it 
       optionally selects the number of clusters per division; (c) 
       the shapes of groups in feature space are not limited to 
       spherical or other regular geometric shapes (thanks 
       to the underlying Isomap algorithm) and (d) the distance measure 
       used for the initial distance matrix can be freely defined. 
       
       Currently, the \code{plot} and \code{identify} methods 
       for class \code{isopam} simply link to the 
       \link[stats]{hclust} object \code{$dendro} resulting 
       from \code{isopam} in case of hierarchical partitioning. 
       The methods work just like \code{plot.hclust} and 
       \code{identify.hclust}.

       The preset distance measure is Bray-Curtis (Odum 1950). 
       Distance measures are passed to 
       \link[vegan:vegdist]{vegdist} in \pkg{vegan}. If vegan does 
       not know the given measure it is passed to 
       \link[proxy:dist]{dist} in \pkg{proxy}. Measures available 
       in vegan are listed in \link[vegan:vegdist]{vegdist}. 
       Measures registered in \pkg{proxy} can be listed with 
       \code{summary(pr_DB)} once \pkg{proxy} is loaded. New 
       measures can be defined and registered as described in 
       \code{?pr_DB}. Isopam can't deal with distance matrices
       as a replacement for the original data matrix because it
       operates on individual descriptors (species).       
     }                                     
     \references{
       Odum, E.P. (1950): Bird populations in the Highlands (North 
       Carolina) plateau in relation to plant succession and avian 
       invasion. \emph{Ecology} \bold{31}: 587--605.
       
       Kaufman, L., Rousseeuw, P.J. (1990): \emph{Finding groups in 
       data}. Wiley.
       
       Schmidtlein, S., \enc{Tichý}{Tichy}, L., Feilhauer, H., Faude, U.
       (2010): A brute force approach to vegetation classification.
       \emph{Journal of Vegetation Science} \bold{21}: 1162--1171.
       
       Tenenbaum, J.B., de Silva, V., Langford, J.C. (2000): A global 
       geometric framework for nonlinear dimensionality reduction. 
       \emph{Science} \bold{290}, 2319--2323.        
     }
     \note{
       For large datasets, Isopam may need too much memory or too much 
       computation time. The optimization procedure (selection 
       of Isomap dimensions and -\emph{k}, optionally selection of cluster 
       numbers) is based on a brute force approach that takes its time 
       with large data sets. Low speed is inherent to the method, so don't 
       complain. If used with data not representing species in plots make 
       sure that the indicator approach is appropriate.
       
       With very small datasets, the indicator based optimization may fail. 
       In such cases consider using \code{filtered = FALSE} instead of 
       the default method.
     }
     \author{
       Sebastian Schmidtlein with contributions from Jason Collison and 
       Lubomir \enc{Tichý}{Tichy}        
       }	
     \seealso{
       \code{\link{isotab}} for a table of descriptor (species) 
       frequency in clusters.
     }
     \examples{
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
     
     }
     \keyword{ cluster }
