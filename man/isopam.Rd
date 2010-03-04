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
     isopam (dat, max.level = FALSE, fixed.number = FALSE, 
            opt.number = TRUE, max.number = 6, stopat = c(1,7), 
            filtered = TRUE, thresh = 3.5, distance = 'bray', 
            k.max = 100, juice = FALSE) 
     
     \method{identify}{isopam}(x, \dots)
     \method{plot}{isopam}(x, \dots)
     }
     \arguments{
       \item{dat}{data matrix: each row corresponds to an 
         object (typically a plot), each column corresponds 
         to a descriptor (typically a species). All variables must 
         be numeric. Missing values (NAs) are not allowed. At least
         3 rows (plots) are required.}
       \item{max.level}{maximum number of hierarchy levels. Defaults 
         to \code{FALSE} (no maximum number). Note that divisions
         may stop well before this number is reached (see 
         \code{stopat}). Use \code{max.level = 1} for non-hierarchical 
         partitioning (or use \code{fixed.number}).}
       \item{fixed.number}{number of clusters (defaults to \code{FALSE}). 
         If a number is given, non-hierarchical partitioning is performed, 
         \code{opt.number} and \code{max.number} are ignored and 
         \code{max.level} is set to one.}
       \item{opt.number}{if \code{TRUE} (the default) cluster 
         numbers are optimized in the range between 2 and 
         \code{max.number} (slow and thorough). If \code{FALSE} 
         groups are divided into two subgroups (quick and dirty).}
       \item{max.number}{maximum number of clusters per partition.
         Applies to all partitioning steps if 
         \code{opt.number = TRUE}.}
       \item{stopat}{vector with stopping rules for hierarchical
         clustering. Two values define if a partition should be 
         retained in hierarchical clustering: the first determines 
         how many indicators must be present per cluster, the second 
         defines the standardized G-value that must be reached by 
         these indicators. \code{stopat} is not effective at the first 
         hierarchy level or in non-hierarchical partitioning.}                                      
       \item{filtered}{logical. If \code{TRUE} (the deafult), only 
         descriptors (species) exceeding a threshold defined by 
         \code{thresh} are used in the search for a good clustering 
         solution. Their number is multiplied with their mean 
         standardized G-value. The product is used as optimality 
         criterion. If \code{FALSE} all descriptors are used for 
         optimization.}
       \item{thresh}{threshold (standardized G value) for descriptors 
         (species) to be considered in the search for a good 
         clustering solution. Effective with \code{filtered = TRUE}.}
       \item{distance}{distance measure for the distance matrix
         used as a starting point for Isomap. All but the 
         Bray-Curtis and the Jaccard measure are passed to the 
         \code{method} argument in package \pkg{proxy} (see details).} 
       \item{k.max}{maximum Isomap \emph{k}.}
       \item{juice}{logical. If \code{TRUE} input files for Juice are
         generated.}
       \item{x}{an \code{isopam} result object.}
       \item{...}{other arguments to S3 functions \code{plot} and 
         \code{identify} corresponding to \link[stats]{hclust}).}         
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
         the final partitioning steps. These are the name of the 
         parent cluster (0 in case of the first partition), the
         number of subgroups, Isomap dimensions, Isomap \emph{k} 
         used, and the number of indicators reaching or exceeding 
         \code{thresh}.} 
       \item{dendro}{an object of class \code{hclust} representing
         the clustering. Not present with only one level of 
         partitioning.}  
       \item{dat}{data used}         
     }        
     \details{
       Isopam is described in Schmidtlein et al. (2010). 
       It consists of dimensionality reduction (Isomap: Tenenbaum 
       et al. 2000; \code{isomap} in package \pkg{vegan}) and 
       partitioning of the resulting ordination space (PAM: Kaufman 
       & Rousseeuw 1990; \code{pam} in package \pkg{cluster}). 
       The classification is performed either as a hierarchical, 
       divisive method, or as non-hierarchical partitioning.  
       Compared to other clustering methods, it has the following 
       features: (a) it optimizes partitions for the performance of 
       group indicators (typically species); (b) it optionally 
       selects the number of clusters per division; (c) the shapes 
       of groups in feature space are not limited to spherical or 
       other regular geometric shapes (thanks to the underlying Isomap 
       algorithm) and (d) the distance measure used for the initial 
       distance matrix can be freely defined. 
       
       Currently, the \code{plot} and \code{identify} methods 
       for class \code{isopam} simply link to the 
       \link[stats]{hclust} object \code{$dendro} resulting 
       from \code{isopam} in case of hierarchical partitioning. 
       The methods work just like \code{plot.hclust} and 
       \code{identify.hclust}.

       The preset distance measure is Bray-Curtis 
       (Odum 1950). Bray-Curtis (\code{'bray'}) and Jaccard 
       distances (\code{'jaccard'}) are passed to \code{vegdist} 
       in \pkg{vegan}. All other measures are passed to the 
       \code{method} argument in package \pkg{proxy}. This package
       contains a growing number of relevant measures. Measures 
       registered in \pkg{proxy} can be listed with 
       \code{summary(pr_DB)} once \pkg{proxy} is loaded. New 
       measures can be defined and registered as described in 
       \code{?pr_DB}.       
     }                                     
     \references{
       Odum, E.P. (1950): Bird populations in the Highlands (North 
       Carolina) plateau in relation to plant succession and avian 
       invasion. \emph{Ecology} \bold{31}: 587--605.
       
       Kaufman, L., Rousseeuw, P.J. (1990): \emph{Finding groups in 
       data}. Wiley.
       
       Schmidtlein, S., \enc{Tichý}{Tichy}, L., Feilhauer, H., Faude, U.
       (2010): A brute force approach to vegetation classification.
       in review.
       
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
     ip<-isopam(andechs,fixed.number=3)
     ip$flat 
     }
     \keyword{ cluster }