     \encoding{UTF-8}
     \name{isopam}
     \Rdversion{1.1}
     \alias{isopam}
     \alias{plot.isopam}
     \alias{identify.isopam}
     \title{Isopam (Hierarchical Clustering)}
     \description{
       Performs Isopam which consists of dimensionality reduction
       and partitioning of the resulting feature space. Optimizes 
       clusters and cluster numbers for maximum performance of 
       group indicators. Developed for matrices representing 
       species abundances in plots.        
     }
     \usage{
     isopam(dat, c.num = FALSE, c.max = 10, filtered = TRUE, 
             distance = 'bray', g.min = 3.5, k.max = 100, 
             stopcrit = c(2,7), maxlev = FALSE,
             juice = FALSE) 
     
     \method{identify}{isopam}(x, \dots)
     \method{plot}{isopam}(x, \dots)
     }
     \arguments{
       \item{dat}{data matrix: each row corresponds to an 
         object (typically a plot), each column corresponds 
         to a descriptor (typically a species). All variables must 
         be numeric. Missing values (NAs) are not allowed. At least
         4 rows (plots) are required.}
       \item{c.num}{number of clusters to be computed. If \code{FALSE}
         (the default), cluster numbers are optimized in the range
         between 2 and c.max. If a number is given, non-hierarchical 
         partitioning is performed (\code{maxlev = 1}).}
       \item{c.max}{maximum number of clusters per partition.}
       \item{distance}{distance measure for the distance matrix
         used as a starting point for Isomap. All but the 
         Bray-Curtis and the Jaccard measure are passed to the 
         \code{method} argument in package \pkg{proxy} (see details).} 
       \item{filtered}{logical. If \code{TRUE}, only descriptors 
         (species) exceeding a standardized G-value of g.min are 
         used in the search for the best partition. Their number 
         is multiplied with their mean standardized G-value and 
         the product is used as optimization criterion.}
       \item{g.min}{threshold for descriptors (species) to be 
         considered as indicators during clustering (standardized 
         G-value). Effective with (default) option 
         \code{filtered = TRUE}.}
       \item{k.max}{maximum Isomap \emph{k}.}
       \item{stopcrit}{vector with stopping rules for hierarchical
         clustering. Two values define if a partition should be 
         retained: the first determines how many indicators must be 
         present, the second defines the standardized G-value that 
         must be reached by these indicators. \code{stopcrit} is not 
         used at the first hierarchy level or in non-hierarchical 
         partitioning.}                                      
       \item{maxlev}{maximum number of hierarchy levels. Defaults to
         \code{FALSE} (no maximum number).}
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
         \code{g.min}.} 
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
       Compared to other hierarchical clustering methods, it has 
       the following features: (a) it optimizes partitions for 
       the performance of group indicators (typically species); 
       (b) in this process it selects the number of clusters per 
       division; (c) the shapes of groups in feature space are 
       not limited to spherical or other geometric shapes (thanks 
       to the underlying Isomap algorithm) and (d) the distance 
       measure used for the initial distance matrix can be freely 
       defined. 
       
       In order to examine frequencies of descriptors (species) 
       for an Isopam clustering solution use \link{isotab}. 
                   
       The included dimensionality reduction is not suitable for 
       very small data sets. Therefore, tables or groups with few 
       rows (plots) or columns (species) are not partitioned. 
       For division, a group must consist of at least three 
       members and the partition must be supported by at least
       \code{stopcrit [1]} descriptors (species) reaching a 
       standardized G-value of \code{stopcrit [2]}.
                                           
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
       For large datasets, Isopam may need too much memory (consider using
       \link{isopam.2}) or too much computation time. The optimization 
       procedure (selection of Isomap dimensions and -\emph{k}, selection 
       of cluster numbers) is based on a brute force approach that takes its 
       time with large data sets. Low speed is inherent to the method, so 
       don't complain. If used with data not representing species in plots 
       make sure that the indicator approach is appropriate.
     }
     \author{
       Sebastian Schmidtlein with contributions from Jason Collison.        
       }	
     \seealso{
       \code{\link{isotab}}, \code{\link{isopam.2}}
     }
     \examples{
     ## load data to the current environment
     data(andechs)
     
     ## call isopam with the standard options
     ip<-isopam(andechs)
     
     ## examine cluster hierarchy
     plot(ip)
     
     ## examine frequency table (second 
     ## hierarchy level)
     isotab(ip, 2)
     }
     \keyword{ cluster }
