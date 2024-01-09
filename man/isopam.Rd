     \encoding{UTF-8}
     \name{isopam}
     \Rdversion{1.1}
     \alias{isopam}
     \alias{plot.isopam}
     \alias{identify.isopam}
     \alias{summary.isopam}
     \alias{print.isopam}
     \title{Isopam (Clustering)}
     \description{
       Isopam classification is performed either as a
       hierarchical, divisive method or as non-hierarchical
       partitioning. Isopam is designed for matrices 
       representing species abundances in plots and with a 
       diagnostic species approach in mind. It optimises 
       clusters and cluster numbers for concentration of 
       indicative species in groups. Predefined indicative 
       species and cluster medoids can optionally be added 
       for a semi-supervised classification.
     }
     \usage{
     isopam(dat, c.fix = FALSE, c.max = 6, l.max = FALSE, stopat = c(1,7),
            sieve = TRUE, Gs = 3.5, ind = NULL, centers = NULL, 
            distance = 'bray', k.max = 100, d.max = 7, juice = FALSE, \dots)

     \method{identify}{isopam}(x, \dots)
     \method{plot}{isopam}(x, \dots)
     \method{summary}{isopam}(object, \dots)
     \method{print}{isopam}(x, \dots)
     }
     \arguments{
       \item{dat}{data matrix: each row corresponds to an
         object (typically a plot), each column corresponds
         to a descriptor (typically a species). All variables must
         be numeric. Missing values (NAs) are not allowed. At least
         3 rows (plots) are required.}
       \item{c.fix}{number of clusters (defaults to \code{FALSE}).
         If a number is given, non-hierarchical partitioning is
         performed, \code{c.max} is ignored and \code{l.max} is 
         set to one.}
       \item{c.max}{maximum number of clusters per partition.
         Applies to all splits.}
       \item{l.max}{maximum number of hierarchy levels. Defaults
         to \code{FALSE} (no maximum number). Note that divisions
         may stop well before this number is reached (see
         \code{stopat}). Use \code{l.max = 1} for non-hierarchical
         partitioning (or use \code{c.fix}).}
       \item{stopat}{vector with stopping rules for hierarchical
         clustering. Two values define if a partition should be
         retained in hierarchical clustering: the first determines
         how many indicator species must be present per cluster, 
         the second defines the standardized G-value that must be 
         reached by these indicators. \code{stopat} is not effective 
         at the first hierarchy level or in non-hierarchical 
         partitioning.}
       \item{sieve}{logical. If \code{TRUE} (the deafult), only
         species exceeding a threshold defined by \code{Gs} are 
         used in the search for a good clustering solution. Their 
         number is multiplied with their mean standardized G-value. 
         The product is used as optimality criterion. If \code{FALSE} 
         all species are used for optimization.}
       \item{Gs}{threshold (standardized G value) for species
         to be considered in the search for a good clustering solution. 
         Effective with \code{sieve = TRUE}.}
       \item{ind}{optional vector of column names from \code{dat}
         defining species used as indicators. This turns Isopam 
         in an expert system. Replaces the automated selection of 
         indicators with \code{sieve = TRUE} (\code{ind} overrules 
         \code{sieve}).}
       \item{centers}{optional vector with indices (numeric) or 
         names (character) of observations used as cluster cores 
         (supervised classification).}
	     \item{distance}{name of a dissimilarity index for the 
         distance matrix used as a starting point for Isomap.
         Any distance measure implemented in packages \pkg{vegan}
         (predefined or using a designdist equation) or 
         \pkg{proxy} can be used (see details).}
       \item{k.max}{maximum Isomap \emph{k}.}
       \item{d.max}{maximum number of Isomap dimensions.}
       \item{juice}{logical. If \code{TRUE} input files for Juice are
         generated.} 
       \item{...}{other arguments used by juice or passed to S3 
         functions \code{plot} and \code{identify} (see 
         \link[stats]{dendrogram} and \link[stats]{hclust}).}
       \item{x}{\code{isopam} result object in methods 
         \code{plot}, \code{print} and \code{identify}.}
       \item{object}{\code{isopam} result object in method 
         \code{summary}.}
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
         the partitioning steps. \code{Name}: name of the
         respective parent cluster (0 in case of the first partition);
         \code{Subgroups}: number of subgroups; \code{Isomap.dim}:
         Isomap dimensions used; \code{Isomap.k.min}: minimum
         possible Isomap \emph{k}; \code{Isomap.k}: Isomap
         \emph{k} used; \code{Isomap.k.max}: maximum possible
         Isomap \emph{k}; \code{Ind.N}: number of indicators
         reaching or exceeding \code{Gs}; \code{Ind.Gs}: the average
         standardized G value of these indicators; and
         \code{Global.Gs}: the average standardized G value of all
         descriptors (species).}
       \item{centers_usr}{Cluster centers suggested by user.}
       \item{ind_usr}{Indicators suggested by user.}
       \item{indicators}{Indicators used.}
       \item{dendro}{an object of class \code{hclust} representing
         the clustering (as used by \code{plot}). Not present
         with only one level of partitioning.}
       \item{dat}{data used}
     }
     \details{
       Isopam is described in Schmidtlein et al. (2010).
       It consists of dimensionality reduction (Isomap: Tenenbaum
       et al. 2000; \link[vegan]{isomap} in \pkg{vegan}) and
       partitioning of the resulting ordination space (PAM: Kaufman
       & Rousseeuw 1990; \link[cluster:pam]{pam} in \pkg{cluster}).
       The classification is performed either as a hierarchical, 
       divisive method, or as non-hierarchical partitioning. It has 
       the following features: partitions are optimized for the 
       occurrence of species with high fidelity to groups; it 
       optionally selects the number of clusters per division; 
       the shapes of groups in feature space are not restricted to 
       spherical or other regular geometric shapes (thanks to the 
       underlying Isomap algorithm); the distance measure 
       used for the initial distance matrix can be freely defined.
       
       In semi-supervised mode, clusters are build around the provided
       medoids. Pre-defined indicator species are not as constraining,
       even if preference is given to cluster solutions in which their
       fidelity is maximized. It depends on the data how much they affect 
       the result.
       
       The preset distance measure is Bray-Curtis (Odum 1950).
       Distance measures are passed to \link[vegan:vegdist]{vegdist}
       or to \link[vegan:designdist]{designdist} in \pkg{vegan}.
       If this does not work it is passed to \link[proxy:dist]{dist}
       in \pkg{proxy}. Measures available in vegan are listed in
       \link[vegan:vegdist]{vegdist}. Isopam does not accept distance 
       matrices as a replacement for the original data matrix because 
       it operates on individual descriptors (species).

       Isopam is slow with large data sets. It switches to a slow mode
       when an internally used lookup array for the results of the
       search for an optimal parameterisation (selection of Isomap 
       dimensions and -\emph{k}, optionally selection of cluster numbers) 
       does not fit into RAM.

       \code{plot} creates (and silently returns) an object of class 
       \code{dendrogram} and calls the S3 plot method for that class.
       \code{identify} works just like \code{identify.hclust}. 

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
       With very small datasets, the indicator based optimization may 
       fail. In such cases consider using \code{sieve = FALSE} instead 
       of the default method.
     }
     \author{
       Sebastian Schmidtlein with contributions from Jason Collison and
       Lubomir \enc{Tichý}{Tichy}
       }
     \seealso{
       \code{\link{isotab}} for a table of descriptor (species)
       frequencies in clusters and fidelity measures. There is a plot 
       method associated to \code{\link{isotab}} objects that visualizes 
       species fidelities to clusters.
     }
     \examples{
     ## load data to the current environment
     data(andechs)
     
     ## call isopam with the standard options
     ip <- isopam(andechs)

     ## print function
     ip
     
     ## examine cluster hierarchy
     plot(ip)

     ## retrieve cluster vectors
     clusters <- ip$flat
     clusters
     
     ## same but hierarchical style (available with cluster trees)
     hierarchy <- ip$hier 
     hierarchy

     ## frequency table
     it <- isotab(ip)
     it

     ## plot with species fidelities (equalized phi)
     plot(it)

     ## non-hierarchical partitioning with three clusters
     ip <- isopam(andechs, c.fix = 3)
     ip

     ## limiting the set of species used in cluster search
     ip <- isopam(andechs, ind = c("Car_pan", "Sch_fer"), c.fix = 2)
     ip

     ## supervised mode with fixed cluster medoids
     ip <- isopam(andechs, centers = c("p20", "p22"))
     ip

     }
     \keyword{ cluster }
