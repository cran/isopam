\encoding{UTF-8}
\name{isotab}
\Rdversion{1.1}
\alias{isotab}
\alias{print.isotab}

\title{
Fidelity and frequency of species in clusters
}
\description{
Calculates the fidelity of species to clusters. Returns equalized
phi coefficients of association, an ordered frequency table and 
Fisher`s exact test for the probability of obtaining the observed 
frequencies. Isopam objects as well as other combinations of tables 
and cluster vectors are accepted as input data. An associated plotting 
method visualises how closely individual species are associated with 
clusters.   
}
\usage{
isotab(x, level = NULL, clusters = NULL, phi.min = "isotab", p.max = .05)
\method{print}{isotab}(x, n = NA, \dots)
}
\arguments{
  \item{x}{
   Object either of class \code{isopam} or a dataframe or matrix with rownames
   (plot names) and column names (species names) that is accompanied by a 
   cluster vector (\code{clusters}) with named elements corresponding to the
   sites in \code{x}. Tibbles need a column with plot names (<chr>), while the 
   other columns are of class <dbl> or <int>. In method \code{print}, \code{x} 
   is an object of class \code{isotab}.}
  \item{clusters}{Vector with assignments of clusters to plots, only
   needed if \code{x} is not an isopam object. The names of the elements need
   to be identical to the rownames of \code{x}.}
  \item{level}{
   Level in cluster hierarchy starting with 1 = first division.}
  \item{phi.min}{
   Threshold of equalized \emph{phi} determining which species are listed 
   in the upper part of the table. Applies only to species passing the 
   criterion defined by \code{p.max}. If \code{phi.min = "isopam"} (the 
   default) isotab suggests a value based on the numbers of observations.}  
  \item{p.max}{
   Threshold of Fisher's \emph{p} determining which species are listed in the 
   upper part of the table. Applies only to species passing the criterion 
   defined by \code{phi.min}.}
  \item{n}{number of lines used by \code{print}. If \code{NA} (the default), 
   n is oriented on the number of diagnosic species. Use \code{n = Inf} to print
   all rows.}
  \item{...}{other arguments used by \code{print}.}
}  
\details{
   \code{phi.min} is based on the 'equalized \emph{phi}' value according to
   \enc{Tichý}{Tichy} & \enc{Chitrý}{Chitry} 2006. The threshold proposed if 
   \code{phi.min} is set to "isotab" should be adjusted to local conditions. 
   The significance (Fisher's \emph{p}) refers to the probability that the 
   observed frequency is reached. The test is two-tailed, which means that 
   exceptionally low frequencies can result as highly significant as well as 
   exceptionally high frequencies. This allows positive and negative 
   characterisation of a cluster by species.   
}
\value{
   \item{call}{generating call}
   \item{depth}{Number of levels in the cluster hierarchy from the original
     clustering procedure.}
   \item{level}{Level chosen for isotab.}
   \item{tab}{Ordered species by cluster table with frequencies and their 
     significance. The latter is derived from Fisher's exact test (see 
     \code{fisher_p} and details, \emph{p} <= 0.05: *, \emph{p} <= 0.01: **, 
     \emph{p} <= 0.001: ***).}   
   \item{phi}{Dataframe with equalized phi values (see details).} 
   \item{fisher_p}{Numerical results from Fisher's exact test (see details)}
   \item{n}{Matrix with cluster sizes.}
   \item{thresholds}{\code{phi.min} and \code{p.max} used for table sorting.} 
   \item{typical}{Text with items (often species) typically 
     found in clusters (according to \code{thresholds}).}
   \item{typical_vector}{\code{typical} as a single character vector.}
   \item{sorted_table}{Ordered species by plot table.}
}
\references{
   \enc{Tichý}{Tichy}, L., \enc{Chytrý}{Chitry}, M. (2006): Statistical 
   determination of diagnostic species for site groups of unequal size.
   \emph{Journal of Vegetation Science} \bold{17}: 809-–818.

   Schmidtlein, S., \enc{Tichý}{Tichy}, L., Feilhauer, H., Faude, U.
   (2010): A brute force approach to vegetation classification.
   \emph{Journal of Vegetation Science} \bold{21}: 1162--1171.
}
\author{
   Sebastian Schmidtlein
}

\seealso{
   \code{\link{isopam}}, \code{\link{plot.isotab}}
}
\examples{
   ## load data to the current environment
   data(andechs)
     
   ## call isopam with the standard options
   ip<-isopam(andechs)
    
   ## build table
   it <- isotab(ip)
   it

   ## change phi threshold
   it <- isotab(ip, phi.min = 0.8)

   ## switch cluster level
   it <- isotab(ip, level = 1)
   it
}
