\encoding{UTF-8}
\name{plot.isotab}
\Rdversion{1.1}
\alias{plot.isotab}
\title{
Plot species fidelities to clusters
}
\description{
Function to plot \code{\link{isotab}} results. Based on equalised phi 
values according to \enc{Tichý}{Tichy} & \enc{Chitrý}{Chitry} (2006),
the method visualises how closely how many species are associated with 
clusters.
}
\usage{
\method{plot}{isotab}(x, labels = FALSE, text.size = 15, title = NULL,
       phi.min = "isotab", p.max = "isotab", \dots)
}
\arguments{
  \item{x}{Object of class \code{isotab}.}
  \item{labels}{
   Logical. Whether the bars should be labeled with species names. You may 
   need to enlarge the figure height to accommodate these names (or decrease 
   \code{text.size}).} 
  \item{text.size}{
   Text size }
  \item{title}{
   Optional text string with title}
  \item{phi.min}{
   Threshold of equalized \emph{phi} determining which species are shown. 
   Applies only to species passing the criterion defined by \code{p.max}. 
   If \code{phi.min = "isotab"} (the default) the threshold used by isotab
   is applied. Use \code{phi.min = 0} to remove the filter.}  
  \item{p.max}{
   Threshold of Fisher's \emph{p} determining which species are shown. 
   Applies only to species passing the criterion defined by \code{phi.min}.
   Note that this value relates to frequencies rather than phi. If 
   \code{p.max = "isotab"} (the default) the threshold used by isotab
   is applied. Use \code{p.max = 1} to remove the filter.}
  \item{...}{Other arguments (ignored)}
}  
\details{
   The thresholds are explained in \code{\link{isotab}}.  
}
\value{
   Prints and returns (invisibly) an object of class \code{ggplot}.
}
\references{
   \enc{Tichý}{Tichy}, L., \enc{Chytrý}{Chitry}, M. (2006): Statistical 
   determination of diagnostic species for site groups of unequal size.
   \emph{Journal of Vegetation Science} \bold{17}: 809-–818.
}
\author{
   Sebastian Schmidtlein
}
\seealso{
   \code{\link{isopam}}, \code{\link{isotab}}
}
\examples{
   ## load data to the current environment
   data(andechs)
     
   ## call isopam with the default options
   ip <- isopam(andechs)
       
   ## calculate fidelities
   it <- isotab(ip)
   
   ## plotting
   plot(it)

   ## show species labels
   plot(it, labels = TRUE)

   ## show all species
   plot(it, phi.min = 0)
}
