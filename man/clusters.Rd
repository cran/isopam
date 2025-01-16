\encoding{UTF-8}
\name{clusters}
\Rdversion{1.1}
\alias{clusters}

\title{
Returns clusters of an Isopam clustering
}
\description{
Provides a vector or data frame with cluster affiliations determined using Isopam.  
}
\usage{
clusters(x, level = NULL, k = NULL, style = c("flat", "hierarchical"))
}
\arguments{
  \item{x}{
   Object of class \code{isopam}.}
  \item{level}{
   An integer scalar or vector with the desired cluster level(s). Level numbers start with 1 for the first division.}
  \item{k}{
   An integer scalar or vector with the desired number(s) of groups}  
  \item{style}{
   Whether the labels of the clusters are consecutive ('flat') or nested ('hierarchical', i.e. 1.1, 1.2 etc.).}
}
\value{
   Factor vector or data frame with the cluster assignments.
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
     
   ## call isopam with the standard options
   ip <-isopam(andechs)
    
   ## return clusters
   clusters(ip)

   ## clusters of level 2, with labels reflecting the hierarchy
   clusters(ip, 2, style = "hierarchical")

   ## cluster solution with 3 classes
   clusters(ip, k = 3)
}
