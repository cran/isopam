\encoding{UTF-8}
\name{isotab}
\Rdversion{1.1}
\alias{isotab}
\title{
Ordered frequency table for Isopam clusters
}
\description{
Computes an ordered frequency table based on Isopam clustering results. 
The upper part of the table lists typical descriptors (usually species) with
a significant binding to single clusters (according to customisable thresholds). 
The lower part of the table is ordered by descending overall frequency.  
}
\usage{
isotab(ip, level = 1, phi.min = 'auto', p.max = .05)
}
\arguments{
  \item{ip}{
   object of class \code{isopam}.
}
  \item{level}{
   level in cluster hierarchy starting with 1 = first division.
}
  \item{phi.min}{
   threshold of \emph{phi} determining which descriptors (species) 
   are listed in the upper part of the table. Applies only to 
   descriptors passing the criterion defined by \code{p.max}. If
   \code{phi.min = 'auto'} (the default) isotab suggests a 
   suitable value based on the numbers of clusters, observations, 
   and descriptors.
}  
  \item{p.max}{
   threshold of Fisher's \emph{p} determining which descriptors 
   (species) are listed in the upper part of the table. Applies 
   only to descriptors passing the criterion defined by 
   \code{phi.min}.
}
}
\details{
   \code{phi.min} is based on the standardized \emph{phi} value according to 
   \enc{Chitrý}{Chitry} et al. 2002.  
}
\value{
   \item{tab}{dataframe with ordered frequencies and their significance. 
   The latter is derived from Fisher's exact test (\emph{p} <= 0.05: *, 
   \emph{p} <= 0.01: **, \emph{p} <= 0.001: ***).}
   \item{n}{matrix with cluster sizes.}
   \item{thresholds}{\code{phi.min} and \code{p.max} used.} 
   \item{typical}{dataframe with items (often species) typically found in 
   clusters (according to \code{thresholds}).}
}
\references{
   \enc{Chitrý}{Chitry}, M., \enc{Tichý}{Tichy}, L., Holt, J., 
   \enc{Botta-Dukát}{Botta-Dukat}, Z. (2002): Determination of diagnostic 
   species with statistical fidelity measures. \emph{Journal of Vegetation 
   Science} \bold{13}, 79--90.

   Schmidtlein, S., \enc{Tichý}{Tichy}, L., Feilhauer, H., Faude, U.
   (2010): A brute force approach to vegetation classification.
   Journal of Vegetation Science (in press).
}
\author{
   Sebastian Schmidtlein
}

\seealso{
   \code{\link{isopam}}
}
\examples{
   ## load data to the current environment
   data(andechs)
     
   ## call isopam with the standard options
   ip<-isopam(andechs)
    
   ## build table (uppermost hierarchy level)
   isotab(ip)

   ## build table (lower hierarchy level)
   isotab(ip,2)
}