\encoding{UTF-8}
\name{isotab}
\Rdversion{1.1}
\alias{isotab}
\title{
Frequency table for Isopam clusters
}
\description{
Computes an ordered frequency table based on Isopam clustering results. 
The upper part of the table lists descriptors (usually species) with
a significant binding to clusters. The lower part is ordered by descending 
overall frequency. The threshold between these two parts can be shifted. 
}
\usage{
isotab(ip, level = 1, phi.min = .5, p.max = .05)
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
   descriptors passing the criterion defined by \code{p.max}.
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
   \item{n}{matrix with cluster sizes.}
   \item{tab}{Data frame with ordered frequencies and their significance. 
   The latter is derived from Fisher's exact test (\emph{p} <= 0.05: *, 
   \emph{p} <= 0.01: **, \emph{p} <= 0.001: ***).}
}
\references{
   \enc{Chitrý}{Chitry}, M., \enc{Tichý}{Tichy}, L., Holt, J., 
   \enc{Botta-Dukát}{Botta-Dukat}, Z. (2002): Determination of diagnostic 
   species with statistical fidelity measures. \emph{Journal of Vegetation 
   Science} \bold{13}, 79--90.

   Schmidtlein, S., \enc{Tichý}{Tichy}, L., Feilhauer, H., Faude, U.
   (2010): A brute force approach to vegetation classification.
   in review.
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