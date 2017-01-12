\name{replacementBootstrap}
\alias{replacementBootstrap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Entropy-Based Bootstrap For Dependent Data}
\description{This bootstrap scheme accounts for temporal dependence in a sample by maximizing the relative entropy of the bootstrapped sample with respect to the initial sample. The entropy itself is computed with the lz4 library.}
\usage{
replacementBootstrap(x,seed=NA,replacementFraction=1.0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{vector of sample values.  It can be an xts object.  NAs and Inf values are removed.}
  \item{seed}{random-number generator seed.  If specified, set.seed uses this value.}
  \item{replacementFraction}{by default, set to 1.  Other values are not implemented for the time being.}
}
\details{The algorithm accounts for serial dependence by maximizing the relative entropy of the bootstrapped samples with respect to that of the original distribution.  Entropy values are computed by zipping both sample vectors with lz4 library (see https://github.com/bwlewis/lz4)
}
\value{A vector of bootstrapped sample values.
}
\references{
Sani, Amir, Alessandro Lazaric, and Daniil Ryabko. "The replacement bootstrap for dependent data." 2015 IEEE International Symposium on Information Theory (ISIT). IEEE, 2015.
}
\author{Amir Sani and Damien Challet}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{boot}}, ~~~
}
\examples{
x=rnorm(100)
replacementBootstrap(x)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ manip }% use one of  RShowDoc("KEYWORDS")
\keyword{ ts }% __ONLY ONE__ keyword per line
