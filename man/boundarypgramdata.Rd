\name{boundarypgramdata}

\alias{boundarypgramdata}

\title{
compute data about specific parallelograms in the boundary of a zonohedron
}

\description{
The boundary of a zonohedron is the union of parallelograms,
where some of them may be facets, and some may be tiles
in the standard tiling of more complex facets.
The edges of each parallelogram are given by a pair of distinct simplified generators.
If a zonohedron has \eqn{n} of these generators,
then there are \eqn{n(n{-}1)/2} such pairs.
For each pair of generators,
there are 2 parallelograms which are antipodal to each other.
This function computes data about one parallelogram from this antipodal pair.
The total number of parallelograms is \eqn{n(n{-}1)}.
}

\usage{
boundarypgramdata( x, gndpair, cube=FALSE )
}


\arguments{
\item{x}{
a \bold{zonohedron} object as returned by the constructor \code{zonohedron()}
}

\item{gndpair}{
an Mx2 integer matrix.
Each row of \code{gndpair} must contain a pair of points
in the ground set of the simplified matroid of the zonohedron \code{x}.
The 1st point must be less that the 2nd point.
\cr
\code{gndpair} can also be a numeric vector that can be converted
to such a matrix, by row.

}

\item{cube}{
if \code{TRUE}, then a point of the cube that maps to the center
of the given parallelogram is returned, see \bold{Value}.
}

}



\value{
\code{boundarypgramdata()} returns a \code{data.frame}
with M rows and these columns:

\item{gndpair}{the given \code{gndpair}}

\item{hyperplaneidx}{
the index of the hyperplane in the simplified matroid
of \code{x} that contains \code{gndpair}
}

\item{center}{the center of the parallelogram \emph{relative to the center of
the zonohedron}.
For the antipodal parallelogram, multiply this by -1.
}

\item{transitions}{
the number of transitions in \code{pcube} - a point in the n-cube
that maps to \code{center}.
This is a non-negative even integer.
}

And if \code{cube} is \code{TRUE}, then this column is added:

\item{pcube}{
a point in the n-cube that maps to \code{center}.
For the antipodal parallelogram, subtract this from 1;
the number of transitions is the same.
}

If a row of \code{gndpair} has an invalid pair,
the other columns are filled with \code{NA}s.

In case of global error, the function returns \code{NULL}.
}


\section{WARNING}{
In this version of the package,
when the parallelogram is a tile in a more complex facet,
the function may return incorrect results.
This will be fixed in a future version.
}




\seealso{
\code{\link{zonohedron}()}
}




% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{adaptation}% use one of  RShowDoc("KEYWORDS")

