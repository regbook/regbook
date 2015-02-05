\encoding{utf8}
\name{galton}
\docType{data}
\alias{galton}

\title{galton}

\description{
골턴의 키 자료 \code{galton}.
}

\usage{galton}

\format{이 데이터프레임은 다음과 같은 열로 구성되어 있다.
  \describe{
    \item{sheight}{아들의 키(inch)}
    \item{fheight}{아버지의 키(inch)}
    }
}

\examples{
  with(galton, {
    son <- sheight * 2.54
    father <- fheight * 2.54
    plot(son ~ father)
  })
}

\references{
  Weisberg. 1985.
}

\source{
}

\keyword{datasets}
