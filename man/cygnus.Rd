\encoding{utf8}
\name{cygnus}
\docType{data}
\alias{cygnus}

\title{
  cygnus
}

\description{
  자료 \code{cygnus}는 백조자리 성운 방향에 있는 47개 별들의 밝기와
  표면온도를 조사한 것이다.

}

\usage{cygnus}

\format{이 데이터프레임은 다음과 같은 열로 구성되어 있다.
  \describe{
    \item{light}{빛의 강도. 로그 변환된 값이다.}
    \item{temp}{표면 온도. 로그 변환된 값이다.}
    }
}

\examples{
plot(light ~ temp, cygnus)
}

\references{
  Rousseuw & Leroy, 1987
}

\source{
}

\keyword{datasets}
