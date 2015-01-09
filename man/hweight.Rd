\encoding{utf8}
\name{hweight}
\docType{data}
\alias{hweight}

\title{hweight}

\description{ 기술표준원에서는 우리 국민이 사용하는 의류, 가구, 주택,
  자동차 등 모든 산업제품과 생활공간 설계에 필요한 인체치수자료를
  확보하기 위하여 2003년부터 2년간 한국인의 인체치수 및 인체형상을 측정
  조사하는 제5차 한국인 인체치수조사사업(Size Korea)을 수행하였다. 자료
  \code{hweight}는 이 중에서 표본으로 추출된 20--24세의 남녀들의 키와
  몸무게에 대한 것이다. }

\usage{hweight}

\format{이 데이터프레임은 다음과 같은 열로 구성되어 있다.
  \describe{
    \item{weight}{몸무게(kg)}
    \item{height}{키(cm)}
    \item{gender}{성별. 남성(M), 여성(F)}
    }
}

\examples{
x <- with(hweight, height[gender=="M"])
hist(x, probability=TRUE)
lines(density(x))
}

\references{
  ...
}

\source{ ... }

\keyword{datasets}