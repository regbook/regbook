\encoding{utf8}
\name{steel}
\docType{data}
\alias{steel}

\title{
  steel
}

\description{
  기계장비용 강판을 생산하는 회사에서는 납품하는 강판의 가격이 그 강판을
  사용하는 기계의 가격으로부터 어떠한 영향을 받는지를 조사하기 위하여
  과거 16개월 동안의 자료를 수집하였다.
}

\usage{steel}

\format{이 데이터프레임은 다음과 같은 열로 구성되어 있다.
  \describe{
    \item{steelprice}{강판의 가격(단위: 십만원)}
    \item{machineprice}{기계의 가격(단위: 백만원)}
    }
}

\examples{
fit <- lm(steelprice ~ machineprice, steel)

par(mfrow=c(1,2))
plot(steelprice ~ machineprice, steel)
abline(fit)
plot(resid(fit))
abline(h=0)
}

\references{
}

\source{

}

\keyword{datasets}
