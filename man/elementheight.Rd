\encoding{utf8}
\name{elementheight}
\docType{data}
\alias{elementheight}

\title{
  elementheight
}

\description{ 자료 \code{elementheight}에는 지식경계부 산하
기술표준원에서 실시한 제5차 한국인 인체치수조사사업(2003--2004)에서
발표한 초등학교에 해당한는 6세에서 12세까지의 나이별 키 자료가 나와있다.
각 나이별로 250명 정도의 표본에서 구한 키의 표본평균과 표준편차가 주어져
있는데 나이에 따라 키와 함께 표준편차도 증가하고 있음을 알 수 있다.
}

\usage{elementheight}

\format{이 데이터프레임은 다음과 같은 열로 구성되어 있다.
  \describe{
    \item{meanheight}{키의 표본평균}
    \item{sdheight}{키의 표본표준편차}
    \item{age}{나이}
    }
}

\examples{
fit <- lm(meanheight ~ age, elementheight, weight=1/sdheight^2)
anova(fit)
summary(fit)
}

\references{
}

\source{
}

\keyword{datasets}
