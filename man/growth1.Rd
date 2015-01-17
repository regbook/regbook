\encoding{utf8}
\name{growth1}
\docType{data}
\alias{growth1}

\title{
  growth1
}

\description{
  자료 \code{growth1}는 어느 식물의 성장 속도에 관한 실험 결과이다. 

}

\usage{growth1}

\format{이 데이터프레임은 다음과 같은 열로 구성되어 있다.
  \describe{
    \item{y}{성장치}
    \item{x}{시간}
    }
}

\examples{
fit <- nls(y ~ beta1 * x/ (beta2 + x), data=growth1, start=list(beta1=22, beta2=9), trace=TRUE)
summary(fit)

plot(y ~ x, growth1)
lines(fitted(fit) ~ x, growth1)

plot(resid(fit) ~ fitted(fit))
abline(h=0)

plot(resid(fit) ~ x, growth1)
abline(h=0)
}

\references{

}

\source{
}

\keyword{datasets}
