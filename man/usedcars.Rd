\encoding{utf8}
\name{usedcars}
\docType{data}
\alias{usedcars}

\title{usedcars}

\description{ 중고차 가격에는 어떠한 변수들이 영향을 주는 걸까? 국내의
  중고차 사이트에는 기본적으로 연식, 차종, 색, 도어, 엔진, 배기량, 연료,
  주행거리 등의 정보가 올라와 있다. 자료 \code{usedcars}는 국내 유명
  중고차 사이트에서 2007년에 수집한 것으로 가솔린 엔진의 중형차를
  대상으로 각 중고차에 대한 가격, 연식, 주행거리, 배기량, 변속기종류를
  조사하였다. }

\usage{usedcars}

\format{이 자료는 다음과 같은 열로 구성되어 있다.
  \describe{
    \item{price}{가격(만원)}
    \item{year}{연식(개월)은 현재까지 몇 개월이 지났는지를
      계산하였다.}
    \item{mileage}{주행거리(km)}
    \item{cc}{배기량(cc)}
    \item{automatic}{변속기종류. 수동(0), 자동(1).}
    }
}


\examples{
pairs(usedcars)

fit <- lm(price ~ year + mileage + cc + automatic, usedcars)
anova(fit)
summary(fit)

predict(fit, se.fit=TRUE)
residuals(fit)
rstandard(fit)
cooks.distance(fit)
}

\references{
}
\source{
}