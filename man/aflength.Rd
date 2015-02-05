\encoding{utf8}
\name{aflength}
\docType{data}
\alias{aflength}

\title{aflength}

\description{ 의학적으로 정확하지는 않지만 인체의 특성 중에는 양팔을
  벌린 길이는 키와 같다든지, 목둘레는 허리둘레 또는 팔목굵기와 관계가
  있다든지 하는 흥미로운 사실들이 많이 있다. 자료 \code{aflength}는
  S대학교 학생들 중 무작위로 남녀 각각 16명을 추출하여 조사한 자료로
  발길이와 팔안쪽길이이다. }

\usage{aflength}

\format{이 데이터프레임은 다음과 같은 열로 구성되어 있다.
  \describe{
    \item{foot}{발길이(cm). 양말을 벗은 상태로 측정하였고 오른쪽 발만 측정하였다.}
    \item{forearm}{팔안쪽길이(cm). 손목부터 팔꿈치가 접히는 부분까지의
  길이이다. 오른쪽 팔만 측정하였다.}
    \item{gender}{성별. 남성(M), 여성(F).}
    }
}

\examples{
  plot(foot ~ forearm, aflength)
  aflength.lm <- lm(foot ~ forearm, aflength)
  abline(aflength.lm)
  with(aflength, lines(lowess(forearm, foot), col="red"))
}

\references{
}

\source{
}

\keyword{datasets}
