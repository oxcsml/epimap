gig_lpdf <- function(x,p,a,b) {
  (p * 0.5 * log(a / b)
      - log(2 * besselK(sqrt(a * b), p))
      + (p - 1) * log(x)
      - (a * x + b / x) * 0.5
  )
}

x <- exp(seq(-10,log(10),length.out=200))
plot(x,exp(gig_lpdf(x,2,1,1)))
