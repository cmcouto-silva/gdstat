myseq <- 1:100

vec <- logical(length = 786)

for(i in c(1:786)) {
  from <- seq(1, (length(myseq) - window) + 1, by = i)
  to <- from + (window - 1)
  vec[i] <- to[length(to)] > myseq[length(myseq)]
}

any(vec)

from <- seq(1, (length(myseq) - window) + 1, by = 5)
to <- from + (window - 1)
wid <- seq_along(from)

dt <- data.table(from = from, to = to, wid = wid)

10

window <- 20
slide <- 7


wid[12] * slide - (slide - 1) 

(slide - 1) / slide


f <- wid[5] * slide - (slide - 1) 
t <- f + (window - 1)
f;t




length(to)

length(21:40)
