myseq <- 1:100

vec <- logical(length = 786)

for(i in c(1:786)) {
  from <- seq(1, (length(myseq) - window) + 1, by = i)
  to <- from + (window - 1)
  vec[i] <- to[length(to)] > myseq[length(myseq)]
}

any(vec)

from <- seq(1, (length(myseq) - window) + 1, by = 7)
to <- from + (window - 1)
from;to

exceds <- which(to > myseq[length(myseq)])
from <- from[-exceds]
to <- to[-exceds]
wid <- seq_along(from)

from[5]
to[5]


3 * slide - (slide - 1) 

dt <- data.table(from = from, to = to, wid = wid)

window <- 20
slide <- 7

f <- wid[5] * slide - (slide - 1) 
t <- f + (window - 1)
f;t




length(to)

length(21:40)
