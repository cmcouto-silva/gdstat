refIdx <- function(N, window_size, step) {
  from <- seq(1, N - window_size + 1, by = step)
  to <- from + (window_size - 1)
  data.table(from = from, to = to, wid = seq_along(to))
}
