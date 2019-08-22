program_on_path <- function(program) {
  logi <- ifelse(Sys.which(program) == "", FALSE, TRUE)
  if(!logi) stop(program, " not installed on system path!")
}
