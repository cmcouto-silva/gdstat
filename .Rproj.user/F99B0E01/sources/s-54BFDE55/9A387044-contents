distribplot <- function(ihs.data, file.name, file.type = 'png', 
                        lty = 1, 
                        lwd = 1.5, 
                        col = c("blue", "red"),
                        main = "Genome-wide distribution",
                        xlab = "iHS value",
                        cex.main = 1.5,
                        cex.lab = 1.25,
                        qqplot = TRUE) {
  
  substrRight <- function(x, n) {
    substr(x, nchar(x)-n+1, nchar(x))
  }
  
  filename <- ifelse(test = substrRight(file.name, 4L) == paste0(".", file.type), 
                     yes = paste0(file.name, ".ihsplot"), no = paste0(file.name, ".ihsplot.", file.type))
  
  Cairo::Cairo(file = filename, type = file.type, width = 1200, height = 600, units = "px", pointsize = 14, dpi = 72)
  par(mar = c(5, 5, 4, 2) + 0.1)
  
  plot(density(ihs.data, na.rm = TRUE), main = main, xlab = xlab, col = col[1], 
       lty = lty, lwd = lwd, cex.main = cex.main, cex.lab = cex.lab)
  
  curve(dnorm, col = col[2], add = TRUE, bty='n')
  
  legend("topright", c("Observed", "Gaussian"), bty = "n", 
         col = col, lty = lty, lwd = lwd)
  
  dev.off()
  
  if(qqplot) {
    
    par(mar = c(5, 5, 4, 2) + 0.1)
    filename <- ifelse(test = substrRight(file.name, 4L) == paste0(".", file.type), 
                       yes = paste0(file.name, ".ihsplot"), no = paste0(file.name, ".ihsqqplot.", file.type))
    
    Cairo::Cairo(file = filename, type = file.type, width = 1200, height = 600, units = "px", pointsize = 14, dpi = 72)
    
    qqnorm(ihs.data[!is.na(ihs.data)], cex.main = cex.main, cex.lab = cex.lab, pch = 16, cex = 0.75)
    abline(a = 0, b = 1, lty = 2)
    dev.off()
  }
  
}
