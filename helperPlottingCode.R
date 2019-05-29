
cut <- c(ageD[11], ageD[16], ageD[21])
colours <- c("blue3", "orange2", "gray83")
plot(ylim = c(0, 0.00055), xlim = c(0,12))
l <- 0

for (i in cut) {
  l <- l+1
  
  m <- sapply(levels(factor(i$Age)), function(x) {
    mean(i$Frequency[i$Age==x])
  })
  
  df <- data.frame(Age = as.numeric(names(m)), Frequency = unname(m))
  mod <- nls(Frequency ~ a*Age^b, data = df, start = list(a=1,b=1), control = list(maxiter=1000) )
  points(df,  cex = 1.5, pch=20, col=colours[[l]], add=T)
  a <- mod$m$getPars()[[1]]
  b <- mod$m$getPars()[[2]]
  curve(a*x^b, add = T, lwd = 2, col=colours[[l]])

}