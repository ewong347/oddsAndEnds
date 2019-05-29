## args <- c("Tennessee/robustness/tnBadRunGD", "Tennessee/tnGraph.rds")

args <- commandArgs(trailingOnly = T)

runs <- lapply(list.files(args[1]), function(x) {readRDS(file=paste0(args[1], "/", x))})
g <- readRDS(file=args[2])

gaics <- lapply(runs, function(run){
  sapply(run, function(x) {x$gaic})
})

vars <- lapply(runs, function(run){
  sapply(run, function(x) {var(x$csize)})
})

norm <- lapply(runs, function(run){
  sapply(run, function(x) {x$gaic/log(var(x$csize))})
})

mins <- sapply(gaics, function(x){which(x==min(x))}) 
plot(as.numeric(names(mins)), main = "Robustness Test3, Tennessee", xlab = "Run (runX shows the minimum GAIC obtained after removing x-1 years)", ylab = "Cutoff obtaining minimum GAIC")

mins <- sapply(norm, function(x){which(x==min(x))}) 
plot(as.numeric(names(mins)), main = "Robustness Test4, Tennessee", xlab = "Run (runX shows the minimum GAIC obtained after removing x-1 years)", ylab = "Cutoff obtaining minimum GAIC normalized to log variance of cluster size")


for (i in 1:length(runs)){
  plot(gaics[[i]]/log(vars[[i]]))
}

