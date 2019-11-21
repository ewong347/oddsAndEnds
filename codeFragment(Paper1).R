v1 <- vl[order(vl$ID),]
v2 <- vm[order(vm$ID),]


Times <- as.numeric(names(table(subG$v$Time)))

EdgeFrequency <- sapply(Times, function(t) {
  oldE <- subset(subG$e, t1<=t | t2<=t)
  oldV <- subset(subG$v, Time==t)
  return(nrow(oldE)/nrow(oldV) )
})
plot(Times,EdgeFrequency)


g1 <- readRDS("~/Data/Tennessee/analysis_PRO/tn93TnsubB_nomet_G.rds")
g2 <- readRDS("~/Data/Tennessee/analysis_PRO/tn93TnsubB_met_G.rds")

g1$f <- likData(g1)
g2$f <- likData(g2)

saveRDS(g1, "~/Data/Tennessee/analysis_2cv/tn93TnsubB_nomet_G.rds")
saveRDS(g2, "~/Data/Tennessee/analysis_2cv/tn93TnsubB_met_G.rds")