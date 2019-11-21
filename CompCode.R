saveRDS(res1, "res1.rds")
saveRDS(res2, "res2.rds")

res3 <- readRDS("res1.rds")
res2 <- readRDS("res2.rds")

gaics1 <- sapply(res1, function(i){i$gaic})
gaics2 <- sapply(res2, function(i){i$gaic})

growth1 <- sapply(res1, function(i){mean(i$g)})
growth2 <- sapply(res2, function(i){mean(i$growth)})

ageD1 <- lapply(res1, function(i){i$f})
ageD2 <- lapply(res2, function(i){i$ageD})

mod1 <- sapply(res1, function(i){summary(i$mod)$aic})
mod2 <- sapply(res2, function(i){summary(i$mod)$aic})