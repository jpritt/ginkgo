#!/usr/bin/env Rscript

s = paste("Hello", "world", sep=" ")
saveRDS(s, file=paste("/local1/work/ginkgodev", "/out.txt", sep=""))

w = readRDS(paste("/local1/work/ginkgodev", "/out.txt", sep=""))
print(w)
