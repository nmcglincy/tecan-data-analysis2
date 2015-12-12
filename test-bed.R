setwd("~/Dropbox/tecan-repo/")
list.files()
source("tecan-functions.R")
dat = init.proc(tecanOutput = "chx-test-20151024.csv",
                sorting.var1 = "strain")

gdat = filt.and.norm(dat = dat,
                     bad.wells = c("A1", "B5", "B6", "C1", "C6"),
                     sorting.var1 = "strain",
                     sorting.var2 = "chx_ugml",
                     var.to.norm = "od600")

head(gdat)

source("growth-curve-analysis.R")

param.list = grofit.analysis(dat = gdat,
                             sorting.var1 = "strain",
                             sorting.var2 = "chx_ugml")



