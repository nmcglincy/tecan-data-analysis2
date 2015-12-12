filt.and.norm2  = function(dat, 
                          # bad.wells, 
                          sorting.var1,
                          sorting.var2,
                          var.to.norm) {
  require(plyr)
  require(dplyr)
  
  write.csv(x = dat %>%
              # filter(well %in% setdiff(dat$well, bad.wells)) %>%
              filter(strain == "YEPD") %>%
              group_by_(sorting.var1, sorting.var2) %>%
              summarise(med = median(od600)),
            file = "yepd-data.csv",
            row.names = FALSE)

  good.norm.dat = ldply(Map(function(x, y){mutate(x, norm.val = x[,var.to.norm] - y)},
  							dlply(.data = dat %>% 
										    # filter(well %in% setdiff(dat$well, bad.wells)) %>%
										    filter(strain != "YEPD"), 
								  .variables = c(sorting.var1, sorting.var2)),
							(dat %>%
				                 # filter(well %in% setdiff(dat$well, bad.wells)) %>%
				                 filter(strain == "YEPD") %>%
				                 group_by_(sorting.var1, sorting.var2) %>%
				                 summarise(med = median(od600)) %>%
				                 select(med))$med),
						.id = NULL)
						
  write.csv(x = good.norm.dat,
            file = "filt-norm-tecan-data.csv",
            row.names = FALSE)
  return(good.norm.dat)
}