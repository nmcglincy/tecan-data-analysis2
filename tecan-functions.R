# data.file is character string of the tecan output .csv formatted for number in excel
# record whole plate, but annotate empty wells as "NA" in sample-info.csv
# sample-info.csv must have column named "strain", and be produced by the function.

# Generates a template for sample.info

generate.sample.info.template = function(){
  plate.template = data.frame(plate.letters = rep(LETTERS[1:8], each = 12),
                              plate.numbers = rep(1:12, 8),
                              well = paste(rep(LETTERS[1:8], each = 12), 
					rep(1:12, 8), 
					sep = ""))
  write.csv(plate.template, 
	  file = "plate-template.csv", 
	  row.names = FALSE)
  return(plate.template)
}

# Initial manipulation and examining temperature trend

init.proc = function(tecanOutput, 
		sorting.var1){
	suppressMessages(require(readr))
	suppressMessages(require(plyr))
	suppressMessages(require(dplyr))
	suppressMessages(require(tidyr))
	suppressMessages(require(ggplot2))
	suppressMessages(require(ggthemes))

	data = read_csv(file = tecanOutput,
	              skip = 35,
	              n_max = 96,
	              col_names = TRUE) %>%
	    na.omit()
	names(data)[2] = "Temp.C"
	names(data) = make.names(names = names(data), unique = TRUE, allow_ = TRUE)
	data$Cycle.Nr. = as.numeric(data$Cycle.Nr.)

	data = data %>%
		select(c(1:2, 4, seq(from = 3, to = length(data), by = 2))) %>%
		gather("well", "od600", A1:H12) %>%
		mutate(Time.hrs = Time..ms./3600000)
	
	sample.info = read_csv(file = "sample-info.csv", col_names = TRUE)
	# would dplyr::full_join() work here?
	data = join(data, sample.info, by = "well")

	ggplot(data,
	     aes_string(x = "Time.hrs",
	                y = "od600",
	                colour = sorting.var1)) +
	geom_line(size = 1) +
	facet_wrap( ~ well, ncol = 12) +
	xlab("Time, hrs") +
	ylab("OD, 600 nm") +
	scale_x_continuous(breaks = c(0, 10)) +
	scale_colour_colorblind() +
	guides(colour = guide_legend(override.aes = list(size=10))) +
	theme(panel.border = element_rect(fill = NA, colour = "black"),
	      axis.title.x = element_text(vjust = 0, size = 16),
	      axis.title.y = element_text(vjust = 1, size = 16),
	      axis.text.x = element_text(size=16),
	      axis.text.y  = element_text(size=16),
	      plot.title = element_text(size = 20),
	      legend.text = element_text(size = 16),
	      legend.title = element_text(size = 18),
	      strip.text.x = element_text(size = 16),
	      strip.text.y = element_text(size = 16))
	ggsave("by-well.png") 

	# Analysis of temperature drift
	cycle.data  = data %>% 
		select(Cycle.Nr., Temp.C) %>%
		distinct(Cycle.Nr.)

	png(file = "cycle-temp-plot.png", width = 7, height = 7, units = "in", res = 300)
	plot(x = cycle.data$Cycle.Nr.,
	   y = cycle.data$Temp.C,
	   type = "b",
	   ylim = c(29, 31),
	   xlim = c(0, 96),
	   xlab = "Cycle",
	   ylab = "Temperature, C")
	with(cycle.data, abline(lm(Temp.C ~ Cycle.Nr.), col = "blue"))
	with(cycle.data, abline(line(Cycle.Nr., Temp.C), col = "red"))
	dev.off()

	cycle.temp.stats = capture.output(with(cycle.data, summary(lm(Temp.C ~ Cycle.Nr.))))
	cat(cycle.temp.stats, file = "cycle-temp-stats.txt", sep = "\n", append = TRUE)
	cat("Spearman's correlation coefficient", file = "cycle-temp-stats.txt", sep = "\n", append = TRUE)
	cycle.temp.stats = capture.output(cor(cycle.data$Cycle.Nr. , cycle.data$Temp.C, method = "spearman"))
	cat(cycle.temp.stats, file = "cycle-temp-stats.txt", sep = "\n", append = TRUE)
	cat("Pearson's correlation coefficient", file = "cycle-temp-stats.txt", sep = "\n", append = TRUE)
	cycle.temp.stats = capture.output(cor(cycle.data$Cycle.Nr., cycle.data$Temp.C, method = "pearson"))
	cat(cycle.temp.stats, file = "cycle-temp-stats.txt", sep = "\n", append = TRUE)

	return(data)
}
	
# Filter and correct
# Define object bad.wells

filt.and.norm  = function(dat, 
                          bad.wells,
                          sorting.var1,
                          sorting.var2,
                          var.to.norm){
  suppressMessages(require(plyr))
  suppressMessages(require(dplyr))
  
  write.csv(x = dat %>%
              filter(well %in% setdiff(dat$well, bad.wells)) %>%
              filter(strain == "YEPD") %>%
              group_by_(sorting.var1, sorting.var2) %>%
              summarise(med = median(od600)),
            file = "yepd-data.csv",
            row.names = FALSE)

  good.norm.dat = ldply(Map(function(x, y){mutate(x, norm.val = x[,var.to.norm] - y)},
				dlply(.data = dat %>% 
						filter(well %in% setdiff(dat$well, bad.wells)) %>%
						filter(strain != "YEPD"), 
					.variables = c(sorting.var1, sorting.var2)),
				(dat %>%
					filter(well %in% setdiff(dat$well, bad.wells)) %>%
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

# Function to fit spline on pooled replicates per condition

spline.fit = function(dat,
		sorting.var1,
		sorting.var2){
    	suppressMessages(require(plyr))
	gcl2 = function(foo) {
		suppressMessages(require(grofit))
		spl.fit = smooth.spline(foo$Time.hrs, foo$norm.val, keep.data = FALSE)
		data.frame(Time.hrs = spl.fit$x, abs.fit = spl.fit$y)
	}
	lgdf = dlply(.data = dat,
			.variables = c(sorting.var1, sorting.var2))
	df.fin = merge(ldply(lgdf),
			ldply(lapply(lgdf, gcl2)),
			by = c(".id", "Time.hrs"))
	return(df.fin)
}

# Using package grofit to look at growth curve parameters

grofit.analysis = function(dat,
                           sorting.var1,
                           sorting.var2){
  gca = function(foo) {
    suppressMessages(require(grofit))
    # fit smooth spline
    spl.fit = smooth.spline(foo$Time.hrs, foo$norm.val)
    # calculate first derivative on smooth spline values
    firstDeriv = predict(spl.fit, sort(foo$Time.hrs), deriv = 1)
    # mu.time is the time at which the maximal growth rate happens
    # this is defined at the time at which the first derivative (the speed/ rate of change) is maximal
    mu.time = firstDeriv$x[which.max(firstDeriv$y)]
    # mu: the maximal growth rate in abs units/hour
    mu = firstDeriv$y[which.max(firstDeriv$y)]
    # mu.abs is the od at which the the growth rate is maximal
    mu.abs = spl.fit$y[which.max(firstDeriv$y)]
    # the norm.growth.rate is the maximal rate normalised for the od at which it is happening
    norm.growth.rate = mu/mu.abs
    # lambda: the length of the lag-phase in hrs the x-axis intercept of the 
    # straight line from the log-phase. This seems to be a bit conservative (there is 
    # a clear inflection of the curve lambda - an "acceleration phase"), but is more 
    # stable and better reflects my visual impression of the curves.
    lambda = -(spl.fit$y[which.max(firstDeriv$y)] - (mu * mu.time))/mu
    # A: Absorbance of maximal growth
    A = max(spl.fit$y)
    # AUC: Area under the curve - measure of total growth
    AUC = grofit::low.integrate(spl.fit$x, spl.fit$y)
    # 
    system("mkdir gca-plots")
    jpeg(file = paste("gca-plots/", paste(foo$strain, foo$well, sep = "-"), ".jpeg", sep = ""), 
         width = 7, height = 7, units = "in", res = 300)
    plot(foo$Time.hrs, 
         foo$norm.val, 
         ylim = c(0,1.2),
         xlab = "Time, hrs", 
         ylab = "Absorbance, A.U.")
    lines(spl.fit, 
          col = "red")
    lines(firstDeriv$x, 
          firstDeriv$y/max(firstDeriv$y), 
          col = "blue")
    points(firstDeriv$x, 
           firstDeriv$y/max(firstDeriv$y), 
           col = "blue")
    abline(v = mu.time, 
           col = "darkgreen", 
           lty = 2)
    abline(h = A, 
           col = "darkgreen", 
           lty = 2)
    abline(a = spl.fit$y[which.max(firstDeriv$y)] - (mu * mu.time) , 
           b = mu)
    abline(v = lambda, 
           col = "orange", 
           lty = 2)
    dev.off()
    list(mu =  mu, lambda = lambda, A = A, AUC = AUC, norm.growth.rate = norm.growth.rate)
  }
  suppressMessages(require(plyr))
  suppressMessages(require(dplyr))
  suppressMessages(require(tidyr))
  
  lgca = lapply(dlply(.data = dat,
                      .variables = c(sorting.var1, sorting.var2, "well")),
                gca)
  
  sample.info = read.csv(file = "sample-info.csv",
                         header = TRUE,
                         stringsAsFactors = FALSE)
  
  params = merge(ldply(dlply(.data = sample.info, 
                             .variables = c(sorting.var1, sorting.var2, "well"))),
                 ldply(lapply(lgca, ldply, .id = "param")) %>% spread(param, V1),
                 by = ".id")
  write.csv(x = params,
            file = "grofit-params.csv",
            row.names = FALSE)
  
  params.summ = params %>%
    group_by(strain, chx_ugml) %>%
    summarise(N = length(mu),
              m_mu = mean(mu),
              se_mu = sd(mu)/sqrt(N),
              m_lam = mean(lambda),
              se_lam = sd(lambda)/sqrt(N),
              m_A = mean(A),
              se_A = sd(A)/sqrt(N),
              m_AUC = mean(AUC),
              se_AUC = sd(AUC)/sqrt(N),
              m_ngr = mean(norm.growth.rate),
              se_ngr = sd(norm.growth.rate)/sqrt(N))
  write.csv(x = params.summ,
            file = "grofit-params-summ.csv",
            row.names = FALSE)
  
  return(list(params = params, params.summ = params.summ))
}
