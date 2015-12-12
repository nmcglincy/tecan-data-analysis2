spline.fit = function(dat,
					  sorting.var1,
					  sorting.var2){
    suppressMessages(require(plyr))
	lgdf = dlply(.data = dat,
				 .variables = c(sorting.var1, sorting.var2))
	growth.curve.lite2 = function(foo) {
				   suppressMessages(require(grofit))
				   spl.fit = smooth.spline(foo$Time.hrs, foo$norm.val, keep.data = FALSE)
				   data.frame(Time.hrs = spl.fit$x, abs.fit = spl.fit$y)
				 }
	lspl = lapply(lgdf, growth.curve.lite2)
	spl.df = ldply(lspl)
	# names(spl.df)[2] = "Time.hrs"
	df.fin = merge(ldply(lgdf),
	      		   spl.df,
	      		   by = c(".id", "Time.hrs"))
	return(df.fin)
}

spline.fit.slim = function(dat,
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
