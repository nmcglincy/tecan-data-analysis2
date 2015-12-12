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