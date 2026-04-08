#source functions
#interval.eweek: used for controlling temporal aggregation

#Part 1: directly used in DSC_sector, the only part allowed to use parallel computing

#1. get counterfactual for treated sector (bootstrap, placebo test included)
getCf <- function(id, t.info, t.s, c.s, time.ind, space.ind, no.trap, no.mosq, lambda.ind = 1, bootstrap = FALSE, space.placebo = TRUE, time.placebo = TRUE, time.lag = NULL, lag.T0 = NULL, typeq = NULL, space.agg = NULL, interval.eweek = 2, beta = FALSE, time.weight.ind = FALSE, cores = 1){
  in.main <- isFALSE(bootstrap) && isTRUE(space.placebo) && isTRUE(time.placebo)
  # Prepare for calculate lambdas. Need to know 
  #1. when treatment starts 2. Targer value (GAI) 3. quantile function for each sector each time
  T0 <- t.info[which(t.info[space.ind] == id), time.ind]
  c.s["E2week"] <- floor((c.s[[time.ind]] - T0) / interval.eweek)
  t.s <- t.s %>% filter(!!sym(space.ind) == id)
  t.s["E2week"] <- floor((t.s[[time.ind]] - T0) / interval.eweek)
  c.s <- c.s %>% rename(c(no.trap = no.trap, no.mosq = no.mosq))
  t.s <- t.s %>% rename(c(no.trap = no.trap, no.mosq = no.mosq))
  if (is.null(typeq)) typeq <- 1
  c.s[["GAI"]] <- c.s$no.mosq/c.s$no.trap
  t.s[["GAI"]] <- t.s$no.mosq/t.s$no.trap
  
  # get time representation (treatment start week = 0) and a vector of control sector ID
  e2weeks <- unlist(unique(t.s$E2week)) |> sort() |> as.vector()
  control.s <- unique(c.s$Sector_ID) |> sort.sector()
  
  #get E2week vs e-week table (convert E2week to readable e-week), will be used in plotting line chart
  convertweek <-t.s %>% select("E2week","Eweek_tot","Eyear","Eweek") %>% distinct()
  convertweek["eyew"] <- apply(convertweek %>% select("Eyear","Eweek") , MARGIN = 1, paste, collapse = '-')
  convertweek <- convertweek %>% select("E2week","Eweek_tot","eyew")
  E2week_n <- convertweek %>% filter(E2week < 0) %>% group_by(E2week) %>% summarize(n_wks = n()) %>% ungroup()
  
  if (time.weight.ind){
    valid_wks = E2week_n %>% filter(n_wks == interval.eweek) %>% select(E2week) %>% unlist()
    valid_wks = c(valid_wks, e2weeks[e2weeks >= 0])
    c.s <- c.s[(c.s$E2week |> unlist()) %in% valid_wks,]
    t.s <- t.s[(t.s$E2week |> unlist()) %in% valid_wks,]
    e2weeks <- unname(valid_wks) |> sort() |> as.vector()
    convertweek <-t.s %>% select("E2week","Eweek_tot","Eyear","Eweek") %>% distinct()
    convertweek["eyew"] <- apply(convertweek %>% select("Eyear","Eweek") , MARGIN = 1, paste, collapse = '-')
    convertweek <- convertweek %>% select("E2week","Eweek_tot","eyew")
  }
  
  #outside bootstrap already use parallel computation. DO NOT SET ANY parallel calculation for bootstrap!
  if (bootstrap){
    c.s <-c.s %>% group_by(!!sym(space.ind), E2week) %>% slice_sample(prop = 1, replace = TRUE) %>% ungroup()
    #sampling Njt records (individuals) for bootstrapping WITH replacement
    quantFun <- lapply(e2weeks, getQuantile, 
                       prob_list = seq(0,1,length.out = 1001), 
                       control = c.s,
                       treatment = t.s,
                       time.ind = "E2week",
                       space.ind = space.ind,
                       val.ind = "GAI",
                       q.type = typeq
    )
    names(quantFun) <- e2weeks
    if (lambda.ind == 1){
      lambdas <- lapply(c(min(t.s$E2week):-1),getLambda,
                        control = c.s,
                        treatment = t.s,
                        time.ind = "E2week",
                        space.ind = space.ind,
                        val.ind = "GAI",
                        q.type = typeq,
                        beta = beta)
      f.lambdas <- unlist(lapply(control.s, function(x){
        return (lapply(c(1:length(lambdas)),function(y){return (lambdas[[y]][x])}) |> unlist() |> mean())
      }))
      names(f.lambdas) <- control.s
    } else {
      f.lambdas <- getLambda2(time = c(min(t.s$E2week):-1), control = c.s, treatment = t.s, time.ind = "E2week", space.ind = space.ind, val.ind = "GAI", q.type = typeq, beta = beta)
      names(f.lambdas) <- control.s
    }
    
    cf <- do.call(cbind,lapply(c(1:length(quantFun)),
                               function(x){return(as.matrix(cbind(quantFun[[x]][[2]])) %*% f.lambdas)}
    ))
    t_quant <- do.call(cbind,lapply(c(1:length(quantFun)), function(x){return(quantFun[[x]][[1]])}))
    diff <- t_quant - cf
    colnames(diff) <- e2weeks
    colnames(cf) <- e2weeks
    return (list(diff, f.lambdas, cf))
  }
  
  #--if not call by bootstrap
  if (isTRUE(in.main)) {print("Data preprocessed.")}
  quantFun <- pbmclapply(e2weeks, getQuantile, 
                         prob_list = seq(0,1,length.out = 1001), 
                         control = c.s,
                         treatment = t.s,
                         time.ind = "E2week",
                         space.ind = space.ind,
                         val.ind = "GAI",
                         q.type = typeq,
                         mc.cores = cores,
                         mc.silent = TRUE)
  names(quantFun) <- e2weeks
  if (isTRUE(in.main)) {print("get all quantile function")} #quantile function is used for later conterfactual calculation and difference check
  
  #lambda calculation
  if (lambda.ind == 1){
    if (!is.null(lag.T0)){ #For in-time placebo test
      lag.time = -lag.T0
      lambdas <- mclapply(c(min(t.s$E2week):lag.time),getLambda,
                          control = c.s,
                          treatment = t.s,
                          time.ind = "E2week",
                          space.ind = space.ind,
                          val.ind = "GAI",
                          q.type = typeq,
                          mc.cores = cores,
                          beta = beta)
    } else { #For non in-time placebo test
      lambdas <- mclapply(c(min(t.s$E2week):-1),getLambda,
                          control = c.s,
                          treatment = t.s,
                          time.ind = "E2week",
                          space.ind = space.ind,
                          val.ind = "GAI",
                          q.type = typeq,
                          beta = beta,
                          mc.cores = cores)
    }
    f.lambdas <- unlist(lapply(control.s, function(x){
      return (lapply(c(1:length(lambdas)),function(y){return (lambdas[[y]][x])}) |> unlist() |> mean())
    }))
    
  } else {
    if (!is.null(lag.T0)){ #For in-time placebo test
      lag.time = -lag.T0
      f.lambdas <- getLambda2(time = c(min(t.s$E2week):lag.time), control = c.s, treatment = t.s, time.ind = "E2week", space.ind = space.ind, val.ind = "GAI", q.type = typeq, beta = beta)
    } else { #For non in-time placebo test
      f.lambdas <- getLambda2(time = c(min(t.s$E2week):-1), control = c.s, treatment = t.s, time.ind = "E2week", space.ind = space.ind, val.ind = "GAI", q.type = typeq, beta = beta)
    }
    lambdas = NULL
  }
  names(f.lambdas) <- control.s
  print("get lambdas")
  
  #calculate prediction diff for each quantile and time
  cf <- do.call(cbind,lapply(c(1:length(quantFun)),
                             function(x){return(as.matrix(cbind(quantFun[[x]][[2]])) %*% f.lambdas)}
  ))
  t_quant <- do.call(cbind,lapply(c(1:length(quantFun)), function(x){return(quantFun[[x]][[1]])}))
  diff <- t_quant - cf
  colnames(diff) <- e2weeks
  colnames(cf) <- e2weeks
  colnames(t_quant) <- e2weeks
  if (isTRUE(in.main)) {print("Get prediction diff")}
  
  
  if (time.placebo & !is.null(time.lag)){
    print("---Start in-time placebo test---")
    time.lag = time.lag[time.lag * interval.eweek <= 52]
    placebo.res.time <- lapply(time.lag, getCf, id = id,
                               t.info = treatment_info, 
                               t.s = t.s,
                               c.s = c.s, 
                               time.ind = "Eweek_tot", 
                               space.ind = "Sector_ID", 
                               no.trap = "Number.of.functional.Gravitrap", 
                               no.mosq = "Mean.number.of.F.Ae.albopictus.caught.per.functional.trap",
                               lambda.ind = lambda.ind,
                               bootstrap = FALSE, 
                               space.placebo = FALSE,
                               time.placebo = FALSE,
                               time.lag = NULL,
                               typeq = typeq,
                               interval.eweek = interval.eweek,
                               time.weight.ind = time.weight.ind,
                               cores = cores)
    
    names(placebo.res.time) <- time.lag
    extra <- list()
    for (i in time.lag){
      extra <- c(extra, placebo.res.time[[which(names(placebo.res.time) == i)]][5])
      
    }
    names(extra) <- time.lag
    time.placebo.plots <- mapply(time.placebo.visual,time.lag = time.lag, placebo.f = extra, MoreArgs = list(tid = id, time.range = e2weeks, actual = t_quant, interval.eweek = interval.eweek, save.path = "../Rimages/time_placebo/"))
    print("---End in-time placebo test---")
  } else {
    time.placebo.plots <- NULL
  }
  
  #placebo test part
  if (space.placebo){
    print("---Start in-space placebo test---")
    perfect_distribution <- matrix(0,nrow = nrow(cf), ncol = ncol(cf)) |> as.data.frame()
    rownames(perfect_distribution) <- rownames(cf)
    colnames(perfect_distribution) <- colnames(cf)
    res.d1 <- sapply(c(1:ncol(cf)), function(x) {return (get_d(x,cf,t_quant))})
    res.d2 <- sapply(c(1:ncol(cf)), function(x) {return (get_d(x,cf,perfect_distribution))})
    treat.res <- res.d1/res.d2 #(df2 - df1)/df1
    
    treat.res <- mclapply(c(1:ncol(cf)), function(x) {return (get_IE(x,cf,t_quant))},mc.cores = cores)
    placebo.res <- pbmclapply(control.s, placebo.test, 
                              control = c.s,
                              time.ind = "E2week",
                              space.ind = space.ind,
                              val.ind = "GAI",
                              lambda.ind = lambda.ind,
                              q.type = typeq,
                              load.q.func = quantFun,
                              beta = beta,
                              mc.cores = cores)
    
    names(placebo.res) <- control.s
    treat.res <- unlist(treat.res)
    names(treat.res) <- sort(e2weeks)
    placebo.res[[length(placebo.res)+1]] <- treat.res
    names(placebo.res)[length(placebo.res)] <- id
    print("---End in-space placebo test---")
  } else {
    placebo.res = NULL
  }
  
  #real space placebo
  if (space.placebo){
    real.space.placebo.res <- pbmclapply(control.s, real.space.placebo.test, 
                              control = c.s,
                              time.ind = "E2week",
                              space.ind = "Sector_ID",
                              val.ind = "GAI",
                              lambda.ind = lambda.ind,
                              q.type = typeq,
                              load.q.func = quantFun,
                              beta = beta,
                              mc.cores = cores)
    names(real.space.placebo.res) <- control.s
  } else {
    real.space.placebo.res = NULL
  }
  return (list(diff, f.lambdas, convertweek, placebo.res, cf, t_quant, quantFun, lambdas, time.placebo.plots, real.space.placebo.res))
}

#------------------------------Part 2: auxilary functions, no need to use parallel computing---------------------------------#
#1. change percentage to decimal number
rm.last <- function(string){
  return (round(as.numeric(substr(string,1,(nchar(string) - 1)))/100,3))
}

#1.1 auxilary function for sector_id sorting
sort.sector <- function(sector.vector){
  nums <- sapply(sector.vector, function(x) {
    if (grepl("[A-Za-z]",substring(x,nchar(x)))){
      return (substr(x, 3, nchar(x) - 1))
    } else {
      return (substr(x, 3, nchar(x)))
    }})
  nums <- data.frame(true.id = names(nums), num.id = as.numeric(nums))
  nums <- nums[order(nums$num.id,nums$true.id),]
  return (nums$true.id)
}

#2. get quantile function of a set of treatment and control sectors
getQuantile <- function(time, prob_list, control, treatment, time.ind, space.ind, val.ind, q.type){
  c.temp <- control %>% filter(!!sym(time.ind) == time)
  t.temp <- treatment %>% filter(!!sym(time.ind) == time)
  c.s <- sort.sector(unique(c.temp[[space.ind]]))
  t.s <- unique(t.temp[[space.ind]])
  c.quant <- lapply(c.s,function(x){
    return (stats::quantile((c.temp %>% filter(!!sym(space.ind) == x))[[val.ind]], prob = prob_list, type = q.type))
  })
  c.quant <- do.call(cbind.data.frame, c.quant)
  colnames(c.quant) <- c.s
  t.quant <- stats::quantile(t.temp[[val.ind]], prob = prob_list, type = q.type)
  #rownames(c.quant) <- names(t.quant)
  
  res <- list(t.quant,c.quant)
  names(res) <- c(t.s,'others')
  return (res)
}

#3. Get weight of each control sector (lambda) at a certain time
getLambda <- function(time, control, treatment, time.ind, space.ind, val.ind, q.type, beta){
  #Use m ~ U(0,1) to sample 500 points to value of corresponding quantiles
  if (beta){
    m <- rbeta(500,3,3)
  } else {
    m <- runif(500)
  }
  
  res2 <- getQuantile(time, m, control, treatment, time.ind, space.ind, val.ind, q.type)
  
  #solve QP with constraint interpolation
  treat.vec <- res2[[1]]
  control.df <- res2[[2]]
  
  lambda <- Variable(length(control.df))
  objective <- Minimize(sum((treat.vec - control.df %*% lambda)^2))
  constraints <- list(sum(lambda) == 1, lambda >= 0)
  problem <- Problem(objective, constraints)
  sol <- solve(problem, solver = "ECOS") #one possible solver with no solver error
  res.lambda <- round(sol$getValue(lambda), digits = 6)
  names(res.lambda) <- names(control.df)
  return(res.lambda)
}

#4. calculate 2-Wasserstein distance
get_d <- function(time, df1, df2){
  d1x <- sapply(rownames(df1), rm.last)
  d1y <- df1[,time]
  d1 <- list(as.vector(d1x),as.vector(d1y))
  names(d1) <- c("x","y")
  d2x <- sapply(rownames(df2), rm.last)
  d2y <- df2[,time]
  d2 <- list(as.vector(d2x),as.vector(d2y))
  names(d2) <- c("x","y")
  dist <- dist4den(d1, d2,fctn_type = "quantile")
  return (dist)
}

#4.1 possible IE calculation
get_IE <- function(time, df.1, df.2){
  d1x <- sapply(rownames(df.1), rm.last)
  d1y <- df.1[,time]
  d1x0 <- d1x
  d1y0 <- c(rep(0,length(d1y)))
  d1 <- list(as.vector(d1x),as.vector(d1y))
  names(d1) <- c("x","y")
  d10 <- list(as.vector(d1x0),as.vector(d1y0))
  names(d10) <- c("x","y")
  
  d2x <- sapply(rownames(df.2), rm.last)
  d2y <- df.2[,time]
  d2 <- list(as.vector(d2x),as.vector(d2y))
  names(d2) <- c("x","y")
  d2x0 <- d2x
  d2y0 <- c(rep(0,length(d2y)))
  d20 <- list(as.vector(d2x0),as.vector(d2y0))
  names(d20) <- c("x","y")
  
  dist1 <- dist4den(d1, d10,fctn_type = "quantile")
  dist2 <- dist4den(d2, d20,fctn_type = "quantile")
  
  IE <- (dist2 - dist1)/ dist1
  return (IE)
}

#5. permuration test - Not used in albopictus-related work
placebo.test <- function(f_treatment, min.week, control, time.ind, space.ind, val.ind, q.type, load.q.func = NULL, lambda.ind = 1, beta){
  postal_c <- control %>% filter(!!sym(space.ind) != f_treatment) #rest of 247 control sectors
  postal_p <- control %>% filter(!!sym(space.ind) == f_treatment) #one control sector considered as "treated"
  
  
  # get time representation (treatment start week = 0) and a vector of control sector ID
  e2weeks <- unlist(unique(postal_p$E2week)) |> sort() |> as.vector()
  control.s2 <- unique(postal_c[[space.ind]]) |> sort.sector()
  
  if (lambda.ind == 1){
    lambdas <- lapply(c(min(postal_p[[time.ind]]):-1),function(x){
      return (getLambda(x, postal_c, postal_p, time.ind, space.ind, val.ind, q.type, beta))})
    
    f.lambdas <- unlist(lapply(control.s2, function(x){
      return (mean(unlist(lapply(c(1:length(lambdas)),function(y){
        return (lambdas[[y]][x])
      }))))
    }))
    names(f.lambdas) <- control.s2
  } else{
    f.lambdas <- getLambda2(time = c(min(postal_p[[time.ind]]):-1), control = postal_c, treatment = postal_p, time.ind = time.ind, space.ind = space.ind, val.ind = val.ind, beta = beta, q.type = q.type)
    names(f.lambdas) <- control.s2
  }
  
  if (is.null(load.q.func)){
    load.q.func = getQuantile(time, seq(0,1,length.out = 1001), postal.c, postal.p, time.ind, space.ind, val.ind, q.type)
    total.time.length = length(unlist(unique(postal_p[time.ind])))
    e2weeks <- NULL
    cf.placebo <- do.call(cbind,lapply(c(1:total.time.length),
                                       function(x){as.matrix(do.call(cbind,
                                                                     load.q.func[[x]][[2]][which(names(load.q.func[[x]][[2]]) != f_treatment)])) %*% f.lambdas}
    ))
    t_quant.placebo <- do.call(cbind,lapply(c(1:total.time.length), function(x){return(load.q.func[[x]][[1]])}))
    
  } else {
    total.time.length = length(load.q.func)
    e2weeks <- names(load.q.func)
    cf.placebo <- do.call(cbind,lapply(e2weeks,
                                       function(x){as.matrix(cbind(load.q.func[[x]][[2]][which(names(load.q.func[[which(names(load.q.func) == x)]][[2]]) != f_treatment)])) %*% f.lambdas}
    ))
    t_quant.placebo <- do.call(cbind,lapply(e2weeks, function(x){return(load.q.func[[which(names(load.q.func) == x)]][[2]][[which(names(load.q.func[[which(names(load.q.func) == x)]][[2]]) == f_treatment)]])}))
  }
  
  rownames(cf.placebo) <- rownames(load.q.func[[1]][[2]])
  rownames(t_quant.placebo) <- rownames(load.q.func[[1]][[2]])
  perfect_distribution <- matrix(0,nrow = nrow(cf.placebo), ncol = ncol(cf.placebo)) |> as.data.frame()
  rownames(perfect_distribution) <- rownames(cf.placebo)
  colnames(perfect_distribution) <- colnames(cf.placebo)
  res.d1 <- sapply(c(1:ncol(cf.placebo)), function(x) {return (get_d(x,cf.placebo,t_quant.placebo))})
  res.d2 <- sapply(c(1:ncol(cf.placebo)), function(x) {return (get_d(x,cf.placebo,perfect_distribution))})
  #res.d <- res.d1/res.d2
  
  res.d <- sapply(c(1:ncol(cf.placebo)), function(x) {return (get_IE(x,cf.placebo,t_quant.placebo))})
  
  names(res.d) <- e2weeks
  return (res.d)
}

#6.1 leave-one-out sampling on each sector (space.ind) each 2 eweeks (time.ind) on a dataframe
getlos <- function(id, df, time.ind, space.ind, n){
  df <- df %>% filter(!!sym(space.ind) == id) %>% group_by(!!sym(time.ind)) %>% sample_n(n) %>% ungroup()  
  return (df)
}

#12 in-time placebo test is done by visual inspection. Plots are required.
time.placebo.visual <- function(tid, time.range, time.lag, interval.eweek = 2, save.path = "../Rimages/time_placebo/", ...){
  kwargs <- list(...)
  placebo.time = time.range[time.range >= -time.lag & time.range < 0] |> rev()
  plot.list <- lapply(placebo.time, function(curweek){
    df.tot = data.frame(E_2week = numeric(0), group = character(0), GAI = numeric(0))
    for (name in names(kwargs)){
      df <- data.frame(quantile = seq(0,1,length.out = 1001), group = rep(name,length(1001)))
      temp.vector <- kwargs[[which(names(kwargs) == name)]]
      temp.vector <- temp.vector[,which(colnames(temp.vector) == curweek)]
      df["GAI"] <- temp.vector
      df.tot <- rbind(df.tot,df)
    }
    if (length(names(kwargs)) == 2){
      groups = names(kwargs)
      gai1 <- df.tot %>% filter(group == groups[1])
      gai1.1 <- gai1$GAI |> unlist()
      names(gai1.1) <- gai1$quantile
      gai2 <- df.tot %>% filter(group == groups[2])
      gai2.2 <- gai2$GAI |> unlist()
      names(gai2.2) <- gai2$quantile
      d1x <- names(gai1.1)|>as.numeric()
      d1y <- gai1.1
      d1 <- list(as.vector(d1x),as.vector(d1y))
      names(d1) <- c("x","y")
      d2x <- d1x
      d2y <- gai2.2
      d2 <- list(as.vector(d2x),as.vector(d2y))
      names(d2) <- c("x","y")
      dist <- dist4den(d1, d2,fctn_type = "quantile")
    }
    p <- ggplot(data = df.tot, mapping = aes(x = quantile, y = GAI, group = group, color = group)) +
      geom_line() + 
      scale_x_continuous(breaks = c(0,0.25,0.5,0.75,1), labels = c("0%","25%","50%","75%","100%")) +
      labs(x = "Quantile", y = "GAI value", color = c("actual value","counterfactual value"))+
      theme(legend.text = element_text(size = rel(0.7)), 
            plot.subtitle=element_text(size = rel(0.7)),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line.x = element_line(colour = 'black', linewidth = 0.5, linetype='solid'),
            axis.line.y = element_line(colour = 'black', linewidth = 0.5, linetype='solid'),
            panel.background = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5))+
      coord_cartesian(ylim = c(0,NA)) +
      ggtitle(paste(tid,"time placebo."), subtitle = paste(-curweek*13-12 , " to ",-curweek*13, " wks before true treat. Assume fake 'treat' at ",time.lag * 13," wks before true one. Diff = ",round(dist,3), sep = ''))
    #system(paste("mkdir -p", save.path))
    #ggsave(paste(tid,"_",time.lag,"_",curweek,"_time_placebo.png",sep = ""), plot = p, device = "png",
    #       width = 2100, height = 900, units = "px", path = save.path)
    return (p)
  })    
  return (plot.list)
}


#15 change data stored as e2week into eweek
convertdata <- function(x, conversion, bef.ind, aft.ind){
  time <- as.numeric(names(x))
  ew_tot_min <- min(conversion[which(conversion[bef.ind] == min(time)),][[aft.ind]])
  ew_tot_max <- max(conversion[which(conversion[bef.ind] == max(time)),][[aft.ind]])
  ew_tot_all <- lapply(c(ew_tot_min:ew_tot_max), function(l){
    rank_this_week = conversion[which(conversion[aft.ind] == l),][[bef.ind]] |> as.character()
    return (unname(x[which(time == rank_this_week)]))
  })
  names(ew_tot_all) <- c(ew_tot_min:ew_tot_max)
  return (ew_tot_all)
}

#16 New getLambda function; consider a convex function for all pretreatment period
getLambda2 <- function(time, control, treatment, time.ind, space.ind, val.ind, q.type, beta){
  #Use m ~ U(0,1) to sample 500 points to value of corresponding quantiles
  if(beta) {m <- rbeta(500,3,3)} else {m <- runif(500)}
  
  res2 <- lapply(time, getQuantile, prob_list = m, control = control, treatment = treatment, time.ind = time.ind, space.ind = space.ind, val.ind = val.ind, q.type = q.type)
  #solve QP with constraint interpolation
  treat.vec <- lapply(res2, function(x) {return(x[[1]])})
  control.df <- lapply(res2, function(x) {return(x[[2]])})
  
  lambda <- Variable(ncol(control.df[[1]]))
  optimization.func <- sapply(c(1:length(control.df)), function(x){
    return (sum((treat.vec[[x]] - control.df[[x]] %*% lambda)^2))
  })
  optimization.func <- Reduce('+',optimization.func)
  objective <- Minimize(optimization.func)
  constraints <- list(sum(lambda) == 1, lambda >= 0)
  problem <- Problem(objective, constraints)
  sol <- solve(problem, solver = "ECOS") #one possible solver with no solver error
  res.lambda <- round(sol$getValue(lambda), digits = 6)
  names(res.lambda) <- names(control.df)
  return(res.lambda)
}

#17 trim control sectors to a smaller size, 248 sectors -> around 30~40 with high correlation to treated sector
trim_cs <- function(id, t.info, t.s, c.s, time.ind, space.ind, no.trap, no.mosq, interval.eweek = 13, cores = 1, typeq = quantile_type){
  T0 <- t.info[which(t.info[space.ind] == id), time.ind]
  c.s[["E2week"]] <- floor((c.s[time.ind] - T0) / interval.eweek)
  t.s <- t.s %>% filter(!!sym(space.ind) == id)
  t.s[["E2week"]] <- floor((t.s[time.ind] - T0) / interval.eweek)
  c.s <- c.s %>% rename(c(no.trap = no.trap, no.mosq = no.mosq)) %>% filter(E2week < 0)
  t.s <- t.s %>% rename(c(no.trap = no.trap, no.mosq = no.mosq)) %>% filter(E2week < 0)
  c.s[["GAI"]] <- c.s$no.mosq/c.s$no.trap
  t.s[["GAI"]] <- t.s$no.mosq/t.s$no.trap
  
  e2weeks <- unlist(unique(t.s$E2week)) |> sort() |> as.vector()
  control.s <- unique(c.s$Sector_ID) |> sort.sector()
  quantFun <- pbmclapply(e2weeks, getQuantile, 
                         prob_list = seq(0,1,length.out = 1001), 
                         control = c.s,
                         treatment = t.s,
                         time.ind = "E2week",
                         space.ind = space.ind,
                         val.ind = "GAI",
                         q.type = typeq,
                         mc.cores = cores,
                         mc.silent = TRUE)
  names(quantFun) <- e2weeks
  t.quantile <- do.call(cbind,lapply(c(1:length(quantFun)), function(x){return(quantFun[[x]][[1]])}))
  
  #pretre.avg <- mclapply(control.s, function(c){
  #  temp.c <- c.s %>% filter(!!sym(space.ind) == c)
  #  res.cor <- lapply(c(min(e2weeks):-1), function(t){
  #    control.t <- temp.c %>% filter(E2week == t) %>% select(GAI) |> unlist() |> mean()
  #    treat.t <- t.s %>% filter(E2week == t) %>% select(GAI) |> unlist() |> mean()
  #    return (list(control.t, treat.t))
  #  })
  #  control.t <- sapply(res.cor, function(x) {return(x[[1]])})
  #  treat.t <- sapply(res.cor, function(x) {return(x[[2]])})
  #  correlation <- stats::cor(treat.t, control.t)
  #  return (correlation)
  #}, mc.cores = cores) |> unlist()
  
  pretre.avg <- mclapply(control.s, function(c){
    c.quantile <- sapply(c(1:length(quantFun)), function(x){
      k <- quantFun[[x]][[2]][[c]]
      return(k)
    })
    if(is.null(rownames(c.quantile))) rownames(c.quantile) = rownames(t.quantile)
    return (lapply(c(1:length(e2weeks)),get_d, df1 = t.quantile, df2 = c.quantile) |> unlist() |> mean())
  }, mc.cores = cores)
  names(pretre.avg) <- control.s
  selected.s <- sort(unlist(pretre.avg), decreasing = FALSE)[1:40]
  return (pretre.avg[control.s %in% names(selected.s)])
}
#
square <- function(x){
  return (x^2)
}
#15 in-space placebo test - no need in albopictus-related work
real.space.placebo.test <- function(f_treatment, min.week, control, time.ind, space.ind, val.ind, q.type, load.q.func = NULL, lambda.ind = 1, beta){
  postal_c <- control %>% filter(!!sym(space.ind) != f_treatment) #rest of 247 control sectors
  postal_p <- control %>% filter(!!sym(space.ind) == f_treatment) #one control sector considered as "treated"
  
  
  # get time representation (treatment start week = 0) and a vector of control sector ID
  e2weeks <- unlist(unique(postal_p$E2week)) |> sort() |> as.vector()
  control.s2 <- unique(postal_c[[space.ind]]) |> sort.sector()
  
  if (lambda.ind == 1){
    lambdas <- lapply(c(min(postal_p[[time.ind]]):-1),function(x){
      return (getLambda(x, postal_c, postal_p, time.ind, space.ind, val.ind, q.type, beta))})
    
    f.lambdas <- unlist(lapply(control.s2, function(x){
      return (mean(unlist(lapply(c(1:length(lambdas)),function(y){
        return (lambdas[[y]][x])
      }))))
    }))
    names(f.lambdas) <- control.s2
  } else{
    f.lambdas <- getLambda2(time = c(min(postal_p[[time.ind]]):-1), control = postal_c, treatment = postal_p, time.ind = time.ind, space.ind = space.ind, val.ind = val.ind, beta = beta, q.type = q.type)
    names(f.lambdas) <- control.s2
  }
  
  if (is.null(load.q.func)){
    load.q.func = getQuantile(time, seq(0,1,length.out = 1001), postal.c, postal.p, time.ind, space.ind, val.ind, q.type)
    total.time.length = length(unlist(unique(postal_p[time.ind])))
    e2weeks <- NULL
    cf.placebo <- do.call(cbind,lapply(c(1:total.time.length),
                                       function(x){as.matrix(do.call(cbind,
                                                                     load.q.func[[x]][[2]][which(names(load.q.func[[x]][[2]]) != f_treatment)])) %*% f.lambdas}
    ))
    t_quant.placebo <- do.call(cbind,lapply(c(1:total.time.length), function(x){return(load.q.func[[x]][[1]])}))
    
  } else {
    total.time.length = length(load.q.func)
    e2weeks <- names(load.q.func)
    cf.placebo <- do.call(cbind,lapply(e2weeks,
                                       function(x){as.matrix(cbind(load.q.func[[x]][[2]][which(names(load.q.func[[which(names(load.q.func) == x)]][[2]]) != f_treatment)])) %*% f.lambdas}
    ))
    t_quant.placebo <- do.call(cbind,lapply(e2weeks, function(x){return(load.q.func[[which(names(load.q.func) == x)]][[2]][[which(names(load.q.func[[which(names(load.q.func) == x)]][[2]]) == f_treatment)]])}))
  }
  
  rownames(cf.placebo) <- rownames(load.q.func[[1]][[2]])
  rownames(t_quant.placebo) <- rownames(load.q.func[[1]][[2]])
  colnames(cf.placebo) <- e2weeks
  colnames(t_quant.placebo) <- e2weeks
  return (list(cf.placebo, t_quant.placebo))
}

