#
#getLambda3 has not been modified. (lambda = 3)
#18 getCf: 
getCf.binary <- function(tid, t.info, t.s, c.s, time.ind, space.ind, no.trap, no.mosq, lambda.ind = 1, bootstrap = FALSE, space.placebo = TRUE, time.placebo = TRUE, time.lag = NULL, lag.T0 = NULL, interval.eweek = 2, time.weight.ind = FALSE, cores = 1){
  in.main <- isFALSE(bootstrap) && isTRUE(space.placebo) && isTRUE(time.placebo)
  T0 <- t.info[which(t.info[[space.ind]] == tid), time.ind]
  c.s[["E2week"]] <- floor((c.s[time.ind] - T0) / interval.eweek)
  t.s <- t.s %>% filter(!!sym(space.ind) == tid)
  t.s[["E2week"]] <- floor((t.s[[time.ind]] - T0) / interval.eweek)
  c.s <- c.s %>% rename(c(no.trap = no.trap, no.mosq = no.mosq))
  t.s <- t.s %>% rename(c(no.trap = no.trap, no.mosq = no.mosq))
  c.s[["GAI"]] <- sapply(c.s$no.mosq, function(x){if (x == 0){return (0)} else {return (1)}}) # change to binary classification: either 0 or 1
  t.s[["GAI"]] <- sapply(t.s$no.mosq, function(x){if (x == 0){return (0)} else {return (1)}}) # change to binary classification: either 0 or 1
  
  e2weeks <- unlist(unique(t.s$E2week)) |> sort() |> as.vector()
  control.s <- unique(c.s$Sector_ID) |> sort.sector()
  
  #get E2week vs e-week table (convert E2week to readable e-week), will be used in plotting line chart
  convertweek <- t.s[,c("E2week","Eweek_tot","Eyear","Eweek")] %>% distinct()
  convertweek["eyew"] <- apply(convertweek[,c("Eyear","Eweek")] , MARGIN = 1, paste, collapse = '-')
  convertweek <- convertweek %>% dplyr::select("E2week","Eweek_tot","eyew")
  E2week_n <- convertweek %>% filter(E2week < 0) %>% group_by(E2week) %>% summarize(n_wks = n()) %>% ungroup()

  if (time.weight.ind){
    valid_wks = E2week_n %>% filter(n_wks == interval.eweek) %>% select(E2week) |> unlist() |> unname()
    valid_wks = c(valid_wks, e2weeks[e2weeks >= 0])
    c.s <- c.s[(c.s$E2week |> unlist()) %in% valid_wks,]
    t.s <- t.s[(t.s$E2week |> unlist()) %in% valid_wks,]
    e2weeks <- unname(valid_wks) |> sort() |> as.vector()
  }
  
  
  if (bootstrap){
    cores = 1
    c.s <-c.s %>% group_by(!!sym(space.ind), E2week) %>% slice_sample(prop = 1, replace = TRUE) %>% ungroup() #sampling Njt records (individuals) for bootstrapping WITH replacement
    #c.s <- lapply(control.s,getlos, df = c.s, time.ind = "E2week", space.ind = space.ind, n = 1) |> do.call(what = rbind) |> suppressMessages(anti_join(x = c.s)) 
    #c.s2 <- c.s %>% tibble::rownames_to_column(var = "rowname") %>% group_by(!!sym(space.ind), E2week) %>% 
    #  slice_sample(n = 1) %>% ungroup() %>% select(rowname) %>% unlist() %>% as.numeric()
    #c.s <- c.s[-c.s2,]
    #t.s <- getlos(tid, t.s,"E2week", space.ind) |> suppressMessages(anti_join(x = t.s))
    quantFun <- lapply(e2weeks, getECDF, time.ind = "E2week", space.ind = space.ind, val.ind = "GAI",c.s = c.s, t.s = t.s, control.s = control.s, tid = tid)
    names(quantFun) <- e2weeks
    
    if (lambda.ind == 1){
      lambdas <- mclapply(c(min(t.s$E2week):-1),getLambda3,
                          control = c.s,
                          treatment = t.s,
                          time.ind = "E2week",
                          space.ind = space.ind,
                          val.ind = "GAI",
                          control.s = control.s,
                          tid = tid,
                          lambda.ind = 1,
                          mc.cores = cores)
      f.lambdas <- lambdas |> unlist() |> array(dim = c(length(control.s), abs(min(t.s$E2week)))) |> apply(MARGIN = 1, FUN = mean)
    } else if (lambda.ind == 2){
      f.lambdas <- getLambda3(time = c(min(t.s$E2week):-1),
                              control = c.s,
                              treatment = t.s,
                              time.ind = "E2week",
                              space.ind = space.ind,
                              val.ind = "GAI",
                              control.s = control.s,
                              tid = tid,
                              lambda.ind = 2,
                              cores = 1)
      lambdas = NULL
    } else if (lambda.ind == 3){
      lambdas <- lapply(c(min(t.s$E2week):-1),getLambda3,
                         control = c.s,
                         treatment = t.s,
                         time.ind = "E2week",
                         space.ind = space.ind,
                         val.ind = "GAI",
                         control.s = control.s,
                         tid = tid,
                         lambda.ind = 3,
                         mc.cores = cores)
      f.lambdas <- lambdas |> unlist() |> array(dim = c(length(control.s), abs(min(t.s$E2week)))) |> apply(MARGIN = 1, FUN = mean)
    }
    names(f.lambdas) <- control.s
    return (f.lambdas)
  }
  if(in.main) {print("Data preprocessed.")}
  
  #get ECDF for each sector (as a function, use fn(x) to get y)
  quantFun <- mclapply(e2weeks, getECDF, time.ind = "E2week", space.ind = space.ind, val.ind = "GAI",c.s = c.s, t.s = t.s, control.s = control.s, tid = tid, mc.cores = cores)
  names(quantFun) <- e2weeks
  
  if (lambda.ind == 1){
    if (!is.null(lag.T0)){ #For in-time placebo test
      lag.time = -lag.T0
      lambdas <- mclapply(c(min(t.s$E2week):lag.time),getLambda3,
                          control = c.s,
                          treatment = t.s,
                          time.ind = "E2week",
                          space.ind = space.ind,
                          val.ind = "GAI",
                          control.s = control.s,
                          tid = tid,
                          lambda.ind = 1,
                          mc.cores = cores)
    } else {
      lambdas <- mclapply(c(min(t.s$E2week):-1),getLambda3,
                          control = c.s,
                          treatment = t.s,
                          time.ind = "E2week",
                          space.ind = space.ind,
                          val.ind = "GAI",
                          control.s = control.s,
                          tid = tid,
                          lambda.ind = 1,
                          mc.cores = cores)
    }
    
    f.lambdas <- lambdas |> unlist() |> array(dim = c(length(control.s), abs(min(t.s$E2week)))) |> apply(MARGIN = 1, FUN = mean)
  } else if (lambda.ind == 2) {
    if (!is.null(lag.T0)){ #For in-time placebo test
      lag.time = -lag.T0
      f.lambdas <- getLambda3(time = c(min(t.s$E2week):lag.time),
                              control = c.s,
                              treatment = t.s,
                              time.ind = "E2week",
                              space.ind = space.ind,
                              val.ind = "GAI",
                              control.s = control.s,
                              tid = tid,
                              lambda.ind = 2,
                              cores = cores)
      
    } else {
      f.lambdas <- getLambda3(time = c(min(t.s$E2week):-1),
                              control = c.s,
                              treatment = t.s,
                              time.ind = "E2week",
                              space.ind = space.ind,
                              val.ind = "GAI",
                              control.s = control.s,
                              tid = tid,
                              lambda.ind = 2,
                              cores = cores)
    }
    lambdas = NULL
  } else if (lambda.ind == 3){
    if (!is.null(lag.T0)){ #For in-time placebo test
      lag.time = -lag.T0
      f.lambdas <- getLambda3(time = c(min(t.s$E2week):lag.time),
                            control = c.s,
                            treatment = t.s,
                            time.ind = "E2week",
                            space.ind = space.ind,
                            val.ind = "GAI",
                            control.s = control.s,
                            tid = tid,
                            lambda.ind = 3)
    } else {
      f.lambdas <- getLambda3(time = c(min(t.s$E2week):-1),
                            control = c.s,
                            treatment = t.s,
                            time.ind = "E2week",
                            space.ind = space.ind,
                            val.ind = "GAI",
                            control.s = control.s,
                            tid = tid,
                            lambda.ind = 3)
      }
    }
  names(f.lambdas) <- control.s
  if (in.main) {print("get lambdas")}
  
  cf <- sapply(c(1:length(quantFun)),function(x){return(sapply(quantFun[[x]][[2]], function(y){return(y(0))}) %*% f.lambdas)})
  t_quant <- sapply(c(1:length(quantFun)), function(x){
    k <- quantFun[[x]][[1]](0)
    return(k)
  })
  
  diff <- cf - t_quant #diff on percentage of 0
  names(diff) <- e2weeks
  names(cf) <- e2weeks
  names(t_quant) <- e2weeks
  if (!is.null(lag.T0)) {return (cf)}
  if (in.main) {print("Get prediction diff")}
  
  #----------------Robustness check for the sector-----------------
  #in-time placebo test
  if (time.placebo & !is.null(time.lag)){
    print("---Start in-time placebo test---")
    time.lag = time.lag[time.lag * interval.eweek <= 52]
    placebo.res.time <- lapply(time.lag, getCf.binary, tid = tid,
                               t.info = treatment_info, 
                               t.s = t.s,
                               c.s = c.s, 
                               time.ind = "Eweek_tot", 
                               space.ind = "Sector_ID", 
                               no.trap = "Number.of.functional.Gravitrap", 
                               no.mosq = "Number.of.wild.type.female.Ae..aegypti",
                               lambda.ind = lambda.ind,
                               bootstrap = FALSE, 
                               space.placebo = FALSE,
                               time.placebo = FALSE,
                               time.lag = NULL,
                               interval.eweek = interval.eweek,
                               cores = cores)
    
    names(placebo.res.time) <- time.lag
    time.placebo.plots <- mapply(time.placebo.visual.binary,time.lag = time.lag, placebo.f = placebo.res.time, MoreArgs = list(tid = tid, time.range = e2weeks, actual = t_quant, interval.eweek = interval.eweek, save.path = "../Rimages/time_placebo_binary/"))
    print("---End in-time placebo test---")
  } else {
    time.placebo.plots <- NULL
  }
  
  #Permutation test
  if (space.placebo){
    treat.res <- get_IE.p(cf, t_quant) #treat.res <- get_was(t_quant,cf)
    placebo.res <- pbmclapply(control.s, placebo.test.binary, 
                              control = c.s,
                              time.ind = "E2week",
                              space.ind = space.ind,
                              val.ind = "GAI",
                              lambda.ind = lambda.ind,
                              load.q.func = quantFun,
                              mc.cores = cores)
    names(placebo.res) <- control.s
    treat.res <- unlist(treat.res)
    names(treat.res) <- sort(e2weeks)
    placebo.res[[length(placebo.res)+1]] <- treat.res
    names(placebo.res)[length(placebo.res)] <- tid
    print("---End in-space placebo test---")
  } else {
    placebo.res = NULL
    placebo.res.was = NULL
  }
  
  #real in-space placebo test
  if (space.placebo){
    real.space.placebo.res <- pbmclapply(control.s, real.space.placebo.test.binary, 
                              control = c.s,
                              time.ind = "E2week",
                              space.ind = space.ind,
                              val.ind = "GAI",
                              lambda.ind = lambda.ind,
                              load.q.func = quantFun,
                              mc.cores = cores)
    names(real.space.placebo.res) <- control.s
    res2 <- do.call(rbind, real.space.placebo.res) |> apply(MARGIN = 2, mean, na.rm = TRUE)
    
    e2weeks <- names(res2) |> as.numeric()
    post.weeks <- e2weeks >= 0
    if (max(e2weeks) < 9){
      add_NA <- rep(NA, 10-sum(post.weeks))
      res2 <- c(res2[post.weeks], add_NA)
      names(res2) <- c(0:9)
    }
  } else {
    res2 = NULL
  }

  
  return (list(diff, f.lambdas, convertweek, placebo.res, cf, t_quant, quantFun, placebo.res.was = NULL, time.placebo.plots, res2))
}

#Code for Permutation test
placebo.test.binary <- function(f_treatment, control, time.ind, space.ind, val.ind, lambda.ind, load.q.func){ #no parallel work here
  postal_c <- control %>% filter(!!sym(space.ind) != f_treatment) #rest of 247 control sectors
  postal_p <- control %>% filter(!!sym(space.ind) == f_treatment) #one control sector considered as "treated"
  
  control.s2 <- unique(postal_c[[space.ind]]) |> sort.sector()
  if (lambda.ind == 1){
    lambdas <- lapply(c(min(postal_p$E2week):-1),getLambda3,
                        control = postal_c,
                        treatment = postal_p,
                        time.ind = "E2week",
                        space.ind = space.ind,
                        val.ind = "GAI",
                        control.s = control.s2,
                        tid = f_treatment,
                        lambda.ind = 1)
    f.lambdas <- lambdas |> unlist() |> array(dim = c(length(control.s2), abs(min(postal_p$E2week)))) |> apply(MARGIN = 1, FUN = mean)
  } else if (lambda.ind == 2){
    f.lambdas <- getLambda3(time = c(min(postal_p$E2week):-1),
                            control = postal_c,
                            treatment = postal_p,
                            time.ind = "E2week",
                            space.ind = space.ind,
                            val.ind = "GAI",
                            control.s = control.s2,
                            tid = f_treatment,
                            lambda.ind = 2)
    lambdas = NULL
  } else if (lambda.ind == 3){
    lambdas <- lapply(c(min(postal_p$E2week):-1),getLambda3,
                      control = postal_c,
                      treatment = postal_p,
                      time.ind = "E2week",
                      space.ind = space.ind,
                      val.ind = "GAI",
                      control.s = control.s2,
                      tid = f_treatment,
                      lambda.ind = 3)
    f.lambdas <- lambdas |> unlist() |> array(dim = c(length(control.s2), abs(min(postal_p$E2week)))) |> apply(MARGIN = 1, FUN = mean)
  } 
  if (is.null(load.q.func)){
    load.q.func <- lapply(e2weeks, getECDF, time.ind = "E2week", space.ind = space.ind, val.ind = "GAI",c.s = postal_c, t.s = postal_p, control.s = control.s2, tid = f_treatment)
  }
  names(f.lambdas) <- control.s2
  total.time.length = length(load.q.func)
  e2weeks <- names(load.q.func)
  
  cf.placebo <- sapply(c(1:length(load.q.func)),function(x){return(sapply(load.q.func[[x]][[2]][which(names(load.q.func[[x]][[2]]) != f_treatment)], function(y){return(y(0))}) %*% f.lambdas)})
  t_quant.placebo <- sapply(c(1:length(load.q.func)), function(x){ 
    k <- load.q.func[[x]][[2]][[which(names(load.q.func[[x]][[2]]) == f_treatment)]] 
    return(k(0))
  })
  diff.placebo <- get_IE.p(cf.placebo, t_quant.placebo)
  names(diff.placebo) <- e2weeks
  return (diff.placebo)
}

#Estimate the cumulative density function of the variable (discrete binary variable x: 0/1, so F(x) = 1 for x = 1)
getECDF <- function(time, time.ind, space.ind, val.ind, c.s, t.s, control.s, tid){
  c.temp <- c.s %>% filter(!!sym(time.ind) == time)
  t.temp <- t.s %>% filter(!!sym(time.ind) == time)
  c.quant <- lapply(control.s, function(x){
    return (stats::ecdf((c.temp %>% filter(!!sym(space.ind) == x))[[val.ind]]))
  })
  names(c.quant) <- control.s
  t.quant <- stats::ecdf((t.temp %>% filter(!!sym(space.ind) == tid))[[val.ind]])
  res <- list(t.quant,c.quant)
  names(res) <- c(tid,'others')
  return (res)
}

#Calculate the optimal weight of ECDFs of control sectors to fit the ECDF of treated sector during pre-treatment period
getLambda3 <- function(time, control, treatment, time.ind, space.ind, val.ind, control.s, tid, lambda.ind =1, cores = 1){
  if (lambda.ind == 1){# 1-D Wasserstein distance, each pre-intervention e-week
    res3 <- getECDF(time, time.ind, space.ind, val.ind, control, treatment, control.s, tid)
    treat.fn <- res3[[1]]
    control.fn <- res3[[2]]
    treat.0 <- treat.fn(0)
    control.0 <- sapply(control.fn, function(x){return(x(0))})
    names(control.0) <- control.s
  
    lambda <- Variable(length(control.0))
    control.0 <- array(control.0, dim = c(1,length(control.0)))
    objective <- Minimize(abs(treat.0 - control.0 %*% lambda)) 
    constraints <- list(sum(lambda) == 1, lambda >= 0)
    problem <- Problem(objective, constraints)
    
    sol <- solve(problem, solver = "ECOS")
    #sol <- solve(problem), sometimes the solver get error if not specify the name of solver.
   
    res.lambda <- round(sol$getValue(lambda), digits = 6)[,1]
    names(res.lambda) <- names(control.0)
    
  } else if (lambda.ind == 2) {# 1-D Wasserstein distance, whole pre-intervention period
    res3 <- mclapply(time, getECDF, time.ind = time.ind, space.ind = space.ind, val.ind = val.ind, c.s = control, t.s = treatment, control.s = control.s, tid = tid, mc.cores = min(length(time),cores))
    treat.0 <- sapply(res3, function(x){return (x[[1]](0))})
    control.0 <- lapply(res3, function(x){
      k <- sapply(x[[2]], function(y){return (y(0))})
      return (k)
    }) %>% do.call(what = rbind)
    colnames(control.0) <- control.s
    rownames(control.0) <- time
    
    lambda <- Variable(ncol(control.0))
    objective <- Minimize(sum(abs(treat.0 - control.0 %*% lambda)))
    constraints <- list(sum(lambda) == 1, lambda >= 0)
    problem <- Problem(objective, constraints)
    sol <- solve(problem , solver = "ECOS")
    res.lambda <- round(sol$getValue(lambda), digits = 6)
    names(res.lambda) <- colnames(control.0)
  } else if (lambda.ind == 3) { #JS divergence as objective function, whole pre-intervention period
    res3 <- mclapply(time, getECDF, time.ind = time.ind, space.ind = space.ind, val.ind = val.ind, c.s = control, t.s = treatment, control.s = control.s, tid = tid, mc.cores = min(length(time),cores))
    treat.0 <- sapply(res3, function(x){return (x[[1]](0))})
    control.0 <- lapply(res3, function(x){
      k <- sapply(x[[2]], function(y){return (y(0))})
      return (k)
    }) %>% do.call(what = rbind) %>% array(dim = c(length(time),length(res3[[1]][[2]])))
    colnames(control.0) <- control.s
    rownames(control.0) <- time
  
    lambda <- Variable(ncol(control.0))
    predicted <- control.0 %*% lambda
    kl_divergence <- kl_div(treat.0, (treat.0 + predicted)/2) + kl_div(predicted, (treat.0 + predicted)/2)
    objective <- Minimize(sum(kl_divergence))
    
    constraints <- list(sum(lambda) == 1, lambda >= 0)
    problem <- Problem(objective, constraints)
    sol <- solve(problem, solver = "SCS", abstol = 1e-13,reltol = 1e-13) #solver SCS is used to solve EXP problem
    res.lambda <- round(sol$getValue(lambda), digits = 6)
    names(res.lambda) <- colnames(control.0)
  }
  return (res.lambda)
  
}

#Code for In-time placebo test
time.placebo.visual.binary <- function(tid, time.range, time.lag, interval.eweek = 2, save.path = "../Rimages/time_placebo_binary/", ...){
  kwargs <- list(...)
  placebo.time = time.range[time.range >= -time.lag & time.range < 0]
  plot.list <- lapply(placebo.time, function(curweek){
    df.tot = data.frame(category = character(0), group = character(0), percentage = numeric(0))
    for (name in names(kwargs)){
      df <- data.frame(category = c("0","1"), group = rep(name,2))
      temp.0 <- kwargs[[which(names(kwargs) == name)]]
      temp.0 <- temp.0[which(names(temp.0) == curweek)]
      temp.1 <- 1 - temp.0
      df["percentage"] <- c(temp.0, temp.1)
      df.tot <- rbind(df.tot,df)
    }
    if (length(names(kwargs)) == 2){
      groups = names(kwargs)
      group0.0 <- df.tot %>% filter(category == "0" & group == groups[1])
      group1.0 <- df.tot %>% filter(category == "0" & group == groups[2])
      
      dist <- get_JS(group0.0$percentage, group1.0$percentage)
      print(dist)
    }
    p <- ggplot(data = df.tot, mapping = aes(x = category, y = percentage, fill = group)) +
      geom_bar(position="dodge", stat="identity") + 
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
    return (p)
    #system(paste("mkdir -p", save.path))
    #ggsave(paste(tid,"_",time.lag,"_",curweek,"_time_placebo_binary.png",sep = ""), plot = p, device = "png",
    #       width = 2100, height = 900, units = "px", path = save.path)
  })
  return(plot.list)
}

#Calculate KL-divergence
get_KL <- function(prob1,prob2){
    kl <- mapply(function(x,y){
      if (x == 1){
        return(- x * log(y))
      } else if (x == 0){
        return(log(1-y))
      } else {
        return (-(x * log(y / x) + (1-x) * log((1-y)/(1-x))))
      }
    },prob1, prob2)
    #kl = -(prob1 * log(prob2 / prob1) + (1-prob1) * log((1-prob2)/(1-prob1)))
    names(kl) <- names(prob1)
    return (kl)
}

#Calculate JS divergence
get_JS <- function(prob1,prob2){
    avg_distribution = (prob1 + prob2)/2
    js = 0.5 * get_KL(prob1, avg_distribution) + 0.5 * get_KL(prob2, avg_distribution)
    names(js) <- names(prob1)
    return (js)
}

#Calculate 1-D wasserstein distance
get_was <- function(prob1, prob2){
    p.actual1 <- 1 - prob1
    p.cf1 <- 1 - prob2
    return (abs(p.actual1 - p.cf1))
}

#Calculate intervention efficacy using JS divergence as difference.
get_IE.p <- function(cf.prob, t_quant.prob){
  prob2 <- rep(1, length(cf.prob))
  #IE = ((prob2 - t_quant.prob) - (prob2 - cf.prob))/(prob2 - cf.prob)
  cf.ce <- get_JS(prob2, cf.prob)
  t_quant.ce <- get_JS(prob2, t_quant.prob)
  IE <- (t_quant.ce - cf.ce)/cf.ce
  return (-IE)
}

#Not used - permutation test plot
process.placebo <- function(placebo.res, t.id, time.ind, alpha.num, conversion, before, after, new.xlabel, save.path = "../Rimages/space_placebo/", t.info = treatment_info, interval.eweek = 2){
  e2week <- names(placebo.res[[1]]) |> as.numeric()
  ce_distance <- placebo.res |> unlist() |> matrix(nrow = length(e2week), ncol = length(placebo.res))
  ranking <- ce_distance |> apply(MARGIN = 1, FUN = function(x){return (sum(x >= x[length(x)]))}) 
  p.rank <- ranking / length(placebo.res)
  
  df.tot <- data.frame(diff = numeric(0), E2week = numeric(0), sector = character(0))
  df.tot <- do.call(rbind, lapply(names(placebo.res), function(x, df.tot){
    li <- placebo.res[[which(names(placebo.res) == x)]]
    df <- data.frame(diff = li, E2week = names(li), sector = rep(x, length(li)))
    return (df)
  }))
  rownames(df.tot) <- NULL
  
  ew_tot_min <- min(conversion[which(conversion[before] == min(e2week)),][[after]])
  ew_tot_max <- max(conversion[which(conversion[before] == max(e2week)),][[after]])
  ew.df.tot <- lapply(c(ew_tot_min:ew_tot_max), function(t){
    IE_this_week = conversion[which(conversion[after] == t),][[before]] |> as.character()
    ew.df <- df.tot[which(df.tot[,which(colnames(df.tot) == "E2week")] == IE_this_week),]
    ew.df["Eweek_tot"] <- rep(t, nrow(ew.df))
    return(ew.df)
  }) |> do.call(what = rbind)
  
  x_break = seq(from = ew_tot_min, to = ew_tot_max, by = 26)
  if (!(183 %in% x_break)) {x_break = append(x_break,183)}
  x_label = sapply(x_break, function(x){return(conversion[which(conversion[after] == x),][["eyew"]])})
  treatment_start = t.info[t.info["Sector_ID"] == t.id, "Eweek_tot"]
  manual_colors <- c("treatment_start" = "green")
  
  
  placebo.plot = ggplot() + 
    geom_line(data = ew.df.tot, aes(x = Eweek_tot, y = diff, group = sector), alpha = 0.2) + 
    geom_line(data = ew.df.tot %>% filter(sector == t.id), aes(x = Eweek_tot, y = diff), color = 'red') +
    geom_vline(aes(xintercept = treatment_start, color= "treatment_start"), linetype = "longdash", alpha = 0.5) + 
    scale_color_manual(name = "", values = manual_colors) +
    labs(x = "Epidemiological week", y ="IE", color = "Legend") +
    scale_x_continuous(breaks = x_break, labels = x_label) +
    coord_cartesian(ylim = c(NA,1)) +
    ggtitle(paste("Sector",t.id,"in-space placebo test"), subtitle = paste("Use JS divergence in IE calculation. ", interval.eweek,"wks combined"))
  
  #system(paste("mkdir -p",save.path))
  #ggsave(paste(t.id,"_space_placebo.png",sep = ""), plot = placebo.plot, device = "png",
  #       width = 2100,height = 900, units = "px", path = save.path)
  
  #ranking
  p.rank <- ranking / length(placebo.res)
  names(p.rank) <- e2week
  ranking_ew_tot <- sapply(c(ew_tot_min:ew_tot_max), function(x){
    IE_this_week = conversion[which(conversion[after] == x),][[before]] |> as.character()
    return (p.rank[which(names(p.rank) == IE_this_week)])
  })
  p.rank <- data.frame(eweek = c(ew_tot_min:ew_tot_max))
  p.rank['ranks'] <- ranking_ew_tot
  manual_colors_2 = c("treatment_start" = "green", "0.05" = "lightblue")
  ranking.plot = ggplot() + 
    geom_line(data = p.rank, aes(x = eweek,y = ranks)) +
    geom_vline(aes(xintercept = treatment_start, color= "treatment_start"), linetype = "longdash", alpha = 0.5) + 
    geom_hline(aes(yintercept = 0.05, color= "0.05"), linetype = "solid", alpha = 0.7) +
    scale_color_manual(name = "", values = manual_colors_2) +
    labs(x = "Epidemiological week", y ="p-value", color = "Legend") +
    scale_x_continuous(breaks = x_break, labels = x_label) +
    ggtitle(paste("Sector",t.id,"in-space placebo test"), subtitle = paste("Ranking.",interval.eweek,"wks combined"))
  
  #ggsave(paste(t.id,"_in_space_rank.png",sep = ""), plot = ranking.plot, device = "png",
  #       width = 2100,height = 900, units = "px", path = save.path)
  return(list(placebo.plot, ranking.plot))
}


#select top 40 similar sectors for the treated sector.
trim_cs_binary <- function(tid, t.info, t.s, c.s, time.ind, space.ind, no.trap, no.mosq, interval.eweek = 13, cores = 1, JS = TRUE){
  T0 <- t.info[which(t.info[space.ind] == tid), time.ind]
  c.s[["E2week"]] <- floor((c.s[time.ind] - T0) / interval.eweek)
  t.s <- t.s %>% filter(!!sym(space.ind) == tid)
  t.s[["E2week"]] <- floor((t.s[time.ind] - T0) / interval.eweek)
  c.s <- c.s %>% rename(c(no.trap = no.trap, no.mosq = no.mosq)) %>% filter(E2week < 0)
  t.s <- t.s %>% rename(c(no.trap = no.trap, no.mosq = no.mosq)) %>% filter(E2week < 0)
  c.s[["GAI"]] <- sapply(c.s$no.mosq, function(x){if (x == 0){return (0)} else {return (1)}}) 
  t.s[["GAI"]] <- sapply(t.s$no.mosq, function(x){if (x == 0){return (0)} else {return (1)}}) 
  
  e2weeks <- unlist(unique(t.s$E2week)) |> sort() |> as.vector()
  control.s <- unique(c.s$Sector_ID) |> sort.sector()
  quantFun <- mclapply(e2weeks, getECDF, time.ind = "E2week", space.ind = space.ind, val.ind = "GAI",c.s = c.s, t.s = t.s, control.s = control.s, tid = tid, mc.cores = cores)
  names(quantFun) <- e2weeks
  t.ecdf <- sapply(c(1:length(quantFun)), function(x){
    k <- quantFun[[x]][[1]](0)
    return(k)
  })
  if (JS){
    pretre.avg <- mclapply(control.s, function(c){
    c.ecdf <- sapply(c(1:length(quantFun)), function(x){
      k <- quantFun[[x]][[2]][[c]](0)
      return(k)
    })
    return(get_JS(t.ecdf, c.ecdf) |> unlist() |> mean())
  }, mc.cores = cores) |> unlist()
  } else {
    pretre.avg <- mclapply(control.s, function(c){
      c.ecdf <- sapply(c(1:length(quantFun)), function(x){
        k <- quantFun[[x]][[2]][[c]](0)
        return(k)
      })
      return(get_was(t.ecdf, c.ecdf) |> abs() |> unlist() |> mean())
    }, mc.cores = cores) |> unlist()
  }
  
  names(pretre.avg) <- control.s
  selected.s <- sort(pretre.avg, decreasing = FALSE)[1:40]
  return (pretre.avg[control.s %in% names(selected.s)])
}

#1.0 change percentage to decimal number
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

#1.2
square <- function(x){
  return (x^2)
}

#1.3 resampling records
getlos <- function(id, df, time.ind, space.ind, n){
  df <- df %>% filter(!!sym(space.ind) == id) %>% group_by(!!sym(time.ind)) %>% sample_n(n)
  return (df)
}

#1.4 convert timestamp from "E2week" (relative week number based on treatment start time of the sector) to "Eweek_tot"(a unified week number which "2019-08" = week 01)
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

#1.5 a safe way to calculate log of x.
safe_log <- function(x, tol = 1e-16){
  log.x <- log(x + rep(tol, length(x)))
  return (log.x)
}

#real space placebo test
real.space.placebo.test.binary <- function(f_treatment, control, time.ind, space.ind, val.ind, lambda.ind, load.q.func){ #no parallel work here
  postal_c <- control %>% filter(!!sym(space.ind) != f_treatment) #rest of 247 control sectors
  postal_p <- control %>% filter(!!sym(space.ind) == f_treatment) #one control sector considered as "treated"
  
  control.s2 <- unique(postal_c[[space.ind]]) |> sort.sector()
  if (lambda.ind == 1){
    lambdas <- lapply(c(min(postal_p$E2week):-1),getLambda3,
                      control = postal_c,
                      treatment = postal_p,
                      time.ind = "E2week",
                      space.ind = space.ind,
                      val.ind = "GAI",
                      control.s = control.s2,
                      tid = f_treatment,
                      lambda.ind = 1)
    f.lambdas <- lambdas |> unlist() |> array(dim = c(length(control.s2), abs(min(postal_p$E2week)))) |> apply(MARGIN = 1, FUN = mean)
  } else if (lambda.ind == 2){
    f.lambdas <- getLambda3(time = c(min(postal_p$E2week):-1),
                            control = postal_c,
                            treatment = postal_p,
                            time.ind = "E2week",
                            space.ind = space.ind,
                            val.ind = "GAI",
                            control.s = control.s2,
                            tid = f_treatment,
                            lambda.ind = 2)
    lambdas = NULL
  } else if (lambda.ind == 3){
    lambdas <- lapply(c(min(postal_p$E2week):-1),getLambda3,
                      control = postal_c,
                      treatment = postal_p,
                      time.ind = "E2week",
                      space.ind = space.ind,
                      val.ind = "GAI",
                      control.s = control.s2,
                      tid = f_treatment,
                      lambda.ind = 3)
    f.lambdas <- lambdas |> unlist() |> array(dim = c(length(control.s2), abs(min(postal_p$E2week)))) |> apply(MARGIN = 1, FUN = mean)
  } 
  if (is.null(load.q.func)){
    load.q.func <- lapply(e2weeks, getECDF, time.ind = "E2week", space.ind = space.ind, val.ind = "GAI",c.s = postal_c, t.s = postal_p, control.s = control.s2, tid = f_treatment)
  }
  names(f.lambdas) <- control.s2
  total.time.length = length(load.q.func)
  e2weeks <- names(load.q.func)
  
  cf.placebo <- sapply(c(1:length(load.q.func)),function(x){return(sapply(load.q.func[[x]][[2]][which(names(load.q.func[[x]][[2]]) != f_treatment)], function(y){return(y(0))}) %*% f.lambdas)})
  t_quant.placebo <- sapply(c(1:length(load.q.func)), function(x){ 
    k <- load.q.func[[x]][[2]][[which(names(load.q.func[[x]][[2]]) == f_treatment)]] 
    return(k(0))
  })
  diff.placebo <- get_JS(cf.placebo, t_quant.placebo)
  names(diff.placebo) <- e2weeks
  return (diff.placebo)
}
