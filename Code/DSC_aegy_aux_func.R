#source functions
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
  
  #-------------Three robustness check------------#
  # in-time placebo test
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
                               no.mosq = "Number.of.wild.type.female.Ae..aegypti",
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
  
  # permutation test
  if (space.placebo){
    print("---Start in-space placebo test---")
    treat.res <- mclapply(c(1:ncol(cf)), function(x) {return (get_IE(x,cf,t_quant))},mc.cores = cores) #(df2 - df1)/df1
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
  
  # in-space placebo
  if(space.placebo){
    real.space.placebo <- pbmclapply(control.s, real.space.placebo.test, 
                              control = c.s,
                              time.ind = "E2week",
                              space.ind = "Sector_ID",
                              val.ind = "GAI",
                              lambda.ind = lambda.ind,
                              q.type = typeq,
                              load.q.func = quantFun,
                              beta = beta,
                              mc.cores = cores)
    names(real.space.placebo) <- control.s
  } else {
    real.space.placebo = NULL
    }
    
  
  
  return (list(diff, f.lambdas, convertweek, placebo.res, cf, t_quant, quantFun, lambdas, time.placebo.plots,real.space.placebo))
}


#Visual inspection based on actual CF and bootstrap CF. For less time, directly use the result from bootstrap.
#Plot graph for treatment sectors during pretreatment period and do visual inspection on the graphs
vi <- function(CIres, CIres.b, tid, B, save.path = "../Rimages/VI/", interval.eweek = 2, cores = 1){
  actual <- CIres[[6]]
  cf.0 <- CIres[[5]]
  cf.b <- mclapply(CIres.b, function(x){return (abs(x[[3]] - cf.0))}, mc.cores = cores, mc.silent = TRUE)
  cf.b.arr <- array(data = unlist(cf.b), dim = c(1001, ncol(cf.0), B)) # row 1001 (0-1), column eweek, list B
  cf.b.arr <- apply(cf.b.arr,c(1,2),function(x){return (sort(x, decreasing = FALSE) [ceiling(B*0.95)])})
  
  rownames(cf.b.arr) <- rownames(cf.b[[1]])
  colnames(cf.b.arr) <- colnames(cf.b[[1]])
  
  time.conversion <- CIres[[3]]
  
  preprocessedVI <- list(cf.b.arr, cf.0, time.conversion, actual)
  pbmclapply(colnames(cf.0), plotVI, tid = tid, preprocessed_CI = preprocessedVI,
             time.ind = "E2week", save.path = save.path, interval.eweek = interval.eweek, mc.cores = cores)
  return (NULL)
}

vi3 <- function(CIres, CIres.b, tid, B, save.path = "../Rimages/VI3/", interval.eweek = 2, cores = 1){ 
  #Get maximum diff for each quantile per sector per time interval and then take 5% largest value (0.95)
  #Therefore the CI width are the same inside one image
  actual <- CIres[[6]]
  cf.0 <- CIres[[5]]
  quantFun <- CIres[[7]]
  cf.b <- mclapply(CIres.b, function(temp.lambda){
    cf <- do.call(cbind,lapply(c(1:length(quantFun)),
                               function(y){return(as.matrix(cbind(quantFun[[y]][[2]])) %*% temp.lambda)}
    ))
    return (abs(cf - cf.0))},
    mc.cores = cores, mc.silent = TRUE)
  cf.b.arr <- array(data = unlist(cf.b), dim = c(1001, ncol(cf.0), B)) # row 1001 (0-1), column eweek, list B
  cf.b.arr <- apply(cf.b.arr,c(1,2),function(x){return (sort(x, decreasing = FALSE) [ceiling(B*1)])}) |> apply(c(2), function(x){return (sort(x, decreasing = FALSE) [ceiling(1001*0.95)])})
  
  names(cf.b.arr) <- colnames(actual)
  
  time.conversion <- CIres[[3]]
  
  preprocessedVI <- list(cf.b.arr, cf.0, time.conversion, actual)
  pbmclapply(colnames(actual), plotVI2, tid = tid, preprocessed_CI = preprocessedVI,
             time.ind = "E2week", save.path = save.path, interval.eweek = interval.eweek, mc.cores = 20)
  return (NULL)
}

plotVI <- function(tid, preprocessed_CI, time, time.ind, interval.eweek, save.path = "../Rimages/VI/"){
  widthdf <- preprocessed_CI[[1]]
  syndf <- preprocessed_CI[[2]]
  time.conversion <- preprocessed_CI[[3]]
  actualdf <- preprocessed_CI[[4]]
  
  distance <- get_d(which(colnames(syndf) == time), df1 = syndf, df2 = actualdf)
  width = widthdf[,which(colnames(widthdf) == time)] #get CF GAI variance for a e-2week
  synval = syndf[,which(colnames(syndf) == time)] #get CF GAI for a e-2week
  actual = actualdf[,which(colnames(actualdf) == time)] #get actual GAI for a e-2week
  
  width <- data.frame(quantile = sapply(names(width),rm.last), GAI = unname(width)) #change percentage to decimal
  synval <- data.frame(quantile = sapply(names(synval), rm.last), GAI = unname(synval)) #change percentage to decimal
  actual <- data.frame(quantile = sapply(names(actual), rm.last), GAI = unname(actual))
  
  diff95 <- width$GAI
  time.min <- min(time.conversion[which(time.conversion[[time.ind]] == time),"Eweek_tot"])
  time.max <- max(time.conversion[which(time.conversion[[time.ind]] == time),"Eweek_tot"])
  
  p = ggplot() +    
    geom_line(data = synval, aes(x = quantile, y = GAI, color = "synthetic")) + 
    geom_line(data = actual, aes(x = quantile, y = GAI, color = "actual")) +
    geom_ribbon(data = synval, aes(x = quantile, ymin = GAI-diff95, ymax = GAI+diff95),alpha = 0.2) + # geom_ribbon function is used to add confidence interval 
    labs(x = "percentile", y = "GAI") +
    coord_cartesian(ylim = c(0,NA)) +
    scale_color_manual("",breaks = c("synthetic","actual"), values = c("synthetic" = 'yellow',"actual" = "blue"))+
    ggtitle(paste("Sector",tid,"at e-week",toString(time.min),'-',toString(time.max)), subtitle = paste("Distance = ",round(distance,5),'. ', interval.eweek,"wks combined."))
  
  #save the image of IE plot
  system(paste("mkdir -p",save.path))
  ggsave(paste(tid,"_VI_",time,".png",sep = ""), plot = p, device = "png",
         width = 2100, height = 900, units = "px", path = save.path)
  
  return (NULL)
}

plotVI2 <- function(tid, preprocessed_CI, time, time.ind, interval.eweek, save.path = "../Rimages/VI3/"){
  widthdf <- preprocessed_CI[[1]]
  syndf <- preprocessed_CI[[2]]
  time.conversion <- preprocessed_CI[[3]]
  actualdf <- preprocessed_CI[[4]]
  
  distance <- get_d(which(colnames(syndf) == time), df1 = syndf, df2 = actualdf)
  width = widthdf[which(names(widthdf) == time)] #get CF GAI variance for a e-2week
  synval = syndf[,which(colnames(syndf) == time)] #get CF GAI for a e-2week
  actual = actualdf[,which(colnames(actualdf) == time)] #get actual GAI for a e-2week
  
  synval <- data.frame(quantile = sapply(names(synval), rm.last), GAI = unname(synval)) #change percentage to decimal
  actual <- data.frame(quantile = sapply(names(actual), rm.last), GAI = unname(actual))
  
  diff95 <- width
  time.min <- min(time.conversion[which(time.conversion[[time.ind]] == time),"Eweek_tot"])
  time.max <- max(time.conversion[which(time.conversion[[time.ind]] == time),"Eweek_tot"])
  p = ggplot() +    
    geom_line(data = synval, aes(x = quantile, y = GAI, color = "synthetic")) + 
    geom_line(data = actual, aes(x = quantile, y = GAI, color = "actual")) +
    geom_ribbon(data = synval, aes(x = quantile, ymin = GAI-diff95, ymax = GAI+diff95),alpha = 0.2) + # geom_ribbon function is used to add confidence interval 
    labs(x = "percentile", y = "GAI") +
    coord_cartesian(ylim = c(0,NA)) +
    scale_color_manual("",breaks = c("synthetic","actual"), values = c("synthetic" = 'yellow',"actual" = "blue"))+
    ggtitle(paste("Sector",tid,"at e-week",toString(time.min),'-',toString(time.max)), subtitle = paste("Distance = ",round(distance,5),'. ', interval.eweek,"wks combined."))
  
  #save the image of IE plot
  system(paste("mkdir -p",save.path))
  ggsave(paste(tid,"_VI_",time,".png",sep = ""), plot = p, device = "png",
         width = 2100, height = 900, units = "px", path = save.path)
  
  return (NULL)
}

#7.2 get all plots from start to end
getplots <- function(tid, CIres, CIres.b, time.ind, save.path = "../Rimages/Diff/", interval.eweek = 2, cores = 1){
  preprocessed_CI <- preprocessCI(CIres, CIres.b, time.ind)
  record_start = preprocessed_CI[[5]]
  treatment_length = preprocessed_CI[[4]]
  pbmclapply((record_start:treatment_length),plotCI, tid = tid, preprocessed_CI = preprocessed_CI, time.ind = time.ind, interval.eweek = interval.eweek, save.path = save.path ,mc.cores = cores)
  return (NULL)
}

#8 get in-space placebo test related info. only 1 core used
process.placebo <- function(placebo.res, t.id, time.ind, alpha.num, conversion, before, after, new.xlabel, save.path = "../Rimages/space_placebo/", t.info = treatment_info, interval.eweek = 2){
  placebo.df <- lapply(names(placebo.res), build.placebo.df, placebo.li = placebo.res, conversion = conversion, before = before, after = after, new.xlabel = new.xlabel) 
  placebo.performance <- lapply(placebo.res, function(x){return(x[names(x) < 0])}) |> lapply(mean) |> unlist()
  mean.var <- var(placebo.performance)
  treat.performance <- placebo.performance[names(placebo.performance) == t.id]
  use.sector <- names(placebo.performance)[which(placebo.performance <= (as.numeric(treat.performance) + 1.96*mean.var))]
  x_break <- placebo.df[[1]][[2]]
  if (!(183 %in% x_break)) {x_break = append(x_break,183)}
  x_label <- placebo.df[[1]][[3]]
  ew_tot_min <- placebo.df[[1]][[4]]
  ew_tot_max <- placebo.df[[1]][[5]]
  placebo.df <- lapply(placebo.df, function(x){return(x[[1]])})|> do.call(what = rbind)
  treatment_start = t.info[t.info["Sector_ID"] == t.id, "Eweek_tot"]
  manual_colors <- c("treatment_start" = "green")
  
  min_value = floor(-max(placebo.df$Diff)*2)/2
  min_treat_value = floor((-max(placebo.df %>% filter(Sector_ID == t.id) %>% select(Diff)))*2)/2
  min_value = min(max(-1,min_value),min_treat_value)
  
  placebo.plot = ggplot() + 
    geom_line(data = placebo.df, aes(x = eweek,y = -Diff, group = Sector_ID), alpha = alpha.num) + 
    geom_line(data = placebo.df %>% filter(!!sym("Sector_ID") == t.id), aes(x = eweek, y = -Diff), color = 'red') +
    geom_vline(aes(xintercept = treatment_start, color= "Event start"),linewidth = 0.4, linetype = "longdash", alpha = 0.3) + 
    scale_color_manual(name = "", values = manual_colors) +
    labs(x = "Epidemiological week", y ="Intervention Efficacy", color = "Legend") +
    scale_x_continuous(breaks = x_break, labels = x_label) +
    scale_y_continuous(breaks = seq(min_value, 1, by = 0.5), labels = sapply(seq(min_value, 1, by = 0.5) ,function(x){return(paste(as.character(x*100),'%', sep = ''))}), limits = c(min_value,1)) +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5),
          legend.title = element_text(size = 8), 
          legend.text = element_text(size=6), 
          plot.subtitle=element_text(size=10),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line.x = element_line(colour = 'black', linewidth = 0.5, linetype='solid'),
          axis.line.y = element_line(colour = 'black', linewidth = 0.5, linetype='solid'),
          panel.background = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5)) +
    ggtitle(paste("Sector",t.id,"in-space placebo test"), subtitle = paste("Surveillance data combined every",interval.eweek,"weeks."))
  
  #system(paste("mkdir -p",save.path))
  #ggsave(paste(t.id,"_space_placebo.png",sep = ""), plot = placebo.plot, device = "png",
  #       width = 2100,height = 900, units = "px", path = save.path)
  
  #Some control sectors are not used because of bad performance during pre-treatment period
  no.of.sectors <- 0 #no need this.
  if (no.of.sectors > 0){
    placebo.plot3 = ggplot() +
    geom_line(data = placebo.df %>% filter(!!sym("Sector_ID") %in% use.sector), aes(x = eweek,y = Diff, group = Sector_ID), alpha = alpha.num) + 
    geom_line(data = placebo.df %>% filter(!!sym("Sector_ID") == t.id), aes(x = eweek, y = Diff), color = 'red') +
    geom_vline(aes(xintercept = treatment_start, color= "treatment_start"), linetype = "longdash", alpha = 0.5) + 
    scale_color_manual(name = "", values = manual_colors) +
    labs(x = "Epidemiological week", y ="Diff", color = "Legend") +
    scale_x_continuous(breaks = x_break, labels = x_label) +
    ggtitle(paste("Sector",t.id,"in-space placebo test 2"), subtitle = paste("Use 2-Wasserstein to get diff.",no.of.sectors,"control sectors used.", interval.eweek,"wks combined"))
  
    ggsave(paste(t.id,"_space_placebo_2.png",sep = ""), plot = placebo.plot3, device = "png",
         width = 2100,height = 900, units = "px", path = save.path) 
  }
  
  
  placebo.df2 <- data.frame(placebo.res)
  target.id = which(colnames(placebo.df2) == t.id)
  ranking <- apply(placebo.df2, 1, 
                   function(x){
                     target = x[target.id]
                     d = sum(x <= target)   #IE is from -1 to Inf, smaller, better.
                     return (d/length(x))
                   })
  names(ranking) <- rownames(placebo.df2)
  
  ranking_ew_tot <- sapply(c(ew_tot_min:ew_tot_max), function(x){
    IE_this_week = conversion[which(conversion[after] == x),][[before]] |> as.character()
    return (ranking[which(names(ranking) == IE_this_week)])
  })
  ranking <- data.frame(eweek = c(ew_tot_min:ew_tot_max))
  ranking['ranks'] <- ranking_ew_tot
  percent_pass <- ranking[ranking$eweek >= treatment_start,]
  percent_pass <- (sum(percent_pass$ranks <= 0.05)/ nrow(percent_pass)) |> round(digits = 4)
  manual_colors_2 = c("treatment_start" = "green", "0.05" = "lightblue")
  ranking.plot = ggplot() + 
    geom_line(data = ranking, aes(x = eweek,y = ranks)) +
    geom_vline(aes(xintercept = treatment_start, color= "treatment_start"), linetype = "longdash", alpha = 0.5) + 
    geom_hline(aes(yintercept = 0.05, color= "0.05"), linetype = "solid", alpha = 0.7) +
    scale_color_manual(name = "", values = manual_colors_2) +
    labs(x = "Epidemiological week", y ="p-value", color = "Legend") +
    scale_x_continuous(breaks = x_break, labels = x_label) +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5),
          legend.title = element_text(size = 8), 
          legend.text = element_text(size=6), 
          plot.subtitle=element_text(size=10),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line.x = element_line(colour = 'black', linewidth = 0.5, linetype='solid'),
          axis.line.y = element_line(colour = 'black', linewidth = 0.5, linetype='solid'),
          panel.background = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5))+
    ggtitle(paste("Sector",t.id,"in-space placebo test"), subtitle = paste(percent_pass*100,"% pass. ",interval.eweek,"wks combined"))
  
  #ggsave(paste(t.id,"_in_space_rank.png",sep = ""), plot = ranking.plot, device = "png",
  #       width = 2100,height = 900, units = "px", path = save.path)
  
  if (no.of.sectors > 0){
    placebo.df2 <- placebo.df2[,which(colnames(placebo.df2) %in% use.sector)]
    target.id = which(colnames(placebo.df2) == t.id)
    ranking2 <- apply(placebo.df2, 1, 
                      function(x){
                        target = x[target.id]
                        d = sum(x >= target)
                        return (d/length(x))
                      })
    names(ranking2) <- rownames(placebo.df2)
    
    ranking_ew_tot2 <- sapply(c(ew_tot_min:ew_tot_max), function(x){
      IE_this_week = conversion[which(conversion[after] == x),][[before]] |> as.character()
      return (ranking2[which(names(ranking2) == IE_this_week)])
    })
    ranking2 <- data.frame(eweek = c(ew_tot_min:ew_tot_max))
    ranking2['ranks'] <- ranking_ew_tot2
    manual_colors_2 = c("treatment_start" = "green", "0.05" = "lightblue")
    ranking2.plot = ggplot() + 
      geom_line(data = ranking2, aes(x = eweek,y = ranks)) +
      geom_vline(aes(xintercept = treatment_start, color= "treatment_start"), linetype = "longdash", alpha = 0.5) + 
      geom_hline(aes(yintercept = 0.05, color= "0.05"), linetype = "solid", alpha = 0.7) +
      scale_color_manual(name = "", values = manual_colors_2) +
      labs(x = "Epidemiological week", y ="p-value", color = "Legend") +
      scale_x_continuous(breaks = x_break, labels = x_label) +
      ggtitle(paste("Sector",t.id,"in-space placebo test 2"), subtitle = paste("Ranking with" ,no.of.sectors,"control sectors.", interval.eweek,"wks combined"))
    
    ggsave(paste(t.id,"_in_space_rank2.png",sep = ""), plot = ranking2.plot, device = "png",
           width = 2100,height = 900, units = "px", path = save.path)
  } else {
    ranking2 = NULL
  }
  
  return (list(placebo.plot, ranking.plot))
}


#11 IE permutation testing, using 2-Wasserstein dist IE
IE_permutation <- function(tid, cid, t.info, t.s, c.s, time.ind, space.ind, no.trap, no.mosq, lambda.ind = 1, quantFun.pre = NULL,typeq = quantile_type, interval.eweek = 2, time.weight.ind = FALSE, beta = FALSE){
  all.sector <- c(cid, tid)
  t.s <- t.s %>% filter(!!sym(space.ind) == tid)
  all.s <- rbind(c.s,t.s)
  tid.p <- sample(all.sector, 1)
  c.s <- all.s %>% filter(!!sym(space.ind) != tid.p)
  t.s <- all.s %>% filter(!!sym(space.ind) == tid.p)
  T0.IE <- t.info[which(t.info[space.ind] == tid), time.ind]
  remove(all.s)
  remove(all.sector)
  res.p <- getCf.p(id = tid.p,
                   t.info = t.info,
                   t.s = t.s,
                   c.s = c.s, 
                   time.ind = time.ind, 
                   space.ind = space.ind, 
                   no.trap = no.trap,
                   no.mosq = no.mosq, 
                   p.T0 = T0.IE,
                   quantFun.pre = quantFun.pre,
                   typeq = typeq,
                   space.agg = NULL,
                   lambda.ind = lambda.ind,
                   interval.eweek = interval.eweek,
                   time.weight.ind = time.weight.ind,
                   beta = beta
  )
  cf.p <- res.p[[1]]
  t_quant.p <- res.p[[2]]
  e2weeks <- as.numeric(colnames(cf.p))
  distance.p <- sapply(c(1:length(e2weeks)), get_d, df1 = cf.p, df2 = t_quant.p)
  names(distance.p) <- e2weeks
  dist.avg.bef <- sapply(which(e2weeks < 0), get_d, df1 = cf.p, df2 = t_quant.p) |> square() |> mean() |> sqrt() #!!!!!!
  dist.avg.aft <- sapply(which(e2weeks >= 0), get_d, df1 = cf.p, df2 = t_quant.p) |> square() |> mean() |> sqrt()
  Rj <- dist.avg.aft/dist.avg.bef
  IE.p <- sapply(which(e2weeks >= 0), get_IE, df.1 = cf.p, df.2 = t_quant.p)
  names(IE.p) <- e2weeks[which(e2weeks>=0)]
  return (list(IE.p, Rj, distance.p))
}

#13.1 calculate IE based on averages and its CI
average.IE <- function(actual_qf, cf_qf, cf.b_qf,time_conversion, previous.ind, current.ind, cores = 1){
  actual.mean <- average.qf(actual_qf, time_conversion, previous.ind, current.ind)
  cf.mean <- average.qf(cf_qf,time_conversion, previous.ind, current.ind)
  cf.b.mean <- mclapply(cf.b_qf, average.qf, time_conversion = time_conversion, previous.ind = previous.ind, current.ind = current.ind, mc.cores = cores)
  
  diff.mean <- actual.mean - cf.mean
  diff.b.mean <-  mclapply(cf.b.mean, function(x){return (actual.mean - x)})
  diff.95.mean <- array(unlist(diff.b.mean), dim = c(length(cf.mean),length(cf.b.mean)))
  diff.975.mean <- apply(diff.95.mean, MARGIN = 1, 
                        function(x){
                          if (is.integer(B*0.975)) return (sort(x, decreasing = FALSE) [ceiling(B*0.975)])
                          else {
                            upper <- ceiling(B*0.975)
                            lower <- floor(B*0.975)
                            val <- sort(x, decreasing = FALSE)[upper]*abs(upper - B*0.975) + sort(x, decreasing = FALSE)[lower]* abs(B*0.975 - lower)
                            return (val)
                            }
                          }
                        )
  
  diff.025.mean <- apply(diff.95.mean, MARGIN = 1, 
                         function(x){
                           if (is.integer(B*0.025)) return (sort(x, decreasing = FALSE) [ceiling(B*0.025)])
                           else {
                             upper <- ceiling(B*0.025)
                             lower <- floor(B*0.025)
                             val <- sort(x, decreasing = FALSE)[upper]*abs(upper - B*0.025) + sort(x, decreasing = FALSE)[lower]* abs(B*0.025 - lower)
                             return (val)
                           }
                         }
  )
  
  ie.mean <- (actual.mean - cf.mean)/cf.mean
  ie.b.mean <- mclapply(cf.b.mean, function(x){return ((actual.mean - x)/x)})
  ie.95.mean <- array(unlist(ie.b.mean), dim = c(length(cf.mean),length(cf.b.mean)))
  #ie.95.mean <- apply(ie.95.mean, MARGIN = 1, function(x){return (sort(x, decreasing = FALSE) [ceiling(B*0.95)])})
  ie.975.mean <- apply(ie.95.mean, MARGIN = 1, 
                         function(x){
                           if (is.integer(B*0.975)) return (sort(x, decreasing = FALSE) [ceiling(B*0.975)])
                           else {
                             upper <- ceiling(B*0.975)
                             lower <- floor(B*0.975)
                             val <- sort(x, decreasing = FALSE)[upper]*abs(upper - B*0.975) + sort(x, decreasing = FALSE)[lower]* abs(B*0.975 - lower)
                             return (val)
                           }
                         }
  )
  
  ie.025.mean <- apply(ie.95.mean, MARGIN = 1, 
                         function(x){
                           if (is.integer(B*0.025)) return (sort(x, decreasing = FALSE) [ceiling(B*0.025)])
                           else {
                             upper <- ceiling(B*0.025)
                             lower <- floor(B*0.025)
                             val <- sort(x, decreasing = FALSE)[upper]*abs(upper - B*0.025) + sort(x, decreasing = FALSE)[lower]* abs(B*0.025 - lower)
                             return (val)
                           }
                         }
  )
  
  
  start.week <- time_conversion[time_conversion$E2week == 0,"Eweek_tot"] |> unlist() |> min()
  combine.value <- c(diff.mean[as.numeric(names(diff.mean)) < start.week], ie.mean[as.numeric(names(diff.mean)) >= start.week])
  combine.975 <- c(diff.975.mean[as.numeric(names(diff.mean)) < start.week], ie.975.mean[as.numeric(names(diff.mean)) >= start.week])
  combine.025 <- c(diff.025.mean[as.numeric(names(diff.mean)) < start.week], ie.025.mean[as.numeric(names(diff.mean)) >= start.week])
  
  avg.IE.table <- data.frame(actual = actual.mean, 
                             conterfactual= cf.mean, 
                             row.names = names(actual.mean))
  avg.IE.table$ie.lb <- combine.value + (combine.value - combine.975)
  avg.IE.table$ie.ub <- combine.value + (combine.value - combine.025)
  avg.IE.table$ie <- combine.value
  avg.IE.table$Eweek <-as.numeric(rownames(avg.IE.table))
  return (avg.IE.table)
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

#5. placebo test
placebo.test <- function(f_treatment, min.week, control, time.ind, space.ind, val.ind, q.type, load.q.func = NULL, lambda.ind = 1, beta){
  postal_c <- control %>% filter(!!sym(space.ind) != f_treatment) #rest of control sectors
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
  res.d <- sapply(c(1:ncol(cf.placebo)), function(x) {return (get_IE(x,cf.placebo,t_quant.placebo))})
  names(res.d) <- e2weeks
  return (res.d)
}

#6.1 leave-one-out sampling on each sector (space.ind) each 2 eweeks (time.ind) on a dataframe
getlos <- function(id, df, time.ind, space.ind, n){
  df <- df %>% filter(!!sym(space.ind) == id) %>% group_by(!!sym(time.ind)) %>% sample_n(n) %>% ungroup()  
  return (df)
}

#7.1 preprocessing data from getCI function for plotting
preprocessCI <- function(CIres, CIres.b, time.ind){
  actual.diff <- CIres[[1]]
  t_quant <- CIres[[6]]
  quantFun <- CIres[[7]]
  b.lambda <- CIres.b
  cf.b <- lapply(b.lambda, function(x){
    k <- lapply(c(1:length(quantFun)), function(y){
      return (as.matrix(cbind(quantFun[[y]][[2]])) %*% x)}
    )
    return (k)}) %>% lapply(function(x){
      return (array(unlist(x), dim = c(length(x[[1]]),length(x))))
    })
  bootstrap.diff <- lapply(cf.b, function(x) {return (t_quant - x)})
  b.diff <- lapply(bootstrap.diff, function(x){return (abs(actual.diff - x))})
  #Get conversion of eweek and e2week
  week.conversion <- CIres[[3]]
  treatment_length <- max(week.conversion[time.ind])
  record_start <- min(week.conversion[time.ind])
  #Get 95% largest abs diff across bootstrap
  b.diff.arr <- array(data = unlist(b.diff), dim = c(1001,(treatment_length - record_start + 1), B))
  b.diff.arr <- apply(b.diff.arr,c(1,2),function(x){return (sort(x, decreasing = TRUE) [1])})
  b.diff.arr <- apply(b.diff.arr,c(2), function(x){return(sort(x, decreasing = FALSE)[ceiling(0.95*1001)])})
  names(b.diff.arr) <- c(record_start:treatment_length)
  #b.diff.arr is a 1D matrix containing max abs diff over each quantile (row) and e2week (column)
  return (list(b.diff.arr, actual.diff, week.conversion,treatment_length,record_start))
}

#7 plot confidential band (quantile, GAI)
plotCI <- function(tid, preprocessed_CI, time, time.ind, interval.eweek = 2, save.path = "../Rimages/diff/"){
  width <- preprocessed_CI[[1]] #diff of diff
  syndf <- preprocessed_CI[[2]] # real diff
  time.conversion <- preprocessed_CI[[3]]
  width = width[which(names(width) == time)] #get diff variance for a e-2week
  synval = syndf[,which(colnames(syndf) == time)] #get diff for a e-2week
  synval <- data.frame(quantile = sapply(names(synval), rm.last), GAI = unname(synval)) #change percentage to decimal
  diff95 <- width
  time.min <- min(time.conversion[which(time.conversion[[time.ind]] == time),"Eweek_tot"])
  time.max <- max(time.conversion[which(time.conversion[[time.ind]] == time),"Eweek_tot"])
  if (time >= 0){
    lb.df <- synval
    lb.df["lb"] <- -lb.df$GAI-diff95
    lb.df["previous"] <- lead(lb.df$quantile)
    lb.df <- lb.df %>% filter(lb >0)
    lb.df["later"] <- lead(lb.df$quantile)
    lb.df["continuity"] <- mapply(all.equal,lb.df$later, lb.df$previous) != TRUE
    lb.df["breakpoint"] <- rev(cumsum(rev(lb.df$continuity)))
    lb.df <- lb.df %>% filter(breakpoint == 0)
    if (nrow(lb.df) > 0){
      lb.percent <- min(lb.df$quantile) 
    } else {
      lb.percent = NULL
    }
  }
    
  if (time >= 0){
    if (!is.null(lb.percent)){
      p = ggplot(synval, aes(quantile, -GAI)) +    
      geom_line() + 
      geom_ribbon(aes(ymax = -GAI+diff95, ymin = -GAI-diff95),alpha = 0.2) + # geom_ribbon function is used to add confidence interval 
      labs(x = "quantile", y = "GAI Diff") +
      geom_hline(yintercept=0, linetype="dashed", color = "red")+
      ggtitle(paste("Sector",tid,"at e-week",toString(time.min),'-',toString(time.max)), subtitle = paste("GAI difference. After",lb.percent,"has significant diff. ",interval.eweek,"wks combined"))

    } else {
      p = ggplot(synval, aes(quantile, -GAI)) +    
        geom_line() + 
        geom_ribbon(aes(ymax = -GAI+diff95, ymin = -GAI-diff95),alpha = 0.2) + # geom_ribbon function is used to add confidence interval 
        labs(x = "quantile", y = "GAI Diff") +
        geom_hline(yintercept=0, linetype="dashed", color = "red")+
        ggtitle(paste("Sector",tid,"at e-week",toString(time.min),'-',toString(time.max)), subtitle = paste("GAI difference. No significant diff. ",interval.eweek,"wks combined"))
    }
        
  } else {
    p = ggplot(synval, aes(quantile, -GAI)) +    
      geom_line() + 
      geom_ribbon(aes(ymax = -GAI+diff95, ymin = -GAI-diff95),alpha = 0.2) + # geom_ribbon function is used to add confidence interval 
      labs(x = "quantile", y = "GAI Diff") +
      geom_hline(yintercept=0, linetype="dashed", color = "red")+
      ggtitle(paste("Sector",tid,"at e-week",toString(time.min),'-',toString(time.max)), subtitle = paste("GAI difference between actual and conterfactual.",interval.eweek,"wks combined"))
    
  }
  #save the image of IE plot
  system(paste("mkdir -p",save.path))
  ggsave(paste(tid,"_",time,".png",sep = ""), plot = p, device = "png",
         width = 2100, height = 900, units = "px", path = save.path)
  return (NULL)
}


#8.1 convert in-space placebo lists to dataframe recognizable by ggplot function
build.placebo.df <- function(placebo.li, conversion, sectorID, before, after, new.xlabel){
  time <- as.numeric(names(placebo.li[[1]]))
  placebo.diff <- placebo.li[[sectorID]]
  ew_tot_min <- min(conversion[which(conversion[before] == min(time)),][[after]])
  ew_tot_max <- max(conversion[which(conversion[before] == max(time)),][[after]])
  diff_ew_tot <- sapply(c(ew_tot_min:ew_tot_max), function(x){
    IE_this_week = conversion[which(conversion[after] == x),][[before]] |> as.character()
    return (placebo.diff[which(names(placebo.diff) == IE_this_week)])
  })
  temp <- data.frame(eweek = c(ew_tot_min:ew_tot_max))
  temp['Sector_ID'] <- rep(sectorID, nrow(temp))
  temp['Diff'] <- diff_ew_tot
  
  x_break = seq(from = ew_tot_min, to = ew_tot_max, by = 26)
  if (!(183 %in% x_break)) {x_break = append(x_break,183)}
  x_label = sapply(x_break, function(x){return(conversion[which(conversion[after] == x),][[new.xlabel]])})
  return (list(temp,x_break,x_label,ew_tot_min, ew_tot_max))
}



#9 plot IE
plot_IE <- function(time, IE, IE.width, tid, conversion = NULL, tid.indicator = NULL, before = NULL, after = NULL, new.xlabel = NULL, interval.eweek = 2, save.path = "../Rimages/IE"){
  ew_tot_min <- min(conversion[which(conversion[before] == min(time)),][[after]])
  ew_tot_max <- max(conversion[which(conversion[before] == max(time)),][[after]])
  IE_ew_tot <- sapply(c(ew_tot_min:ew_tot_max), function(x){
    IE_this_week = conversion[which(conversion[after] == x),][[before]] |> as.character()
    return (IE[which(names(IE) == IE_this_week)])
  })
  x_break = seq(from = ew_tot_min, to = ew_tot_max, by = 13)
  if (!(183 %in% x_break)) {x_break = append(x_break,183)}
  x_label = sapply(x_break, function(x){return(conversion[which(conversion[after] == x),][[new.xlabel]])})
  IE.df <- data.frame(eweek = c(ew_tot_min:ew_tot_max), ie = -IE_ew_tot) #get positive number for positive treatment
  
  ie.width <- sapply(c(ew_tot_min:ew_tot_max), function(x){
    IE_this_week = conversion[which(conversion[after] == x),][[before]] |> as.character()
    return (IE.width[which(names(IE.width) == IE_this_week)])
  })
  treatment.start.ew = min(conversion[which(conversion[before] == 0),][[after]])
  
  p = ggplot(IE.df, aes(eweek, ie)) +    
    geom_line() +
    geom_ribbon(data = IE.df, aes(x = eweek, ymin = ie-ie.width, ymax = ie+ie.width),alpha = 0.2)+
    labs(x = "epidemiological week", y = "Intervention Efficacy") +
    geom_hline(yintercept = 1, linetype = "solid", color = "yellow", alpha = 0.7) +
    geom_vline(xintercept = treatment.start.ew, linetype="dashed", color = "red", alpha = 0.5)+
    coord_cartesian(ylim = c(-2,1.5)) +
    scale_x_continuous(breaks = x_break, labels = x_label) +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5)) +
    ggtitle(paste(tid,"intervention efficacy trend"), subtitle = paste("calculated by change ratio of 2-Wasserstein distance,", interval.eweek,"wks combined"))
  
  system(paste("mkdir -p",save.path))
  ggsave(paste(tid,"_IE.png",sep = ""), plot = p, device = "png",
         width = 2100, height = 900, units = "px", path = save.path)
  
  
  return (NULL)
}

#10 draw a plot to show GAI trend in quantile -- currently no use
combine.df <- function(tid, time1, ...){
  kwargs <- list(...)
  df.tot = data.frame(quantile = numeric(0), group = character(0), GAI = numeric(0))
  for (name in names(kwargs)){
    df <- data.frame(quantile = seq(0,1,length.out = 1001), group = rep(name,1001))
    temp.df <- kwargs[[which(names(kwargs) == name)]]
    df["GAI"] <- temp.df[,which(colnames(temp.df) == time1)]
    df.tot <- rbind(df.tot,df)
  }
  
  p <- ggplot(data = df.tot, mapping = aes(x = quantile, y = GAI, group = group, color = group)) +geom_line()
  
  return (p)
}


#11.1 getCf function special for IE permutation test
getCf.p <- function(id, t.info, t.s, c.s, time.ind, space.ind, no.trap, no.mosq, p.T0, lambda.ind = 1, quantFun.pre = NULL, typeq = NULL, space.agg = NULL, interval.eweek = 2, beta = FALSE, time.weight.ind = FALSE){
  T0 <- p.T0
  c.s["E2week"] <- floor((c.s[[time.ind]] - T0) / interval.eweek)
  t.s <- t.s %>% filter(!!sym(space.ind) == id)
  t.s["E2week"] <- floor((t.s[[time.ind]] - T0) / interval.eweek)
  if (!is.null(space.agg)){
    c.s <- c.s %>% group_by(E2week,!!sym(space.agg),!!sym(space.ind),!!sym(time.ind),Eyear,Eweek) %>% summarise(no.trap = sum(!!sym(no.trap)), no.mosq = sum(!!sym(no.mosq)))
    t.s <- t.s %>% group_by(E2week,!!sym(space.agg),!!sym(space.ind),!!sym(time.ind),Eyear,Eweek) %>% summarise(no.trap = sum(!!sym(no.trap)), no.mosq = sum(!!sym(no.mosq)))
    if (is.null(typeq)) typeq <- 7
  } else {
    c.s <- c.s %>% rename(c(no.trap = no.trap, no.mosq = no.mosq))
    t.s <- t.s %>% rename(c(no.trap = no.trap, no.mosq = no.mosq))
    if (is.null(typeq)) typeq <- 1
  }
  c.s["GAI"] <- c.s$no.mosq/c.s$no.trap
  t.s["GAI"] <- t.s$no.mosq/t.s$no.trap
  
  # get time representation (treatment start week = 0) and a vector of control sector ID
  e2weeks <- unlist(unique(t.s$E2week)) |> sort() |> as.vector()
  control.s <- unique(c.s$Sector_ID) |> sort.sector()
  
  #get E2week vs e-week table (convert E2week to readable e-week), will be used in plotting line chart
  convertweek <-t.s %>% select("E2week","Eweek_tot","Eyear","Eweek") %>% distinct()
  convertweek["eyew"] <- apply(convertweek %>% select("Eyear","Eweek") , MARGIN = 1, paste, collapse = '-')
  convertweek <- convertweek %>% select("E2week","Eweek_tot","eyew")
  E2week_n <- convertweek %>% filter(E2week < 0) %>% group_by(E2week) %>% summarize(n_wks = n()) %>% ungroup()
  
  if (time.weight.ind){
    valid_wks = E2week_n %>% filter(n_wks == interval.eweek) %>% select(E2week) %>% unlist() %>% unname()
    valid_wks = c(valid_wks, e2weeks[e2weeks >= 0])
    c.s <- c.s %>% filter(E2week %in% valid_wks)
    t.s <- t.s %>% filter(E2week %in% valid_wks)
    e2weeks <- unname(valid_wks) |> sort() |> as.vector()
  }
  
  if (is.null(quantFun.pre)){
    quantFun <- lapply(e2weeks, getQuantile, 
                       prob_list = seq(0,1,length.out = 1001), 
                       control = c.s,
                       treatment = t.s,
                       time.ind = "E2week",
                       space.ind = space.ind,
                       val.ind = "GAI",
                       q.type = typeq
    )
  } else {
    quantFun <- lapply(quantFun.pre, function(x){
      control.pre <- x[[2]]
      treat.pre <- x[1]
      t..id <- names(treat.pre)
      quantile.name <- names(treat.pre[[1]])
      treat.pre <- list(unname(treat.pre[[1]]))
      names(treat.pre) <- c(t..id)
      tot <- c(control.pre, treat.pre)
      treat.post <- tot[which(names(tot) == id)] |> unname()|> unlist() |> as.double()
      names(treat.post) <- quantile.name
      
      control.post <- tot[which(names(tot)!= id)] 
      sorting <- sort.sector(names(control.post))
      control.post <- control.post[sorting] |> as.data.frame()
      res.post <- list(treat.post,control.post)
      names(res.post) <- c(id,'others')
      return (res.post)
    })
  }
  
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
      return (mean(unlist(lapply(c(1:length(lambdas)),function(y){
        return (lambdas[[y]][x])
      }))))
    }))
  } else {
    f.lambdas <- getLambda2(time = c(min(t.s$E2week):-1), control = c.s, treatment = t.s, time.ind = "E2week", space.ind = space.ind, val.ind = "GAI", q.type = typeq, beta = beta)
  }
  
  
  names(f.lambdas) <- control.s
  
  #calculate diff
  cf <- do.call(cbind,lapply(c(1:length(quantFun)),
                             function(x){return(as.matrix(cbind(quantFun[[x]][[2]])) %*% f.lambdas)}
  ))
  t_quant <- do.call(cbind,lapply(c(1:length(quantFun)), function(x){return(quantFun[[x]][[1]])}))
  colnames(t_quant) <- e2weeks
  colnames(cf) <- e2weeks
  rownames(cf) <- rownames(t_quant)
  return (list(cf,t_quant))
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

#13 calculate average value of cf and calculate IE using formula in Somya's paper:
#calculate IE using average values of CF and actual one and calculate % change

average.qf <- function(qf, time_conversion, previous.ind, current.ind){
  row_names <- rownames(qf) |> sapply(rm.last)
  e2week.qf <- colnames(qf) |> sapply(as.numeric) |> unname()
  means <- apply(qf, MARGIN = 2, function(x){
    rand_num <- runif(2000) %>% round(digit = 3)
    sampled <- sapply(rand_num, function(i){return (x[which(row_names == i)])})
    return (mean(sampled))
  })
  names(means) <- e2week.qf
  ew.tot.min <- min(time_conversion[time_conversion[previous.ind] == min(e2week.qf),][[current.ind]])
  ew.tot.max <- max(time_conversion[time_conversion[previous.ind] == max(e2week.qf),][[current.ind]])
  means.ew.tot <- sapply(c(ew.tot.min:ew.tot.max),function(x){
    IE.this.week = time_conversion[time_conversion[current.ind] == x,][[previous.ind]] |> as.character()
    return (means[which(names(means) == IE.this.week)])
  })
  names(means.ew.tot) <- c(ew.tot.min:ew.tot.max)
  return (means.ew.tot)
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
#15 real space placebo test
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
