source("DSC_binary_aux_func.R")

sector_DSC <- function(treatment_id, B, lambda.ind = 1, interval.eweek = 2, use.cores = 20, time.lag = c(13,26), time.weight.ind = FALSE, trim = FALSE){
  curtime = Sys.time()
  if (trim){
    selected.s <- trim_cs_binary(tid = treatment_id, 
                                 t.info = treatment_info, 
                                 t.s = treat_surveillance,
                                 c.s = control_surveillance, 
                                 time.ind = "Eweek_tot", 
                                 space.ind = "Sector_ID", 
                                 no.trap = "Number.of.functional.Gravitrap", 
                                 no.mosq = "Number.of.wild.type.female.Ae..aegypti",
                                 interval.eweek = interval.eweek,
                                 cores = use.cores)
    control_sectors <- names(selected.s)
    control_surveillance <- control_surveillance %>% filter(Sector_ID %in% names(selected.s))
    print("Trimmed to 40 control sectors.")
    }
  
  #keep control sector earliest time to the same as treatment sector. remove unused records beyond trial period.
  #(e.g. FL350 has the first record on 2019-36, so records before 2019-36 in control sectors are removed)
  min_eweek <- treat_surveillance %>% filter(Sector_ID == treatment_id) %>% select(Eweek_tot) %>% min()
  if (min_eweek > 8){
    control_surveillance <- control_surveillance %>% filter(Eweek_tot >= min_eweek)
  }
  
  res <- getCf.binary(tid = treatment_id, 
                      t.info = treatment_info, 
                      t.s = treat_surveillance,
                      c.s = control_surveillance, 
                      time.ind = "Eweek_tot", 
                      space.ind = "Sector_ID", 
                      no.trap = "Number.of.functional.Gravitrap", 
                      no.mosq = "Number.of.wild.type.female.Ae..aegypti",
                      lambda.ind = lambda.ind,
                      bootstrap = FALSE, 
                      space.placebo = TRUE, 
                      time.placebo = TRUE,
                      time.lag = time.lag,
                      lag.T0 = NULL,
                      interval.eweek = interval.eweek,
                      time.weight.ind = time.weight.ind,
                      cores = use.cores)
  
  system("clear")
  time1 = Sys.time()
  print(paste("DSC get lambda and placebo test finish. Time cost:", round(difftime(time1,curtime, unit = "mins"),2),"mins"))
  
  #---------------------Bootstrapping---------------------------#
  print("Start bootstrapping")
  t_quant <- res[[6]]
  quantFun <- res[[7]]
  res.b <- pbmclapply(rep(treatment_id,B),getCf.binary, 
                      t.info = treatment_info, 
                      t.s = treat_surveillance,
                      c.s = control_surveillance, 
                      time.ind = "Eweek_tot", 
                      space.ind = "Sector_ID", 
                      no.trap = "Number.of.functional.Gravitrap", 
                      no.mosq = "Number.of.wild.type.female.Ae..aegypti", 
                      lambda.ind = lambda.ind,
                      bootstrap = TRUE, 
                      space.placebo = FALSE,
                      time.placebo = FALSE,
                      time.lag = NULL,
                      lag.T0 = NULL,
                      cores = 1,
                      interval.eweek = interval.eweek,
                      time.weight.ind = time.weight.ind,
                      mc.cores = use.cores)
  
  res.b.lambda <- res.b
  res.b.cf <- lapply(res.b.lambda, function(z) {return (sapply(c(1:length(quantFun)),function(x){return(sapply(quantFun[[x]][[2]], function(y){return(y(0))}) %*% z)}))}) |> unlist() |> array(dim = c(length(quantFun),B))
  #res.b.diff <- lapply(res.b.cf, function(x){return (t_quant - x)})
  system('clear')
  time2 = Sys.time()
  print(paste(treatment_id,"bootstrap finished. Time cost:", round(difftime(time2,curtime, unit = "mins"),2), "mins"))
  
  #-------------save permutation test results-----------------#
  convertweek <- res[[3]]
  placebo <- res[[4]]
  placebo_table <- process.placebo(placebo.res = placebo,
                                   t.id = treatment_id,
                                   time.ind ="E2week",
                                   alpha.num = 0.1,
                                   save.path = "../Rimages/space_placebo/",
                                   t.info = treatment_info,
                                   conversion = convertweek,
                                   before = "E2week",
                                   after = "Eweek_tot",
                                   new.xlabel = "eyew",
                                   interval.eweek = interval.eweek
  )
  system("clear")
  time3 = Sys.time()

  #--------------IE calculation from ratio change to baseline by JS divergence --------------------#
  cf <- res[[5]]
  t_quant <- res[[6]]
  e2weeks <- as.numeric(names(cf))
  
  IEs <- get_IE.p(cf, t_quant)
  names(IEs) <- e2weeks
  IEs.b <- mclapply(res.b.cf, function(x){
    ies.b <- get_IE.p(x, t_quant)
    names(ies.b) <- e2weeks
    return (ies.b)
  }, mc.cores = use.cores)
  
  #IEs.b.diff <- array(unlist(IEs.b), dim = c(length(IEs.b[[1]]),length(IEs.b))) |> abs() |> apply(MARGIN = 1, function(x){return (sort(x, decreasing = FALSE) [ceiling(B*0.95)])})
  IEs.diff <- mclapply(IEs.b, function(x){
    ies.diff <- abs(x-IEs)
    return (ies.diff)
  }, mc.cores = use.cores)
  IEs.diff.arr <- array(unlist(IEs.diff), dim = c(length(IEs),B))
  IEs.diff.arr <- apply(IEs.diff.arr, MARGIN = 1, function(x){return (sort(x, decreasing = FALSE) [ceiling(B*0.95)])})
  names(IEs.diff.arr) = names(IEs)
  
  plot_IE(time = e2weeks, IE = IEs, IE.width = IEs.diff.arr, tid = treatment_id, conversion = convertweek,tid.indicator = "Sector_ID", before = "E2week", after = "Eweek_tot", new.xlabel = "eyew", interval.eweek = interval.eweek)
  IEs.was.eweek <- IEs |> convertdata(conversion = convertweek, bef.ind = "E2week", aft.ind = "Eweek_tot")
  
  system("clear")
  time4 = Sys.time()
  print(paste(treatment_id,"IE calculation and plots finished. Time cost:", round(difftime(time4, curtime, unit = "mins"),2), "mins"))

  
  #-------------1-Wasserstein distance for pre-treatment period---------------#
  e2weeks <- as.numeric(e2weeks)
  pre_distance = pbmclapply(res.b.cf, function(x){return(get_was(prob1 = x, prob2 = t_quant))}, mc.cores = use.cores)
  pre_distance <- pre_distance |> unlist() |> array(dim = c(sum(e2weeks<0),B))
  
  pre_dist_treat <- get_was(prob1 = cf[e2weeks < 0], prob2 = t_quant[e2weeks < 0])
  pre_dist_diff <- apply(pre_distance, MARGIN = 2,  function(x) {return(x - pre_dist_treat)}) |> apply(MARGIN = 1, function(x) return(sort(abs(x), decreasing = FALSE) [ceiling(B*0.95)]))
  pre_dist_df <- data.frame(treat = pre_dist_treat, lb = pre_dist_treat - pre_dist_diff, ub = pre_dist_treat + pre_dist_diff, Eweek = e2weeks[e2weeks < 0])
  pre_dist_plot <- ggplot(data = pre_dist_df, aes(x = Eweek)) +
    geom_line(aes(y = treat))+
    #geom_ribbon(aes(ymin = lb, ymax = ub),alpha = 0.2, linetype=0) +
    coord_cartesian(ylim = c(0,NA)) +
    scale_x_continuous(breaks = e2weeks[e2weeks < 0], labels = e2weeks[e2weeks < 0])+
    labs(x = "Aggregated e-week", y = "Difference") +
    ggtitle(paste("Sector",treatment_id,"pre-treatment period differnce"), subtitle = paste("Difference by 1-wasserstein distance. Score = ",round(mean(pre_dist_treat),5)))
  
  ggsave(paste(treatment_id,"_distance_1was.png", sep = ""), plot = pre_dist_plot, device = "png",width = 2100, height = 900, units = "px", path = "../Rimages/")
  
  pre_distance2 = pbmclapply(res.b.cf, function(x){return(get_JS(prob1 = x, prob2 = t_quant))}, mc.cores = use.cores)
  pre_distance2 <- pre_distance2 |> unlist() |> array(dim = c(sum(e2weeks<0),B))
  
  pre_dist_treat2 <- get_JS(prob1 = cf[e2weeks < 0], prob2 = t_quant[e2weeks < 0])
  pre_dist_diff2 <- (pre_distance2 - pre_dist_treat2) |> apply(MARGIN = 1, function(x) return(sort(abs(x), decreasing = FALSE) [ceiling(B*0.95)]))
  pre_dist_df2 <- data.frame(treat = pre_dist_treat2, lb = pre_dist_treat2 - pre_dist_diff2, ub = pre_dist_treat2 + pre_dist_diff2, Eweek = e2weeks[e2weeks < 0])
  pre_dist_plot2 <- ggplot(data = pre_dist_df2, aes(x = Eweek)) +
    geom_line(aes(y = treat))+
    #geom_ribbon(aes(ymin = lb, ymax = ub),alpha = 0.2, linetype=0) +
    coord_cartesian(ylim = c(0,NA)) +
    scale_x_continuous(breaks = e2weeks[e2weeks < 0], labels = e2weeks[e2weeks < 0])+
    labs(x = "Aggregated e-week", y = "Difference") +
    ggtitle(paste("Sector",treatment_id,"pre-treatment period differnce"), subtitle = paste("Difference by JS divergence. Score = ",round(mean(pre_dist_treat),5)))
  
  ggsave(paste(treatment_id,"_distance_JS.png", sep = ""), plot = pre_dist_plot2, device = "png",width = 2100, height = 900, units = "px", path = "../Rimages/")
  
  
  time6 = Sys.time()
  #-------------Ending--------------#
  f.lambdas <- res[[2]]
  #save environmental variables
  if (lambda.ind == 1){
    system("mkdir -p ./RData")
    save(list = ls(all.names = TRUE), file = paste("./RData/",treatment_id,"_data.RData", sep = ''))
  } else if (lambda.ind == 2){
    system("mkdir -p ./RData")
    save(list = ls(all.names = TRUE), file = paste("./RData/",treatment_id,"_data_2.RData", sep = ''))
  } else if (lambda.ind == 3){
    system("mkdir -p ./RData")
    save(list = ls(all.names = TRUE), file = paste("./RData/",treatment_id,"_data_3.RData", sep = ''))
  }
  
  #save three robustness check data of the sector.
  saveRDS(placebo_table, file = paste("./permutation/",treatment_id,"_permutation_plot.RDS", sep = ''))
  saveRDS(res[[9]], file = paste("./timeplacebo/",treatment_id,"_time_placebo_plot.RDS", sep = ''))
  saveRDS(res[[10]], file = paste("./spaceplacebo/",treatment_id,"_space_placebo_plot.RDS", sep = ''))
}