#-----------------------Load library and auxilary functions--------------------#
source("DSC_aegy_aux_func.R")

sector_DSC <- function(treatment_id, quantile_type, B, P,lambda.ind = 1, interval.eweek = 2, beta.ind = FALSE, use.cores = 20, space.agg = NULL, time.lag = c(13,26), time.weight.ind = FALSE, trim = FALSE, treat_surveillance = treat_surveillance, control_surveillance = control_surveillance){
  # Spatial aggregation of surveillance data by summation.
  if (!is.null(space.agg)){
    treat_surveillance <- treat_surveillance %>% group_by(!!sym(space.agg),Sector_ID,Eweek_tot,Eyear,Eweek) %>% summarise(Number.of.functional.Gravitrap = sum(Number.of.functional.Gravitrap), Number.of.wild.type.female.Ae..aegypti = sum(Number.of.wild.type.female.Ae..aegypti)) |> ungroup()
    control_surveillance <- control_surveillance %>% group_by(!!sym(space.agg),Sector_ID,Eweek_tot,Eyear,Eweek) %>% summarise(Number.of.functional.Gravitrap = sum(Number.of.functional.Gravitrap), Number.of.wild.type.female.Ae..aegypti = sum(Number.of.wild.type.female.Ae..aegypti)) |> ungroup()
    space.agg = NULL
  }
  
  curtime = Sys.time()
  treatment_id = treatment_id
  
  # Reduce the number of control sectors to 40. Only the most correlated 40 sectors will be selected.
  if (trim){
    selected.s <- trim_cs(id = treatment_id, 
                        t.info = treatment_info, 
                        t.s = treat_surveillance,
                        c.s = control_surveillance, 
                        time.ind = "Eweek_tot", 
                        space.ind = "Sector_ID", 
                        no.trap = "Number.of.functional.Gravitrap", 
                        no.mosq = "Number.of.wild.type.female.Ae..aegypti",
                        interval.eweek = interval.eweek,
                        cores = use.cores,
                        typeq = quantile_type)
    control_sectors <- names(selected.s)
    control_surveillance <- control_surveillance %>% filter(Sector_ID %in% names(selected.s))
    print("--------------trimmed to 40 control sectors---------------")
  }
  
  #remove surveillance records of control sectors earlier than the first record of treatment sector (e.g. FL350 has the first record on 2019-36)
  min_eweek <- treat_surveillance %>% filter(Sector_ID == treatment_id) %>% select(Eweek_tot) %>% min()
  if (min_eweek > 8){
    control_surveillance <- control_surveillance %>% filter(Eweek_tot >= min_eweek)
  }
  
  #calculate counterfactual outcomes, do placebo tests and permutation tests
  res <- getCf(id = treatment_id, 
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
               typeq = quantile_type, 
               space.agg = space.agg,
               interval.eweek = interval.eweek,
               beta = beta.ind,
               time.weight.ind = time.weight.ind,
               cores = use.cores)
  
  system("clear")
  time1 = Sys.time()
  print(paste("DSC get lambda and placebo test finish. Time cost:", round(difftime(time1,curtime, unit = "mins"),2),"mins"))
  
  #---------------------Bootstrapping---------------------------#
  print("Start bootstrapping")
  quantFun <- res[[7]]
  cf <- res[[5]]
  res.b <- pbmclapply(rep(treatment_id,B),getCf, 
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
                      typeq = quantile_type,
                      space.agg = space.agg,
                      cores = 1,
                      interval.eweek = interval.eweek,
                      beta = beta.ind,
                      time.weight.ind = time.weight.ind,
                      mc.cores = use.cores)
  res.b.lambda <- lapply(res.b, function(x) {return (x[[2]])})
  res.b.cf <- pbmclapply(res.b.lambda, function(z) {
    res.temp <- sapply(c(1:length(quantFun)),function(x){return(as.matrix(cbind(quantFun[[x]][[2]])) %*% z)})
    colnames(res.temp) <- colnames(cf)
    rownames(res.temp) <- rownames(cf)
    return (res.temp)
    },  mc.cores = use.cores)

  system('clear')
  time2 = Sys.time()
  print(paste(treatment_id,"bootstrap finished. Time cost:", round(difftime(time2,curtime, unit = "mins"),2), "mins"))

  #--------------IE calculation--------------------#
  #Here calculate intervention effectiveness of Wolbachia IIT-SIT in the intervened sector for every post-intervention time interval.
  t_quant <- res[[6]]
  e2weeks <- as.numeric(colnames(cf))
  
  IEs <- sapply(seq(1:length(e2weeks)), get_IE, df.1 = cf, df.2 = t_quant)
  names(IEs) <- e2weeks
  IEs.b <- mclapply(res.b.cf, function(x){
    colnames(x) = colnames(t_quant)
    rownames(x) = rownames(t_quant)
    ies.b <- sapply(seq(1:length(e2weeks)), get_IE, df.1 = x, df.2 = t_quant)
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
  
  IEs.was.eweek <- IEs |> convertdata(conversion = convertweek, bef.ind = "E2week", aft.ind = "Eweek_tot")
  
  system("clear")
  time4 = Sys.time()
  print(paste(treatment_id,"IE calculation and plots finished. Time cost:", round(difftime(time4, curtime, unit = "mins"),2), "mins"))
  
  #-------------Ending--------------#
  system("clear")
  #return following:
  #1.treatment id
  #2. final lambda used in giving weight
  #3. bootstrap result
  #4. IE permutation test p-value table
  #5. IE calculated by wasserstein dist
  f.lambdas <- res[[2]]
  if (lambda.ind == 1){
    system("mkdir -p ./RData")
    save(list = ls(all.names = TRUE), file = paste("./RData/",treatment_id,"_data.RData", sep = ''))
  } else {
    system("mkdir -p ./RData")
    save(list = ls(all.names = TRUE), file = paste("./RData/",treatment_id,"_data_2.RData", sep = ''))
  }
  in_space_placebo <- res[[10]]
  system("mkdir -p ./RDSData")
  saveRDS(res[[9]], paste("./timeplacebo/",treatment_id,"_time_placebo_plots.RDS", sep = ''))
  #saveRDS(placebo_table, paste("./permutation/",treatment_id,"_permutation_plots.RDS", sep = ''))
  saveRDS(in_space_placebo,paste("./spaceplacebo/",treatment_id,"_space_placebo_plots.RDS", sep = ''))
  return (c(treatment_id, f.lambdas, res.b, placebo_table, IEs.was.eweek))
}