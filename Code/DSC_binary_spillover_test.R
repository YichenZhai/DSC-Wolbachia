#spillover test
#Load necessary libraries
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(CVXR, lib.loc = "/home/yichen.zhai/Rlibs"))
suppressPackageStartupMessages(library(osqp, lib.loc = "/home/yichen.zhai/Rlibs"))
suppressPackageStartupMessages(library(pbmcapply, lib.loc = "/home/yichen.zhai/Rlibs"))
suppressPackageStartupMessages(library(ECOSolveR, lib.loc = "/home/yichen.zhai/Rlibs"))
suppressPackageStartupMessages(library(frechet, lib.loc = "/home/yichen.zhai/Rlibs"))
suppressPackageStartupMessages(library(sendmailR, lib.loc = "/home/yichen.zhai/Rlibs"))
suppressPackageStartupMessages(library(msm, lib.loc = "/home/yichen.zhai/Rlibs"))
#suppressPackageStartupMessages(library(CEoptim, lib.loc = "/home/yichen.zhai/Rlibs"))
suppressPackageStartupMessages(library(ggplot2))
load("DSC_data3.RData")
source("DSC_classification.R")
options(dplyr.summarise.inform = FALSE)

#Get control sectors that are adjacent to intervened sectors.
get_control_sector <- function(tid, df, treat_status, control_status){
  
  bef.release.info = df %>% filter(Sector_ID == tid & type == control_status)
  aft.release.info = df %>% filter(Sector_ID == tid & type == treat_status)
  
  aft_end_eweek <- aft.release.info$end_on
  bef_start_eweek <- bef.release.info$Eweek_tot
  
  temp_df <- df %>% filter(type == control_status & Eweek_tot <= bef_start_eweek & end_on >= aft_end_eweek)
  
  return (temp_df$Sector_ID |> unlist())
}

treat_control <- function(treat_status, control_status, control_status2 = NULL){
  if (is.null(control_status2)){
    k <- brief.release.info %>% filter(type == control_status & next_stage == treat_status & (end_on - Eweek_tot + 1) >= 52)
    res <- lapply(k$Sector_ID, get_control_sector, df = brief.release.info, treat_status = treat_status, control_status = control_status)
    names(res) <- k$Sector_ID
    return (res)
  } else {
    k <- brief.release.info %>% filter(type == control_status & next_stage == treat_status & (end_on - Eweek_tot + 1) >= 52)
    res <- lapply(k$Sector_ID, get_control_sector, df = brief.release.info, treat_status = treat_status, control_status = control_status2)
    names(res) <- k$Sector_ID
    return (res)
  }
 
}
res <- treat_control("NA.","NN")
treatment_sectors = names(res)

B = 500
treatment_ids = treatment_sectors

wks = 13 #Integer. Number of weeks aggregated
space.agg = "Postal" #Char. NULL: trap-level data; "Postal": postal-level data aggregated by the same postal code.
l.ind = 2 #Integer. 1: calculate optimal weight by averaging optimal weights of pre-intervention time intervals;
          #         2: calculate optimal weight by optimizing the weight of the whole pre-intervention period.
beta.ind = TRUE #Boolean. TRUE: beta distribution random sampling; FALSE: uniform distribution random sampling in weight optimization.
q_type = 7 #Integer. 
time.weight.ind = TRUE # Boolean. TRUE: remove pre-treatment interval shorter than "wks"; FALSE: does not remove
trim.ind = TRUE # Boolean. TRUE: trim to 40 control sectors; FALSE: consider all control sectors

time.lag = c(26,52) / wks # half a year and a year before in-time placebo test
total_surveillance <- rbind.data.frame(treat_surveillance, control_surveillance)
treatment_info <- brief.release.info %>% filter(Sector_ID %in% treatment_ids & type == "NA.") %>% as.data.frame()

while (length(treatment_ids) != 0){
  x = treatment_ids[1]
  #Calculate : 1. start date of the sector being affected by spillover effect
  #            2. end date of the sector being affected by spillover effect (be treated or trial period end)
  #Do:         1. keep only the records when the sector is affected by spillover effect and not affected by any effect.
  #            2. check 1-year pre-treatment period
  #            3. check the record is consistent
  time.start <- brief.release.info %>% filter(Sector_ID == x & type == "NN") %>% select(Eweek_tot) %>% unlist() %>% unname() #1. start date of the sector being affected by spillover effect
  time.end  <- brief.release.info %>% filter(Sector_ID == x & type == "NA.") %>% select(end_on) %>% unlist() %>% unname() #2. end date of the sector being affected by spillover effect (be treated or trial period end)
  treat_surveillance = total_surveillance %>% filter(Eweek_tot >= time.start & Eweek_tot <= time.end & Sector_ID == x) #keep only the records when the sector is affected by spillover effect and not affected by any effect.
  treat_weeks <- treat_surveillance$Eweek_tot |> unique() |> unlist() |> length()
  record_start <- treat_surveillance$Eweek_tot |> unique() |> unlist() |> min()
  NN_end <- brief.release.info %>% filter(Sector_ID == x & type == "NN") %>% select(end_on) %>% unlist() %>% unname() #End date of pre-treatment period
  if (treat_weeks == (time.end - time.start + 1) | (record_start <= (NN_end - 51))){ #check 1-year pre-treatment period and check the record consistency
    time.start = record_start
    controls <- res[names(res) == x] |> unlist() |> unname()
    #keep records during the "trial"
    control_surveillance = total_surveillance %>% filter(Eweek_tot >= time.start & Eweek_tot <= time.end & Sector_ID %in% controls)
    control_selection <- mclapply(controls, function(c){
      temp_cs <- control_surveillance %>% filter(Sector_ID == c) %>% select(Eweek_tot) %>% unlist() %>% unique() %>% length()
      return (temp_cs == treat_weeks)
    }, mc.cores = 20) |> unlist()
    if (length(controls) != sum(control_selection)){
      controls <- controls[control_selection]
      control_surveillance = control_surveillance %>% filter(Sector_ID %in% controls)
    }
    print(paste("-------",x,"-------"))
    p <- sector_DSC(treatment_id = x, B = B, P=P, lambda.ind = l.ind, interval.eweek = wks, use.cores = 20, time.lag = time.lag, time.weight.ind = time.weight.ind, trim = trim.ind)
    system("clear")
  }
  treatment_ids <- treatment_ids[-1]
}

