#-----------------------Load library and auxilary functions--------------------#
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(CVXR, lib.loc = "/home/yichen.zhai/Rlibs"))
#suppressPackageStartupMessages(library(scs, lib.loc = "/home/yichen.zhai/Rlibs"))
suppressPackageStartupMessages(library(osqp, lib.loc = "/home/yichen.zhai/Rlibs"))
suppressPackageStartupMessages(library(pbmcapply, lib.loc = "/home/yichen.zhai/Rlibs"))
suppressPackageStartupMessages(library(ECOSolveR, lib.loc = "/home/yichen.zhai/Rlibs"))
suppressPackageStartupMessages(library(frechet, lib.loc = "/home/yichen.zhai/Rlibs"))
suppressPackageStartupMessages(library(msm, lib.loc = "/home/yichen.zhai/Rlibs"))
suppressPackageStartupMessages(library(ggplot2))
source("DSC_binary_sector.R")
options(dplyr.summarise.inform = FALSE)
#-------------------------------------------------------------------------------------#

# temporal aggregation = 13 E-weeks, spatial aggregation = Postal - Sector level, beta

#----------------------(For quicker running) load environment data -------------------#
load("DSC_data3.RData")

B = 500 #times of resampling for bootstrapping.
treatment_ids = treatment_sectors #a vector of sector ID chars

wks = 13 #Integer. Number of weeks aggregated
l.ind = 2 #Integer. 1: calculate optimal weight by averaging optimal weights of pre-intervention time intervals;
          #         2: calculate optimal weight by optimizing the weight of the whole pre-intervention period.
trim.ind = TRUE # Boolean. TRUE: trim to 40 control sectors; FALSE: consider all control sectors
time.weight.ind = TRUE # Boolean. TRUE: remove pre-treatment interval shorter than "wks"; FALSE: does not remove

time.lag = c(26,52) / wks # 0.5 year and 1 year before in-time placebo test

while (length(treatment_ids) != 0){ #Apply DSC on discrete variable to analyze the proportion change of traps with no mosquito detected.
  x = treatment_ids[1]
  print(paste("-------",x,"-------"))
  p <- sector_DSC(treatment_id = x, B= B, lambda.ind = l.ind, interval.eweek = wks, time.lag = time.lag, use.cores = 20, trim = trim.ind)
  system("clear")
  treatment_ids <- treatment_ids[-1]
}


