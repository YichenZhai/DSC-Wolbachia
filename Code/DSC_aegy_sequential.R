#-----------------------Load library and auxilary functions--------------------#
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(CVXR, lib.loc = "/home/yichen.zhai/Rlibs"))
suppressPackageStartupMessages(library(osqp, lib.loc = "/home/yichen.zhai/Rlibs"))
suppressPackageStartupMessages(library(pbmcapply, lib.loc = "/home/yichen.zhai/Rlibs"))
suppressPackageStartupMessages(library(ECOSolveR, lib.loc = "/home/yichen.zhai/Rlibs"))
suppressPackageStartupMessages(library(frechet, lib.loc = "/home/yichen.zhai/Rlibs"))
suppressPackageStartupMessages(library(sendmailR, lib.loc = "/home/yichen.zhai/Rlibs"))
suppressPackageStartupMessages(library(msm, lib.loc = "/home/yichen.zhai/Rlibs"))
suppressPackageStartupMessages(library(ggplot2))
source("DSC_aegy_sector.R")
options(dplyr.summarise.inform = FALSE)
#-------------------------------------------------------------------------------------#

# temporal aggregation = 13 E-weeks, spatial aggregation = Postal - Sector level, beta

#----------------------(For quicker running) load environment data -------------------#
load("DSC_data3.RData")
B = 500 # number of bootstrap repetition
P = 1000 # not in used

treatment_ids = treatment_sectors

wks = 13 #Integer. Number of weeks aggregated
space.agg = "Postal" #Char. NULL: trap-level data; "Postal": postal-level data aggregated by the same postal code.
l.ind = 2 #Integer. 1: calculate optimal weight by averaging optimal weights of pre-intervention time intervals;
          #         2: calculate optimal weight by optimizing the weight of the whole pre-intervention period.
beta.ind = TRUE #Boolean. TRUE: beta distribution random sampling; FALSE: uniform distribution random sampling in weight optimization.
q_type = 7 #Integer. 
time.weight.ind = TRUE # Boolean. TRUE: remove pre-treatment interval shorter than "wks"; FALSE: does not remove
trim.ind = TRUE # Boolean. TRUE: trim to 40 control sectors; FALSE: consider all control sectors

time.lag = c(26,52) / wks # 0.5 year and 1 year before in-time placebo test
treat_surveillance_original = treat_surveillance
control_surveillance_original = control_surveillance

while (length(treatment_ids) != 0){ #Apply DSC to analyze intervention effectiveness of the vector control technique.
    x = treatment_ids[1]
    print(paste("-------",x,"-------"))
    p <- sector_DSC(treatment_id = x, quantile_type = q_type, B= B, P= P,lambda.ind = l.ind, interval.eweek = wks, space.agg = space.agg, beta.ind = beta.ind, time.lag = time.lag, use.cores = 20, time.weight.ind = time.weight.ind, trim = trim.ind, treat_surveillance = treat_surveillance_original, control_surveillance = control_surveillance_original)
    system("clear")
    treatment_ids <- treatment_ids[-1]
}


