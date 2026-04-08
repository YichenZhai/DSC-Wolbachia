#Figure 2
#Load necessary libraries
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(pbmcapply, lib.loc = "/home/yichen.zhai/Rlibs"))
suppressPackageStartupMessages(library(frechet, lib.loc = "/home/yichen.zhai/Rlibs"))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(RColorBrewer))
options(dplyr.summarise.inform = FALSE)
source("DSC_aegy_aux_func.R")
load("DSC_data3.RData")

#----Figure 2A/2C
#Set working directory to the database
setwd("RData/")

#Expected to obtain only env vars of interested sectors calculated in DSCM.
m = list.files(pattern = "*.RData") 

#Extract only IE and raw difference among different quantiles and time.
res <- pbmclapply(m, function(x){
  load(x)
  actual.diff <- -res[[1]][,e2weeks >= 0]
  cf <- res[[5]][,e2weeks >= 0] |> array(dim = c(1001, sum(e2weeks >=0)), dimnames = list(rownames(res[[5]]), e2weeks[e2weeks>=0]))
  IE.quantile <- mapply(function(x,y){
    if (y == 0 & x == 0){ return(0)}
    if (y == 0 & x < 0) {return(-10)}
    return(x/y)
  }, actual.diff,cf) |> array(dim = c(1001,sum(e2weeks>=0)), dimnames = list(rownames(cf),colnames(cf)))
  
  actual.diff <- actual.diff |> array(dim = c(1001,sum(e2weeks>=0)), dimnames = list(rownames(res[[1]]),e2weeks[e2weeks>=0]))
  eweeks <- colnames(actual.diff) |> as.numeric()
  if (max(eweeks) < 9){
    NA.follow <- 9 - max(eweeks)
    df <- array(rep(NA, nrow(actual.diff) * NA.follow), dim = c(nrow(actual.diff),NA.follow)) |> data.frame()
    colnames(df) <- c((max(eweeks)+1):9)
    rownames(df) <- rownames(actual.diff)
    actual.diff <- cbind(actual.diff,df)
  }
  
  if (max(eweeks) < 9){
    NA.follow <- 9 - max(eweeks)
    df <- array(rep(NA, nrow(IE.quantile) * NA.follow), dim = c(nrow(IE.quantile),NA.follow)) |> data.frame()
    colnames(df) <- c((max(eweeks)+1):9)
    rownames(df) <- rownames(IE.quantile)
    IE.quantile <- cbind(IE.quantile,df)
  }
  
  return (list(actual.diff,IE.quantile))
}, mc.cores = 30)

#Number of sectors in interest, please manually change this if the folder contains other files.
sector_num = length(m)

#convert to 3-D array (1001 quantiles, 10 post-treatment intervals and number of sector in interest)
diffs = lapply(res, function(x){return(x[[1]])}) |> unlist() |> array(dim = c(1001,10,sector_num)) 
IEs = lapply(res, function(x){return(x[[2]])}) |> unlist() |> array(dim = c(1001,10,sector_num))

#calculate the average IE for each time and quantile
avg_diff = apply(diffs, MARGIN = c(1,2), mean, na.rm = TRUE)
avg_IE = apply(IEs, MARGIN = c(1,2), mean, na.rm = TRUE)

#Get pointwise confidence interval (CI) of mean IE/raw diff for each time and quantile by bootstrapping
boot_res = pbmclapply(c(1:1000), function(b){
  boots = sample(c(1:sector_num),size = sector_num, replace = TRUE)
  newdiffs = diffs[,,boots]
  newIE = IEs[,,boots]
  
  newdiffs_mean = apply(newdiffs, MARGIN = c(1,2), mean, na.rm = TRUE)
  newIE_mean = apply(newIE, MARGIN = c(1,2), mean, na.rm = TRUE)
  
  return(list(newdiffs_mean, newIE_mean))
}, mc.cores = 14)

boot_diff = lapply(boot_res, function(x){return(x[[1]])}) |> unlist() |> array(dim = c(1001,10,1000))
boot_diff_ub = apply(boot_diff, MARGIN = c(1,2), quantile, probs = 0.975, na.rm = TRUE)
boot_diff_lb = apply(boot_diff, MARGIN = c(1,2), quantile, probs = 0.025, na.rm = TRUE)

boot_IE = lapply(boot_res, function(x){return(x[[2]])}) |> unlist() |> array(dim = c(1001,10,1000))
boot_IE_ub = apply(boot_IE, MARGIN = c(1,2), quantile, probs = 0.975, na.rm = TRUE)
boot_IE_lb = apply(boot_IE, MARGIN = c(1,2), quantile, probs = 0.025, na.rm = TRUE)


# IE dataframe. Columns: mean IE, upper bond of pointwise CI, lower bond of pointwise CI
IE_avg <- as.data.frame(IE_avg)
colnames(IE_avg) <- c("1-13","14-26","27-39","40-52","53-65","66-78","79-91","92-104","105-117","118-124")
IE_avg$number = seq(0,1,by = 0.001)
total_IE_agg = pivot_longer(IE_avg, cols = !number, names_to = "wks_expand", values_to = "IE")

boot_IE_ub <- as.data.frame(boot_IE_ub)
colnames(boot_IE_ub) <- c("1-13","14-26","27-39","40-52","53-65","66-78","79-91","92-104","105-117","118-124")
boot_IE_ub$number = seq(0,1,by = 0.001)
df_IE_ub = pivot_longer(boot_IE_ub, cols = !number, names_to = "wks_expand", values_to = "ub")

boot_IE_lb <- as.data.frame(boot_IE_lb)
colnames(boot_IE_lb) <- c("1-13","14-26","27-39","40-52","53-65","66-78","79-91","92-104","105-117","118-124")
boot_IE_lb$number = seq(0,1,by = 0.001)
df_IE_lb = pivot_longer(boot_IE_lb, cols = !number, names_to = "wks_expand", values_to = "lb")

total_IE_agg = total_IE_agg %>% inner_join(df_IE_lb,by = join_by(number, wks_expand)) %>% inner_join(df_IE_ub, by = join_by(number,wks_expand))

total_IE_agg = total_IE_agg[order(total_IE_agg$wks_expand, total_IE_agg$number),]

color_palette = brewer.pal(10,"Spectral")

give_expand_color <- function(x){
  m = data.frame(period = c("1-13","14-26","27-39","40-52","53-65","66-78","79-91","92-104","105-117","118-124"),
                 color = color_palette)
  return(m$color[m$period == x])
}

total_IE_agg$ind = mapply(function(lb,ub,wks_expand){
  return(ifelse(lb >= 0 , give_expand_color(wks_expand), "#D3D3D3"))
}, total_IE_agg$lb,total_IE_agg$ub,total_IE_agg$wks_expand) |> factor(levels = c(color_palette, "#D3D3D3"))

total_IE_agg <- total_IE_agg %>% group_by(color_group = with(rle(ind |> as.character()), rep(seq_along(values), lengths))) %>% ungroup()

common_levels = setdiff(levels(total_IE_agg$ind) , "#D3D3D3")
#Plot the graph
Figure2A = ggplot(data = total_IE_agg , aes(x = number, y = IE, group = wks_expand, color = wks_expand)) +
  coord_cartesian(ylim = c(-11,1)) +
  
  geom_ribbon(data = total_IE_agg %>% filter(ind != "#D3D3D3"), aes(x = number, ymin = lb, ymax = ub, fill = wks_expand, group = color_group), alpha = 0.2, color = NA) +
  geom_ribbon(data = total_IE_agg %>% filter(ind == "#D3D3D3"), aes(x = number, ymin = lb, ymax = ub, group = color_group), fill = "#D3D3D3", alpha = 0.1, color = NA) +
  geom_line() +
  scale_color_manual(breaks = c("1-13","14-26","27-39","40-52","53-65","66-78","79-91","92-104","105-117","118-124"),
                     values = color_palette,
                     labels = c("1-13","14-26","27-39","40-52","53-65","66-78","79-91","92-104","105-117","118-124")) +
  scale_fill_manual(breaks = c("1-13","14-26","27-39","40-52","53-65","66-78","79-91","92-104","105-117","118-124"),
                    values = color_palette,
                    labels = c("1-13","14-26","27-39","40-52","53-65","66-78","79-91","92-104","105-117","118-124")) +
  labs(y = "IE,aegypti(%)", x = "GAI Quantile", fill = "Event time (week)", color = "Event time (week)") +
  scale_y_continuous(breaks = seq(-11,1,by = 2), labels = seq(-11,1,by = 2)*100) +
  scale_x_continuous(breaks = seq(0,1,by = .25), labels = sapply(seq(0,1,by = .25), function(x){return(paste(x*100,"%", sep = ''))})) +
  geom_hline(yintercept = 0, linetype = 3, linewidth = 0.3) +
  theme(legend.title = element_text(size = rel(0.8), face = "bold"),
        legend.key.spacing.y = unit(0.1, "cm"),
        legend.key.size = unit(0.1,'cm'),
        plot.margin = margin(0,2,0,2),
        legend.background = element_blank(),
        #legend.position.inside = c(0.4,0.75),
        legend.text = element_text(size = 6),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.line.x = element_line(color = 'black', linewidth = 0.5, linetype = 'solid'),
        axis.line.y = element_line(color = 'black', linewidth = 0.5, linetype = 'solid'),
        panel.background = element_blank(),
        panel.border = element_rect(color = 'black',fill = NA, linewidth = 0.5))

#-------------Figure 2C, raw reduction in aegypti
#diff dataframe
diff_avg <- as.data.frame(diff_avg)
colnames(diff_avg) <- c("1-13","14-26","27-39","40-52","53-65","66-78","79-91","92-104","105-117","118-124")
diff_avg$number = seq(0,1,by = 0.001)
total_diff_agg = pivot_longer(diff_avg, cols = !number, names_to = "wks_expand", values_to = "diff")

boot_diff_ub <- as.data.frame(boot_diff_ub)
colnames(boot_diff_ub) <- c("1-13","14-26","27-39","40-52","53-65","66-78","79-91","92-104","105-117","118-124")
boot_diff_ub$number = seq(0,1,by = 0.001)
df_IE_ub = pivot_longer(boot_diff_ub, cols = !number, names_to = "wks_expand", values_to = "ub")

boot_diff_lb <- as.data.frame(boot_diff_lb)
colnames(boot_diff_lb) <- c("1-13","14-26","27-39","40-52","53-65","66-78","79-91","92-104","105-117","118-124")
boot_diff_lb$number = seq(0,1,by = 0.001)
df_IE_lb = pivot_longer(boot_diff_lb, cols = !number, names_to = "wks_expand", values_to = "lb")

total_diff_agg = total_diff_agg %>% inner_join(df_IE_lb,by = join_by(number, wks_expand)) %>% inner_join(df_IE_ub, by = join_by(number,wks_expand))

total_diff_agg = total_diff_agg[order(total_diff_agg$wks_expand, total_diff_agg$number),]

color_palette = brewer.pal(10,"Spectral")

give_expand_color <- function(x){
  m = data.frame(period = c("1-13","14-26","27-39","40-52","53-65","66-78","79-91","92-104","105-117","118-124"),
                 color = color_palette)
  return(m$color[m$period == x])
}

total_diff_agg$ind = mapply(function(lb,ub,wks_expand){
  return(ifelse(lb >= 0 | ub <= 0, give_expand_color(wks_expand), "#D3D3D3"))
}, total_diff_agg$lb,total_diff_agg$ub,total_diff_agg$wks_expand) |> factor(levels = c(color_palette, "#D3D3D3"))

total_diff_agg <- total_diff_agg %>% group_by(color_group = with(rle(ind |> as.character()), rep(seq_along(values), lengths))) %>% ungroup()

common_levels = setdiff(levels(total_diff_agg$ind) , "#D3D3D3")
Figure2C = ggplot(data = total_diff_agg , aes(x = number, y = diff, group = wks_expand, color = wks_expand)) +
  geom_line() +
  coord_cartesian(ylim = c(-0.75,0.25)) +
  scale_color_manual(breaks = c("1-13","14-26","27-39","40-52","53-65","66-78","79-91","92-104","105-117","118-124"),
                     values = color_palette,
                     labels = c("1-13","14-26","27-39","40-52","53-65","66-78","79-91","92-104","105-117","118-124")) +
  geom_ribbon(data = total_diff_agg %>% filter(ind != "#D3D3D3"), aes(x = number, ymin = lb, ymax = ub, fill = wks_expand, group = color_group), alpha = 0.2, color = NA) +
  geom_ribbon(data = total_diff_agg %>% filter(ind == "#D3D3D3"), aes(x = number, ymin = lb, ymax = ub, group = color_group), fill = "#D3D3D3", alpha = 0.1, color = NA) +
  scale_fill_manual(breaks = c("1-13","14-26","27-39","40-52","53-65","66-78","79-91","92-104","105-117","118-124"),
                    values = color_palette,
                    labels = c("1-13","14-26","27-39","40-52","53-65","66-78","79-91","92-104","105-117","118-124")) +
  labs(y = "Raw reduction in aegypti", x = "GAI Quantile", fill = "Event time (week)", color = "Event time (week)") +
  scale_y_continuous(breaks = seq(0,1.2,by = .3), labels = c("0.00","0.30","0.60","0.90","1.20")) +
  coord_cartesian(ylim = c(-0.1,1.2))+
  scale_x_continuous(breaks = seq(0,1,by = .25), labels = sapply(seq(0,1,by = .25), function(x){return(paste(x*100,"%", sep = ''))})) +
  geom_hline(yintercept = 0, linetype = 3, linewidth = 0.3) +
  theme(legend.title = element_text(size = rel(0.8), face = "bold"),
        legend.key.spacing.y = unit(0.1, "cm"),
        legend.key.size = unit(0.1,'cm'),
        plot.margin = margin(0,2,0,2),
        legend.background = element_blank(),
        legend.position= c(0.4,0.75),
        legend.text = element_text(size = 6),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.line.x = element_line(color = 'black', linewidth = 0.5, linetype = 'solid'),
        axis.line.y = element_line(color = 'black', linewidth = 0.5, linetype = 'solid'),
        axis.title.x = element_text(size = rel(0.9)),
        panel.background = element_blank(),
        panel.border = element_rect(color = 'black',fill = NA, linewidth = 0.5)) +
  guides(color = guide_legend(ncol = 5, byrow = FALSE, override.aes = list(linetype = 1,size = 5, alpha = 0.2)))


#----Figure 2B/2D is just Figure 2A/2C but running with spillover effect related data
# 1. Change the working directory to the folder that store environmental variables of "Spillover" sectors.
# 2. Change the size of the last dimension (if needed). There should be 91 treated sectors and 40 sectors adjacent to treated one in consideration
# 3. Change image parameters

#------SI Figure 1A
treatment.sectors = treatment_sectors
B = 1000
b <- function(df){
  sample_df <- df %>% t() %>% as.data.frame() %>% slice_sample(prop = 1, replace = TRUE) %>% t() %>% apply(MARGIN = 1, mean)
  return(sample_df)
}

f <- function(tid){
  load(paste("RData/",tid,"_data_2.RData", sep = ''))
  #print(paste("---------",tid,'start ---------'))
  actual.diff <- -res[[1]][,e2weeks >= 0] 
  actual.diff <- actual.diff |> array(dim = c(1001,sum(e2weeks>=0)), dimnames = list(rownames(res[[1]]),e2weeks[e2weeks>=0]))
  
  actual.diff.wks <- apply(actual.diff, MARGIN = 1, convertdata, conversion = convertweek, bef.ind = "E2week", aft.ind = "Eweek_tot")
  actual.diff.wks <-  do.call(rbind,actual.diff.wks)
  
  eweeks <- colnames(actual.diff.wks) |> as.numeric()
  if (min(eweeks) > 60){
    NA.lead <- min(eweeks) - 60
    df <- array(rep(NA, nrow(actual.diff.wks) * NA.lead), dim = c(nrow(actual.diff.wks),NA.lead)) |> data.frame()
    colnames(df) <- c(60:(min(eweeks) -1))
    rownames(df) <- rownames(actual.diff.wks)
    actual.diff.wks <- cbind(df, actual.diff.wks)
  }
  
  eweeks <- colnames(actual.diff) |> as.numeric()
  #max_complete_week <- convertweek %>% filter(E2week >0) %>% group_by(E2week) %>% summarize(week_length = n()) %>% ungroup() %>% filter(week_length == 13) %>% select(E2week) %>% max(0)
  #eweeks <- max_complete_week
  if (max(eweeks) < 9){
    NA.follow <- 9 - max(eweeks)
    df <- array(rep(NA, nrow(actual.diff) * NA.follow), dim = c(nrow(actual.diff),NA.follow)) |> data.frame()
    colnames(df) <- c((max(eweeks)+1):9)
    rownames(df) <- rownames(actual.diff)
    actual.diff <- cbind(actual.diff,df)
  }
  
  #print(paste("---------",tid,'end ---------'))
  return (list(actual.diff.wks,actual.diff))
}

df.tot <- pbmclapply(treatment.sectors, f, mc.cores = 30)
names(df.tot) <- treatment.sectors
df.tot.wks <- lapply(df.tot, function(x){return(x[[1]])})
df.tot2 <- array(df.tot.wks |> unlist(), dim = c(1001,183-60+1,91), dimnames = list(rownames(df.tot.wks[[1]]), colnames(df.tot.wks[[1]]), treatment_sectors))
total_diff_wks <- apply(df.tot2, MARGIN = c(1,2), mean,na.rm = TRUE) |> as.data.frame()
sector.num <- apply(df.tot2, MARGIN = c(1,2), function(x){return(sum(!is.na(x)))}) |> apply(MARGIN = 2,mean)
total_diff_wks['quantile'] <- rownames(total_diff_wks)
total_diff_wks.longer <- pivot_longer(total_diff_wks, !quantile, names_to = "wks", names_transform = list("wks" = as.numeric), values_to = "diff") 
total_diff_wks.longer["number"] <- sapply(total_diff_wks.longer$quantile, function(x){return(substring(x,1,nchar(x)-1)|> as.numeric())})
total_diff_wks.longer <- total_diff_wks.longer[order(total_diff_wks.longer$wks, total_diff_wks.longer$number),]

lb = c()
ub = c()
for (wk in c(60:183)){
  print(paste("---------- week",wk,"start----------"))
  temp_sample <- df.tot2[,which(as.numeric(dimnames(df.tot2)[[2]]) == wk),]
  temp_sample <- temp_sample[,apply(temp_sample, MARGIN=2, function(x){return(!all(is.na(x)))})]
  actual_value <- temp_sample |> apply(MARGIN = 1, mean)
  bootstrap_value <- pbmclapply(replicate(1000, temp_sample, simplify = FALSE), b, mc.cores = 25)
  b.array <- bootstrap_value |> unlist() |> array(dim = c(1001,B))
  b.025 <- b.array |> apply(MARGIN = 1, function(x) {return(sort(x, decreasing = FALSE)[B*0.025])})
  b.975 <- b.array |> apply(MARGIN = 1, function(x) {return(sort(x, decreasing = FALSE)[B*0.975])})
  
  lb.now = 2 * actual_value - b.975
  ub.now = 2 * actual_value - b.025
  
  lb = c(lb,lb.now)
  ub = c(ub,ub.now)
  
}
total_diff_wks.longer["lb"] = lb
total_diff_wks.longer["ub"] = ub

classify <- function(x, number.li, label.li, label.li2){ # change old axis label to new axis label
  pointer = 1
  while(pointer < length(number.li)){
    pointer.end <- pointer + 1
    start <- number.li[pointer]
    end <- number.li[pointer.end]
    start.label <- label.li[pointer]
    end.label  <- label.li2[pointer]
    if (start <= x & x < end){
      return(paste(start.label,"-",end.label,sep = ''))
    }
    pointer = pointer + 1
  }
  return(paste(start.label,"-",end.label,sep = ''))
}
total_diff_quantile_group = sapply(total_diff_wks.longer$number, classify, number.li = c(0,25,50,75,100), 
                                   label.li = c("0%","25%","50%","75%"), label.li2 = c("24.9%","49.9%","74.9%","100%"))

total_diff_wks_group_26wk = sapply(total_diff_wks.longer$wks, classify, number.li = c(60,79,106,132,158,183), 
                                   label.li = c("EW8","EW27","EW1","EW27","EW1"),
                                   label.li2 = c("EW26,2020","EW53,2020","EW26,2021","EW52,2021","EW26,2022"))

total_diff_wks.longer["quantile_group"] = total_diff_quantile_group
total_diff_wks.longer["wks_group"] = total_diff_wks_group_26wk

total_diff_wks.longer$wks_group <- factor(total_diff_wks.longer$wks_group, levels = c("EW8-EW26,2020","EW27-EW53,2020","EW1-EW26,2021","EW27-EW52,2021","EW1-EW26,2022"))
diff.plot.box <- ggplot(data = total_diff_wks.longer, aes(x = quantile_group, y = diff, fill = wks_group)) + 
  coord_cartesian(ylim = c(0,0.8)) + 
  geom_boxplot(lwd = 0.3, fatten = 0.3,outlier.shape = NA) +
  labs(x = "GAI Quantile", y = "Raw Reduction in aegypti")+
  scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8)) +
  theme(legend.position = "right",
        legend.title = element_text(size=rel(0.45), face="bold"),
        legend.key.spacing.y = unit(0.2, 'cm'),
        legend.text = element_text(size=6), 
        plot.subtitle = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(colour = 'black', linewidth = 0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', linewidth = 0.5, linetype='solid'),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5))+
  geom_hline(yintercept = 0, linetype = 3, linewidth = 0.3) +
  labs(fill="Week interval") +
  guides(fill = guide_legend(ncol=5,byrow=TRUE),
         color = guide_legend(ncol=5,byrow=TRUE)) 

SI1A <- diff.plot.box

#------SI Figure 1B, graph-related parameters can be adjusted.
f <- function(tid){
  load(paste("RData/",tid,"_data_2.RData", sep = ''))
  actual.diff <- -res[[1]][,e2weeks >= 0]
  cf <- res[[5]][,e2weeks >= 0] |> array(dim = c(1001, sum(e2weeks >=0)), dimnames = list(rownames(res[[5]]), e2weeks[e2weeks>=0]))
  IE.quantile <- mapply(function(x,y){
    if (y == 0 & x == 0){ return(0)}
    if (y == 0 & x < 0) {return(-10)}
    return(x/y)
  }, actual.diff,cf) |> array(dim = c(1001,sum(e2weeks>=0)), dimnames = list(rownames(cf),colnames(cf)))
  
  IE.quantile.wks <- apply(IE.quantile, MARGIN = 1, convertdata, conversion = convertweek, bef.ind = "E2week", aft.ind = "Eweek_tot")
  IE.quantile.wks <-  do.call(rbind,IE.quantile.wks)
  
  eweeks <- colnames(IE.quantile.wks) |> as.numeric()
  if (min(eweeks) > 60){
    NA.lead <- min(eweeks) - 60
    df <- array(rep(NA, nrow(IE.quantile.wks) * NA.lead), dim = c(nrow(IE.quantile.wks),NA.lead)) |> data.frame()
    colnames(df) <- c(60:(min(eweeks) -1))
    rownames(df) <- rownames(IE.quantile.wks)
    IE.quantile.wks <- cbind(df, IE.quantile.wks)
  }
  
  eweeks <- colnames(IE.quantile) |> as.numeric()
  
  if (max(eweeks) < 9){
    NA.follow <- 9 - max(eweeks)
    df <- array(rep(NA, nrow(IE.quantile) * NA.follow), dim = c(nrow(IE.quantile),NA.follow)) |> data.frame()
    colnames(df) <- c((max(eweeks)+1):9)
    rownames(df) <- rownames(IE.quantile)
    IE.quantile <- cbind(IE.quantile,df)
  }
  
  #print(paste("---------",tid,'end ---------'))
  return (list(IE.quantile.wks,IE.quantile))
}

df.tot <- pbmclapply(treatment.sectors, f, mc.cores = 30)
names(df.tot) <- treatment.sectors
df.tot.wks <- lapply(df.tot, function(x){return(x[[1]])})
df.tot2 <- array(df.tot.wks |> unlist(), dim = c(1001,183-60+1,91), dimnames = list(rownames(df.tot.wks[[1]]), colnames(df.tot.wks[[1]]), treatment_sectors))
total_IE_wks <- apply(df.tot2, MARGIN = c(1,2), mean,na.rm = TRUE) |> as.data.frame()
sector.num <- apply(df.tot2, MARGIN = c(1,2), function(x){return(sum(!is.na(x)))}) |> apply(MARGIN = 2,mean)
total_IE_wks['quantile'] <- rownames(total_IE_wks)
total_IE_wks.longer <- pivot_longer(total_IE_wks, !quantile, names_to = "wks", names_transform = list("wks" = as.numeric), values_to = "IE") 
total_IE_wks.longer["number"] <- sapply(total_IE_wks.longer$quantile, function(x){return(substring(x,1,nchar(x)-1)|> as.numeric())})
total_IE_wks.longer <- total_IE_wks.longer[order(total_IE_wks.longer$wks, total_IE_wks.longer$number),]

lb = c()
ub = c()
for (wk in c(60:183)){
  print(paste("---------- week",wk,"start----------"))
  temp_sample <- df.tot2[,which(as.numeric(dimnames(df.tot2)[[2]]) == wk),]
  temp_sample <- temp_sample[,apply(temp_sample, MARGIN=2, function(x){return(!all(is.na(x)))})]
  actual_value <- temp_sample |> apply(MARGIN = 1, mean)
  bootstrap_value <- pbmclapply(replicate(1000, temp_sample, simplify = FALSE), b, mc.cores = 25)
  b.array <- bootstrap_value |> unlist() |> array(dim = c(1001,B))
  b.025 <- b.array |> apply(MARGIN = 1, function(x) {return(sort(x, decreasing = FALSE)[B*0.025])})
  b.975 <- b.array |> apply(MARGIN = 1, function(x) {return(sort(x, decreasing = FALSE)[B*0.975])})
  
  lb.now = 2 * actual_value - b.975
  ub.now = 2 * actual_value - b.025
  
  lb = c(lb,lb.now)
  ub = c(ub,ub.now)
  
}
total_IE_wks.longer["lb"] = lb
total_IE_wks.longer["ub"] = ub

total_IE_quantile_group = sapply(total_IE_wks.longer$number, classify, number.li = c(0,25,50,75,100), 
                                 label.li = c("0%","25%","50%","75%"), label.li2 = c("24.9%","49.9%","74.9%","100%"))
total_IE_wks_group = sapply(total_IE_wks.longer$wks, classify, number.li = c(60,79,106,132,158,183), 
                            label.li = c("EW8","EW27","EW1","EW27","EW1"),
                            label.li2 = c("EW26,2020","EW53,2020","EW26,2021","EW52,2021","EW26,2022"))

total_IE_wks.longer["quantile_group"] = total_IE_quantile_group
total_IE_wks.longer["wks_group"] = total_IE_wks_group

bad_names <- total_IE_wks.longer$wks_group[total_IE_wks.longer$lb < 0] |> unique()

total_IE_wks.longer$wks_group <- factor(total_IE_wks.longer$wks_group, levels = c("EW8-EW26,2020","EW27-EW53,2020","EW1-EW26,2021","EW27-EW52,2021","EW1-EW26,2022"))

group_names <- levels(total_IE_wks.longer$wks_group)
iqr.low = -5
ie.plot.box <- ggplot(data = total_IE_wks.longer, aes(x = quantile_group, y = IE, fill = wks_group)) + 
  geom_boxplot(lwd = 0.3, fatten = 0.3,outlier.shape = NA) +
  labs(x = NULL, y = "IE, aegypti(%)")+
  scale_y_continuous(breaks = seq(iqr.low * 100, 100, by = 25)/100, labels = seq(iqr.low * 100, 100, by = 25) |> as.character()) +
  coord_cartesian(ylim=c(iqr.low,1)) +
  labs(color = NULL, fill = NULL) +
  theme(plot.subtitle=element_text(size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(colour = 'black', linewidth = 0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', linewidth = 0.5, linetype='solid'),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5)) +
  labs(fill=NA) +
  geom_hline(yintercept = 0, linetype = 3, linewidth = 0.3) 
#ggtitle(paste("Intervention Efficacy for Ae. aegypti"), subtitle = "Calendar time. Treatment effect on 4 quantile ranges")


s <- ggplot_build(ie.plot.box)
label_x <- s$data[[1]]$x
x_need_arrow <- label_x[s$data[[1]]$ymin < iqr.low] |> as.numeric()
y_canvas_min <- s$layout$panel_params[[1]]$y.range[1]
t <- ie.plot.box$data$wks_group |> levels()
groups <- rep(t,length(x_need_arrow))[s$data[[1]]$ymin < iqr.low]
arrow_df <- data.frame(X = x_need_arrow, Y = rep(iqr.low,length(x_need_arrow)),Yend = rep(y_canvas_min,length(x_need_arrow)), wks_group = groups)
ie.plot.box <- ie.plot.box + geom_segment(data = arrow_df, arrow = arrow(length=unit(.2, 'cm'), type = "closed"), aes(x = X, y = Y, yend = Yend),linewidth = 0.6)

