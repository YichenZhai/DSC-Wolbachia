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
treatment.sectors <- treatment_sectors
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

#---- Figure 3
#Step 1: Figure 3 map IE data preprocessing
#Aggregate based on calendar time interval

tot.res <- pbmclapply(treatment_sectors, function(id){
  load(paste(id,"_data_2.RData", sep = ''))
  
  #from Eweek 79 to 105 -> 2020 2nd half
  #from Eweek 106 to 131 -> 2021 1st half
  #from Eweek 132 to 157 -> 2021 2nd half
  #from Eweek 158 to 183 -> 2022 1st half
  all_IEs <- rep(NA, 183-79+1)
  names(all_IEs) <- seq(79,183)
  
  post_wks <- convertweek %>% filter(E2week >= 0) %>% select(Eweek_tot) %>% unlist()
  get_IEs <- IEs.was.eweek[as.numeric(names(IEs.was.eweek)) >= 79 & as.numeric(names(IEs.was.eweek)) <= 183] |> unlist()
  get_IEs <- get_IEs[names(get_IEs) %in% post_wks]
  if (max(post_wks) < 79){
    return(list(NA,NA,NA,NA))
  }
  min_eweek <- names(get_IEs) |> as.numeric() |> min()
  max_eweek <- names(get_IEs) |> as.numeric() |> max()
  ind1 <- which(names(all_IEs) |> as.numeric() == min_eweek)
  ind2 <- which(names(all_IEs) |> as.numeric() == max_eweek)
  
  all_IEs[ind1:ind2] = -get_IEs                         
  yr1 <- mean(all_IEs[1:27], na.rm = TRUE)
  yr2 <- mean(all_IEs[28:53], na.rm = TRUE)
  yr3 <- mean(all_IEs[54:79], na.rm = TRUE)
  yr4 <- mean(all_IEs[80:105], na.rm = TRUE)
  print(paste(id," end", sep = ''))
  return(list(yr1,yr2,yr3,yr4))
},mc.cores = 31)

tot.res.df <- unlist(tot.res) |> matrix(nrow = length(treatment_sectors), ncol = 4)
tot.res.df$Sector_ID <- treatment_sectors

colnames(IEs) <- c("2020","2021_1","2021_2","2022","Sector_ID")
IEs_longer <- pivot_longer(data = IEs, cols = !Sector_ID, names_to = "Event week", values_to = "IE") 

IEs_longer$IE2 <- sapply(IEs_longer$IE, function(x){ #number "2" is a symbol of no data in study sectors
  if (is.null(x) | is.na(x)){
    return(2)
  } else {
    return (x)
  }
})

# Step 2: Combine geographical data to the corresponding sectors
setwd("../mapCodeFile/")
possible_group <- IEs_longer$`Event week` |> unique()
map_info <- st_read("original_sector_boundary/Sectorboundary_7Nov2022.shp")
map_info <- st_set_crs(map_info, 3414)
map_info <- st_transform(map_info, 4326)
IEs_longer2 <- inner_join(IEs_longer, map_info, by = "Sector_ID")
if ("2020" %in% possible_group){
  for (i in 1:length(possible_group)){
    group = possible_group[i]
    id = c("2020","2021_1","2021_2","2022")[i]
    #id = c(0:4)[i]
    data <- IEs_longer2 %>% filter(`Event week` == group) %>% select(!`Event week`)
    st_write(data, paste("sector_map_IE_shp/binary_IE_",id,".shp",sep = ''), append=FALSE)
  }
  
} else {
  for (i in 1:length(possible_group)){
    group = possible_group[i]
    #id = c("2020","2021_1","2021_2","2022")[i]
    id = c(0:4)[i]
    data <- IEs_longer2 %>% filter(`Event week` == group) %>% select(!`Event week`)
    st_write(data, paste("sector_map_IE_shp/albo_treat_IE_",id,".shp",sep = ''), append=FALSE)
  }
}
# End of data preparation for this map. These two steps connect the IEs of treated sectors to the location of treated sectors by constructing a dataframe.
# The remaining step will not be shown. 
# A possible procedure could be: 1. using QGIS to plot the map and assign different colors to different values of IE
#                                2. cropping the map to leave the areas in interest

#-----Figure 4
#-----Figure 4-1
# get post-intervention IE for treated sectors (stored in "treatment_sectors"), and calculate average mosquito abundance of each treated sector during pre-treatment
res <- pbmclapply(treatment_sectors,function(id){
  load(paste("RData/",id,"_data_2.RData", sep = ''))
  
  post.IEs <- IEs[as.numeric(names(IEs)) >= 0] #get post-intervention IE for sector
  
  alldist <- apply(t_quant, MARGIN = 2, function(x){ #get mosquito abundance for each time interval
    xdf <- data.frame(x = seq(0,1,by = 0.001), y = x)
    perfect_distribution <- data.frame(x = seq(0,1,by = 0.001), y = rep(0,1001))
    dist <- dist4den(xdf,perfect_distribution, fctn_type = "quantile")
    return(dist)
  })
  
  alldist <- alldist[as.numeric(colnames(t_quant)) <= -1] |> mean() #only consider pre-treatment period and get the historical mosquito abundance
  
  if (length(post.IEs) < 10){ #we expect there are totally 10 post-intervention intervals, but not all sectors had experienced such a long treatment (staggered adoption)
    addNA <- rep(NA, 10-length(post.IEs))
    post.IEs <- c(post.IEs, addNA)
    names(post.IEs) <- c(0:9)
    
  }
  
  return(c(post.IEs, alldist))
  
}, mc.cores = 31)

allres <- do.call(rbind,res) |> as.data.frame()

allres_longer <- pivot_longer(allres, !c(V11), names_to = "wks", values_to = "IE")
allres_longer2 <- allres_longer[c(which(!is.na(allres_longer$IE))),] %>% group_by(wks) %>% mutate(r  = rank(V11, na.last = TRUE), nrow = n()) %>% ungroup()
allres_longer2$rank <- (allres_longer2$r - 1) / (allres_longer2$nrow - 1)
res <- allres_longer2
res$IE <- -res$IE
res$r <- NULL
res$nrow <- NULL
colnames(res) <- c("dist_ori","groupname","IE","dist")

k <- ggplot(data = res) +
  geom_point(aes(x = dist,y = IE, group = groupname, color = groupname)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(colour = 'black', linewidth = 0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', linewidth = 0.5, linetype='solid'),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5),
        legend.title = element_text(size = rel(0.8), face = "bold"),
        legend.margin = margin(0,0,0,0),
        legend.key.spacing = unit(0.1,"cm"),
        legend.text = element_text(size = rel(0.8)),
        plot.margin = margin(0,0,0,2)) +
  geom_smooth(method = mgcv::gam, formula = y~x, mapping = aes(x = dist, y = IE, colour = groupname), se = FALSE, linewidth = 0.3) +
  labs(y = "IE (%)", x = "Sector GAI quantile", color = "Event time (week)") +
  scale_x_continuous(breaks = seq(0,1,by = 0.25), labels = sapply(seq(0,1,by = 0.25), function(x){
    return(paste(x*100,"%",sep = ''))
  })) +
  scale_y_continuous(breaks = seq(0,1,by = 0.25), labels = seq(0,1,by = 0.25) * 100) +
  scale_color_discrete(breaks = c(0:9),labels = c("1-13","14-26","27-39","40-52","53-65","66-78","79-91","92-104","105-117","118-124"
  ))

Figure4A <- k
#saveRDS(Figure4A, "Figure4_1.RDS")

#-----Figure 4-2
# Step 1: Aggregate IE based on townships 

aegypti_ie_13 <- pbmclapply(treatment_sectors, function(id) {
  load(paste(id,"_data_2.RData",sep = ''))
  post.IE <- IEs[as.numeric(names(IEs)) >= 0] 
  if (length(post.IE) < 10){
    addNA <- rep(NA, 10-length(post.IE))
    post.IE <- c(post.IE, addNA)
    names(post.IE) <- c(0:9)
  }
  names(post.IE) <- c("1-13","14-26","27-39","40-52","53-65","66-78","79-91","92-104","105-117","118-124")
  return (post.IE)
}, mc.cores = 31)
cont_IE_agg_df <- do.call(rbind, aegypti_ie_13) |> as.data.frame()
cont_IE_agg_df$Sector_ID <- treatment_sectors

a <- readxl::read_excel("../raw_data/Release_matrix_2016_2022ew26.xlsx")
b <- a %>% select(Sector_ID, StudyArea) %>% unique()
join_df <- left_join(cont_IE_agg_df,b,by = "Sector_ID", 
                     copy = FALSE)
cdf <- join_df %>% group_by(StudyArea) %>% summarize(a = mean(`1-13`, na.rm = TRUE),
                                                     b = mean(`14-26`, na.rm = TRUE),
                                                     c = mean(`27-39`, na.rm = TRUE),
                                                     d = mean(`40-52`, na.rm = TRUE),
                                                     e = mean(`53-65`, na.rm = TRUE),
                                                     f = mean(`66-78`, na.rm = TRUE),
                                                     g = mean(`79-91`, na.rm = TRUE),
                                                     h = mean(`92-104`, na.rm = TRUE),
                                                     i = mean(`105-117`, na.rm = TRUE),
                                                     j = mean(`118-124`, na.rm = TRUE))


# Step 2: Calculate average IE of the treated sectors based on the township.
B = 1000

b <- function(df){
  res <- apply(df, MARGIN = 2, function(x){
    na.value <- x %>% is.na()
    x <- x[!na.value]
    return(sample(x, size = length(x), replace = TRUE))
  })
  return(res)
}

lb <- c()
ub <- c()
studyarea <- join_df$StudyArea |> unique()
for (i in studyarea){
  total_cont_IE <- cdf[cdf$StudyArea == i,2:10]
  cont_IE_agg_df = join_df %>% filter(StudyArea == i)
  cont_IE_agg_df = cont_IE_agg_df[,1:9] #last post-intervention interval is removed, since it contains weeks fewer than 13 weeks.
  bootstrap_A <- lapply(replicate(1000, cont_IE_agg_df[1:9], simplify = FALSE), b) #bootstrap to get pointwise confidence intervals

  bootstrap_b2 <- lapply(bootstrap_A, function(x){return(sapply(x,mean))}) #calculate average IE
  bootstrap_b3 <- bootstrap_b2 |> unlist() |> array(dim = c(9, 1000))
  b.025.b <- bootstrap_b3 |> apply(MARGIN = 1, function(x) {return(sort(x, decreasing = FALSE)[B*0.025])}) #lower bound
  b.975.b <- bootstrap_b3 |> apply(MARGIN = 1, function(x) {return(sort(x, decreasing = FALSE)[B*0.975])}) #upper bound
  
  lb.b =  b.025.b
  ub.b =  b.975.b
  
  lb <- c(lb, lb.b)
  ub <- c(ub, ub.b)
}

# Step 3: Generate the plot
ribbon_df <- cbind(unlist(lb),unlist(ub)) |> as.data.frame()
colnames(ribbon_df) <- c("lb","ub")

ribbon_df$StudyArea <- c(rep("Bukit Batok",9),rep("Choa Chu Kang",9), rep("Tampines",9), rep("Yishun",9))
ribbon_df$wks <- rep(c("a","b","c","d","e","f","g","h","i"),4)

x_break2 = c("a","b","c","d","e","f","g","h","i")
cdf <- cdf[1:10] %>% pivot_longer(cols = !StudyArea, names_to = "wks", values_to = "IE")
ie.plot3 <- ggplot() + 
  geom_line(data = cdf,aes(x = wks,y = -IE, group = StudyArea,  color = StudyArea)) + 
  geom_ribbon(data = ribbon_df, aes(x = wks, ymin = -lb, ymax = -ub, group = StudyArea, fill = StudyArea), alpha = 0.2)+
  geom_hline(yintercept = 1, linetype = "dashed", color = "darkgreen", alpha = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkgrey", alpha = 1) +
  labs(x = "post-event week interval", y = "IE(%)")+
  scale_x_discrete(breaks = x_break2, labels = c("1-13","14-26","27-39","40-52","53-65","66-78","79-91","92-104","105-117"))+
  scale_y_continuous(breaks = c(-6,-4,-2,0,1), labels= c(-6,-4,-2,0,1)*100) +
  coord_cartesian(ylim = c(-7.5,1)) +
  theme(#axis.text.x = element_text(angle = 45, vjust = 0.7, hjust=0.5), 
    plot.subtitle = element_text(size = rel(0.5)),
    legend.title = element_blank(),
    axis.line.x = element_line(colour = 'black', linewidth = 0.5, linetype='solid'),
    axis.line.y = element_line(colour = 'black', linewidth = 0.5, linetype='solid'),
    panel.background = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5),
    axis.title.x = element_text(margin = margin(0,0,0,0)),
    legend.margin = margin(0, 0, 0, 0),
    legend.position = c(0.8,0.3),
    legend.background = element_blank()) 

Figure4B <- ie.plot3

#------ Figure 5
f <- function(tid){
  load(paste("RData/",tid,"_data_2.RData", sep = ''))
  #print(paste("---------",tid,'start ---------'))
  diffs <- -res[[1]]
  diffs.original <- diffs
  diffs <- diffs.original[e2weeks >= 0]
  diffs.wks <- convertdata(diffs, conversion = convertweek, bef.ind = "E2week", aft.ind = "Eweek_tot")
  
  eweeks <- names(diffs.wks) |> as.numeric()
  if (min(eweeks) > 60){
    NA.lead <- min(eweeks) - 60
    df <- rep(NA, NA.lead)
    names(df) <- c(60:(min(eweeks) -1))
    diffs.wks <- c(df, diffs.wks)
  }
  
  eweeks <- names(diffs) |> as.numeric()
  
  if (max(eweeks) < 9){
    NA.follow <- 9 - max(eweeks)
    df <- rep(NA, NA.follow)
    names(df) <- c((max(eweeks)+1):9)
    diffs <- c(diffs,df)
  }
  return (list(diffs.wks,diffs))
}

df.tot <- pbmclapply(treatment.sectors, f, mc.cores = 30)
names(df.tot) <- treatment.sectors

#Figure 5A
df.tot.agg <- lapply(df.tot, function(x){return(x[[2]])})

df.agg.tot2 <- array(df.tot.agg |> unlist(), dim = c(9-0+1,91), dimnames = list(names(df.tot.agg[[1]]), treatment_sectors))
total_diff_agg <- df.agg.tot2 |> as.data.frame()

total_diff_agg['wks'] <- rownames(total_diff_agg)
total_diff_agg.longer <- pivot_longer(total_diff_agg, !wks, names_to = "Sector_ID", values_to = "diff") 


total_diff_wks_group2 = sapply(total_diff_agg.longer$wks, classify, number.li = c(0,4,8,9), 
                               label.li = c("1","53","105"),
                               label.li2 = c("52","104","124"))
total_diff_wks_group3 = sapply(total_diff_agg.longer$wks, classify, number.li = c(0:10), 
                               label.li = c("1","14","27","40","53","66","79","92","105","118"),
                               label.li2 = c("13","26","39","52","65","78","91","104","117","124"))

total_diff_agg.longer["wks_group"] = total_diff_wks_group2
total_diff_agg.longer["wks_expand"] = total_diff_wks_group3

total_diff_agg.longer$wks_expand <- factor(total_diff_agg.longer$wks_expand, levels = c("1-13","14-26","27-39","40-52","53-65","66-78","79-91","92-104","105-117","118-124"))
diff.plot2 <- ggplot(data = total_diff_agg.longer, aes(x = wks_expand,y = diff, group = wks_expand, fill = wks_expand)) + 
  geom_boxplot(lwd = 0.3, fatten = 0.3,outlier.shape = NA) +
  labs(x = "Event time past (week)", y = "Raw proportion increment")+
  scale_y_continuous(breaks = c(-0.25,0,0.25,0.5,0.75,1), labels = c("-25%","0%", "25%","50%","75%","100%")) +
  coord_cartesian(ylim = c(-0.25,0.75)) +
  guides(fill = "none", color ="none")+
  theme(legend.title = element_text(size = 8), 
        legend.text = element_text(size=6), 
        plot.subtitle=element_text(size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(colour = 'black', linewidth = 0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', linewidth = 0.5, linetype='solid'),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5)) +
  geom_hline(yintercept = 0, linetype = 3, linewidth = 0.3) +
  ggtitle(paste("Raw proportion increment on 0 female Ae. Aegypti mosquito trap"), subtitle = "Event time. Treatment sectors considered.")

Figure5A <- diff.plot2

#Figure 5B
df.tot.wks <- lapply(df.tot, function(x){return(x[[1]])})
df.tot.agg <- lapply(df.tot, function(x){return(x[[2]])})
df.tot2 <- array(df.tot.wks |> unlist(), dim = c(183-60+1,91), dimnames = list(names(df.tot.wks[[1]]), treatment_sectors))
total_diff_wks <- as.data.frame(df.tot2)
total_diff_wks['wks'] <- rownames(total_diff_wks) |> as.numeric()
total_diff_wks.longer <- pivot_longer(total_diff_wks,!wks, names_to = "Sector_ID", values_to = "diff") 

total_diff_wks_group = sapply(total_diff_wks.longer$wks, classify, number.li = c(60,79,106,132,158,183), 
                              label.li = c("2020EW8","2020EW27","2021EW1","2021EW27","2022EW1"),
                              label.li2 = c("2020EW26","2020EW53","2021EW26","2021EW52","2022EW26"))

total_diff_wks.longer["wks_group"] = total_diff_wks_group

total_diff_wks.longer$wks_group <- factor(total_diff_wks.longer$wks_group, levels = c("2020EW8-2020EW26","2020EW27-2020EW53","2021EW1-2021EW26","2021EW27-2021EW52","2022EW1-2022EW26"))
diff.plot.box <- ggplot(data = total_diff_wks.longer, aes(x = wks_group, y = diff, fill = wks_group)) + 
  geom_boxplot(lwd = 0.3, fatten = 0.3,outlier.shape = NA) +
  labs(x = "Calendar time interval", y = "Raw proportion increment")+
  scale_y_continuous(breaks = c(-0.25,0,0.25,0.5,0.75,1), labels = c("-25%","0%", "25%","50%","75%","100%")) +
  coord_cartesian(ylim=c(-0.25,0.75)) +
  theme(legend.title = element_text(size = 8), 
        legend.text = element_text(size=6), 
        plot.subtitle=element_text(size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(colour = 'black', linewidth = 0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', linewidth = 0.5, linetype='solid'),
        panel.background = element_blank(),
        axis.text.x = element_text(size = rel(0.7)),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5)) +
  geom_hline(yintercept = 0, linetype = 3, linewidth = 0.3) +
  guides(fill = "none") +
  ggtitle(paste("Raw proportion increment on 0 female Ae.aegypti mosquito trap"), subtitle = "Calendar time. Treatment sectors considered.")
Figure5B <- diff.plot.box

#Figure 5C
#Use the spillover effect dataset!
f2 <- function(tid){
  load(paste("RData/",tid,"_data_2.RData", sep = ''))
  #print(paste("---------",tid,'start ---------'))
  diffs <- -res[[1]][e2weeks >= 0]
  #cf <- cf[e2weeks >= 0] |> array(dim = c(1001, sum(e2weeks >=0)), dimnames = list(rownames(res[[5]]), e2weeks[e2weeks>=0]))
  diffs.original <- diffs
  diffs.wks <- convertdata(diffs, conversion = convertweek, bef.ind = "E2week", aft.ind = "Eweek_tot")
  
  eweeks <- names(diffs.wks) |> as.numeric()
  if (min(eweeks) > 60){
    NA.lead <- min(eweeks) - 60
    df <- rep(NA, NA.lead)
    names(df) <- c(60:(min(eweeks) -1))
    diffs.wks <- c(df, diffs.wks)
  }
  if (max(eweeks) < 183){
    NA.follow <- 183-max(eweeks)
    df <- rep(NA, NA.follow)
    names(df) <- c((max(eweeks) +1):183)
    diffs.wks <- c(diffs.wks,df)
  }
  
  eweeks <- names(diffs) |> as.numeric()
  
  if (max(eweeks) < 9){
    NA.follow <- 9 - max(eweeks)
    df <- rep(NA, NA.follow)
    names(df) <- c((max(eweeks)+1):9)
    diffs <- c(diffs,df)
  }
  
  #print(paste("---------",tid,'end ---------'))
  return (list(diffs.wks,diffs))
}

df.tot <- pbmclapply(treatment.sectors, f2, mc.cores = 30)
names(df.tot) <- treatment.sectors

df.tot.agg <- lapply(df.tot, function(x){return(x[[2]])})

df.agg.tot2 <- array(df.tot.agg |> unlist(), dim = c(9-0+1,length(treatment.sectors)), dimnames = list(names(df.tot.agg[[1]]), treatment_sectors))
total_diff_agg <- df.agg.tot2 |> as.data.frame()

total_diff_agg['wks'] <- rownames(total_diff_agg)
total_diff_agg.longer <- pivot_longer(total_diff_agg, !wks, names_to = "Sector_ID", values_to = "diff") 


total_diff_wks_group2 = sapply(total_diff_agg.longer$wks, classify, number.li = c(0,4,8,9), 
                               label.li = c("1","53","105"),
                               label.li2 = c("52","104","124"))
total_diff_wks_group3 = sapply(total_diff_agg.longer$wks, classify, number.li = c(0:10), 
                               label.li = c("1","14","27","40","53","66","79","92","105","118"),
                               label.li2 = c("13","26","39","52","65","78","91","104","117","124"))

total_diff_agg.longer["wks_group"] = total_diff_wks_group2
total_diff_agg.longer["wks_expand"] = total_diff_wks_group3

total_diff_agg.longer$wks_expand <- factor(total_diff_agg.longer$wks_expand, levels = c("1-13","14-26","27-39","40-52","53-65","66-78","79-91","92-104","105-117","118-124"))
ie.plot2 <- ggplot(data = total_diff_agg.longer, aes(x = wks_expand,y = diff, group = wks_expand, fill = wks_expand)) + 
  geom_boxplot(lwd = 0.3, fatten = 0.3,outlier.shape = NA) +
  labs(x = "Event time past (week)", y = "Raw proportion increment")+
  scale_y_continuous(breaks = c(-0.2,0,0.2,0.4), labels = c("-20%","0%", "20%","40%")) +
  coord_cartesian(ylim = c(-0.2,0.4)) +
  guides(fill = "none", color ="none")+
  theme(legend.title = element_text(size = 8), 
        legend.text = element_text(size=6), 
        plot.subtitle=element_text(size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(colour = 'black', linewidth = 0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', linewidth = 0.5, linetype='solid'),
        panel.background = element_blank(),
        panel.border = element_rect(color = 'black', fill = NA, linewidth = 0.5)) +
  geom_hline(yintercept = 0, linetype = 3, linewidth = 0.3) +
  ggtitle(paste("Percentage increase for 0 female Ae.aegypti mosquito trap"), subtitle = "Event time. Spillover effect sectors.")

Figure5C <- ie.plot2

#Figure 5D
df.tot.wks <- lapply(df.tot, function(x){return(x[[1]])})
df.tot2 <- array(df.tot.wks |> unlist(), dim = c(183-60+1,length(treatment.sectors)), dimnames = list(names(df.tot.wks[[1]]), treatment_sectors))
total_diff_wks <- as.data.frame(df.tot2)
total_diff_wks['wks'] <- rownames(total_diff_wks) |> as.numeric()
total_diff_wks.longer <- pivot_longer(total_diff_wks,!wks, names_to = "Sector_ID", values_to = "diff") 

total_diff_wks_group = sapply(total_diff_wks.longer$wks, classify, number.li = c(60,79,106,132,158,183), 
                              label.li = c("EW8","EW27","EW1","EW27","EW1"),
                              label.li2 = c("EW26,2020","EW53,2020","EW26,2021","EW52,2021","EW26,2022"))


total_diff_wks.longer["wks_group"] = total_diff_wks_group

total_diff_wks.longer$wks_group <- factor(total_diff_wks.longer$wks_group, levels = c("EW8-EW26,2020","EW27-EW53,2020","EW1-EW26,2021","EW27-EW52,2021","EW1-EW26,2022"))
diff.plot.box <- ggplot(data = total_diff_wks.longer, aes(x = wks_group, y = diff, fill = wks_group)) + 
  geom_boxplot(lwd = 0.3, fatten = 0.3,outlier.shape = NA) +
  labs(x = "Calendar time interval", y = "Raw proportion increment")+
  scale_y_continuous(breaks = c(-0.2,0,0.2,0.4), labels = c("-20%","0%", "20%","40%")) +
  coord_cartesian(ylim=c(-0.2,0.4)) +
  theme(legend.title = element_text(size = 8), 
        legend.text = element_text(size=6), 
        plot.subtitle=element_text(size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(colour = 'black', linewidth = 0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', linewidth = 0.5, linetype='solid'),
        panel.background = element_blank(),
        axis.text.x = element_text(size = rel(0.7)),
        panel.border = element_rect(color = 'black', fill = NA, linewidth = 0.5)) +
  geom_hline(yintercept = 0, linetype = 3, linewidth = 0.3) +
  guides(fill = "none") +
  ggtitle(paste("Percentage increase for 0 female Ae.aegypti mosquito trap"), subtitle = "Calendar time. Spillover effect sectors.")

Figure5D <- diff.plot.box

#------ (Optional) Combine Figure 5A,B,C,D to one figure. Adjust parameters for a better visualization.

a <- Figure5A
b <- Figure5B
e <- Figure5C
f <- Figure5D

modified <- lapply(list(a,b,e,f), function(x){
  x + theme(axis.text.y = element_text(size = rel(1)), 
            axis.title.y = element_text(size = rel(1)),
            axis.title.x = element_text(size = rel(1)),
            axis.text.x = element_text(size = rel(1)),
            plot.margin = margin(0, 0, 0, 0)) + 
    labs(title = NULL, subtitle = NULL, y = NULL) 
})
a <- modified[[1]]
b <- modified[[2]]
e <- modified[[3]]
f <- modified[[4]]


modified11 <- lapply(list(a,b), function(x){
  x + theme(axis.title.x = element_blank())
})

a <- modified11[[1]]
b <- modified11[[2]]


e <- e + labs(x = "Time since intervention (weeks)")

b <- b + scale_x_discrete(labels = c("EW8-EW26, 2020","EW27-EW53, 2020","EW1-EW26, 2021","EW27-EW52, 2021","EW1-EW26, 2022")) + theme(axis.text.x = element_text(size = 5))
f <- f + scale_x_discrete(labels = c("EW8-EW26, 2020","EW27-EW53, 2020","EW1-EW26, 2021","EW27-EW52, 2021","EW1-EW26, 2022")) + theme(axis.text.x = element_text(size = 5))
a <- a + theme(axis.text.x = element_text(size = 6))
e <- e + theme(axis.text.x = element_text(size = 6))

Figure5 <- wrap_plots(
  a,b,
  e,f,
  ncol = 2, nrow = 2
) + plot_annotation(tag_levels = list(c("A. Increase in percentage of zero traps (direct)", "B. Increase in percentage of zero traps (direct)",
                                        #"C. Proportion relative change (%, direct)", "D. Proportion relative change (%, direct)",
                                        "C. Increase in percentage of zero traps (spillover)", "D. Increase in percentage of zero traps (spillover)"
                                        #"G. Proportion relative change (%, spillover)", "H. Proportion relative change (%, spillover)"
))) +
  plot_layout(guides = "collect") & theme(panel.spacing = unit(0, "cm"), 
                                          plot.margin = margin(3, 0, 3, 0),
                                          plot.tag.position = c(0.5,0.93),
                                          plot.tag = element_text(size = 7),
                                          axis.title.x = element_text(size = 7),
                                          axis.title.y = element_text(size = 7),
                                          #axis.text.x  = element_text(size = 6),
                                          axis.text.y  = element_text(size = 6),
                                          plot.title   = element_text(size = 7, face = "bold")
  ) # Remove axis.~ arguments to get larger texts

#ggsave("./Figure5_300dpi.tiff",Figure5, units  ="in", width = 7.5, height = 3, dpi = 300)


#------SI Figure 1A
treatment.sectors = treatment_sectors
B = 1000
b <- function(df){
  sample_df <- df %>% t() %>% as.data.frame() %>% slice_sample(prop = 1, replace = TRUE) %>% t() %>% apply(MARGIN = 1, mean)
  return(sample_df)
}

f <- function(tid){
  load(paste("RData/",tid,"_data_2.RData", sep = ''))
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
  print(paste("week",wk,"start"))
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


s <- ggplot_build(ie.plot.box)
label_x <- s$data[[1]]$x
x_need_arrow <- label_x[s$data[[1]]$ymin < iqr.low] |> as.numeric()
y_canvas_min <- s$layout$panel_params[[1]]$y.range[1]
t <- ie.plot.box$data$wks_group |> levels()
groups <- rep(t,length(x_need_arrow))[s$data[[1]]$ymin < iqr.low]
arrow_df <- data.frame(X = x_need_arrow, Y = rep(iqr.low,length(x_need_arrow)),Yend = rep(y_canvas_min,length(x_need_arrow)), wks_group = groups)
ie.plot.box <- ie.plot.box + geom_segment(data = arrow_df, arrow = arrow(length=unit(.2, 'cm'), type = "closed"), aes(x = X, y = Y, yend = Yend),linewidth = 0.6)

#-----------SI Figure 2 ----------
#SI Figure 2A
#Step 1: aggregate IEs of treated sectors for all each event time interval into one dataframe
aegypti_ie_13 <- pbmclapply(treatment.sectors, function(id){
  load(paste(id,"_data_2.RData", sep = ''))
  post.IE <- IEs[as.numeric(names(IEs)) >= 0]
  if (length(post.IE)<10){
    addNA <- rep(NA, 10-length(post.IE))
    post.IE <- c(post.Ie, addNA)
    names(post.IE) <- c(0:9)
  }
  names(post.IE) <- c("1-13","14-26","27-39","40-52","53-65","66-78","79-91","92-104","105-117","118-124")
  return(post.IE)
}, mc.cores = 31)

cont_IE_agg_df <- do.call(rbind, aegypti_ie_13) |> as.data.frame()

#Step 2: Based on the status (buffer/core) of the treated sector, generate two dataframes, one for core sectors and one for buffer sectors.
direct_C <- brief.release.info %>% filter(sector_type == "treat" & type %in% c("NN","NA.") & next_stage %in% c("C","CNN"))
B2C <- brief.release.info %>% filter(sector_type == "treat" & type %in% c("NN","NA.") & next_stage %in% c("B"))
B_interval <- brief.release.info %>% filter(sector_type == "treat" & type == "B")
prod(B_interval$Sector_ID == B2C$Sector_ID) == 1
B_intervals <- B_interval$end_on - B_interval$Eweek_tot + 1
B2C["interval"] <- B_intervals
B2C["B_agg_week"] <- floor(B_intervals/13) + c(B_interval$next_stage == "B") * c(B_intervals %% 13 > 0) #for those final status is Buffer site, add the last aggegated week data

#Buffer areas may convert to core areas, but core areas never convert to buffer areas
#So we first generate a dataframe for buffer areas.
B.df <- apply(B2C, MARGIN = 1, function(x){
  sector_id = x["Sector_ID"]
  agg_week = x["B_agg_week"]
  k <- cont_IE_agg_df[cont_IE_agg_df$Sector_ID == sector_id,c(1:agg_week)] |> unlist()
  if (length(k) < 10){
    NA.following <- rep(NA,10-length(k))
    k <- c(k, NA.following)
  }
  names(k) <- c(0:9)
  return (k)
}) |> as.data.frame() |> t() |> as.data.frame()

rownames(B.df) <- B2C$Sector_ID

#Second we generate a dataframe for core areas. This part of the core sectors is those sectors converted from the buffer area.
C.df1 <- apply(B2C, MARGIN = 1, function(x){
  sector_id = x["Sector_ID"]
  agg_week = x["B_agg_week"] |> as.numeric() + 1
  if (as.numeric(x["interval"]) %%13 > 0){
    agg_week = agg_week + 1
  }
  if (agg_week > 9){
    return (rep(NA,10))
  }
  k <- cont_IE_agg_df[cont_IE_agg_df$Sector_ID == sector_id,c(agg_week:10)] |> unlist()
  NA.lead <- rep(NA, 10-length(k))
  k <- c(NA.lead,k)
  names(k) <- c(0:9)
  return (k)
}) |> as.data.frame() |> t()
rownames(C.df1) <- B2C$Sector_ID
colnames(C.df1) <- c(0:9)

#This part of the core sectors is the remaining sectors, which are always the core sectors.
C.df2 <- cont_IE_agg_df[!(cont_IE_agg_df$Sector_ID %in% B2C$Sector_ID),]
rownames(C.df2) <- C.df2$Sector_ID
colnames(C.df2) <- c(0:9)

#Combine two 'core sector' dataframes together
C.df <- rbind(C.df1|>as.data.frame(), C.df2[1:10])

total_cont_IE_B <- apply(B.df, MARGIN = 2, mean,na.rm = TRUE)
total_cont_IE_C <- apply(C.df, MARGIN = 2, mean,na.rm = TRUE)

#Bootstrap CI
B = 1000

b <- function(df){
  res <- apply(df, MARGIN = 2, function(x){
    na.value <- x %>% is.na()
    x <- x[!na.value]
    return(sample(x, size = length(x), replace = TRUE))
  })
  return(res)
}

if (Sys.info()['sysname'] %in% c("Linux","Unix")){
  bootstrap_B <- pbmclapply(replicate(1000, B.df , simplify = FALSE), b, mc.cores = 25)
  bootstrap_C <- pbmclapply(replicate(1000, C.df, simplify = FALSE), b, mc.cores = 25)
} else {
  bootstrap_B <- lapply(replicate(1000, B.df, simplify = FALSE), b) 
  bootstrap_C <- lapply(replicate(1000, C.df, simplify = FALSE), b)
}

#For buffer
bootstrap_b2 <- lapply(bootstrap_B, function(x){return(sapply(x,mean))})
bootstrap_b3 <- bootstrap_b2 |> unlist() |> array(dim = c(10, 1000))
b.025.b <- bootstrap_b3 |> apply(MARGIN = 1, function(x) {return(sort(x, decreasing = FALSE)[B*0.025])})
b.975.b <- bootstrap_b3 |> apply(MARGIN = 1, function(x) {return(sort(x, decreasing = FALSE)[B*0.975])})

lb.b = 2 * total_cont_IE_B - b.975.b
ub.b = 2 * total_cont_IE_B - b.025.b

#For core
bootstrap_c4 <- lapply(bootstrap_C, function(x){return(sapply(x,mean))})
bootstrap_c5 <- bootstrap_c4 |> unlist() |> array(dim = c(10, 1000))
b.025.c <- bootstrap_c5 |> apply(MARGIN = 1, function(x) {return(sort(x, decreasing = FALSE)[B*0.025])})
b.975.c <- bootstrap_c5 |> apply(MARGIN = 1, function(x) {return(sort(x, decreasing = FALSE)[B*0.975])})

lb.c = 2 * total_cont_IE_C - b.975.c
ub.c = 2 * total_cont_IE_C - b.025.c


graph_df_B <- data.frame(weeks = names(total_cont_IE_B) |> as.numeric(), 
                         value = total_cont_IE_B, 
                         lb = lb.b, 
                         ub = ub.b,
                         group = "Buffer area")
graph_df_C <- data.frame(weeks = names(total_cont_IE_C) |> as.numeric(), 
                         value = total_cont_IE_C, 
                         lb = lb.c, 
                         ub = ub.c,
                         group = "Core area")

graph_df_BC <- rbind(graph_df_B,graph_df_C)

x_break2 = seq(from = 0, to = 9, by = 1)
ie.plot3 <- ggplot(data = graph_df_BC, aes(x = weeks,y = -value, fill = group, color = group)) + geom_line(aes()) + 
  geom_ribbon(aes(ymin = -lb, ymax = -ub, color = NULL),alpha = 0.2)+
  geom_hline(yintercept = 1, linetype = "dashed", color = "darkgreen", alpha = 1) +
  labs(x = "post-event week interval", y = "IE(%)")+
  scale_x_continuous(breaks = x_break2, labels = c("1-13","14-26","27-39","40-52","53-65","66-78","79-91","92-104","105-117","118-124"))+
  scale_y_continuous(breaks = seq(0,1,by = 0.25), labels = seq(0,1,by = 0.25)*100) +
  coord_cartesian(ylim = c(0,1)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.7, hjust=0.5), 
        plot.subtitle = element_text(size = rel(0.5)),
        legend.title = element_blank(),
        axis.line.x = element_line(colour = 'black', linewidth = 0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', linewidth = 0.5, linetype='solid'),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5),
        axis.title.x = element_text(margin = margin(0,0,0,0)),
        legend.margin = margin(0, 0, 0, 0),
        legend.position = c(0.8,0.3),
        legend.background = element_blank()) 

SIFigure2A <- ie.plot3

#Figure 2B ------- The procedure to get Figure 2B is similar to that to get Figure 2A.
#Step 1 :aggregate IEs of all treated sectors for each calendar-time interval into one dataframe
aegypti_ie_cal <- pbmclapply(treatment_sectors, function(id) {
  load(paste(id,"_data_2.RData",sep = ''))
  treatment_start <- convertweek %>% filter(E2week == 0) %>% select(Eweek_tot) %>% min()

  post.IE <- IEs.was.eweek[as.numeric(names(IEs.was.eweek)) >= 60 & as.numeric(names(IEs.was.eweek)) >= treatment_start] |> unlist()
  if(length(post.IE) == 0){
    res <- rep(NA,5)
    names(res) <- c("EW8-26, 2020", "EW27-53, 2020","EW1-53, 2021","EW27-52, 2021","EW1-26, 2022")
    return (res)
  }
  
  if (min(as.numeric(names(post.IE))) > 60) {
    addNA <- rep(NA, min(as.numeric(names(post.IE))) - 60)
    post.IE <- c(addNA,post.IE)
    names(post.IE) <- c(60:max(as.numeric(names(post.IE)), na.rm = TRUE))
  }
  
  if (length(post.IE) < 183-60 + 1){
    addNA <- rep(NA, 183-60+1-length(post.IE))
    post.IE <- c(post.IE, addNA)
    names(post.IE) <- c(60:183)
  }
  
  mean0 <- mean(post.IE[1:19], na.rm = TRUE)
  
  post.IE <- post.IE[20:length(post.IE)]
  mean1 <- mean(post.IE[1:27], na.rm = TRUE)
  mean2 <- mean(post.IE[28:53], na.rm = TRUE)
  mean3 <- mean(post.IE[54:79], na.rm = TRUE)
  mean4 <- mean(post.IE[80:105], na.rm = TRUE)
  
  res <- c(mean0, mean1,mean2,mean3,mean4)
  names(res) <- c("EW8-26, 2020","EW27-53, 2020","EW1-53, 2021","EW27-52, 2021","EW1-26, 2022")
  return (res)
}, mc.cores = 31)

cont_IE_agg_df <- do.call(rbind, aegypti_ie_cal) |> as.data.frame()
cont_IE_agg_df$Sector_ID <- treatment_sectors

#Step 2 Based on the status (buffer/core) of the treated sector, generate two dataframes, one for core sectors and one for buffer sectors.
direct_C <- brief.release.info %>% filter(sector_type == "treat" & type %in% c("NN","NA.") & next_stage %in% c("C","CNN"))
B2C <- brief.release.info %>% filter(sector_type == "treat" & type %in% c("NN","NA.") & next_stage %in% c("B"))
B_interval <- brief.release.info %>% filter(sector_type == "treat" & type == "B")
prod(B_interval$Sector_ID == B2C$Sector_ID) == 1

# || EW8-26, 2020  || EW27-53, 2020  || EW1-26, 2021 ||  EW27-52, 2021 || EW1-26, 2022 ||
# || 60 - 78       || 79 - 105       || 106 - 131    ||  132 - 157     || 158 - 183    ||
#Step 2-1: For sectors in buffer areas
B.df <- apply(B2C, MARGIN = 1, function(x){
  sector_id = x["Sector_ID"]
  end_on = B_interval[B_interval$Sector_ID == sector_id,]$end_on
  
  if (end_on >=60 & end_on < 78){
    return(rep(NA,5))
  }
  if (end_on >=78 & end_on < 105){
    return(c(cont_IE_agg_df[cont_IE_agg_df$Sector_ID == sector_id,1],rep(NA,4)))
  }
  if (end_on >=105 & end_on < 131){
    return(c(cont_IE_agg_df[cont_IE_agg_df$Sector_ID == sector_id,c(1:2)] |> unlist(),rep(NA,3)))
  }
  if (end_on >=131 & end_on < 157){
    return(c(cont_IE_agg_df[cont_IE_agg_df$Sector_ID == sector_id,c(1:3)]|> unlist(),rep(NA,2)))
  }
  if (end_on >=157 & end_on < 183){
    return(c(cont_IE_agg_df[cont_IE_agg_df$Sector_ID == sector_id,c(1:4)]|> unlist(),rep(NA,1)))
  }
  if (end_on == 183){
    return(cont_IE_agg_df[cont_IE_agg_df$Sector_ID == sector_id,c(1:5)]|> unlist())
  }
}) |> as.data.frame() |> t()

rownames(B.df) <- B2C$Sector_ID

#Step 2-2: For sectors in core areas
C.df1 <- apply(B2C, MARGIN = 1, function(x){
  sector_id = x["Sector_ID"]
  end_on = B_interval[B_interval$Sector_ID == sector_id,]$end_on
  IEs <- cont_IE_agg_df[cont_IE_agg_df$Sector_ID == sector_id,c(1:5)]
  if (end_on < 60){
    return(IEs)
  }
  if (end_on >=60 & end_on <= 78){
    return(c(NA,IEs[2:5]))
  }
  if (end_on >78 & end_on <= 105){
    return(c(rep(NA,2),IEs[3:5]))
  }
  if (end_on >105 & end_on <= 131){
    return(c(rep(NA,3),IEs[4:5]))
  }
  if (end_on >131 & end_on <= 157){
    return(c(rep(NA,4),IEs[5]))
  }
  if (end_on >157 & end_on <= 183){
    return(rep(NA,5))
  }
  
}) |> as.data.frame() |> t()
rownames(C.df1) <- B2C$Sector_ID
colnames(C.df1) <- c(0:4)

C.df2 <- cont_IE_agg_df[!(cont_IE_agg_df$Sector_ID %in% B2C$Sector_ID),]

rownames(C.df2) <- C.df2$Sector_ID
colnames(C.df2) <- c(0:4)

C.df <- rbind(C.df1|>as.data.frame(), C.df2[1:5])

total_cont_IE_B <- apply(B.df, MARGIN = 2, function(x){return(x |> as.numeric() |> mean(na.rm = TRUE))})
total_cont_IE_C <- apply(C.df, MARGIN = 2, function(x){return(x |> as.numeric() |> mean(na.rm = TRUE))})

#Step 3: Use bootstrap to estimate the pointwise confidence intervals
#Bootstrap CI
B = 1000

b <- function(df){
  res <- apply(df, MARGIN = 2, function(x){
    na.value <- x %>% is.na()
    x <- x[!na.value]
    return(sample(x, size = length(x), replace = TRUE))
  })
  return(res)
}

if (Sys.info()['sysname'] %in% c("Linux","Unix")){
  bootstrap_B <- pbmclapply(replicate(1000, B.df , simplify = FALSE), b, mc.cores = 25)
  bootstrap_C <- pbmclapply(replicate(1000, C.df, simplify = FALSE), b, mc.cores = 25)
} else {
  bootstrap_B <- lapply(replicate(1000, B.df, simplify = FALSE), b) 
  bootstrap_C <- lapply(replicate(1000, C.df, simplify = FALSE), b)
}

#For buffer areas
bootstrap_b2 <- lapply(bootstrap_B, function(x){return(sapply(x,mean))})
bootstrap_b3 <- bootstrap_b2 |> unlist() |> array(dim = c(5, 1000))
b.025.b <- bootstrap_b3 |> apply(MARGIN = 1, function(x) {return(sort(x, decreasing = FALSE)[B*0.025])})
b.975.b <- bootstrap_b3 |> apply(MARGIN = 1, function(x) {return(sort(x, decreasing = FALSE)[B*0.975])})

lb.b = 2 * total_cont_IE_B - b.975.b
ub.b = 2 * total_cont_IE_B - b.025.b

#For core areas
bootstrap_c4 <- lapply(bootstrap_C, function(x){return(sapply(x,function(x){return(x |> unlist() |> mean())}))})
bootstrap_c5 <- bootstrap_c4 |> unlist() |> array(dim = c(5, 1000))
b.025.c <- bootstrap_c5 |> apply(MARGIN = 1, function(x) {return(sort(x, decreasing = FALSE)[B*0.025])})
b.975.c <- bootstrap_c5 |> apply(MARGIN = 1, function(x) {return(sort(x, decreasing = FALSE)[B*0.975])})

lb.c = 2 * total_cont_IE_C - b.975.c
ub.c = 2 * total_cont_IE_C - b.025.c

#Step 4: Plot the graph (SI Figure 2B)
graph_df_B <- data.frame(weeks = c(0:4) |> as.numeric(), 
                         value = total_cont_IE_B, 
                         lb = lb.b, 
                         ub = ub.b,
                         group = "Buffer area")
graph_df_C <- data.frame(weeks = c(0:4) |> as.numeric(), 
                         value = total_cont_IE_C, 
                         lb = lb.c, 
                         ub = ub.c,
                         group = "Core area")

graph_df_BC <- rbind(graph_df_B,graph_df_C)

x_break2 = seq(from = 0, to = 4, by = 1)
ie.plot3 <- ggplot(data = graph_df_BC, aes(x = weeks,y = -value, fill = group, color = group)) + geom_line(aes()) + 
  geom_ribbon(aes(ymin = -lb, ymax = -ub, color = NULL),alpha = 0.2)+
  geom_hline(yintercept = 1, linetype = "dashed", color = "darkgreen", alpha = 1) +
  labs(x = "Calendar week interval", y = "IE(%)")+
  scale_x_continuous(breaks = x_break2, labels = c("EW8-26, 2020","EW27-53, 2020","EW1-26, 2021","EW27-52, 2021","EW1-26, 2022"))+
  scale_y_continuous(breaks = seq(0,1,by = 0.25), labels = seq(0,1,by = 0.25)*100) +
  coord_cartesian(ylim = c(0,1)) +
  theme(plot.subtitle = element_text(size = rel(0.5)),
        legend.title = element_blank(),
        axis.line.x = element_line(colour = 'black', linewidth = 0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', linewidth = 0.5, linetype='solid'),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5),
        axis.title.x = element_text(margin = margin(0,0,0,0)),
        legend.margin = margin(0, 0, 0, 0),
        legend.position = c(0.8,0.3),
        legend.background = element_blank()) 

SIFigure2B <- ie.plot3