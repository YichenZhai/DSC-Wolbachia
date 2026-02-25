# The databases used in this study are the property of the National Environment Agency, Singapore.
# Permission to access mosquito abundance data should be obtained from the National Environment Agency, Singapore.

#1. Load necessary libraries
library(dplyr)
library(readxl)
#2. Set working directory that contains mosquito surveillance dataset.
setwd("~/Desktop/gravitrap/")

#3. Load surveillance dataset and merge all surveillance data from 2019 to 2024
raw_data <- data.frame()
for (i in c(2019:2024)){
  for (j in c("Q1","Q2","Q3","Q4")){
    tryCatch({
      df <- readxl::read_excel(paste(as.character(i),j,"weekly_trap_albo.xlsx",sep = "_"))
      if (length(colnames(raw_data))!=0){
        if (all.equal(colnames(df),colnames(raw_data)) != TRUE){
          colnames(df) <- colnames(raw_data)
        }
      }
      raw_data <- rbind(raw_data, df)
      print(paste(as.character(i),j,"weekly_trap_albo.xlsx finished",sep = "_"))
    },error = function(e){
      print(paste(as.character(i),j,"weekly_trap_albo.xlsx NOT FOUND",e,sep = "_"))
    },finally = {})
    
  }
}
remove(df)

#4. Data Preprocessing
  # Keep only necessary rows (2019 EW8- 2022 EW26) and columns (related to GAI, sector and time)
raw_data <- raw_data[,which(colnames(raw_data) %in% c("Eyear","Emonth","Eweek","Sector_ID","Postal","Number of wild-type female Ae. aegypti","Number of functional Gravitrap","Mean number of F.Ae.albopictus caught per functional trap"))]
raw_data <- raw_data %>% filter((Eyear == 2019 & Eweek >= 8) | (Eyear > 2019 & Eyear < 2022)|(Eyear == 2022 & Eweek <= 26))

  # Get the treatment sectors and control sectors
treatment_info <- readxl::read_excel("Release_matrix_2016_2022ew26.xlsx")

  # Get treatment e-year and e-week for treated sectors
treatment_sectors <- as.data.frame(treatment_info) %>% filter(ReleaseOrNot_revised == 1)
treatment_sectors <- aggregate(treatment_sectors, by = list(treatment_sectors$Sector_ID,treatment_sectors$StudyArea), FUN = min)[,c("Sector_ID","StudyArea","Eyear")] 
treatment_sectors["Eweek"] <- unlist(lapply(treatment_sectors$Sector_ID,function(x){
  y = treatment_sectors[treatment_sectors$Sector_ID == x, "Eyear"]
  return((treatment_info %>% filter(Sector_ID == x & Eyear == y & ReleaseOrNot_revised == 1))$Eweek %>% min())
}))

  # Remove those treated before 2020 EW1
treatment_sectors <- treatment_sectors[which(treatment_sectors$Eyear>=2020),]


  # Calculate treatment starting week based on [ 2019 EW1 = Week 1 ]
correction <- raw_data[,c("Eyear","Eweek")] %>% unique() %>% group_by(Eyear) %>% summarize(last_week = max(Eweek))
correction["correct"] <- c(0,cumsum(correction$last_week)[1:3])
treatment_sectors["Eweek_tot"] <- unlist(lapply(treatment_sectors$Eyear, function(y){
  return(correction[correction$Eyear == y,]$correct)
}))+treatment_sectors$Eweek

  # Find first surveillance record for each treatment sector and check whether the sector has at least 1-year pre-treatment period.
treat_surveillance <- raw_data %>% filter(Sector_ID %in% treatment_sectors$Sector_ID)
control_surveillance <- raw_data %>% filter(!Sector_ID %in% treatment_sectors$Sector_ID)
treat_surveillance["Eweek_tot"] <- unlist(lapply(treat_surveillance$Eyear, function(y){
  return(correction[correction$Eyear == y,]$correct)
}))+treat_surveillance$Eweek
treat_record <- treat_surveillance[,c("Sector_ID","Eweek_tot")] %>% unique() %>% group_by(Sector_ID) %>% summarize(min_eweek = min(Eweek_tot))
treatment_sectors <- treatment_sectors |> left_join(treat_record, by = "Sector_ID")
treatment_sectors["EW_diff"] <- treatment_sectors$Eweek_tot - unlist(lapply(treatment_sectors$min_eweek,function(z){return(max(z,8))}))
treatment_info <- treatment_sectors %>% filter(EW_diff >= 52)
treatment_sectors <-treatment_info$Sector_ID
treat_surveillance <- treat_surveillance %>% filter(Sector_ID %in% treatment_sectors)

remove(treat_record)

#Get surveillance record for all possible control sectors
release.sector.info <- read.csv("TASK1_REDONE_4-levels_Sectors.csv")

control.weeks <- control_surveillance %>% group_by(Sector_ID) %>% summarise(nrows = n_distinct(Eyear,Eweek))
control_sectors <- control.weeks[control.weeks$nrows == 183-8+1,"Sector_ID"] |> unlist() |> unname()
control_sectors <- control_sectors[grepl("\\d",control_sectors)]
r <- release.sector.info$Sector_ID |>unique() |>unlist()
r2 <- release.sector.info %>% filter(ReleaseOrNot == 1) %>% select (Sector_ID) |>unique() |>unlist()
r <- r[!r%in% r2]

all_sectors <- c(control_sectors, treatment_sectors)
release.sector.info <- release.sector.info[release.sector.info$Sector_ID |> unlist() %in% all_sectors,]
release.sector.info["Eyear"] <- sapply(release.sector.info$ey.ew, function(x) {return (substr(x,1,4))})
release.sector.info["Eweek"] <- sapply(release.sector.info$ey.ew, function(x) {return (substr(sprintf("%.2f", x),6,7))})
release.sector.info["Eweek_tot"] <- sapply(release.sector.info[["Eweek"]], as.numeric) + sapply(release.sector.info[["Eyear"]], function(y) {return(correction[correction$Eyear == as.numeric(y),"correct"])})|>unlist()


get_release_info <- function(sector_id){
  if (sector_id %in% treatment_sectors) index = "treat" else index = "control"
  temp_df <- release.sector.info %>% filter(Sector_ID == sector_id)
  min_NN <- temp_df %>% filter(NN == 1 & C == 0)
  if (nrow(min_NN) > 0){
    min_NN <- min_NN[min_NN$Eweek_tot == min(min_NN$Eweek_tot),]
    min_NN["type"] <- "NN"
  }
  min_NA <- temp_df %>% filter(NA. == 1)
  if (nrow(min_NA) > 0){
    min_NA <- min_NA[min_NA$Eweek_tot == min(min_NA$Eweek_tot),]
    min_NA["type"] <- "NA."
  }
  min_B <- temp_df %>% filter(B == 1)
  if (nrow(min_B) > 0){
    min_B <- min_B[min_B$Eweek_tot == min(min_B$Eweek_tot),]
    min_B["type"] <- "B"
  }
  min_C <- temp_df %>% filter(C == 1 & NN == 0)
  if (nrow(min_C) > 0){
    min_C <- min_C[min_C$Eweek_tot == min(min_C$Eweek_tot),]
    min_C["type"] <- "C"
  }
  min_CNN <- temp_df %>% filter(C == 1 & NN == 1)
  if (nrow(min_CNN) > 0){
    min_CNN <- min_CNN[min_CNN$Eweek_tot == min(min_CNN$Eweek_tot),]
    min_CNN["type"] <- "CNN"
  }
  combine_df <- rbind(min_NN,min_NA,min_B, min_C, min_CNN)
  combine_df["sector_type"] <- rep(index, nrow(combine_df))
  return (combine_df)
}

brief.release.info <- lapply(all_sectors, get_release_info) %>% do.call(what = rbind)
brief.release.info <- brief.release.info %>% group_by(Sector_ID) %>% mutate(end_on = ifelse(is.na(lead(Eweek_tot, n=1, order_by=Eweek_tot)),184,lead(Eweek_tot, n=1, order_by=Eweek_tot)) - 1) %>% ungroup()
brief.release.info <- brief.release.info %>% group_by(Sector_ID) %>% mutate(next_stage = ifelse(is.na(lead(type, n=1, order_by=Eweek_tot)),type,lead(type, n=1, order_by=Eweek_tot))) %>% ungroup()


control_sectors <- control_sectors[control_sectors %in% r]
never_NA_control <- brief.release.info %>% filter(type == "NN" & end_on == 183) %>% select(Sector_ID) %>% unlist() %>% unique()
control_sectors <- intersect(control_sectors,never_NA_control)
remove(control.weeks)
remove(r)
remove(r2)


control_surveillance <- raw_data %>% filter(Sector_ID %in% control_sectors)
control_surveillance["Eweek_tot"] <- unlist(lapply(control_surveillance$Eyear, function(y){
  return(correction[correction$Eyear == y,]$correct)
}))+control_surveillance$Eweek

  #put spillover sector data into treat_surveillance
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
spillover_sectors = names(res)

spillover_surveillance <- raw_data %>% filter(Sector_ID %in% spillover_sectors)
spillover_surveillance["Eweek_tot"] <- unlist(lapply(spillover_surveillance$Eyear, function(y){
  return(correction[correction$Eyear == y,]$correct)
}))+spillover_surveillance$Eweek

treat_surveillance <- rbind(treat_surveillance, spillover_surveillance)
  # Remove the original dataset to save storage.
remove(raw_data)

#5. Save workspace for later work
save.image(file = "DSC_data3.RData")

#End of data preprocessing.


