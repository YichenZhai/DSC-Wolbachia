#Final paper use the following code for info extraction
#Table 1 in manuscript
library(pbapply)

setwd("~/../Desktop/DSC_redo_all/aegypti/Rfiles/")
library(dplyr)
load("DSC_data3.RData")

setwd("~/../Desktop/gravitrap/")

#Merge all surveillance data together from 2019 to 2024
raw_data <- data.frame()
for (i in c(2019:2022)){
  for (j in c("Q1","Q2","Q3","Q4")){
    tryCatch({
      df <- readxl::read_excel(paste(as.character(i),j,"weekly_trap_albo.xlsx",sep = "_"))
      if (length(colnames(raw_data))!=0){
        if (all.equal(colnames(df),colnames(raw_data)) != TRUE){
          colnames(df) <- colnames(raw_data)
        }
      }
      raw_data <- rbind(raw_data, df)
      print(paste(as.character(i),j,"weekly_trap.csv finished",sep = "_"))
    },error = function(e){
      print(paste(as.character(i),j,"weekly_trap.csv NOT FOUND",e,sep = "_"))
    },finally = {})
    
  }
}

raw_data2 <- raw_data %>% filter((Eyear == 2019 & Eweek >= 8) | Eyear == 2020 | Eyear == 2021 |(Eyear == 2022 & Eweek <=26))
raw_data2["Eweek_tot"] <- unlist(lapply(raw_data2$Eyear, function(y){
  return(correction[correction$Eyear == y,]$correct)
}))+raw_data2$Eweek

#save.image(file = "table1_preprocessing.RData")
load("table1_preprocessing.RData")
#The number of traps will be counted as follow:
#1. Filter the surveillance dataset ("raw_data2") to the specific sort of sectors
#   E.g. "Core": should be directly treated sectors surrounded by treated sectors during at least one week of the interval
#2. Remove duplicated records that shares the same trap characteristics ("Sector_ID","Block","Postal","Level","Unit")
#   Note: One combination of trap characteristics is corresponding to one unique trap
#   e.g. In Sector_ID = FL424,with the same "Block","Postal","Level",different "Unit" number means different traps.
#   Use "group_by" function
#3. Count the number of records left
#   Use "summarize" function

# Interval  (60-78)|(79-105)|(106-131) | (132-157) | (158-183)
# Use Min =   59   |  78    |  105     |   131     |    157   (not inclusive)
# Use Max =   78   |  105   |  131     |   157     |    183   (inclusive)
mins = c(59,78,105,131,157)
maxs = c(78,105,131,157,183)

#(1) All Directly treated sectors. Split the records by "c" core and "b" buffer
raw_data3 <- raw_data2 %>% filter(Sector_ID %in% treatment_sectors)
buffer_df <- raw_data3[0,]
core_df <-raw_data3[0,]

raw_data4 <- left_join(raw_data3, brief.release.info, by = join_by(Sector_ID,Eweek_tot >= Eweek_tot, Eweek_tot <= end_on)) %>%
  rename(Eweek_tot = Eweek_tot.x)

#1. Buffer
buffer_records <- raw_data4 %>% filter(type == "B")

  #Count the number of traps in Buffer area
sapply(c(1:5), function(x){
  a <- mins[x]
  b <- maxs[x]
  return (buffer_records %>% filter(Eweek_tot > a & Eweek_tot <= b) %>% group_by(Sector_ID,Block,Postal,Level,Unit) %>% summarize(a = n()) %>% ungroup() %>% nrow()
)
})
buffer_records %>% filter(Eweek_tot > 59 & Eweek_tot <= 183) %>% group_by(Sector_ID,Block,Postal,Level,Unit) %>% summarize(a = n()) %>% ungroup() %>% nrow()

  #Count the number of sectors in Buffer area.
sapply(c(1:5), function(x){
  a <- mins[x]
  b <- maxs[x]
  return (buffer_records %>% filter(Eweek_tot > a & Eweek_tot <= b) %>% group_by(Sector_ID) %>% summarize(a = n()) %>% ungroup() %>% nrow()
  )
})

buffer_records %>% filter(Eweek_tot > 59 & Eweek_tot <= 183) %>% group_by(Sector_ID) %>% summarize(a = n()) %>% ungroup() %>% nrow()

remove(buffer_records)

#2. Core
core_records <- raw_data4 %>% filter(type %in% c("C","CNN"))
sapply(c(1:5), function(x){
  a <- mins[x]
  b <- maxs[x]
  return (core_records %>% filter(Eweek_tot > a & Eweek_tot <= b) %>% group_by(Sector_ID,Block,Postal,Level,Unit) %>% summarize(a = n()) %>% ungroup() %>% nrow()
  )
})

#Number of sectors
sapply(c(1:5), function(x){
  a <- mins[x]
  b <- maxs[x]
  return (core_records %>% filter(Eweek_tot > a & Eweek_tot <= b) %>% group_by(Sector_ID) %>% summarize(a = n()) %>% ungroup() %>% nrow()
  )
})

#Total
core_records %>% filter(Eweek_tot > 59 & Eweek_tot <= 183) %>% group_by(Sector_ID,Block,Postal,Level,Unit) %>% summarize(a = n()) %>% ungroup() %>% nrow()
core_records %>% filter(Eweek_tot > 59 & Eweek_tot <= 183) %>% group_by(Sector_ID) %>% summarize(a = n()) %>% ungroup() %>% nrow()

remove(core_records)

#3. Direct (should be smaller than the sum of core and buffer, because some sectors change their status from Buffer to Core and cause double counts)
direct_records <- raw_data4 %>% filter(type %in% c("B","C","CNN"))
sapply(c(1:5), function(x){
  a <- mins[x]
  b <- maxs[x]
  return (direct_records %>% filter(Eweek_tot > a & Eweek_tot <= b) %>% group_by(Sector_ID,Block,Postal,Level,Unit) %>% summarize(a = n()) %>% ungroup() %>% nrow()
  )
}) 
direct_records %>% filter(Eweek_tot > 59 & Eweek_tot <= 183) %>% group_by(Sector_ID,Block,Postal,Level,Unit) %>% summarize(a = n()) %>% ungroup() %>% nrow()

#Number of trap in direct treated sectors
sapply(c(1:5), function(x){
  a <- mins[x]
  b <- maxs[x]
  return (direct_records %>% filter(Eweek_tot > a & Eweek_tot <= b) %>% group_by(Sector_ID) %>% summarize(a = n()) %>% ungroup() %>% nrow()
  )
}) 
direct_records %>% filter(Eweek_tot > 59 & Eweek_tot <= 183) %>% group_by(Sector_ID) %>% summarize(a = n()) %>% ungroup() %>% nrow()


remove(direct_records)

remove(raw_data3)

#4. Sectors in spillover area
raw_data3 <- raw_data2 %>% filter(Sector_ID %in% spillover_sectors)
raw_data4 <- left_join(raw_data3, brief.release.info, by = join_by(Sector_ID,Eweek_tot >= Eweek_tot, Eweek_tot <= end_on)) %>%
  rename(Eweek_tot = Eweek_tot.x)

spillover_records <-  raw_data4 %>% filter(type == "NA.")
sapply(c(1:5), function(x){
  a <- mins[x]
  b <- maxs[x]
  return (spillover_records %>% filter(Eweek_tot > a & Eweek_tot <= b) %>% group_by(Sector_ID,Block,Postal,Level,Unit) %>% summarize(a = n()) %>% ungroup() %>% nrow()
  )
}) 
spillover_records %>% filter(Eweek_tot > 59 & Eweek_tot <= 183) %>% group_by(Sector_ID,Block,Postal,Level,Unit) %>% summarize(a = n()) %>% ungroup() %>% nrow()
 
  #Number of sectors
sapply(c(1:5), function(x){
  a <- mins[x]
  b <- maxs[x]
  return (spillover_records %>% filter(Eweek_tot > a & Eweek_tot <= b) %>% group_by(Sector_ID) %>% summarize(a = n()) %>% ungroup() %>% nrow()
  )
}) 
spillover_records %>% filter(Eweek_tot > 59 & Eweek_tot <= 183) %>% group_by(Sector_ID) %>% summarize(a = n()) %>% ungroup() %>% nrow()

remove(spillover_records)
remove(raw_data3)
remove(raw_data4)
#4. Control
raw_data3 <- raw_data2 %>% filter(Sector_ID %in% control_sectors)
control_records <- raw_data3

sapply(c(1:5), function(x){
  a <- mins[x]
  b <- maxs[x]
  return (control_records %>% filter(Eweek_tot > a & Eweek_tot <= b) %>% group_by(Sector_ID,Block,Postal,Level,Unit) %>% summarize(a = n()) %>% ungroup() %>% nrow()
  )
}) 
control_records %>% filter(Eweek_tot > 59 & Eweek_tot <= 183) %>% group_by(Sector_ID,Block,Postal,Level,Unit) %>% summarize(a = n()) %>% ungroup() %>% nrow()

  #For the number of sectors
sapply(c(1:5), function(x){
  a <- mins[x]
  b <- maxs[x]
  return (control_records %>% filter(Eweek_tot > a & Eweek_tot <= b) %>% group_by(Sector_ID) %>% summarize(a = n()) %>% ungroup() %>% nrow()
  )
}) 
control_records %>% filter(Eweek_tot > 59 & Eweek_tot <= 183) %>% group_by(Sector_ID) %>% summarize(a = n()) %>% ungroup() %>% nrow()

remove(control_records)
remove(raw_data3)

#5. Calculate the total number of traps in analysis (Consider the whole study, from 2019 to 2022)
raw_data3 <- raw_data2 %>% filter(Sector_ID %in% c(control_sectors, treatment_sectors,spillover_sectors))
raw_data3 %>% filter(Eweek_tot > 7 & Eweek_tot <= 183) %>% group_by(Sector_ID,Block,Postal,Level,Unit) %>% summarize(a = n()) %>% ungroup() %>% nrow()
