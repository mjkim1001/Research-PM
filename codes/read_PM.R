library(readxl)
library(lubridate)
xdata <- read_excel("data/2020년 1월.xlsx")
tdata <- tibble(code = xdata$측정소코드, date = substring(xdata$측정일시,1,8), PM25 = xdata$PM25)
tdata = tdata %>% group_by(code, date)%>% summarise(PM25 = mean(PM25,na.rm=T)) %>% mutate(date = as.Date(date,"%Y%m%d")) %>% ungroup()

##################################################3
#Read Data fucntion

read.data<-function(tdata, years){
  for(year in years)
  for(file in list.files(path= paste("./data/",year,sep = ""))){
    if(endsWith(file,".xlsx")){
      xdata<-read_excel(paste("data/",year,"/",file,sep = ""))
    }else{next}
    print(sprintf("Opening %s",file))
    temp <- tibble(code = xdata$측정소코드, date = substring(xdata$측정일시,1,8), PM25 = xdata$PM25)
    temp = temp %>% group_by(code, date)%>% summarise(PM25 = mean(PM25,na.rm=T)) %>% mutate(date = as.Date(date,"%Y%m%d")) %>% ungroup()
    
    tdata = rbind(tdata,temp)
  }
  return(tdata)
}
tdata = read.data(tdata, c(2019,2018))

saveRDS(tdata, "dataset_daily.rds")

#6hr
tdata <- tibble(code = xdata$측정소코드, date_org = xdata$측정일시, PM25 = xdata$PM25)
tdata = tdata %>% mutate(date = substring(date_org,1,8)) %>% mutate(date_idx = as.integer((as.numeric(substring(date_org,9,10))-1)/6)+1)
tdata = tdata %>% group_by(code, date,date_idx)%>% summarise(PM25 = mean(PM25,na.rm=T)) %>% mutate(date = as.Date(date,"%Y%m%d")) %>% ungroup()
##################################################3
#Read Data fucntion

read.data<-function(tdata, years){
  for(year in years)
    for(file in list.files(path= paste("./data/",year,sep = ""))){
      if(endsWith(file,".xlsx")){
        xdata<-read_excel(paste("data/",year,"/",file,sep = ""))
      }else{next}
      print(sprintf("Opening %s",file))
      temp <- tibble(code = xdata$측정소코드, date_org = xdata$측정일시, PM25 = xdata$PM25)
      temp = temp %>% mutate(date = substring(date_org,1,8)) %>% mutate(date_idx = as.integer((as.numeric(substring(date_org,9,10))-1)/6)+1)
      temp = temp %>% group_by(code, date,date_idx)%>% summarise(PM25 = mean(PM25,na.rm=T)) %>% mutate(date = as.Date(date,"%Y%m%d")) %>% ungroup()
      
      tdata = rbind(tdata,temp)
    }
  return(tdata)
}
tdata = read.data(tdata, 2018:2021)
saveRDS(tdata, "dataset_6hr.rds")
