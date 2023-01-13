library(fields)
library(lubridate)
library(quantreg)
library(ggplot2)
library(ggmap)
library(grid)
library(gridExtra)
###X : averaged anomalies (data for the analysis)
# mu_st : matrix for mean_effect (X = Observed_data - mean_effect), the average for each specific place and each time of the year by pooling the five years together 
#true_pm : true_particulate_matter values(Observed data after averaged by hr hours)
avg_data<-function(hr=4, remove_place = FALSE){
  if(hr==1){return()}
  X = PM25_all$dust
  
  ### erase place with too many NAs
  na_length = c()
  for(i in 1:dim(X)[2]){
    na_length[i] = sum(is.na(X[,i]))
  }
  if(remove_place){
    na.rm = which(na_length>3000) # the treshold(3000) can be changed.
    X = X[,-na.rm]
    print(na.rm)
  }
  ### average
  for(i in 1:(nrow(X)/hr)){
    X[(hr*(i-1)+1), ] = apply(X[((hr*(i-1)+1):(hr*i)),],2, mean, na.rm=T)
  }
  X = (X[hr*(1:(nrow(X)/hr))-hr+1,])
  true_pm = X
  X.Rdate = strptime(rownames(X), "%Y%m%d%H")
  
  ### substract the mean effect
  mu_st<- matrix(nrow=length(unique(substr(rownames(X),5,10))),ncol=ncol(X))
  rownames(mu_st)<-unique(substr(rownames(X),5,10))
  colnames(mu_st) <- PM25_all$place.name
  if(remove_place){
    colnames(mu_st) <- PM25_all$place.name[-na.rm]
  }
  for(day in unique(substr(rownames(X),5,10))){
    day.idx = which(substr(rownames(X),5,10)==day)
    if(length(day.idx)==1){
      mu_st[day, ] = X[day.idx,]
    }else{
      mu_st[day, ]= apply(X[day.idx,], 2, mean, na.rm=TRUE)
    }
  }
  for(i in 1:ncol(mu_st)){
    print(which(is.na(mu_st[,i])))
  }
  #smoothed by moving average
  for(i in 1:dim(mu_st)[2]){
    library(pracma)
    mu_st[,i] = movavg(mu_st[,i],11 ,type = "s")
  }
  #calculate anomaly
  for(day in rownames(mu_st)){
    day.idx = which(substr(rownames(X),5,10)==day)
    X[day.idx,] = X[day.idx, ] -rep(1,length(day.idx))%*%t(mu_st[day,])
  }
  #missing imputatation
  for(i in 1:dim(X)[2]){
    n <- which(is.na(X[,i]))
    idx = 1
    while( idx <= length(n)){
      run = 1
      while(idx<length(n) & n[idx]+1 == n[idx+1]){
        idx = idx+1
        run = run+1
      }
      idx= idx+1
      start = n[idx-run]
      end = n[idx-1]
      if(run <= 6){ #연속한 NA가 6개 이하 -> linearly interpolate
        if(start != 1 & end != dim(X)[1]){
          X[(start:end),i] = X[(start-1), i] + (X[(end+1), i] - X[(start-1), i])/(run+1) * (1:run) 
        }else if(start ==1 & run ==1){
          X[start,i] = X[start+1, i]
        }
        #print(X[(start-1), i] + (X[(end+1), i] - X[(start-1), i])/(run+1) * (1:run) )
      }
    }
  }
  naidx=c() #else: REMOVE ALL
  for(i in 1:dim(X)[2]){
    #print(length(which(is.na(X[,i]))))
    naidx = c(naidx,which(is.na(X[,i])))
  }
  print(length(unique(naidx)))
  X = X[-unique(naidx), ]
  return(list(X = X, mu_st = mu_st, true_pm = true_pm))
}

avg4 = avg_data(hr=4)
X = avg4$X
X.mean = apply(X,2,mean)
X.std = apply(X, 2, var)
Te = nrow(X)
N = ncol(X)
X = (X-rep(1,Te)%*%t(X.mean))/rep(1,Te)%*%t(sqrt(X.std))
mu_st = avg4$mu_st
true_pm = avg4$true_pm
rm(avg4)
X.Rdate = strptime(rownames(X), "%Y%m%d%H") #vector of times of Rdate-type
#save(X,file = 'X_all.RDS')
save(dust_predicted.list, mons,arima.result.list, X, X.Rdate,mu_st, qfm.binary.list, twtot.list,train.idx.list,   file = "4hr spring_manydays.rda")
save(mons,loc,tau.vec,days, X, X.Rdate, mu_st, train.idx.list, file = "run4spring.rda")
##############################################################
##############################################################
######                    QFM fitting                   ######
##############################################################
##############################################################
# see doQR.r for the codes for find_r and do_qfm
X = readRDS('X_all.RDS')
r.vec = find_r(X)
pred.quant.i =do_qfm(X, r.vec, AB=0, save=1,fname="QFM.RDS") 
pred.quant.sort <- array(unlist(lapply(pred.quant.i, function(x){x$F %*% t(x$L)})), dim=c(nrow(X),ncol(X),length(tau.vec)))
pred.quant.sort <- aperm(apply(pred.quant.sort, c(1,2), sort), c(2,3,1))
#save(pred.quant.i,file="predquanti_rvec.rda")

xk <- -2.5+seq(0.02,8,0.02)
n.xk <- length(xk)
tau.vec <- seq(from=0.05, to=0.95, by=0.05)
#extreme 안쓰는 경우 여기서 끝
##############################################################
##############################################################
######             Cross Validation Result              ######
##############################################################
##############################################################
#extreme 쓰는 경우 - Hill estimator 찾기 위해 cross-validation
tau.vec.ex <- c(0.96, 0.97, 0.98, 0.99, 0.995, 0.999, 0.9995)
tau.vec.bind = c(tau.vec,tau.vec.ex)
gamma.cal <- function(data, k=floor(c*length(data)^(1/3))){
  data.sorted <- sort(data, dec=TRUE)
  return(mean(log(data.sorted[1:k]))-log(data.sorted[k]))
}

set.seed(2020)
samp.idx = sample(1:5,length(X),replace = TRUE) #for 5-fold cross-validation (to find Hill estimator)
table(samp.idx)
index.validation.list=list()
validation.list=list()
X.fill.list = list()
#####5개 경우 별로 validation set, X.fill 등을 만들고 list에 저장한다.
#X.fill은 X에서 validation set 제거 후 아래 방법에 따라 제거된 값을 다시 채워넣어준 값.
for(s.ind in 1:5){
  X = readRDS('X_all.RDS')
  index.validation=which(samp.idx==s.ind)
  validation = X[index.validation]
  X[index.validation] = NA
  index.validation.list[[s.ind]] = index.validation
  validation.list[[s.ind]] = validation
  #랜덤하게 validation set을 정했고, 이 값들을 채워주어야 함
  #missing imputation1 - as did before
  for(i in 1:dim(X)[2]){
    n <- which(is.na(X[,i]))
    idx = 1
    while( idx <= length(n)){
      run = 1
      while(idx<length(n) & n[idx]+1 == n[idx+1]){
        idx = idx+1
        run = run+1
      }
      idx= idx+1
      start = n[idx-run]
      end = n[idx-1]
      if(run <= 6){ #연속한 NA가 6개 이하 -> linearly interpolate 4시간은 3, 6시간은 4
        if(start != 1 & end != dim(X)[1]){
          X[(start:end),i] = X[(start-1), i] + (X[(end+1), i] - X[(start-1), i])/(run+1) * (1:run) 
        }else if(start ==1 & run ==1){
          X[start,i] = X[start+1, i]
        }
        #print(X[(start-1), i] + (X[(end+1), i] - X[(start-1), i])/(run+1) * (1:run) )
      }
    }
  }
  #else : missing imputation2 - use the average value of adjacent places & times
  loc = cbind(PM25_all$lon,PM25_all$lat) 
  globalMean = mean(X[which(!is.infinite(X))], na.rm=TRUE)
  for(j in 1:ncol(X)){ # for each location
    n <- which(is.na(X[,j]))
    dist.j <- drop(rdist.earth(x1=matrix(loc[j,],nrow=1),x2=loc,miles=FALSE)) # compute the distance between the j-th location and all other locations
    loc.nei <- which(dist.j<=10) # find its neighbors (within the specified radius of 5km)
    for(i in n){ # for each hour
      time.i <- i+c(-3:4) # define the day for i-th hour (i.e., temporal neighborhood)
      time.i <- time.i[time.i>=1 & time.i<=nrow(X)] # make sure the day contains valid indices
      tmp.mean <- mean(X[time.i,loc.nei], na.rm=T) #####
      if(is.infinite(tmp.mean)){
        tmp.mean = globalMean
      }
      X[i,j] <- tmp.mean
    }
  }
  X.fill = X
  X[index.validation] = NA
  X.fill.list[[s.ind]] = X.fill
}

do_validation<-function(s.ind,twtot.mat,c.vec){
  X.fill = X.fill.list[[s.ind]]
  validation = validation.list[[s.ind]]
  Te = nrow(X.fill)
  N = ncol(X.fill)
  poo.time <- matrix(rep(1:Te, N), ncol=N)
  poo.loc <- matrix(rep(1:N, Te), ncol=N, byrow=T) 
  index.time.validation.list = list()
  index.loc.validation.list = list()
  index.time.validation <- poo.time[index.validation.list[[s.ind]]]
  index.loc.validation <- poo.loc[index.validation.list[[s.ind]]]
  index.time.validation.list[[s.ind]]=index.time.validation
  index.loc.validation.list[[s.ind]]=index.loc.validation
  
  n.validation <- length(validation)
  ##Hill estimator를 찾는 validation
  prediction.list = list()
  pred.quant.i =do_qfm(X.fill, r.vec, AB=0, save=1,s.ind)
#  return()
  pred.quant.sort <- array(unlist(lapply(pred.quant.i, function(x){x$F %*% t(x$L)})), dim=c(nrow(X.fill),ncol(X.fill),length(tau.vec)))
  pred.quant.sort <- aperm(apply(pred.quant.sort, c(1,2), sort), c(2,3,1))
  
  twtot = rep(0,length(c.vec))
  for(c in c.vec){
    gamma.hat <- apply(X.fill, 2, gamma.cal, k=floor(c*Te^(1/3)))
    pred.quant.ex <- array(dim=c(Te,N,length(tau.vec.ex)))
    for(i in 1:length(tau.vec.ex)){
      tau.ex <- tau.vec.ex[i]
      pred.quant.ex[,,i] <- (((1-tau.vec[length(tau.vec)])/(1-tau.ex))^(matrix(rep(gamma.hat, Te), ncol=N, byrow=TRUE))*pred.quant.sort[,,length(tau.vec)])*(pred.quant.sort[,,length(tau.vec)]>=0) + pred.quant.sort[,,length(tau.vec)]*(pred.quant.sort[,,length(tau.vec)]<0)
    }
    tau.vec.bind <- c(tau.vec, tau.vec.ex)
    
    # idx.index.validation : which index belongs to the current cluster among validation indices? 
    # For the members of idx.index.validation, inversion of estimated quantile process is executed. 
    # prediction.forsave is the final prediction object at the current cluster.
    idx.index.invalidation <- NULL
    prediction <- matrix(nrow=n.validation, ncol=n.xk) #xk 별로 누적분포함수 값
    pred.quant.final <- matrix(nrow=n.validation, ncol=length(tau.vec)+length(tau.vec.ex))#tau.vec 별로 quantile 값
    for(i in 1:N){
      idx.index <- which(i==index.loc.validation)
      idx.index.invalidation <- c(idx.index.invalidation, idx.index)
      for(k in idx.index){
        poo <- index.time.validation[k]
        for(j in 1:length(tau.vec)){
          pred.quant.final[k, j] <- pred.quant.sort[poo, i, j]
        }
        for(j in 1:length(tau.vec.ex)){
          pred.quant.final[k, length(tau.vec)+j] <- pred.quant.ex[poo, i, j]
        }
        prediction[k,] <- sapply(xk, cdf.est, pred.quant=pred.quant.final[k,], tau.vec=tau.vec.bind)
      }
    }
    for(i in 1:dim(prediction)[1]){
      twtot[which(c==c.vec)] = twtot[which(c==c.vec)] + sum((prediction[i,]-(validation.list[[s.ind]][i] <= xk))^2 *pnorm((xk-1.5)))/100
    }
    twtot[which(c==c.vec)] = twtot[which(c==c.vec)]/dim(prediction)[1]
    prediction.list[[which(c == c.vec)]] = prediction
    #saveRDS(prediction.forsave, file=filename)
  }
  twtot.mat[s.ind,]=twtot
  filename <- paste("prediction(sind=",s.ind,"hr=4).rda", sep="")
  save(r.vec, pred.quant.i, prediction.list,  file = filename)
  return(twtot.mat)
}

c.vec = 1+seq(0.1,1,0.2)
twtot.mat = matrix(nrow=5,ncol=length(c.vec))
for(i in 1:5){
  twtot.mat = do_validation(i,twtot.mat,c.vec)}
which.min(apply(twtot.mat,2,mean))

colnames(twtot.mat)=c.vec
twtot.save = twtot.mat


save(dust_predicted.list, mons,arima.result.list, X, X.Rdate,mu_st, qfm.binary.list, twtot.list,train.idx.list,   file = "4hr spring_manydays.rda")
save(mons,loc,tau.vec,days, X, X.Rdate, mu_st, train.idx.list, file = "run4spring.rda")
save(mu_st, pred.quant.i,pred.quant.sort,pred.quant.ex,pred.quant.final,prediction, r.vec,twtot.mat,X.Rdate, file = "all_place_4hr.rda")


##############################################################
##############################################################
######                     Result                      ######
##############################################################
##############################################################
par(family="AppleGothic")
#Figure 1 in poster
#퀀타일 추정값 ##figure 1 in poster
pred.quant.final = pred.quant.sort
par(mfrow=c(2,1))
par(mar=c(3.5,3.5,2,1))#space around figure!!
par(oma=c(1,1,1,1)) #outer margin
par(mgp=c(2.5,1,0)) #축과의 거리
idx = which(year(X.Rdate)==2019) 
for(place in c(77)){#station 77 : Gwanak
  cl <- rainbow(3)
  plot(X[idx,place],type = 'l', xaxt="n", xlab="time", ylab = "Anomaly",main=sprintf("Gwanak-gu"))
  axis(1,at =1:length(X[idx,1]) , labels=substr(rownames(X[idx,]),3,8))
  for (tau in c(1,3)) {
    lines(pred.quant.final[idx,place,c(1,10,19)[tau]], col =  c("blue","green","red")[tau],type="l")
  }
  legend("topright", legend = tau.vec.bind[c(1,19)], col = c("blue","red"), lwd = 1, cex = 0.5, title = "quantile")
  #lines(X[idx,place],type="l")
  my.quant=rep(0,Te)
  for(i in idx){
    my.quant[i] = cdf.est(X[i,place], pred.quant=pred.quant.final[i,place,1:19], tau.vec)
  }
  plot(my.quant[idx],type='l', xaxt="n", xlab="time", ylab = "Probability")
  abline(h=0.95, col = "Red")
  axis(1,at =1:length(X[idx,1]) , labels=substr(rownames(X[idx,]),3,8))
  #Sys.sleep(10)
}

#Figure 2 in poster
###factor loading 지도에 나타내고 싶을 때
#############################
library(ggmap)
library(ggplot2)
library(raster)
library(rgeos)
library(maptools)
library(rgdal)
#################map그림

#시군구 정보 홈페이지에서 다운로드
maps <- shapefile("./SIG_201703/TL_SCCO_SIG.shp")
maps <- spTransform(maps, CRSobj = CRS('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')) 
maps$SIG_KOR_NM <- iconv(maps$SIG_KOR_NM, from = "CP949", to = "UTF-8", sub = NA, mark = TRUE, toRaw = FALSE)
new_map <- fortify(maps, region = 'SIG_CD')
new_map$id <- as.numeric(new_map$id)
seoul_map <- new_map[new_map$id <= 11740,] 
order =c(2,1,3,5,4,7,6, 8,10,12,13,14,16,17,19,20,21,23,22,24,25,18,9,15,11)
maps$SIG_KOR_NM[order]
for(i in 1:length(seoul_map$id)){
  seoul_map$id[ i]=(KPM25_seoul$place.code[which(seoul_map$id[i] == as.numeric(maps[[1]][order]))] )
}

seoul.idx=which(PM25_all$place.code %in% KPM25_seoul$place.code)
p=list()
for( ind in c(1,19)){
  for(lambda in 1:2){
    factor_95 = pred.quant.i[[ind]]
    loading = data.frame(id = KPM25_seoul$place.code, loadings = (ifelse(lambda==2 & ind==19,-1,1)*factor_95$L[seoul.idx,lambda]))
    P_merge <- merge(seoul_map, loading, by='id')
    p[[lambda+2*(which(ind == c(1,19))-1)]] = ggplot() + geom_polygon(data = P_merge, aes(x=long, y=lat, group=group, fill =loadings))+
      xlab(paste(tau.vec[ind]*100, "% quantile, ",lambda,ifelse(lambda==1,"st","nd")," factor loading", sep=""))+
      theme(#axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),axis.text.y=element_blank(),
        axis.ticks.x=element_blank(),axis.ticks.y=element_blank())
    #+ guides(x=FALSE)
  }
}
grid.arrange(grobs=p,ncol=2)
for(i in 10*1:3){
  plot(X[,i],type='l')
  plot(X_all[,i],type='l')
}
#factor score
par(mfrow=c(4,3))
par(mar=c(3.5,3.5,1,1))#space around figure!!
par(oma=c(3,2,1,1.5)) #outer margin
par(mgp=c(2.5,1,0)) #축과의 거리
idx = which(year(X.Rdate)==2019 )
for(ind in c(1,7,13,19)){
  for(lambda in 1:3){
    factor_95 = pred.quant.i[[ind]]
    plot(factor_95$F[idx,lambda],type = 'l', xaxt="n",xlab=paste(tau.vec[ind]*100, "% quantile, ",lambda,ifelse(lambda==1,"st","nd")," factor score", sep=""),ylab="factor")
    axis(1,at =1:length(X.Rdate[idx]) , labels=substr(rownames(X)[idx],3,8))
  }
}

####EVA 2019 Springer 논문 그림
####draw Figure 6 
#factor: 각 tau별 Fhat and Lhat
#prediction: 각 validation point별 xk별 분포(F(-10 ~ 30))
library(lubridate)

Nex<-function(u){
  ntemp <- rep(0,Te) 
  for(i in 1:nrow(q)){
    ntemp[i] =sum(q[i,]>u)
  }
  nmonth = matrix(0,nrow=5,ncol=12)
  for(yr in 2015:2019){
    for(i in 1:12){
      nmonth[(yr-2014),i] = sum(ntemp[which(month(X.Rdate)==i & year(X.Rdate)==yr )])
    }
  }
  nmonth = as.data.frame(nmonth)
  rownames(nmonth)=2015:2019
  colnames(nmonth)=1:12
  nmonth$row = 2015:2019
  return(nmonth)
}  

par(mfrow=c(3,1))
for(u in c(2,3,4)){
  nmonth = Nex(u)
  library(reshape2)
  new_vv <- melt(nmonth, id.vars="row")
  bxpdat <- as.data.frame(boxplot(value~variable, data=new_vv, xlab = "month", ylab = "counts")[c("out", "group")])
  title(sprintf("Median counts of extreme PM2.5 (u = %d) for each month",u))
  #Get matches with boxplot data
  text_guide <- do.call(rbind, apply(bxpdat, 1, 
                                     function(x) new_vv[new_vv$value==x[1]&new_vv$variable==x[2], ]))
  #Add labels
  with(text_guide, text(x=as.numeric(variable)+0.2, y=value, labels=row))
}

####plot Figure 7
####Spatial distribution of predicted 95% quantiles at selected moments

library(maps)
par(mfrow=c(1,2))
daylist = c(as.integer(median(which(month(X.Rdate)==10 & year(X.Rdate)==2015)))+73,
            as.integer(median(which(month(X.Rdate)==12 & year(X.Rdate)==2017)))+179, 
            as.integer(median(which(month(X.Rdate)==1 & year(X.Rdate)==2019)))-47, 
            as.integer(median(which(month(X.Rdate)==3 & year(X.Rdate)==2019)))-236)
for(day in daylist){
  library(maps)
  map("world", xlim=c(126.35, 129.45), ylim=c(33.4, 38), fill=T, col="gray") #xlim: 위도 #ylim: 경도
  quilt.plot(loc[,1], loc[,2], q[day,], add=T)
  #text(loc[,1], loc[,2], labels=PM25_all$place.name, cex=1, pos=3,col=1)
  title(paste(X.Rdate[day], " predicted")) 
  
  map("world", xlim=c(126.35, 129.45), ylim=c(33.4, 38), fill=T, col="gray")
  quilt.plot(loc[,1], loc[,2], q[day,], add=T)
  #text(loc[,1], loc[,2], labels=PM25_all$place.name, cex=1, pos=3,col=1)
  title(paste(X.Rdate[day], " real"))
}
day = as.integer(median(which(month(X.Rdate)==10 & year(X.Rdate)==2015)))+73


####Spatial Extreme
###Figure 5
seoul.idx=which(PM25_all$place.code %in% KPM25_seoul$place.code)
Nsp<-function(u,tau){
  q = pred.quant.final[,,tau]
  ntemp <- rep(0,N) 
  for(i in 1:ncol(q)){
    ntemp[i] =sum(q[,i]>u)
  }
  return(ntemp)
}
p=list()
for(ind in c(19,26)){
  u=4
  nsp = data.frame(id = KPM25_seoul$place.code, counts = log(1+Nsp(u,ind)[seoul.idx],base=10),text=1:25)
  P_merge <- merge(seoul_map, nsp, by='id')
  p[[which(ind==c(19,26))]] = ggplot() + geom_polygon(data = P_merge, aes(x=long, y=lat, group=group, fill =counts)) +
    theme(axis.title.x=element_blank(),axis.title.y=element_blank(),
          axis.text.x=element_blank(),axis.text.y=element_blank(),
          axis.ticks.x=element_blank(),axis.ticks.y=element_blank())+# geom_text(aes(x=KPM25_seoul$lon, y=KPM25_seoul$lat, label=1:25,color = "red", hjust=-0.5))
    geom_text(aes(x=KPM25_seoul$lon[c(9,13)]-0.01, y=KPM25_seoul$lat[c(9,13)]+0.01, label=c("B","A")),color = "red", hjust=-0.5)
  #+ ggtitle(paste("Log Count of Extremes" when u = ", u, ", tau = ",tau.vec.bind[ind])) #+scale_fill_distiller(palette = "Set3")
}
grid.arrange(grobs=p,ncol=2)
which(KPM25_seoul$place.code[13] == PM25_all$place.code)
par(mfrow=c(1,2))
for(u in c(10,15,20)){
  map("world", xlim=c(126.75, 127.15), ylim=c(37.4, 37.7), fill=T, col="gray")
  quilt.plot(loc[,1], loc[,2], log(1+Nsp(u),base=10), add=T,nx=20,ny=20)
  text(loc[,1], loc[,2], labels=PM25_all$place.name, cex=1, pos=3,col=1)
  title(paste("Log Count of Extremes when u = ", u))
  
}


