
##5-10 cm as example
library(sp)
library(GA)
library(lattice)
library(rgdal)
library(raster)
library(maptools)
library(shapefiles)
library(corrplot)
library(gstat)
library(raster)
library(automap)
library(MASS)
library(soiltexture)
library(aqp)
library(plyr)
library(ithir)
library(hydroGOF)
library(ggplot2)
library(ggsci)
library(reshape2)
library(ncdf4) 
library(maptools) 
library(terra)
library(randomForest)
library(caret)
library(reticulate)
use_condaenv("r-reticulate")
py_config()
sklearn <- import("sklearn")
np <- import("numpy")
config <- py_config()
config$numpy

time_cal<-Sys.time()

df<-read.csv("SMAR_data_new.csv")
para<-read.csv("indicator_all.csv")
para_dat<-read.csv("para_spatial.csv")

site_info<-data.frame(Site=para$Site)
site_info$Site<-as.factor(site_info$Site)
site_name<-levels(site_info$Site)


path<-"~ECdata_1km/EC_no_Time"
ECfile_1km_no_time <- dir(path, pattern = ".tif$", full.names = TRUE)
for(i in 1 : length(ECfile_1km_no_time)){
  assign(paste("ec", i, sep = ""),raster(ECfile_1km_no_time[i]))
}
ec_raster_1km<-ec1
for(j in 2: length(ECfile_1km_no_time)){
  ec_raster_1km<-stack(ec_raster_1km,get(paste("ec",j,sep="")))
}

ec_raster_1km$Sand1<-ec_raster_1km$Sand1/10
ec_raster_1km$Sand2<-ec_raster_1km$Sand2/10
ec_raster_1km$Sand3<-ec_raster_1km$Sand3/10
ec_raster_1km$Silt1<-ec_raster_1km$Silt1/10
ec_raster_1km$Silt2<-ec_raster_1km$Silt2/10
ec_raster_1km$Silt3<-ec_raster_1km$Silt3/10
ec_raster_1km$Clay1<-ec_raster_1km$Clay1/10
ec_raster_1km$Clay2<-ec_raster_1km$Clay2/10
ec_raster_1km$Clay3<-ec_raster_1km$Clay3/10

path.ST<-"~soil texture"
ST.file <- dir(path.ST, pattern = ".tif$", full.names = TRUE)
ST_0_5<-readAll(raster(ST.file [1]))
ST_5_10<-readAll(raster(ST.file [2]))
ST_10_20<-readAll(raster(ST.file [3]))
ST_20_40<-readAll(raster(ST.file [4]))

n_0_5<-ST_0_5
names(n_0_5)<-"n_0_5"
n_0_5@data@values[which(n_0_5@data@values==0)]<-1
for(texture.count in 1:12){
n_0_5@data@values[which(n_0_5@data@values==texture.count)]<-para_dat[para_dat$id==texture.count,]$n
}
n_0_5<-n_0_5/1

n_5_10<-ST_5_10
names(n_5_10)<-"n_5_10"
n_5_10@data@values[which(n_5_10@data@values==0)]<-1
for(texture.count in 1:12){
n_5_10@data@values[which(n_5_10@data@values==texture.count)]<-para_dat[para_dat$id==texture.count,]$n
}
n_5_10<-n_5_10/1

n_10_20<-ST_10_20
names(n_10_20)<-"n_10_20"
n_10_20@data@values[which(n_10_20@data@values==0)]<-1
for(texture.count in 1:12){
n_10_20@data@values[which(n_10_20@data@values==texture.count)]<-para_dat[para_dat$id==texture.count,]$n
}
n_10_20<-n_10_20/1

n_20_40<-ST_20_40
names(n_20_40)<-"n_20_40"
n_20_40@data@values[which(n_20_40@data@values==0)]<-1
for(texture.count in 1:12){
n_20_40@data@values[which(n_20_40@data@values==texture.count)]<-para_dat[para_dat$id==texture.count,]$n
}
n_20_40<-n_20_40/1

ss=1
point<-subset(df,Day==ss)
point<-cbind(point,para)
if("NST25"%in%point$Site)
{point <- point[point$Site != "NST25", ]}

trainX<-as.matrix(point[,11:32])
trainY<-as.matrix(point[,40:41])
#colnames(trainY)

MultiOutputRegressor <- sklearn$multioutput$MultiOutputRegressor
LinearRegression <- sklearn$linear_model$LinearRegression
base_regressor <- LinearRegression()
multi_output_model <- MultiOutputRegressor(base_regressor)
multi_output_model$fit(trainX, trainY)


X_test<- as.matrix(getValues(ec_raster_1km))
non_na_indices <- complete.cases(X_test)
X_test_non_na <- X_test[non_na_indices, ]
y_pred_non_na <- multi_output_model$predict(X_test_non_na)

print("predict done")
y_pred <- matrix(NA, nrow = nrow(X_test), ncol = ncol(y_pred_non_na))
y_pred[non_na_indices, ] <- y_pred_non_na


ncol_raster<-ec_raster_1km@ncols
nrow_raster<-ec_raster_1km@nrows
colnames(y_pred)<-colnames(trainY)

sc1_sp_matrix <- matrix(y_pred[, 1], nrow = nrow_raster, ncol = ncol_raster, byrow = TRUE)
sc1_sp <- raster(sc1_sp_matrix)
extent(sc1_sp) <- extent(ec_raster_1km)
crs(sc1_sp) <- crs(ec_raster_1km)

sw2_sp_matrix <- matrix(y_pred[, 2], nrow = nrow_raster, ncol = ncol_raster, byrow = TRUE)
sw2_sp <- raster(sw2_sp_matrix)
extent(sw2_sp) <- extent(ec_raster_1km)
crs(sw2_sp) <- crs(ec_raster_1km)
##we have already obtaind sc1 and sw2 (spatial prediction) 


##===RK
trainY2 <- multi_output_model$predict(trainX)
resid_sc1<-trainY[,1]-trainY2[,1]
resid_sw2<-trainY[,2]-trainY2[,2]

coordinates(point)=~x+y
proj4string(point)<- proj4string(ec_raster_1km)
point@data <- cbind(point@data, resid_1 = resid_sc1,resid_2 = resid_sw2)
formMod_1 <- resid_1 ~ 1
formMod_2 <- resid_2 ~ 1
rstPixDF <- as(ec_raster_1km[[1]], "SpatialPixelsDataFrame")


  variog_1 <- autofitVariogram(formMod_1, point, model  = c("Exp","Gau","Sph"))
  variog_2 <- autofitVariogram(formMod_2, point,model  = c("Exp","Gau","Sph"))

RK_1 <- krige(formula =  formMod_1 ,
              locations = point, 
              model = variog_1$var_model,
              newdata = rstPixDF )
RK_2 <- krige(formula =  formMod_2 ,
              locations = point, 
              model =variog_2$var_model,
              newdata = rstPixDF )

resid_sc1_pred <- as(RK_1, "RasterLayer")
resid_sw2_pred <- as(RK_2, "RasterLayer")


sc1_sp<-sc1_sp+resid_sc1_pred
sw2_sp<-sw2_sp +resid_sw2_pred


for(ss in 1: 122){#=================================================

#path.SM.1km<-"~RF"
path.SM.1km<-"~GLM"
SM.file <- dir(path.SM.1km, pattern = ".tif$", full.names = TRUE)
SM_0_5<-raster(SM.file[ss])

s1<-SM_0_5/n_0_5
Zr1<-abs(0-5)
Zr2<-abs(5-10)
n1<-n_0_5 
n2<-n_5_10 
sc1<-sc1_sp 
sw2<-sw2_sp

v.temp1<-1-SM_0_5/0.55
v.temp1@data@values[which(v.temp1@data@values<0)]<-0.01
v.temp2<--log(v.temp1)
V<-v.temp2/0.976 
b=(n1*Zr1)/((1-sw2)*n2*Zr2) 
a_temp<-1/((1-sw2)*n2*Zr2)
a=V*a_temp 

#==========================================================================
setwd("~GLMRK")
writeRaster(sc1,"sc1_5_10.tif",overwrite=T)
writeRaster(sw2,"sw2_5_10.tif",overwrite=T)
writeRaster(b,"b_5_10.tif",overwrite=T)

a_plot_name<-paste("a_5_10_",substring(names(SM_0_5),8),".tif",sep="")
writeRaster(a,a_plot_name,overwrite=T)
#==========================================================================


if(ss==1){S2_SMAR<-stack(sc1,sc1)}

judge.s1sc<-sc1
judge.s1sc@data@values<-0
judge.s1sc<-judge.s1sc/1
judge.temp<-(s1-sc1)/1
judge.s1sc@data@values[which(judge.temp@data@values>=0)]<-1
judge.s1sc<-judge.s1sc/1

S2_SMAR[[2]]=sw2+((S2_SMAR[[1]]-sw2)*exp(-a))+(judge.s1sc*(1-sw2)*b*(s1-sc1))
S2_SMAR[[1]]<-S2_SMAR[[2]]

sm_SMAR<-S2_SMAR[[2]]*n2
output_deep_SM<-paste("SM_5_10_",substring(names(SM_0_5),8),".tif",sep="")

#writeRaster(sm_SMAR,output_deep_SM,overwrite=T)
print(paste("The", ss, "day is done!----", round((ss/122)*100,2),"%"))
}
Sys.time()-time_cal



