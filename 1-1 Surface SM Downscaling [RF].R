
library(ncdf4) 
library(raster) 
library(rgdal) 
library(ggplot2) 
library(maptools) 
library(terra)
library(randomForest)
library(caret)
library(raster)

df_all_set<-data.frame()
setwd("~SMAP")
path.sm<-""~SMAP""
SM_file_36km <- dir(path.sm, pattern = ".tif$", full.names = TRUE)
time_cal<-Sys.time()

for(ss in 1:length(SM_file_36km)){

#Read soil moisture 36km data
setwd("~SMAP")
path.sm<-""~SMAP""
SM_file_36km <- dir(path.sm, pattern = ".tif$", full.names = TRUE)

#--Read and merge 36km environment variables (non-time changes)
path<-"~EC/ECdata_36km/EC_no_Time"
ECfile_36km_no_time <- dir(path, pattern = ".tif$", full.names = TRUE)
for(i in 1 : length(ECfile_36km_no_time)){
  assign(paste("ec", i, sep = ""),raster(ECfile_36km_no_time[i]))
}
ec_raster_36km<-ec1
for(j in 2: length(ECfile_36km_no_time)){
  ec_raster_36km<-stack(ec_raster_36km,get(paste("ec",j,sep="")))
}

#Merge 36km other daily environment variables
path.lstd.36km<-"~ECdata_36km/EC_with_Time/LST_D"
lstd.file <- dir(path.lstd.36km, pattern = ".tif$", full.names = TRUE)
LST_D<-raster(lstd.file[ss])
LST_D@data@names<-"LST_D"
ec_raster_36km<-stack(ec_raster_36km,LST_D)

path.lstn.36km<-"~ECdata_36km/EC_with_Time/LST_N"
lstn.file <- dir(path.lstn.36km, pattern = ".tif$", full.names = TRUE)
LST_N<-raster(lstn.file[ss])
LST_N@data@names<-"LST_N"
ec_raster_36km<-stack(ec_raster_36km,LST_N)

path.ndvi.36km<-"~ECdata_36km/EC_with_Time/NDVI"
ndvi.file <- dir(path.ndvi.36km, pattern = ".tif$", full.names = TRUE)
NDVI<-raster(ndvi.file[ss])
NDVI@data@names<-"NDVI"
ec_raster_36km<-stack(ec_raster_36km,NDVI)

path.prcp.36km<-"~ECdata_36km/EC_with_Time/Prcp"
prcp.file <- dir(path.prcp.36km, pattern = ".tif$", full.names = TRUE)
Prcp<-raster(prcp.file[ss])
Prcp@data@names<-"Prcp"
ec_raster_36km<-stack(ec_raster_36km,Prcp)


#Data Integration
SMdata<-readAll(raster(SM_file_36km[ss]))
df<-data.frame(SM=SMdata@data@values)
name.ec36km<-names(ec_raster_36km)

for(ec36.count in 1:length(name.ec36km)){
temp.ec36<-readAll(ec_raster_36km[[name.ec36km[ec36.count]]])
temp.ec36.values<-temp.ec36@data@values
df<-data.frame(df,temp.ec36.values)
colnames(df)[(ec36.count+1)]<-name.ec36km[ec36.count]
}

df<-na.omit(df)
#print("data frame build complete")
name.ec36km<-names(ec_raster_36km)
name.ec36km<-c("SM",name.ec36km)
colnames(df)<-name.ec36km

Model<-train( df[,2:ncol(df)],df[,'SM'],method='rf',ntree=1000,importance=T,
                  trControl = trainControl(method = "none"))

#--Read and merge 1km environment variables (non-time-varying)
path<-"~ECdata_1km/EC_no_Time"
ECfile_1km_no_time <- dir(path, pattern = ".tif$", full.names = TRUE)
for(i in 1 : length(ECfile_1km_no_time)){
  assign(paste("ec", i, sep = ""),raster(ECfile_1km_no_time[i]))
}
ec_raster_1km<-ec1
for(j in 2: length(ECfile_1km_no_time)){
  ec_raster_1km<-stack(ec_raster_1km,get(paste("ec",j,sep="")))
}
#Merge 1km other daily environment variables
path.lstd.1km<-"~ECdata_1km/EC_with_Time/LST_D"
lstd.file <- dir(path.lstd.1km, pattern = ".tif$", full.names = TRUE)
LST_D<-raster(lstd.file[ss])
LST_D@data@names<-"LST_D"
ec_raster_1km<-stack(ec_raster_1km,LST_D)

path.lstn.1km<-"~ECdata_1km/EC_with_Time/LST_N"
lstn.file <- dir(path.lstn.1km, pattern = ".tif$", full.names = TRUE)
LST_N<-raster(lstn.file[ss])
LST_N@data@names<-"LST_N"
ec_raster_1km<-stack(ec_raster_1km,LST_N)

path.ndvi.1km<-"~ECdata_1km/EC_with_Time/NDVI"
ndvi.file <- dir(path.ndvi.1km, pattern = ".tif$", full.names = TRUE)
NDVI<-raster(ndvi.file[ss])
NDVI@data@names<-"NDVI"
ec_raster_1km<-stack(ec_raster_1km,NDVI)

path.prcp.1km<-"~ECdata_1km/EC_with_Time/Prcp"
prcp.file <- dir(path.prcp.1km, pattern = ".tif$", full.names = TRUE)
Prcp<-raster(prcp.file[ss])
Prcp@data@names<-"Prcp"
ec_raster_1km<-stack(ec_raster_1km,Prcp)

Pre <- predict(ec_raster_1km, Model)

#writeRaster(Pre,paste("SM_1km_",substring(SM_file_36km[ss],54),sep=""))
print(paste("Mission complete for ",as.character(ss),
        " (",as.character(round(ss/length(SM_file_36km)*100,1)),"%)",sep=""))
}
Sys.time()-time_cal


