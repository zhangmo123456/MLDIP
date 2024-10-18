##----------------------------------------------------------------##
##--This code is based on the SMAR model to predict the root zone soil moisture --##
##--Use the real data of the Qinghai-Tibet Plateau for experiments --##
##--All steps, including data extraction, modeling, and statistical calculations, are completed in R --##
##----------------------------------------------------------------##
#This version:
# uses a lookup table for soil porosity
# uses 0-5cm (surface) to predict 5-10/10-20/20-40cm (deep)
# uses soil moisture [expanded data set] for root zone prediction (2016.06.01-2016.09.30)
# adds site calibration of genetic algorithm-GA method
# continues to improve the GA genetic algorithm, only adjusting the
# two parameters sc1 and sw2, and taking the b parameter as a function of sw2
# extracts the calibration parameters of 72 sites (sc1, sw2), optimized the drawing program V2
# uses 0-5cm (surface layer), predicts 5-10cm, 10-20cm, 20-40cm (deep layer), and 
# extracts the calibration parameters of each layer

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
ubrmse<-function(a,b){
(sum(((a-mean(a))-(b-mean(b)))^2)/length(a))^0.5
}

#Import soil moisture data/site information
setwd("~SMAR")
#SM_data<-read.csv("data_Tibet_Obs.csv")
#SM_data$Site<-as.factor(SM_data$Site)
#Import parameters (according to soil texture type)
para_dat<-read.csv("initial_para.csv")
df_all<-data.frame()
##According to the data of each day, extract 1000m covariates

for(ss in 1:122){

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

path.smrf.1km<-"~RF"
smrf.file <- dir(path.smrf.1km, pattern = ".tif$", full.names = TRUE)
SMrf<-raster(smrf.file[ss])
SMrf@data@names<-"SMrf"
ec_raster_1km<-stack(ec_raster_1km,SMrf)

path.smglm.1km<-"~GLM"
smglm.file <- dir(path.smglm.1km, pattern = ".tif$", full.names = TRUE)
SMglm<-raster(smglm.file[ss])
SMglm@data@names<-"SMglm"
ec_raster_1km<-stack(ec_raster_1km,SMglm)


#Processing data set dimensions and coefficient units (mainly soil particle size)
ec_raster_1km$Sand1<-ec_raster_1km$Sand1/10
ec_raster_1km$Sand2<-ec_raster_1km$Sand2/10
ec_raster_1km$Sand3<-ec_raster_1km$Sand3/10
ec_raster_1km$Silt1<-ec_raster_1km$Silt1/10
ec_raster_1km$Silt2<-ec_raster_1km$Silt2/10
ec_raster_1km$Silt3<-ec_raster_1km$Silt3/10
ec_raster_1km$Clay1<-ec_raster_1km$Clay1/10
ec_raster_1km$Clay2<-ec_raster_1km$Clay2/10
ec_raster_1km$Clay3<-ec_raster_1km$Clay3/10

site_info<-read.csv("Site_Tibet_Obs.csv")
lat_y<-site_info$lat
lon_x<-site_info$lon
xydata<-data.frame(x=lon_x,y=lat_y)
#site_info$Site<-as.factor(site_info$Site)
coordinates(site_info)=~lon+lat
proj4string(site_info)<- proj4string(ec_raster_1km)
point.ec<-extract(ec_raster_1km,site_info)

site_name<-site_info$Site
SM_data<-read.csv("data_Tibet_Obs.csv")
#SM_data$Site<-as.factor(SM_data$Site)
SM_data_temp<-subset(SM_data,Day==ss)
SM_data_temp_all<-cbind(SM_data_temp, xydata,  point.ec)
rownames(SM_data_temp_all)<-NULL
df_all<-rbind(df_all,SM_data_temp_all)
#df_all$Site<-as.factor(df_all$Site)
print(ss)
}
write.csv(df_all,"SMAR_data.csv") 

##------5-10cm
setwd("~SMAR")
data<-read.csv("SMAR_data.csv")
para_dat<-read.csv("initial_para.csv")
data<-na.omit(data)
data$Site<-as.factor(data$Site)
site_name<-levels(data$Site)

df_vali<-data.frame()
data_all<-data.frame()
for(i in 1:length(site_name)){

df<-subset(data,Site==site_name[i])
sand_0_5<-df$Sand1[1]
sand_5_10<-df$Sand2[1]
sand_10_20<-df$Sand2[1]
sand_20_40<-df$Sand3[1]

silt_0_5<-df$Silt1[1]
silt_5_10<-df$Silt2[1]
silt_10_20<-df$Silt2[1]
silt_20_40<-df$Silt3[1]

clay_0_5<-df$Clay1[1]
clay_5_10<-df$Clay2[1]
clay_10_20<-df$Clay2[1]
clay_20_40<-df$Clay3[1]


#calculate soil texture for each layer (0-5 5-10 10-20 20-40cm)
soil.psf_0_5<-data.frame("CLAY"=clay_0_5, "SILT"=silt_0_5,"SAND"=sand_0_5)
soil.tex_0_5<-TT.points.in.classes( tri.data=soil.psf_0_5, class.sys="USDA.TT",PiC.type="t",tri.sum.tst=F )
soil.psf_5_10<-data.frame("CLAY"=clay_5_10, "SILT"=silt_5_10,"SAND"=sand_5_10)
soil.tex_5_10<-TT.points.in.classes( tri.data=soil.psf_5_10, class.sys="USDA.TT",PiC.type="t",tri.sum.tst=F )
soil.psf_10_20<-data.frame("CLAY"=clay_10_20, "SILT"=silt_10_20,"SAND"=sand_10_20)
soil.tex_10_20<-TT.points.in.classes( tri.data=soil.psf_10_20, class.sys="USDA.TT",PiC.type="t",tri.sum.tst=F )
soil.psf_20_40<-data.frame("CLAY"=clay_20_40, "SILT"=silt_20_40,"SAND"=sand_20_40)
soil.tex_20_40<-TT.points.in.classes( tri.data=soil.psf_20_40, class.sys="USDA.TT",PiC.type="t",tri.sum.tst=F )


#Experimental design: Based on the surface soil moisture, 
#the soil moisture of each lower layer is gradually predicted based on SMAR

obs_sm_0_5<-df$SM_5cm
obs_sm_5_10<-df$SM_10cm
obs_sm_10_20<-df$SM_20cm
obs_sm_20_40<-df$SM_40cm

#Divided into 2 layers, upper layer and lower layer, the upper layer 
#is the known layer, the lower layer is the predicted layer
#This version uses 0-5cm to predict 5-10cm

por_0_5<-para_dat[which(para_dat$texture==soil.tex_0_5),]$n
por_5_10<-para_dat[which(para_dat$texture==soil.tex_5_10),]$n
por_10_20<-para_dat[which(para_dat$texture==soil.tex_10_20),]$n
por_20_40<-para_dat[which(para_dat$texture==soil.tex_20_40),]$n

#Upper layer 0-5cm, lower layer 5-10cm
obs_sm_0_5<-df$SM_5cm
s1<-obs_sm_0_5/as.numeric(por_0_5)
Zr1<-abs(0-5)
Zr2<-abs(5-10)
n1<-por_0_5
n2<-por_5_10
sc1<-para_dat[which(para_dat$texture==soil.tex_0_5),]$sc
sw2<-para_dat[which(para_dat$texture==soil.tex_5_10),]$sw
V=-log(1-obs_sm_0_5/0.55)/0.976

b=as.numeric((n1*Zr1)/((1-sw2)*n2*Zr2))
a_temp<-1/((1-sw2)*n2*Zr2)
a=V*as.numeric(a_temp)

S2_SMAR=rep(0,length(s1))
S2_SMAR[1]=sc1
#t=2
for(t in 2:length(s1)){
judge.s1sc<-0
if(s1[t]>=sc1){judge.s1sc<-1}
S2_SMAR[t]=sw2+((S2_SMAR[t-1]-sw2)*exp(-a[t]))+(judge.s1sc*(1-sw2)*b*(s1[t]-sc1))
}
sm_SMAR<-S2_SMAR*n2
day<-df$Day

#plot version 3
df.ggplot<-data.frame(obs_sm_0_5=obs_sm_0_5,obs_sm_5_10=obs_sm_5_10,
                                       sm_SMAR=sm_SMAR,day=day)
mycolor<-pal_cosmic("signature_substitutions")(6)

##Use GA algorithm to calibrate the model
#After overall simulation, calibrate to minimize the average rmse
#Only adjust 2 of the four parameters: [sc1, sw2]
#According to the mechanism, the predicted layers are close, SC1 should be greater than SW2
#Therefore, this limitation should be considered

observed <- obs_sm_5_10
ini_sc1<-sc1
ini_sw2<-sw2

fitness_function <- function(params) {
sc1_calib <- params[1]
sw2_calib <- params[2]

if(sc1_calib>sw2_calib){
s1_calib=s1#Surface soil moisture saturation
n2_calib=n2#Porosity of the lower soil
Zr1_calib<-abs(0-5)
Zr2_calib<-abs(5-10)
n1_calib<-por_0_5
b_calib=as.numeric((n1_calib*Zr1_calib)/((1-sw2_calib)*n2_calib*Zr2_calib))
a_calib=a
S2_SMAR=rep(0,length(s1))
S2_SMAR[1]=s1_calib[1]#The initial lower layer prediction value is represented by the surface layer sc (i.e. sc1)
days_length<-length(a)
for(t in 2:days_length){
judge.s1sc<-0
if(s1_calib[t]>=sc1_calib){judge.s1sc<-1}
S2_SMAR[t]=sw2_calib+((S2_SMAR[t-1]-sw2_calib)*exp(-a_calib[t]))+
                      (judge.s1sc*(1-sw2_calib)*b_calib*(s1_calib[t]-sc1_calib))
}
model_predicted<-S2_SMAR*n2_calib
rmse <- rmse(observed,model_predicted)
return(-rmse)
}
else(return(-9999))
}
#Optimize the parameters of the initial value range

lower_bound <- c( sc1_calib = ini_sc1*0.5, sw2_calib =  ini_sw2*0.5)
upper_bound <-c( sc1_calib = min(ini_sc1*1.5,1), sw2_calib = min(ini_sw2*1.5,1))

result <- ga(type = "real-valued", fitness = fitness_function, 
                     lower = lower_bound, upper = upper_bound,
                      popSize =50,maxiter = 1000,optim=T,monitor = F)
result@solution
result@fitnessValue

params_pre<-result@solution
params_pre<-params_pre[1,]
sc1_calib <- params_pre[1]
sw2_calib <- params_pre[2]
Zr1_calib<-abs(0-5)
Zr2_calib<-abs(5-10)
n1_calib<-por_0_5
s1_calib=s1#Surface soil moisture saturation
n2_calib=n2#Porosity of the lower soil

#b_calib=as.numeric((n1_calib*Zr1_calib)/((1-sw2_calib)*n2_calib*Zr2_calib))
#a_temp_calib<-1/((1-sw2_calib)*n2_calib*Zr2_calib)
#a_calib=V*as.numeric(a_temp_calib)
b_calib=as.numeric((n1_calib*Zr1_calib)/((1-sw2_calib)*n2_calib*Zr2_calib))
a_calib<-a


S2_SMAR=rep(0,length(s1))
S2_SMAR[1]=sc1#The initial lower layer prediction value is represented by the surface layer sc (i.e. sc1)
days_length<-length(a)
for(t in 2:days_length){
judge.s1sc<-0
if(s1_calib[t]>=sc1_calib){judge.s1sc<-1}
S2_SMAR[t]=sw2_calib+((S2_SMAR[t-1]-sw2_calib)*exp(-a_calib[t]))+
                      (judge.s1sc*(1-sw2_calib)*b_calib*(s1_calib[t]-sc1_calib))
}
model_predicted<-S2_SMAR*n2_calib

df$SMAR_5_10<-sm_SMAR
df$SMAR_cb_5_10<-model_predicted
data_all<-rbind(data_all,df)

##--------------------------------------------------------------------------------------------------------
df.ggplot<-data.frame(obs_sm_0_5=obs_sm_0_5,obs_sm_5_10=obs_sm_5_10,
                                       sm_SMAR=sm_SMAR,day=day)
df.ggplot2<-data.frame(day=day,obs_sm_0_5=obs_sm_0_5,obs_sm_5_10=obs_sm_5_10,
                                       sm_SMAR_ori=sm_SMAR,sm_SMAR_cali=model_predicted)

dd.plot.new<-melt(df.ggplot2,id = "day")
max.value<-max(dd.plot.new$value)+0.1
min.value<-min(dd.plot.new$value)-0.05
mylabels = c("Observation (0-5 cm)", "Observation (5-10 cm)", "SMAR_ori (5-10 cm)","SMAR_cb (5-10 cm)")

plot_title_new<-paste(site_name[i]," at 5-10 cm (a=",round(mean(a_calib),2),",b=",round(b_calib,2),
                                      ",sc1=",round(sc1_calib,2),",sw2=",round(sw2_calib,2),")",sep="")
  new_plot<-ggplot(data=dd.plot.new, aes(x=day, y=value, colour=variable,linetype = variable)) +
  geom_line(size=1.1)+
  ylim(0,max.value)+
  scale_x_continuous(breaks = seq(day[1],day[length(day)], by=13)) +
  #scale_y_continuous(breaks = seq(0,max.value, by=0.1)) +
  labs(title=plot_title_new,y =expression("Soil moisture (m"^3*"/m"^3*")"), x = "Day") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
  panel.grid.major.y = element_line(color = "gray",size = 0.5,linetype = 2),
   panel.grid.major.x = element_blank(),
        legend.position = c(0.99,0.99),
        legend.key = element_rect(colour = NA, fill = NA),
        legend.justification=c(1, 1))+
  theme(axis.text.x = element_text( color="black",  size=18),axis.title.x=element_text(size=18),
          axis.text.y = element_text( color="black",  size=18),axis.title.y=element_text(size=18))+
  theme(legend.title = element_blank(),
              plot.title = element_text(size=18),
              legend.key.height = unit(0.7, 'cm'), 
              legend.key.width = unit(1, 'cm'), 
              legend.text = element_text(size=15)) +
  scale_linetype_manual(values=c("dashed", "solid", "solid", "solid"),labels=mylabels) +
  scale_color_manual(values=c("gray91", "black", "FireBrick1", "green4"), labels=mylabels)


observed<-observed[-1:(-10)]
model_predicted<-model_predicted[-1:(-10)]
sm_SMAR<-sm_SMAR[-1:(-10)]

me_after<-round(me(observed,model_predicted),3)
me_before<-round(me(observed,sm_SMAR),3)
rmse_after<-round(rmse(observed,model_predicted),3)
rmse_before<-round(rmse(observed,sm_SMAR),3)
bR2_after<-round(br2(observed,model_predicted),3)
bR2_before<-round(br2(observed,sm_SMAR),3)
ubrmse_after<-round(ubrmse(observed,model_predicted),3)
ubrmse_before<-round(ubrmse(observed,sm_SMAR),3)
nse_after<-round(nse(observed,model_predicted),3)
nse_before<-round(nse(observed,sm_SMAR),3)
  
metric_label = paste(" Before Calibration: \n"," ME (5-10 cm) = ",me_before,"\n"," RMSE (5-10 cm) = ",rmse_before,"\n", 
                                    " bR2 (5-10 cm) = ",bR2_before,"\n"," ubRMSE (5-10 cm) = ",ubrmse_before,"\n", 
                                    " After Calibration: ","\n"," ME (5-10 cm) = ",me_after,"\n"," RMSE (5-10 cm) = ",rmse_after,"\n",
                                    " bR2 (5-10 cm) = ",bR2_after,"\n"," ubRMSE (5-10 cm) = ",ubrmse_after,sep="")
metric_label2 = paste(" ME = ",me_before," / ",me_after,"\n",
                                      " RMSE = ",rmse_before," / ",rmse_after,"\n",
                                      " ubRMSE = ",ubrmse_before," / ",ubrmse_after,"\n",
                                      " R² = ",bR2_before," / ",bR2_after,"\n",
                                      " NSE = ",nse_before," / ",nse_after,sep="")
new_plot2<-new_plot+
  annotate("text", x = -Inf, y = Inf, colour = "black", size=5.5,label = metric_label2 ,hjust = -0.05,vjust = 1.1)
output_plot_name<-paste(site_name[i],"_Depth_5_10.jpg",sep="")
#ggsave(output_plot_name,plot = new_plot2,dpi = 600,height = 6.5,width = 8)

para_summary_temp<-c(site_name[i],sc1,sw2,round(sc1_calib,3),round(sw2_calib,3),
                                           me_before,rmse_before,ubrmse_before,nse_before,bR2_before,
                                           me_after,rmse_after,ubrmse_after,nse_after,bR2_after)

para_summary_temp<-c(site_name[i],sc1,sw2,round(sc1_calib,3),round(sw2_calib,3),
                                           me_before,rmse_before,ubrmse_before,nse_before,bR2_before,
                                           me_after,rmse_after,ubrmse_after,nse_after,bR2_after)
names(para_summary_temp)<-c("Site","Sc1","Sw2","Sc1_cb","Sw2_cb",
                                                         "ME","RMSE","ubRMSE","NSE","R2",
                                                         "ME_cb","RMSE_cb","ubRMSE_cb","NSE_cb","R2_cb")
df_vali<-rbind(df_vali,para_summary_temp)
colnames(df_vali)<-c("Site","Sc1","Sw2","Sc1_cb","Sw2_cb",
                                      "ME","RMSE","ubRMSE","NSE","R2",
                                      "ME_cb","RMSE_cb","ubRMSE_cb","NSE_cb","R2_cb")
print(i)
new_plot2
}
#write.csv(df_vali,"indicator_5_10.csv")
#write.csv(data_all,"data_all_5_10.csv")

##10-20cm

setwd("~SMAR")
data<-read.csv("SMAR_data.csv")
para_dat<-read.csv("initial_para.csv")
data<-na.omit(data)
data$Site<-as.factor(data$Site)
site_name<-levels(data$Site)

df_vali<-data.frame()
data_all<-data.frame()
for(i in 1:length(site_name)){

df<-subset(data,Site==site_name[i])
sand_0_5<-df$Sand1[1]
sand_5_10<-df$Sand2[1]
sand_10_20<-df$Sand2[1]
sand_20_40<-df$Sand3[1]

silt_0_5<-df$Silt1[1]
silt_5_10<-df$Silt2[1]
silt_10_20<-df$Silt2[1]
silt_20_40<-df$Silt3[1]

clay_0_5<-df$Clay1[1]
clay_5_10<-df$Clay2[1]
clay_10_20<-df$Clay2[1]
clay_20_40<-df$Clay3[1]


#calculate soil texture for each layer (0-5 5-10 10-20 20-40cm)
soil.psf_0_5<-data.frame("CLAY"=clay_0_5, "SILT"=silt_0_5,"SAND"=sand_0_5)
soil.tex_0_5<-TT.points.in.classes( tri.data=soil.psf_0_5, class.sys="USDA.TT",PiC.type="t",tri.sum.tst=F )
soil.psf_5_10<-data.frame("CLAY"=clay_5_10, "SILT"=silt_5_10,"SAND"=sand_5_10)
soil.tex_5_10<-TT.points.in.classes( tri.data=soil.psf_5_10, class.sys="USDA.TT",PiC.type="t",tri.sum.tst=F )
soil.psf_10_20<-data.frame("CLAY"=clay_10_20, "SILT"=silt_10_20,"SAND"=sand_10_20)
soil.tex_10_20<-TT.points.in.classes( tri.data=soil.psf_10_20, class.sys="USDA.TT",PiC.type="t",tri.sum.tst=F )
soil.psf_20_40<-data.frame("CLAY"=clay_20_40, "SILT"=silt_20_40,"SAND"=sand_20_40)
soil.tex_20_40<-TT.points.in.classes( tri.data=soil.psf_20_40, class.sys="USDA.TT",PiC.type="t",tri.sum.tst=F )

obs_sm_0_5<-df$SM_5cm
obs_sm_5_10<-df$SM_10cm
obs_sm_10_20<-df$SM_20cm
obs_sm_20_40<-df$SM_40cm

por_0_5<-para_dat[which(para_dat$texture==soil.tex_0_5),]$n
por_5_10<-para_dat[which(para_dat$texture==soil.tex_5_10),]$n
por_10_20<-para_dat[which(para_dat$texture==soil.tex_10_20),]$n
por_20_40<-para_dat[which(para_dat$texture==soil.tex_20_40),]$n

obs_sm_0_5<-df$SM_5cm
s1<-obs_sm_0_5/as.numeric(por_0_5)
Zr1<-abs(0-5)
Zr2<-abs(10-20)
n1<-por_0_5
n2<-por_10_20
sc1<-para_dat[which(para_dat$texture==soil.tex_0_5),]$sc
sw2<-para_dat[which(para_dat$texture==soil.tex_10_20),]$sw
V=-log(1-obs_sm_0_5/0.55)/0.976

b=as.numeric((n1*Zr1)/((1-sw2)*n2*Zr2))
a_temp<-1/((1-sw2)*n2*Zr2)
a=V*as.numeric(a_temp)

S2_SMAR=rep(0,length(s1))
S2_SMAR[1]=sc1
#t=2
for(t in 2:length(s1)){
judge.s1sc<-0
if(s1[t]>=sc1){judge.s1sc<-1}
S2_SMAR[t]=sw2+((S2_SMAR[t-1]-sw2)*exp(-a[t]))+(judge.s1sc*(1-sw2)*b*(s1[t]-sc1))
}
sm_SMAR<-S2_SMAR*n2
day<-df$Day

#plot version 3
df.ggplot<-data.frame(obs_sm_0_5=obs_sm_0_5,obs_sm_10_20=obs_sm_10_20,
                                       sm_SMAR=sm_SMAR,day=day)
mycolor<-pal_cosmic("signature_substitutions")(6)

observed <- obs_sm_10_20
ini_sc1<-sc1
ini_sw2<-sw2

fitness_function <- function(params) {
sc1_calib <- params[1]
sw2_calib <- params[2]

if(sc1_calib>sw2_calib){
s1_calib=s1
n2_calib=n2
Zr1_calib<-abs(0-5)
Zr2_calib<-abs(10-20)
n1_calib<-por_0_5
b_calib=as.numeric((n1_calib*Zr1_calib)/((1-sw2_calib)*n2_calib*Zr2_calib))
a_calib=a
S2_SMAR=rep(0,length(s1))
S2_SMAR[1]=s1_calib[1]
days_length<-length(a)
for(t in 2:days_length){
judge.s1sc<-0
if(s1_calib[t]>=sc1_calib){judge.s1sc<-1}
S2_SMAR[t]=sw2_calib+((S2_SMAR[t-1]-sw2_calib)*exp(-a_calib[t]))+
                      (judge.s1sc*(1-sw2_calib)*b_calib*(s1_calib[t]-sc1_calib))
}
model_predicted<-S2_SMAR*n2_calib
rmse <- rmse(observed,model_predicted)
return(-rmse)
}
else(return(-9999))
}

lower_bound <- c( sc1_calib = ini_sc1*0.5, sw2_calib =  ini_sw2*0.5)
upper_bound <-c( sc1_calib = min(ini_sc1*1.5,1), sw2_calib = min(ini_sw2*1.5,1))

result <- ga(type = "real-valued", fitness = fitness_function, 
                     lower = lower_bound, upper = upper_bound,
                      popSize =50,maxiter = 1000,optim=T,monitor = F)
result@solution
result@fitnessValue

params_pre<-result@solution
params_pre<-params_pre[1,]
sc1_calib <- params_pre[1]
sw2_calib <- params_pre[2]
Zr1_calib<-abs(0-5)
Zr2_calib<-abs(10-20)
n1_calib<-por_0_5
s1_calib=s1
n2_calib=n2

#b_calib=as.numeric((n1_calib*Zr1_calib)/((1-sw2_calib)*n2_calib*Zr2_calib))
#a_temp_calib<-1/((1-sw2_calib)*n2_calib*Zr2_calib)
#a_calib=V*as.numeric(a_temp_calib)
b_calib=as.numeric((n1_calib*Zr1_calib)/((1-sw2_calib)*n2_calib*Zr2_calib))
a_calib<-a


S2_SMAR=rep(0,length(s1))
S2_SMAR[1]=sc1
days_length<-length(a)
for(t in 2:days_length){
judge.s1sc<-0
if(s1_calib[t]>=sc1_calib){judge.s1sc<-1}
S2_SMAR[t]=sw2_calib+((S2_SMAR[t-1]-sw2_calib)*exp(-a_calib[t]))+
                      (judge.s1sc*(1-sw2_calib)*b_calib*(s1_calib[t]-sc1_calib))
}
model_predicted<-S2_SMAR*n2_calib

df$SMAR_10_20<-sm_SMAR
df$SMAR_cb_10_20<-model_predicted
data_all<-rbind(data_all,df)

##--------------------------------------------------------------------------------------------------------
df.ggplot<-data.frame(obs_sm_0_5=obs_sm_0_5,obs_sm_10_20=obs_sm_10_20,
                                       sm_SMAR=sm_SMAR,day=day)
df.ggplot2<-data.frame(day=day,obs_sm_0_5=obs_sm_0_5,obs_sm_10_20=obs_sm_10_20,
                                       sm_SMAR_ori=sm_SMAR,sm_SMAR_cali=model_predicted)

dd.plot.new<-melt(df.ggplot2,id = "day")
max.value<-max(dd.plot.new$value)+0.1
min.value<-min(dd.plot.new$value)-0.05
mylabels = c("Observation (0-5 cm)", "Observation (10-20 cm)", "SMAR_ori (10-20 cm)","SMAR_cb (10-20 cm)")

plot_title_new<-paste(site_name[i]," at 10-20 cm (a=",round(mean(a_calib),2),",b=",round(b_calib,2),
                                      ",sc1=",round(sc1_calib,2),",sw2=",round(sw2_calib,2),")",sep="")
  new_plot<-ggplot(data=dd.plot.new, aes(x=day, y=value, colour=variable,linetype = variable)) +
  geom_line(size=1.1)+
  ylim(0,max.value)+
  scale_x_continuous(breaks = seq(day[1],day[length(day)], by=13)) +
  #scale_y_continuous(breaks = seq(0,max.value, by=0.1)) +
  labs(title=plot_title_new,y =expression("Soil moisture (m"^3*"/m"^3*")"), x = "Day") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
  panel.grid.major.y = element_line(color = "gray",size = 0.5,linetype = 2),
   panel.grid.major.x = element_blank(),
        legend.position = c(0.99,0.99),
        legend.key = element_rect(colour = NA, fill = NA),
        legend.justification=c(1, 1))+
  theme(axis.text.x = element_text( color="black",  size=18),axis.title.x=element_text(size=18),
          axis.text.y = element_text( color="black",  size=18),axis.title.y=element_text(size=18))+
  theme(legend.title = element_blank(),
              plot.title = element_text(size=18),
              legend.key.height = unit(0.7, 'cm'), 
              legend.key.width = unit(1, 'cm'), 
              legend.text = element_text(size=15)) +
  scale_linetype_manual(values=c("dashed", "solid", "solid", "solid"),labels=mylabels) +
  scale_color_manual(values=c("gray91", "black", "FireBrick1", "green4"), labels=mylabels)


observed<-observed[-1:(-10)]
model_predicted<-model_predicted[-1:(-10)]
sm_SMAR<-sm_SMAR[-1:(-10)]

me_after<-round(me(observed,model_predicted),3)
me_before<-round(me(observed,sm_SMAR),3)
rmse_after<-round(rmse(observed,model_predicted),3)
rmse_before<-round(rmse(observed,sm_SMAR),3)
bR2_after<-round(br2(observed,model_predicted),3)
bR2_before<-round(br2(observed,sm_SMAR),3)
ubrmse_after<-round(ubrmse(observed,model_predicted),3)
ubrmse_before<-round(ubrmse(observed,sm_SMAR),3)
nse_after<-round(nse(observed,model_predicted),3)
nse_before<-round(nse(observed,sm_SMAR),3)
  
metric_label = paste(" Before Calibration: \n"," ME (10-20 cm) = ",me_before,"\n"," RMSE (10-20 cm) = ",rmse_before,"\n", 
                                    " bR2 (5-10 cm) = ",bR2_before,"\n"," ubRMSE (10-20 cm) = ",ubrmse_before,"\n", 
                                    " After Calibration: ","\n"," ME (10-20 cm) = ",me_after,"\n"," RMSE (10-20 cm) = ",rmse_after,"\n",
                                    " bR2 (10-20 cm) = ",bR2_after,"\n"," ubRMSE (10-20 cm) = ",ubrmse_after,sep="")
metric_label2 = paste(" ME = ",me_before," / ",me_after,"\n",
                                      " RMSE = ",rmse_before," / ",rmse_after,"\n",
                                      " ubRMSE = ",ubrmse_before," / ",ubrmse_after,"\n",
                                      " R² = ",bR2_before," / ",bR2_after,"\n",
                                      " NSE = ",nse_before," / ",nse_after,sep="")
new_plot2<-new_plot+
  annotate("text", x = -Inf, y = Inf, colour = "black", size=5.5,label = metric_label2 ,hjust = -0.05,vjust = 1.1)
output_plot_name<-paste(site_name[i],"_Depth_10_20.jpg",sep="")
#ggsave(output_plot_name,plot = new_plot2,dpi = 600,height = 6.5,width = 8)

para_summary_temp<-c(site_name[i],sc1,sw2,round(sc1_calib,3),round(sw2_calib,3),
                                           me_before,rmse_before,ubrmse_before,nse_before,bR2_before,
                                           me_after,rmse_after,ubrmse_after,nse_after,bR2_after)

para_summary_temp<-c(site_name[i],sc1,sw2,round(sc1_calib,3),round(sw2_calib,3),
                                           me_before,rmse_before,ubrmse_before,nse_before,bR2_before,
                                           me_after,rmse_after,ubrmse_after,nse_after,bR2_after)
names(para_summary_temp)<-c("Site","Sc1","Sw2","Sc1_cb","Sw2_cb",
                                                         "ME","RMSE","ubRMSE","NSE","R2",
                                                         "ME_cb","RMSE_cb","ubRMSE_cb","NSE_cb","R2_cb")
df_vali<-rbind(df_vali,para_summary_temp)
colnames(df_vali)<-c("Site","Sc1","Sw2","Sc1_cb","Sw2_cb",
                                      "ME","RMSE","ubRMSE","NSE","R2",
                                      "ME_cb","RMSE_cb","ubRMSE_cb","NSE_cb","R2_cb")
print(i)
new_plot2
}
write.csv(df_vali,"indicator_10_20.csv")
write.csv(data_all,"data_all_10_20.csv")

#20-40cm

setwd("~SMAR")
data<-read.csv("SMAR_data.csv")
para_dat<-read.csv("initial_para.csv")
data<-na.omit(data)
data$Site<-as.factor(data$Site)
site_name<-levels(data$Site)

df_vali<-data.frame()
data_all<-data.frame()
for(i in 1:length(site_name)){

df<-subset(data,Site==site_name[i])
sand_0_5<-df$Sand1[1]
sand_5_10<-df$Sand2[1]
sand_10_20<-df$Sand2[1]
sand_20_40<-df$Sand3[1]

silt_0_5<-df$Silt1[1]
silt_5_10<-df$Silt2[1]
silt_10_20<-df$Silt2[1]
silt_20_40<-df$Silt3[1]

clay_0_5<-df$Clay1[1]
clay_5_10<-df$Clay2[1]
clay_10_20<-df$Clay2[1]
clay_20_40<-df$Clay3[1]


#calculate soil texture for each layer (0-5 5-10 10-20 20-40cm)
soil.psf_0_5<-data.frame("CLAY"=clay_0_5, "SILT"=silt_0_5,"SAND"=sand_0_5)
soil.tex_0_5<-TT.points.in.classes( tri.data=soil.psf_0_5, class.sys="USDA.TT",PiC.type="t",tri.sum.tst=F )
soil.psf_5_10<-data.frame("CLAY"=clay_5_10, "SILT"=silt_5_10,"SAND"=sand_5_10)
soil.tex_5_10<-TT.points.in.classes( tri.data=soil.psf_5_10, class.sys="USDA.TT",PiC.type="t",tri.sum.tst=F )
soil.psf_10_20<-data.frame("CLAY"=clay_10_20, "SILT"=silt_10_20,"SAND"=sand_10_20)
soil.tex_10_20<-TT.points.in.classes( tri.data=soil.psf_10_20, class.sys="USDA.TT",PiC.type="t",tri.sum.tst=F )
soil.psf_20_40<-data.frame("CLAY"=clay_20_40, "SILT"=silt_20_40,"SAND"=sand_20_40)
soil.tex_20_40<-TT.points.in.classes( tri.data=soil.psf_20_40, class.sys="USDA.TT",PiC.type="t",tri.sum.tst=F )


obs_sm_0_5<-df$SM_5cm
obs_sm_5_10<-df$SM_10cm
obs_sm_10_20<-df$SM_20cm
obs_sm_20_40<-df$SM_40cm

por_0_5<-para_dat[which(para_dat$texture==soil.tex_0_5),]$n
por_5_10<-para_dat[which(para_dat$texture==soil.tex_5_10),]$n
por_10_20<-para_dat[which(para_dat$texture==soil.tex_10_20),]$n
por_20_40<-para_dat[which(para_dat$texture==soil.tex_20_40),]$n

obs_sm_0_5<-df$SM_5cm
s1<-obs_sm_0_5/as.numeric(por_0_5)
Zr1<-abs(0-5)
Zr2<-abs(20-40)
n1<-por_0_5
n2<-por_20_40
sc1<-para_dat[which(para_dat$texture==soil.tex_0_5),]$sc
sw2<-para_dat[which(para_dat$texture==soil.tex_20_40),]$sw
V=-log(1-obs_sm_0_5/0.55)/0.976

b=as.numeric((n1*Zr1)/((1-sw2)*n2*Zr2))
a_temp<-1/((1-sw2)*n2*Zr2)
a=V*as.numeric(a_temp)

S2_SMAR=rep(0,length(s1))
S2_SMAR[1]=sc1
#t=2
for(t in 2:length(s1)){
judge.s1sc<-0
if(s1[t]>=sc1){judge.s1sc<-1}
S2_SMAR[t]=sw2+((S2_SMAR[t-1]-sw2)*exp(-a[t]))+(judge.s1sc*(1-sw2)*b*(s1[t]-sc1))
}
sm_SMAR<-S2_SMAR*n2
day<-df$Day

#plot version 3
df.ggplot<-data.frame(obs_sm_0_5=obs_sm_0_5,obs_sm_20_40=obs_sm_20_40,
                                       sm_SMAR=sm_SMAR,day=day)
mycolor<-pal_cosmic("signature_substitutions")(6)

observed <- obs_sm_20_40
ini_sc1<-sc1
ini_sw2<-sw2

fitness_function <- function(params) {
sc1_calib <- params[1]
sw2_calib <- params[2]

if(sc1_calib>sw2_calib){
s1_calib=s1
n2_calib=n2
Zr1_calib<-abs(0-5)
Zr2_calib<-abs(20-40)
n1_calib<-por_0_5
b_calib=as.numeric((n1_calib*Zr1_calib)/((1-sw2_calib)*n2_calib*Zr2_calib))
a_calib=a
S2_SMAR=rep(0,length(s1))
S2_SMAR[1]=s1_calib[1]
days_length<-length(a)
for(t in 2:days_length){
judge.s1sc<-0
if(s1_calib[t]>=sc1_calib){judge.s1sc<-1}
S2_SMAR[t]=sw2_calib+((S2_SMAR[t-1]-sw2_calib)*exp(-a_calib[t]))+
                      (judge.s1sc*(1-sw2_calib)*b_calib*(s1_calib[t]-sc1_calib))
}
model_predicted<-S2_SMAR*n2_calib
rmse <- rmse(observed,model_predicted)
return(-rmse)
}
else(return(-9999))
}

lower_bound <- c( sc1_calib = ini_sc1*0.5, sw2_calib =  ini_sw2*0.5)
upper_bound <-c( sc1_calib = min(ini_sc1*1.5,1), sw2_calib = min(ini_sw2*1.5,1))

result <- ga(type = "real-valued", fitness = fitness_function, 
                     lower = lower_bound, upper = upper_bound,
                      popSize =50,maxiter = 1000,optim=T,monitor = F)
result@solution
result@fitnessValue

params_pre<-result@solution
params_pre<-params_pre[1,]
sc1_calib <- params_pre[1]
sw2_calib <- params_pre[2]
Zr1_calib<-abs(0-5)
Zr2_calib<-abs(20-40)
n1_calib<-por_0_5
s1_calib=s1
n2_calib=n2

#b_calib=as.numeric((n1_calib*Zr1_calib)/((1-sw2_calib)*n2_calib*Zr2_calib))
#a_temp_calib<-1/((1-sw2_calib)*n2_calib*Zr2_calib)
#a_calib=V*as.numeric(a_temp_calib)
b_calib=as.numeric((n1_calib*Zr1_calib)/((1-sw2_calib)*n2_calib*Zr2_calib))
a_calib<-a


S2_SMAR=rep(0,length(s1))
S2_SMAR[1]=sc1
days_length<-length(a)
for(t in 2:days_length){
judge.s1sc<-0
if(s1_calib[t]>=sc1_calib){judge.s1sc<-1}
S2_SMAR[t]=sw2_calib+((S2_SMAR[t-1]-sw2_calib)*exp(-a_calib[t]))+
                      (judge.s1sc*(1-sw2_calib)*b_calib*(s1_calib[t]-sc1_calib))
}
model_predicted<-S2_SMAR*n2_calib

df$SMAR_20_40<-sm_SMAR
df$SMAR_cb_20_40<-model_predicted
data_all<-rbind(data_all,df)

##--------------------------------------------------------------------------------------------------------
df.ggplot<-data.frame(obs_sm_0_5=obs_sm_0_5,obs_sm_20_40=obs_sm_20_40,
                                       sm_SMAR=sm_SMAR,day=day)
df.ggplot2<-data.frame(day=day,obs_sm_0_5=obs_sm_0_5,obs_sm_20_40=obs_sm_20_40,
                                       sm_SMAR_ori=sm_SMAR,sm_SMAR_cali=model_predicted)

dd.plot.new<-melt(df.ggplot2,id = "day")
max.value<-max(dd.plot.new$value)+0.1
min.value<-min(dd.plot.new$value)-0.05
mylabels = c("Observation (0-5 cm)", "Observation (20-40 cm)", "SMAR_ori (20-40 cm)","SMAR_cb (20-40 cm)")

plot_title_new<-paste(site_name[i]," at 20-40 cm (a=",round(mean(a_calib),2),",b=",round(b_calib,2),
                                      ",sc1=",round(sc1_calib,2),",sw2=",round(sw2_calib,2),")",sep="")
  new_plot<-ggplot(data=dd.plot.new, aes(x=day, y=value, colour=variable,linetype = variable)) +
  geom_line(size=1.1)+
  ylim(0,max.value)+
  scale_x_continuous(breaks = seq(day[1],day[length(day)], by=13)) +
  #scale_y_continuous(breaks = seq(0,max.value, by=0.1)) +
  labs(title=plot_title_new,y =expression("Soil moisture (m"^3*"/m"^3*")"), x = "Day") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
  panel.grid.major.y = element_line(color = "gray",size = 0.5,linetype = 2),
   panel.grid.major.x = element_blank(),
        legend.position = c(0.99,0.99),
        legend.key = element_rect(colour = NA, fill = NA),
        legend.justification=c(1, 1))+
  theme(axis.text.x = element_text( color="black",  size=18),axis.title.x=element_text(size=18),
          axis.text.y = element_text( color="black",  size=18),axis.title.y=element_text(size=18))+
  theme(legend.title = element_blank(),
              plot.title = element_text(size=18),
              legend.key.height = unit(0.7, 'cm'), 
              legend.key.width = unit(1, 'cm'), 
              legend.text = element_text(size=15)) +
  scale_linetype_manual(values=c("dashed", "solid", "solid", "solid"),labels=mylabels) +
  scale_color_manual(values=c("gray91", "black", "FireBrick1", "green4"), labels=mylabels)


observed<-observed[-1:(-10)]
model_predicted<-model_predicted[-1:(-10)]
sm_SMAR<-sm_SMAR[-1:(-10)]

me_after<-round(me(observed,model_predicted),3)
me_before<-round(me(observed,sm_SMAR),3)
rmse_after<-round(rmse(observed,model_predicted),3)
rmse_before<-round(rmse(observed,sm_SMAR),3)
bR2_after<-round(br2(observed,model_predicted),3)
bR2_before<-round(br2(observed,sm_SMAR),3)
ubrmse_after<-round(ubrmse(observed,model_predicted),3)
ubrmse_before<-round(ubrmse(observed,sm_SMAR),3)
nse_after<-round(nse(observed,model_predicted),3)
nse_before<-round(nse(observed,sm_SMAR),3)
  

metric_label2 = paste(" ME = ",me_before," / ",me_after,"\n",
                                      " RMSE = ",rmse_before," / ",rmse_after,"\n",
                                      " ubRMSE = ",ubrmse_before," / ",ubrmse_after,"\n",
                                      " R² = ",bR2_before," / ",bR2_after,"\n",
                                      " NSE = ",nse_before," / ",nse_after,sep="")
new_plot2<-new_plot+
  annotate("text", x = -Inf, y = Inf, colour = "black", size=5.5,label = metric_label2 ,hjust = -0.05,vjust = 1.1)
output_plot_name<-paste(site_name[i],"_Depth_20_40.jpg",sep="")
#ggsave(output_plot_name,plot = new_plot2,dpi = 600,height = 4.5,width = 13)

para_summary_temp<-c(site_name[i],sc1,sw2,round(sc1_calib,3),round(sw2_calib,3),
                                           me_before,rmse_before,ubrmse_before,nse_before,bR2_before,
                                           me_after,rmse_after,ubrmse_after,nse_after,bR2_after)

para_summary_temp<-c(site_name[i],sc1,sw2,round(sc1_calib,3),round(sw2_calib,3),
                                           me_before,rmse_before,ubrmse_before,nse_before,bR2_before,
                                           me_after,rmse_after,ubrmse_after,nse_after,bR2_after)
names(para_summary_temp)<-c("Site","Sc1","Sw2","Sc1_cb","Sw2_cb",
                                                         "ME","RMSE","ubRMSE","NSE","R2",
                                                         "ME_cb","RMSE_cb","ubRMSE_cb","NSE_cb","R2_cb")
df_vali<-rbind(df_vali,para_summary_temp)
colnames(df_vali)<-c("Site","Sc1","Sw2","Sc1_cb","Sw2_cb",
                                      "ME","RMSE","ubRMSE","NSE","R2",
                                      "ME_cb","RMSE_cb","ubRMSE_cb","NSE_cb","R2_cb")
print(i)
#new_plot2
}
write.csv(df_vali,"indicator_20_40.csv")
write.csv(data_all,"data_all_20_40.csv")


