#The R code R code used in Zhang et al. (2022), Species distribution model identifies influence of climatic constraints at the leading edge of a native insect outbreak

####Climatic data####
#Ftp://ftp.nrcan.gc.ca/pub/outgoing/North_America_Historical_Grids/geotif
#read all .tif files, ~1.2km resolution
#Data directory: under each category folder (MaxT, MinT, pcp...), a folder for each year (1955,1956...), monthly .tif files under each year folder.
path = "D:/MUN/Data/ClimateData/MaxT"

MaxT01 <- list.files(path = path, recursive=T,pattern = "maxt60_01.tif",full.names = T)
MaxT02 <- list.files(path = path, recursive=T,pattern = "maxt60_02.tif",full.names = T)
MaxT03 <- list.files(path = path, recursive=T,pattern = "maxt60_03.tif",full.names = T)
MaxT04 <- list.files(path = path, recursive=T,pattern = "maxt60_04.tif",full.names = T)
MaxT05 <- list.files(path = path, recursive=T,pattern = "maxt60_05.tif",full.names = T)
MaxT06 <- list.files(path = path, recursive=T,pattern = "maxt60_06.tif",full.names = T)
MaxT07 <- list.files(path = path, recursive=T,pattern = "maxt60_07.tif",full.names = T)
MaxT08 <- list.files(path = path, recursive=T,pattern = "maxt60_08.tif",full.names = T)
MaxT09 <- list.files(path = path, recursive=T,pattern = "maxt60_09.tif",full.names = T)
MaxT10 <- list.files(path = path, recursive=T,pattern = "maxt60_10.tif",full.names = T)
MaxT11 <- list.files(path = path, recursive=T,pattern = "maxt60_11.tif",full.names = T)
MaxT12 <- list.files(path = path, recursive=T,pattern = "maxt60_12.tif",full.names = T)

path = "D:/MUN/Data/ClimateData/MinT"
MinT01 <- list.files(path = path, recursive=T,pattern = "mint60_01.tif",full.names = T)
MinT02 <- list.files(path = path, recursive=T,pattern = "mint60_02.tif",full.names = T)
MinT03 <- list.files(path = path, recursive=T,pattern = "mint60_03.tif",full.names = T)
MinT04 <- list.files(path = path, recursive=T,pattern = "mint60_04.tif",full.names = T)
MinT05 <- list.files(path = path, recursive=T,pattern = "mint60_05.tif",full.names = T)
MinT06 <- list.files(path = path, recursive=T,pattern = "mint60_06.tif",full.names = T)
MinT07 <- list.files(path = path, recursive=T,pattern = "mint60_07.tif",full.names = T)
MinT08 <- list.files(path = path, recursive=T,pattern = "mint60_08.tif",full.names = T)
MinT09 <- list.files(path = path, recursive=T,pattern = "mint60_09.tif",full.names = T)
MinT10 <- list.files(path = path, recursive=T,pattern = "mint60_10.tif",full.names = T)
MinT11 <- list.files(path = path, recursive=T,pattern = "mint60_11.tif",full.names = T)
MinT12 <- list.files(path = path, recursive=T,pattern = "mint60_12.tif",full.names = T)

path = "D:/MUN/Data/ClimateData/pcp"
pcp01 <- list.files(path = path, recursive=T,pattern = "pcp60_01.tif",full.names = T)
pcp02 <- list.files(path = path, recursive=T,pattern = "pcp60_02.tif",full.names = T)
pcp03 <- list.files(path = path, recursive=T,pattern = "pcp60_03.tif",full.names = T)
pcp04 <- list.files(path = path, recursive=T,pattern = "pcp60_04.tif",full.names = T)
pcp05 <- list.files(path = path, recursive=T,pattern = "pcp60_05.tif",full.names = T)
pcp06 <- list.files(path = path, recursive=T,pattern = "pcp60_06.tif",full.names = T)
pcp07 <- list.files(path = path, recursive=T,pattern = "pcp60_07.tif",full.names = T)
pcp08 <- list.files(path = path, recursive=T,pattern = "pcp60_08.tif",full.names = T)
pcp09 <- list.files(path = path, recursive=T,pattern = "pcp60_09.tif",full.names = T)
pcp10 <- list.files(path = path, recursive=T,pattern = "pcp60_10.tif",full.names = T)
pcp11 <- list.files(path = path, recursive=T,pattern = "pcp60_11.tif",full.names = T)
pcp12 <- list.files(path = path, recursive=T,pattern = "pcp60_12.tif",full.names = T)

path = "D:/MUN/Data/ClimateData/sg"
sg01 <- list.files(path = path, recursive=T,pattern = "sg60_01.tif",full.names = T)#Julian date number of start of growing season
sg02 <- list.files(path = path, recursive=T,pattern = "sg60_02.tif",full.names = T)
sg03 <- list.files(path = path, recursive=T,pattern = "sg60_03.tif",full.names = T)
sg04 <- list.files(path = path, recursive=T,pattern = "sg60_04.tif",full.names = T)
sg05 <- list.files(path = path, recursive=T,pattern = "sg60_05.tif",full.names = T)
sg06 <- list.files(path = path, recursive=T,pattern = "sg60_06.tif",full.names = T)
sg07 <- list.files(path = path, recursive=T,pattern = "sg60_07.tif",full.names = T)
sg08 <- list.files(path = path, recursive=T,pattern = "sg60_08.tif",full.names = T)#gdd above base temp for period 1(3 mo prior to start of gs)
sg09 <- list.files(path = path, recursive=T,pattern = "sg60_09.tif",full.names = T)#gdd above base temp for period 1(first six weeks of growing season)
sg10 <- list.files(path = path, recursive=T,pattern = "sg60_10.tif",full.names = T)
sg11 <- list.files(path = path, recursive=T,pattern = "sg60_11.tif",full.names = T)
sg13 <- list.files(path = path, recursive=T,pattern = "sg60_13.tif",full.names = T)
sg14 <- list.files(path = path, recursive=T,pattern = "sg60_14.tif",full.names = T)
sg15 <- list.files(path = path, recursive=T,pattern = "sg60_15.tif",full.names = T)
sg16 <- list.files(path = path, recursive=T,pattern = "sg60_16.tif",full.names = T)

#stack MaxT in May and Jun (05,06) and year 1972-1982
MaxT <- stack(c(MaxT05[3:13],MaxT06[3:13])) #choose May and Jun (05,06) and year 1972-1982 (3:13)

#call in NL boundary and make a 20km buffer, transform projection to the raster
#!!don't transfer raster to vector as it need to recalculte
bd <- "D:/MUN/Data/Scratch/NL_Boundary.shp"
NL_BD <- shapefile(bd)

#combine the 1196 features in NL_BD to one feaure, it also removes data, for speed
NL_BD <- gUnion(NL_BD, NL_BD)
NL_BD_21N <- spTransform(NL_BD, projection(NLproj))

NL_BD_BUFF <- buffer(NL_BD_21N, width=20000) #TAKE TIME! cannot buffer if converting to raster projection first, raster CRS not projected

#change projetion of the boundary polygon to raster's
NL_BD_BUFF_proj <- spTransform(NL_BD_BUFF, projection(MaxT))

#crop the raster stack by the boundary
MaxT_NL <- crop(MaxT,NL_BD_BUFF_proj)

#calculate the mean of raster stack, write to disk
MaxTm <- mean(MaxT_NL)

#change raster projection, after cropping to save memory
MaxT.proj <- projectRaster(MaxTm, crs = NLproj)  

#create a blank raster with the same extent as the NL border layer
EmtRaster = raster(extent(MaxT.proj))
res(EmtRaster) <- c(2000,2000)
projection(EmtRaster) <- "+proj=utm +zone=21 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"

#change raster resolution to 2km
MaxT.agrt <- resample(MaxT.proj, EmtRaster, method='bilinear')

MaxT0506 <- writeRaster(MaxT.agrt,'MaxT0506.tif', overwrite=F)


#Repetitive: stack MaxT in Jun and Jul (06,07) and year 1972-1982
MaxT <- stack(c(MaxT06[3:13],MaxT07[3:13])) 
MaxT_NL <- crop(MaxT,NL_BD_BUFF_proj)
MaxTm <- mean(MaxT_NL)
MaxT.proj <- projectRaster(MaxTm, crs = NLproj) 
MaxT.agrt <- resample(MaxT.proj, EmtRaster, method='bilinear')
MaxT0607 <- writeRaster(MaxT.agrt,'MaxT0607.tif', overwrite=F)

#Repetitive: stack MaxT in Jul and Aug (07,08) and year 1972-1982
MaxT <- stack(c(MaxT07[3:13],MaxT08[3:13])) 
MaxT_NL <- crop(MaxT,NL_BD_BUFF_proj)
MaxTm <- mean(MaxT_NL)
MaxT.proj <- projectRaster(MaxTm, crs = NLproj) 
MaxT.agrt <- resample(MaxT.proj, EmtRaster, method='bilinear')
MaxT0708 <- writeRaster(MaxT.agrt,'MaxT0708.tif', overwrite=F)

##pcp
#Repetitive: stack pcp in May and Jun (05,06) and year 1972-1982
pcp <- stack(c(pcp05[3:13],pcp06[3:13])) 
pcp_NL <- crop(pcp,NL_BD_BUFF_proj)
pcpm <- mean(pcp_NL)
pcp.proj <- projectRaster(pcpm, crs = NLproj) 
pcp.agrt <- resample(pcp.proj, EmtRaster, method='bilinear')
pcp0506 <- writeRaster(pcp.agrt,'pcp0506.tif', overwrite=F)

#Repetitive: stack pcp in Jun and Jul (06,07) and year 1972-1982
pcp <- stack(c(pcp06[3:13],pcp07[3:13])) 
pcp_NL <- crop(pcp,NL_BD_BUFF_proj)
pcpm <- mean(pcp_NL)
pcp.proj <- projectRaster(pcpm, crs = NLproj) 
pcp.agrt <- resample(pcp.proj, EmtRaster, method='bilinear')
pcp0607 <- writeRaster(pcp.agrt,'pcp0607.tif', overwrite=F)

#Repetitive: stack pcp in Jul and Aug (07,08) and year 1972-1982
pcp <- stack(c(pcp07[3:13],pcp08[3:13])) 
pcp_NL <- crop(pcp,NL_BD_BUFF_proj)
pcpm <- mean(pcp_NL)
pcp.proj <- projectRaster(pcpm, crs = NLproj) 
pcp.agrt <- resample(pcp.proj, EmtRaster, method='bilinear')
pcp0708 <- writeRaster(pcp.agrt,'pcp0708.tif', overwrite=F)

par(mfrow=c(1,3))
plot(pcp0506)
plot(pcp0607)
plot(pcp0708)
dev.off()

##growing season
#Repetitive: Julian date number of start of growing season in year 1972-1982
sg <- stack(c(sg01[3:13])) 
sg_NL <- crop(sg,NL_BD_BUFF_proj)
sgm <- mean(sg_NL)
sg.proj <- projectRaster(sgm, crs = NLproj) 
sg.agrt <- resample(sg.proj, EmtRaster, method='bilinear')
sg_date <- writeRaster(sg.agrt,'sg_date.tif', overwrite=F)

#Repetitive: gdd above base temp for period 1(3 mo prior to start of gs) in year 1972-1982
sg <- stack(c(sg08[3:13])) 
sg_NL <- crop(sg,NL_BD_BUFF_proj)
sgm <- mean(sg_NL)
sg.proj <- projectRaster(sgm, crs = NLproj) 
sg.agrt <- resample(sg.proj, EmtRaster, method='bilinear')
sg_gdd_p1 <- writeRaster(sg.agrt,'sg__gdd_p1.tif', overwrite=F)

#Repetitive: gdd above base temp for period 2(first 6wk after the start of gs) in year 1972-1982
sg <- stack(c(sg09[3:13])) 
sg_NL <- crop(sg,NL_BD_BUFF_proj)
sgm <- mean(sg_NL)
sg.proj <- projectRaster(sgm, crs = NLproj) 
sg.agrt <- resample(sg.proj, EmtRaster, method='bilinear')
sg_gdd_p2 <- writeRaster(sg.agrt,'sg__gdd_p2.tif', overwrite=F)

par(mfrow=c(1,2))
plot(sg_gdd_p1)
plot(sg_gdd_p2)

#### DEM ------------------------
#90m resolution (DEM from SRTM90, directly in R): 
#'alt' stands for altitude (elevation); the data were aggregated from SRTM 90 m resolution data between -60 and 60 latitude
elevation <- getData('alt', country='CAN')
plot(elevation)

#drawExt <- drawExtent()  #use two click on screen to reduce map area to ~NL 
#elevation.NL <- crop(elevation, drawExt)
NL_BD_BUFF_elev <- spTransform(NL_BD_BUFF, projection(elevation)) #use original projection
plot(elevation.NL)

#NLproj <- "+proj=longlat +datum=NAD83 +no_defs"
NLproj <- "+proj=utm +zone=21 +ellps=GRS80 +datum=NAD83 +units=m +no_defs" #EPSG:26921, NAD83 / UTM zone 21N
elevation.proj <- projectRaster(elevation.NL, crs = NLproj)

#crop the raster stack by the boundary
NL_BD_BUFF_21N <- spTransform(NL_BD_BUFF, projection(NLproj))
elevation.crop <- crop(elevation.proj,NL_BD_BUFF_21N)

#change resolution to 2km
elevation.agrt <- resample(elevation.crop, EmtRaster, method='bilinear')


#calculate slope and aspect
slope <- terrain(elevation.agrt, opt='slope')
aspect <- terrain(elevation.agrt, opt='aspect')

elevationRaster <- writeRaster(elevation.agrt,'elevationRaster.tif', overwrite=F)
slopeRaster <- writeRaster(slope,'slopeRaster.tif', overwrite=F)
aspectRaster <- writeRaster(aspect,'aspectRaster.tif', overwrite=F)

par(mfrow=c(1,3))
plot(elevation.agrt)
plot(slope)
plot(aspect)


#### Road ------------------------
library(rgeos)

#https://geo.statcan.gc.ca/nrn_rrn/nl/nrn_rrn_nl_SHAPE.zip
Rd <- "D:/MUN/Data/Road/NRN_NL_7_0_ROADSEG.shp"
NL_Rd <- shapefile(Rd)

NLproj <- "+proj=utm +zone=21 +ellps=GRS80 +datum=NAD83 +units=m +no_defs" 

NL_Rd_sub <- subset(NL_Rd, ROADCLASS==c('Freeway','Collector','Expressway / Highway'), drop = T)
NL_Rd_proj <- spTransform(NL_Rd_sub, projection(NLproj))

NLrd <- crop(NL_Rd_proj, extent(NL_BD_BUFF_21N))

#merge all line features into one feature
NLrd <- gUnion(NLrd, NLrd)

#create a blank raster with the same extent as the NL border layer
RDraster = raster(extent(EmtRaster))
#res(RDraster) <- c(0.01666666, 0.01666666) #change the resolution to needed. 0.016666 is from MaxT01(~2km resolution, 60 arcsec)
res(RDraster) <- c(2000,2000)
#RDraster <- projectRaster(RDraster, crs = NLproj)
projection(RDraster) <- "+proj=utm +zone=21 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"

#compute distances from each feature in NLrd to the raster points (!!!take long time at fine resolution)
dd = gDistance(NLrd, as(RDraster,"SpatialPoints"), byid=TRUE) 
#creates a matrix with a column for each feature in NLrd. To get the nearest distance to any feature, apply min over rows
RDraster[] = apply(dd,1,min)

plot(RDraster)
plot(NLrd, add=T)
plot(NL_BD_21N, add=T)

RoadRaster <- writeRaster(RDraster,'DistanceRoadRaster.tif', overwrite=F)


#### River ------------------------
#use Atlas of Canada source at the scale of 1:5 000 000, from nrcan_rncan/vector/canvec/shp/Hydro (https://open.canada.ca/data/en/dataset/9d96e8c9-22fe-4ad2-b5e8-94a6991b744b)

River <- "D:/MUN/Data/River_canvec/canvec_5M_CA_Hydro_shp/canvec_5M_CA_Hydro/watercourse_1.shp"
NL_River <- shapefile(River)

NL_River_sub <- subset(NL_River, pd_en==c('Newfoundland and Labrador'), drop = T)
NL_River_proj <- spTransform(NL_River_sub, projection(NLproj))

NLriver <- crop(NL_River_proj, extent(NL_BD_BUFF_21N))

#merge all line features into one feature
NLriver <- gUnion(NLriver, NLriver)

#create a blank raster with the same extent as the NL border layer
Riverraster = raster(extent(EmtRaster))
res(Riverraster) <- c(2000,2000)
#RDraster <- projectRaster(RDraster, crs = NLproj)
projection(Riverraster) <- "+proj=utm +zone=21 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"

#compute distances from each feature in NLriver to the raster points (take long time at fine resolution)
dd = gDistance(NLriver, as(Riverraster,"SpatialPoints"), byid=TRUE) 
#creates a matrix with a column for each feature in NLriver. To get the nearest distance to any feature, apply min over rows
Riverraster[] = apply(dd,1,min)

##MinT
#Repetitive: stack MinT in May and Jun (05,06) and year 1972-1982
MinT <- stack(c(MinT05[3:13],MinT06[3:13])) 
MinT_NL <- crop(MinT,NL_BD_BUFF_proj)
MinTm <- mean(MinT_NL)
MinT.proj <- projectRaster(MinTm, crs = NLproj) 
MinT.agrt <- resample(MinT.proj, EmtRaster, method='bilinear')
MinT0506 <- writeRaster(MinT.agrt,'MinT0506.tif', overwrite=F)

#Repetitive: stack MinT in Jun and Jul (06,07) and year 1972-1982
MinT <- stack(c(MinT06[3:13],MinT07[3:13])) 
MinT_NL <- crop(MinT,NL_BD_BUFF_proj)
MinTm <- mean(MinT_NL)
MinT.proj <- projectRaster(MinTm, crs = NLproj) 
MinT.agrt <- resample(MinT.proj, EmtRaster, method='bilinear')
MinT0607 <- writeRaster(MinT.agrt,'MinT0607.tif', overwrite=F)

#Repetitive: stack MinT in Jul and Aug (07,08) and year 1972-1982
MinT <- stack(c(MinT07[3:13],MinT08[3:13])) 
MinT_NL <- crop(MinT,NL_BD_BUFF_proj)
MinTm <- mean(MinT_NL)
MinT.proj <- projectRaster(MinTm, crs = NLproj) 
MinT.agrt <- resample(MinT.proj, EmtRaster, method='bilinear')
MinT0708 <- writeRaster(MinT.agrt,'MinT0708.tif', overwrite=F)


####species layer####
#forest tree species data were obtained from the Department of Fisheries Forestry and Agriculture of Newfoundland and Labrador subjecting to data share agreement.
#forest tree species layers were processed in ArcGIS


#---------------------------------------------------------
####Run MaxEnt ####
SpFineRes <- raster("D:/RinD/rasters/SpRasterNA.tif")
extent((SpFineRes))
#read all predictor variables
path <- file.path("D:/RinD/rasters")
files <- list.files(path, full.names=TRUE )

# stack the raster files together
predictors <- stack(files)
predictors
names(predictors) # names of the rasters
plot(predictors,c(1:18))


####----------add SBW defoliation layer
#defoliation dataset was obtained from Natural Resoures Canada under data sharing agreement
Defo <- shapefile("D:/MUN/Data/SBW_org.shp")
Defo <- spTransform(Defo, projection(NLproj))
Defo@data$yes <- 1

#create an empty raster for rasterization
Defo_emtR <- raster(extent(slopeRaster))

res(Defo_emtR) <- c(2000,2000)
Defo_emtR[] <- 0
projection(Defo_emtR) <- "+proj=utm +zone=21 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"

DefoR <- rasterize(Defo,Defo_emtR, 'yes', fun=min, background = 0)


#read defoliation point layer
De_pt <- shapefile("D:/MUN/Data/Scratch/De_pt.shp")
# first rows
head(De_pt)
# remove duplicate rows
De_pt <- De_pt[!duplicated(De_pt$TARGET_FID),]

####remove presence points in NA tree species raster
De_pt_sel <- crop(De_pt,SpFineRes)
dim(De_pt_sel) # all 3090 points are within FineRes Species layer

Dept_intm <- which(extract(!is.na(SpFineRes), De_pt) == 1) # ==1 are any De_pt points within the SpFineRes layer
De_pt_sel <- De_pt[Dept_intm,] 
dim(De_pt_sel) #3090-126=2964


# extract the environmental values at each
# of the occurrence locations.
presvals <- extract(predictors, De_pt)

## ----  Background data for train/test set ------------------------------------
set.seed(10)
backg <- randomPoints(SpFineRes, n=14820, extf = 1.25)# create 5x presence random background points in the species raster
#backg <- backg[1:10000] # because of many NAs in species layer (only 8492 valued cells in species), randomPoints generate less then request, use subset to get requested number
colnames(backg) = c('lon', 'lat')
group <- kfold(backg, 5)
backg_train <- backg[group != 1, ]
backg_test <- backg[group == 1, ]

##check Multicollinearity
sdmdata <- predictors
# correlation structure
library(corrplot)
library(PerformanceAnalytics)
varCor <- cor(na.omit(sdmdata[,-c(1,3,4,14,15,23)]))
corrplot(varCor)
chart.Correlation(varCor, histogram=F, pch=19)


# clustering with dendrogram
allDistNew <- abs(as.dist(cor(sdmdata[,-c(1,3,4,14,15,23)])))
allClusNew <- hclust(1 - allDistNew)
plot(allClusNew, hang=-1)


# Variance Inflation Factor
library(usdm)
vif(sdmdata[,-c(1,3,4,14,15,23)])
vifstep(sdmdata[,-c(1,3,4,14,15,23)])
vifcor(sdmdata[,-c(1,3,4,14,15,23)], th=0.8) 

#remove the 5 predictors after VIF
predictors_sel <- dropLayer(predictors,c(7,8,9,15,16))

####evalute model with ENMeval package
#need presence and background in dataframe with just two coord columes
pres_df <- De_pt_sel@coords
colnames(pres_df) <- c("lon","lat")

predictors_sel$SpRasterNA <- raster::as.factor(predictors_sel$SpRasterNA)

e.mx <- ENMevaluate(occs = pres_df, envs = predictors_sel, bg = backg,
                    algorithm = 'maxent.jar', partitions = 'block',
                    tune.args = list(fc=c("L","LH","LQ","LQH","LQP","LQHP","LQHPT"), rm = seq(0.5,5,by=1)),
                    parallel = TRUE#make it faster
                    #categoricals = predictors_sel$SpRasterNA,
)

write.csv(e.mx@results,"D:/RinD/MaxEnt_evalution.csv",row.names = T)



# We first run the null simulations with 100 iterations to get a reasonable null distribution 
# for comparisons with the empirical values
mod.null <- ENMnulls(e.mx, mod.settings = list(fc = "LQ", rm = 4.5), no.iter = 100)
null.results(mod.null) %>% head()
null.emp.results(mod.null)

#run optimal model
filePath <- "D:/RinD/MaxEnt_out/MaxEnt_out22"
mx.LQHPT35 <- maxent(predictors_sel, # env data as a raster stack
                     De_pt_sel, # presence data
                     a=backg,#if use pre-selected bg
                     #nbg=1847, #number of background points
                     factors='SpRasterNA',
                     args=c("responsecurves",
                            #"-q", # turn off quadratic features
                            #"-p", # turn off product features
                            #"nothreshold", # turn off threshold)
                            "betamultiplier=3.5",
                            c("-m", 10000)), # max iterations = 10K
                     path=filePath) # where to save all the output



mx <- mx_LQHPT35
plot(mx)

#get data from each response curve
rs <- response(mx.LQHPT35.jn, var= "DistanceRiverRaster")
rs
plot(rs,type="l",ylim=c(0,1),col=2)



#plot bF and bS
par(mfrow=c(1,2))
r <- predictors[["SpRasterNA"]]
plot(r==1)
points(De_pt_sel, pch=21, col=rgb(0,0,0,0.5), cex=0.2)
plot(r==2)
points(De_pt, pch=21, col=rgb(0,0,0,0.5), cex=0.2)



# evaluate the model using the test data
e <- evaluate(pres_test, backg_test, mx, predictors_sel)
e

# predict back to geography
mxPred <- predict(mx, predictors_sel)
plot(mxPred, col=rgb.tables(1000))
# Predict in 'raw' format and save raster of prediction
filePath <- "D:/RinD"
mxPred <- predict(mx, predictors_sel, args=c("outputformat=raw"),
                  filename=paste0(filePath, '/maxent_predictionRAW.tif'), overwrite=TRUE)
mxPred <- predict(mx, predictors_sel, args=c("outputformat=cloglog"),
                  filename=paste0(filePath, '/maxent_predictioncloglog.tif'),overwrite=TRUE)

plot(mxPred)

#check probability in bS and bF
#call in cloglog tif
predRast <- raster(paste0(filePath, '/maxent_predictioncloglog.tif'))

bFR <- SpRasterNA==1
bSR <- SpRasterNA==2

bFpred1 <- calc(bFR, fun=function(x){x[x!=1] <- NA; return(x)})#set value 0 to NA
bFpred <- mask(predRast,bFpred1) # use bF raster to mask prediction layer

bSpred1 <- calc(bSR, fun=function(x){x[x!=1] <- NA; return(x)})#set value 0 to NA
bSpred <- mask(predRast,bSpred1) # use bS raster to mask prediction layer


bF20 <- calc(bFpred, fun=function(x){x[x<=0.2] <- 1; return(x)}) #s1
bF20 <- calc(bFpred, fun=function(x){x[x>0.6 & x<=0.8] <- 1; return(x)})#s2
bF20 <- calc(bFpred, fun=function(x){x[x>0.8 & x!=0] <- 1; return(x)}) #s3
bF20 <- calc(bF20, fun=function(x){x[x!=1] <- NA; return(x)}) #run after one from s1,2,or 3
cellStats(bF20, 'sum') #0-20:6377; 20-40:1486; 40-60:1040; 60-80:1155 80-100:427;

bS20 <- calc(bSpred, fun=function(x){x[x<=0.2 & x!=0] <- 1; return(x)})
bS20 <- calc(bSpred, fun=function(x){x[x>0.6 & x<=0.8 & x!=0] <- 1; return(x)})
bS20 <- calc(bSpred, fun=function(x){x[x>0.8 & x!=0] <- 1; return(x)})
bS20 <- calc(bS20, fun=function(x){x[x!=1] <- 0; return(x)})
cellStats(bS20, 'sum') #0-20:2390; 20-40:1541; 40-60:1474; 60-69:1676 40-60:644;

#map bF bS prediction map
library(RColorBrewer)
display.brewer.all()

cuts=c(0,0.2,0.4,0.6,0.8,1) #set breaks
pal <- colorRampPalette(c("gray","red"))

par(mfrow=c(1,2), mar=c(3,5,2,5))
plot(bFpred, breaks=cuts, col = pal(6)) #plot with defined breaks
plot(NL_BD_21N, border="grey", lwd=0.01, add=T)
plot(bSpred, breaks=cuts, col = pal(6)) 
plot(NL_BD_21N, border="grey", lwd=0.01, add=T)