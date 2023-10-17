library(raster)
library(usdm)
library(biomod2)

Poly<-shapefile("F:\\Rusts data\\climate\\Iran.shp")
plot(Poly)
bio3<-raster("F:\\Rusts data\\climate\\Current_clim\\bio3.asc")
bio4<-raster("F:\\Rusts data\\climate\\Current_clim\\bio4.asc")
bio10<-raster("F:\\Rusts data\\climate\\Current_clim\\bio10.asc")
bio12<-raster("F:\\Rusts data\\climate\\Current_clim\\bio12.asc")
bio13<-raster("F:\\Rusts data\\climate\\Current_clim\\bio13.asc")
bio14<-raster("F:\\Rusts data\\climate\\Current_clim\\bio15.asc")
bio15<-raster("F:\\Rusts data\\climate\\Current_clim\\bio14.asc")
bio16<-raster("F:\\Rusts data\\climate\\Current_clim\\bio16.asc")
bio18<-raster("F:\\Rusts data\\climate\\Current_clim\\bio18.asc")
Prec1<-raster("F:\\Rusts data\\climate\\Current_clim\\Prec1.asc")
Prec2<-raster("F:\\Rusts data\\climate\\Current_clim\\Prec2.asc")
Prec3<-raster("F:\\Rusts data\\climate\\Current_clim\\Prec3.asc")
Prec4<-raster("F:\\Rusts data\\climate\\Current_clim\\Prec4.asc")
Prec5<-raster("F:\\Rusts data\\climate\\Current_clim\\Prec5.asc")
Prec11<-raster("F:\\Rusts data\\climate\\Current_clim\\Prec11.asc")
Prec12<-raster("F:\\Rusts data\\climate\\Current_clim\\Prec12.asc")
Tmax3<-raster("F:\\Rusts data\\climate\\Current_clim\\Tmax3.asc")
Tmax4<-raster("F:\\Rusts data\\climate\\Current_clim\\Tmax4.asc")
Tmax5<-raster("F:\\Rusts data\\climate\\Current_clim\\Tmax5.asc")
Tmax6<-raster("F:\\Rusts data\\climate\\Current_clim\\Tmax6.asc")
Tmax9<-raster("F:\\Rusts data\\climate\\Current_clim\\Tmax9.asc")
Tmin3<-raster("F:\\Rusts data\\climate\\Current_clim\\Tmin3.asc")
Tmin5<-raster("F:\\Rusts data\\climate\\Current_clim\\Tmin5.asc")
Tmin6<-raster("F:\\Rusts data\\climate\\Current_clim\\Tmin6.asc")
Tmin7<-raster("F:\\Rusts data\\climate\\Current_clim\\Tmin7.asc")

setwd("C:\\Users\\SONY\\Desktop\\Ayoub Paper")
Wheat<-raster("Wheat_cover.tif")
Wheat<-resample(Wheat, Tmin7)
extent(Wheat)<-extent(Tmin7)
Wheat <- raster(vals=values(Wheat),ext=extent(Tmin7),crs=crs(Tmin7),
              nrows=dim(Tmin7)[1],ncols=dim(Tmin7)[2])

bio.current<-stack(bio3,bio4,bio10,bio12,bio13,
                   bio14,bio15,bio16,bio18,Prec1,
                   Prec2,Prec3,Prec4,Prec5,Prec11,Prec12, 
                   Tmax3,Tmax4,Tmax5,Tmax6,Tmax9,Tmin3,
                   Tmin5,Tmin6,Tmin7, Wheat)




yellow_rust<-read.csv("F:\\Rusts data\\yellow_rust.csv")

Dataset<-yellow_rust
names(Dataset)
dim(Dataset)
yellow_rust<-Dataset [1:254,2:3]
m2 <- cbind(1,1:254)
colnames(m2, do.NULL = FALSE)
colnames(m2) <- c("yellow","ID")
colnames(m2) <- colnames(m2, prefix = "Sub_")
colnames(m2)
SP<-data.frame(m2)
MergedDataset <- merge(SP, yellow_rust, all=TRUE, by="row.names")
rownames(MergedDataset) <- MergedDataset$Row.names
MergedDataset$Row.names <- NULL
names(MergedDataset)
LatLong <- MergedDataset [,3:4] # coordinates of points
Expl.Var <- bio.current
Resp.Var <- MergedDataset$yellow# species occurences

myBiomodData<-BIOMOD_FormatingData(resp.var=MergedDataset$yellow,
                                   expl.var=bio.current,
                                   resp.xy = LatLong,
                                   resp.name = "yellow",
                                   eval.resp.var = NULL,
                                   eval.expl.var = NULL,
                                   eval.resp.xy = NULL,
                                   PA.nb.rep = 1,
                                   PA.nb.absences = 1000,
                                   na.rm = TRUE)



myBiomodOptions <- BIOMOD_ModelingOptions()

myBiomodModelOut <- BIOMOD_Modeling(
  myBiomodData,
  models = c('GBM','ANN','RF','GLM'),
  models.options =myBiomodOptions,
  NbRunEval=1,
  DataSplit=80,
  Yweights=NULL,
  VarImport=1,
  models.eval.meth = c('TSS','ROC'),
  SaveObj = TRUE,
  rescal.all.models = TRUE,
  do.full.models = FALSE,
  modeling.id = paste("yellow","FirstModeling",sep=""))## print model evaluation


MODELSS2<-(myBiomodModelEval <- get_evaluations(myBiomodModelOut))

MODELSS2
resultsmodel<-as.data.frame(MODELSS2)

write.csv(resultsmodel, "resultsmodel.csv")


variables_importance<-get_variables_importance(myBiomodModelOut)


write.csv(variables_importance, "variables_importance.csv")

models_scores_graph(myBiomodModelOut, by = "models" ,
                    metrics = c("ROC","TSS"))


myRFs <- BIOMOD_LoadModels(myBiomodModelOut, models='RF')

myRespPlot2D <- response.plot2(models  = myRFs ,
                               Data = get_formal_data(myBiomodModelOut,'expl.var'), 
                               show.variables= get_formal_data(myBiomodModelOut,'expl.var.names'),
                               do.bivariate = FALSE,
                               fixed.var.metric = 'median',
                               col = c("blue"),
                               save.file="yes",
                               legend = FALSE,
                               ImageSize=1000,
                               data_species = get_formal_data(myBiomodModelOut,'resp.var'))



myBiomodProjection_RF <- BIOMOD_Projection(
  modeling.output = myBiomodModelOut,
  new.env = bio.current,
  proj.name = 'current_RF',
  selected.models = grep('_RF', get_built_models(
    myBiomodModelOut), value=TRUE),
  compress = FALSE,
  build.clamping.mask = FALSE)

plot(myBiomodProjection_RF)



myBiomodProjection_GBM <- BIOMOD_Projection(
  modeling.output = myBiomodModelOut,
  new.env = bio.current,
  proj.name = 'current_GBM',
  selected.models = grep('_GBM', get_built_models(
    myBiomodModelOut), value=TRUE),
  compress = FALSE,
  build.clamping.mask = FALSE)

plot(myBiomodProjection_GBM)


myBiomodProjection_GLM <- BIOMOD_Projection(
  modeling.output = myBiomodModelOut,
  new.env = bio.current,
  proj.name = 'current_GLM',
  selected.models = grep('_GLM', get_built_models(
    myBiomodModelOut), value=TRUE),
  compress = FALSE,
  build.clamping.mask = FALSE)

plot(myBiomodProjection_GLM)

myBiomodProjection_ANN <- BIOMOD_Projection(
  modeling.output = myBiomodModelOut,
  new.env = bio.current,
  proj.name = 'current_ANN',
  selected.models = grep('_ANN', get_built_models(
    myBiomodModelOut), value=TRUE),
  compress = FALSE,
  build.clamping.mask = FALSE)

plot(myBiomodProjection_ANN)

########################################################################

####ensemble 

myBiomodEM <- BIOMOD_EnsembleModeling(
  modeling.output = myBiomodModelOut,
  chosen.models = 'all',
  eval.metric = c('TSS'),
  eval.metric.quality.threshold = c(0.55),
  prob.mean = T,
  prob.cv = T,
  prob.ci = T,
  prob.ci.alpha = 0.05,
  prob.median = T,
  committee.averaging = T,
  prob.mean.weight = T,
  prob.mean.weight.decay = 'proportional' )


myBiomodProj <- BIOMOD_Projection(
  modeling.output = myBiomodModelOut,
  new.env = bio.current,
  proj.name = 'All',
  selected.models = 'all',
  binary.meth = NULL,
  compress = 'xz',
  clamping.mask = F,
  output.format = '.grd')




myBiomodEF <- BIOMOD_EnsembleForecasting(
  projection.output = myBiomodProj ,
  EM.output = myBiomodEM )


plot(myBiomodEF)

proj_All_yellow_ensemble
Ens<-raster("C:\\Users\\SONY\\Documents\\yellow\\proj_All\\proj_All_yellow_ensemble.grd")

plot(Ens)
writeRaster(Ens, filename="Ens.asc", format="ascii", overwrite=TRUE)   
Ensemble<-raster("C:\\Rusts data\\Ens.asc")
plot(Ensemble)

##############################################################################
##########ssp126########################

bio3<-raster("F:\\Rusts data\\climate\\RCP2.6_2050\\bio3.asc")
bio4<-raster("F:\\Rusts data\\climate\\RCP2.6_2050\\bio4.asc")
bio10<-raster("F:\\Rusts data\\climate\\RCP2.6_2050\\bio10.asc")
bio12<-raster("F:\\Rusts data\\climate\\RCP2.6_2050\\bio12.asc")
bio13<-raster("F:\\Rusts data\\climate\\RCP2.6_2050\\bio13.asc")
bio14<-raster("F:\\Rusts data\\climate\\RCP2.6_2050\\bio15.asc")
bio15<-raster("F:\\Rusts data\\climate\\RCP2.6_2050\\bio14.asc")
bio16<-raster("F:\\Rusts data\\climate\\RCP2.6_2050\\bio16.asc")
bio18<-raster("F:\\Rusts data\\climate\\RCP2.6_2050\\bio18.asc")
Prec1<-raster("F:\\Rusts data\\climate\\RCP2.6_2050\\Prec1.asc")
Prec2<-raster("F:\\Rusts data\\climate\\RCP2.6_2050\\Prec2.asc")
Prec3<-raster("F:\\Rusts data\\climate\\RCP2.6_2050\\Prec3.asc")
Prec4<-raster("F:\\Rusts data\\climate\\RCP2.6_2050\\Prec4.asc")
Prec5<-raster("F:\\Rusts data\\climate\\RCP2.6_2050\\Prec5.asc")
Prec11<-raster("F:\\Rusts data\\climate\\RCP2.6_2050\\Prec11.asc")
Prec12<-raster("F:\\Rusts data\\climate\\RCP2.6_2050\\Prec12.asc")
Tmax3<-raster("F:\\Rusts data\\climate\\RCP2.6_2050\\Tmax3.asc")
Tmax4<-raster("F:\\Rusts data\\climate\\RCP2.6_2050\\Tmax4.asc")
Tmax5<-raster("F:\\Rusts data\\climate\\RCP2.6_2050\\Tmax5.asc")
Tmax6<-raster("F:\\Rusts data\\climate\\RCP2.6_2050\\Tmax6.asc")
Tmax9<-raster("F:\\Rusts data\\climate\\RCP2.6_2050\\Tmax9.asc")
Tmin3<-raster("F:\\Rusts data\\climate\\RCP2.6_2050\\Tmin3.asc")
Tmin5<-raster("F:\\Rusts data\\climate\\RCP2.6_2050\\Tmin5.asc")
Tmin6<-raster("F:\\Rusts data\\climate\\RCP2.6_2050\\Tmin6.asc")
Tmin7<-raster("F:\\Rusts data\\climate\\RCP2.6_2050\\Tmin7.asc")


bio.RCP2.6_2040<-stack(bio3,bio4,bio10,bio12,bio13,
                       bio14,bio15,bio16,bio18,Prec1,
                       Prec2,Prec3,Prec4,Prec5,Prec11,Prec12, 
                       Tmax3,Tmax4,Tmax5,Tmax6,Tmax9,Tmin3,
                       Tmin5,Tmin6,Tmin7, Wheat)


myBiomodProj <- BIOMOD_Projection(
  modeling.output = myBiomodModelOut,
  new.env = bio.RCP2.6_2040,
  proj.name = 'RCP2.6_2050',
  selected.models = 'all',
  binary.meth = NULL,
  compress = 'xz',
  clamping.mask = F,
  output.format = '.grd')


myBiomodEF <- BIOMOD_EnsembleForecasting(
  projection.output = myBiomodProj,
  EM.output = myBiomodEM)



###########################################################################

bio3<-raster("F:\\Rusts data\\climate\\RCP8.5_2050\\bio3.asc")
bio4<-raster("F:\\Rusts data\\climate\\RCP8.5_2050\\bio4.asc")
bio10<-raster("F:\\Rusts data\\climate\\RCP8.5_2050\\bio10.asc")
bio12<-raster("F:\\Rusts data\\climate\\RCP8.5_2050\\bio12.asc")
bio13<-raster("F:\\Rusts data\\climate\\RCP8.5_2050\\bio13.asc")
bio14<-raster("F:\\Rusts data\\climate\\RCP8.5_2050\\bio15.asc")
bio15<-raster("F:\\Rusts data\\climate\\RCP8.5_2050\\bio14.asc")
bio16<-raster("F:\\Rusts data\\climate\\RCP8.5_2050\\bio16.asc")
bio18<-raster("F:\\Rusts data\\climate\\RCP8.5_2050\\bio18.asc")
Prec1<-raster("F:\\Rusts data\\climate\\RCP8.5_2050\\Prec1.asc")
Prec2<-raster("F:\\Rusts data\\climate\\RCP8.5_2050\\Prec2.asc")
Prec3<-raster("F:\\Rusts data\\climate\\RCP8.5_2050\\Prec3.asc")
Prec4<-raster("F:\\Rusts data\\climate\\RCP8.5_2050\\Prec4.asc")
Prec5<-raster("F:\\Rusts data\\climate\\RCP8.5_2050\\Prec5.asc")
Prec11<-raster("F:\\Rusts data\\climate\\RCP8.5_2050\\Prec11.asc")
Prec12<-raster("F:\\Rusts data\\climate\\RCP8.5_2050\\Prec12.asc")
Tmax3<-raster("F:\\Rusts data\\climate\\RCP8.5_2050\\Tmax3.asc")
Tmax4<-raster("F:\\Rusts data\\climate\\RCP8.5_2050\\Tmax4.asc")
Tmax5<-raster("F:\\Rusts data\\climate\\RCP8.5_2050\\Tmax5.asc")
Tmax6<-raster("F:\\Rusts data\\climate\\RCP8.5_2050\\Tmax6.asc")
Tmax9<-raster("F:\\Rusts data\\climate\\RCP8.5_2050\\Tmax9.asc")
Tmin3<-raster("F:\\Rusts data\\climate\\RCP8.5_2050\\Tmin3.asc")
Tmin5<-raster("F:\\Rusts data\\climate\\RCP8.5_2050\\Tmin5.asc")
Tmin6<-raster("F:\\Rusts data\\climate\\RCP8.5_2050\\Tmin6.asc")
Tmin7<-raster("F:\\Rusts data\\climate\\RCP8.5_2050\\Tmin7.asc")

bio.RCP8.5_2050<-stack(bio3,bio4,bio10,bio12,bio13,
                       bio14,bio15,bio16,bio18,Prec1,
                       Prec2,Prec3,Prec4,Prec5,Prec11,Prec12, 
                       Tmax3,Tmax4,Tmax5,Tmax6,Tmax9,Tmin3,
                       Tmin5,Tmin6,Tmin7, Wheat)


myBiomodProj <- BIOMOD_Projection(
  modeling.output = myBiomodModelOut,
  new.env = bio.RCP8.5_2050,
  proj.name = 'RCP8.5_2050',
  selected.models = 'all',
  binary.meth = NULL,
  compress = 'xz',
  clamping.mask = F,
  output.format = '.grd')


myBiomodEF <- BIOMOD_EnsembleForecasting(
  projection.output = myBiomodProj,
  EM.output = myBiomodEM)


#########################################################################
##########################################################################
###########################################################################

RF_map<- raster("C:\\Users\\SONY\\Desktop\\Ayoub Paper\\yellow\\proj_current_RF\\proj_current_RF_yellow.grd")
GLM_map<- raster("C:\\Users\\SONY\\Desktop\\Ayoub Paper\\yellow\\proj_current_GLM\\proj_current_GLM_yellow.grd")
GBM_map<- raster("C:\\Users\\SONY\\Desktop\\Ayoub Paper\\yellow\\proj_current_GBM\\proj_current_GBM_yellow.grd")
ANN_map<- raster("C:\\Users\\SONY\\Desktop\\Ayoub Paper\\yellow\\proj_current_ANN\\proj_current_ANN_yellow.grd")
Ensemble<- raster("C:\\Users\\SONY\\Desktop\\Ayoub Paper\\yellow\\proj_All\\proj_All_yellow_ensemble.grd")

plot(RF_map)
plot(GLM_map)
plot(GBM_map)
plot(ANN_map)
plot(Ensemble)

setwd("C:\\Users\\SONY\\Desktop\\Dr. Bakhshi\\Result_Bmediterrana")
writeRaster(GLM_map, "GLM_map.asc")
writeRaster(MARS_map, "MARS_map.asc")
writeRaster(RF_map, "RF_map.asc")
writeRaster(GBM_map, "GBM_map.asc")
writeRaster(Ensemble_map, "Ensemble_map.asc")

library(tidyr)
library(ggplot2)
rasdf <- as.data.frame(RF_map,xy=TRUE)%>%drop_na()
head(rasdf)

RF<-ggplot()+
  geom_raster(aes(x=x,y=y,fill=yellow_PA1_RUN1_RF),data=rasdf)+
  scale_fill_viridis_c('Suitablity',direction = -1)+
  coord_sf(expand=c(0,0))+
  labs(x='Longitude',y='Latitude',
       title="RF")+
  theme_void()

RF

rasdf <- as.data.frame(GLM_map,xy=TRUE)%>%drop_na()
head(rasdf)

GLM<-ggplot()+
  geom_raster(aes(x=x,y=y,fill=yellow_PA1_RUN1_GLM),data=rasdf)+
  scale_fill_viridis_c('Suitablity',direction = -1)+
  coord_sf(expand=c(0,0))+
  labs(x='Longitude',y='Latitude',
       title="GLM")+
  theme_void()

GLM


rasdf <- as.data.frame(GBM_map,xy=TRUE)%>%drop_na()
head(rasdf)

GBM<-ggplot()+
  geom_raster(aes(x=x,y=y,fill=yellow_PA1_RUN1_GBM),data=rasdf)+
  scale_fill_viridis_c('Suitablity',direction = -1)+
  coord_sf(expand=c(0,0))+
  labs(x='Longitude',y='Latitude',
       title="GBM")+
  theme_void()

GBM


rasdf <- as.data.frame(ANN_map,xy=TRUE)%>%drop_na()
head(rasdf)

ANN<-ggplot()+
  geom_raster(aes(x=x,y=y,fill=yellow_PA1_RUN1_ANN),data=rasdf)+
  scale_fill_viridis_c('Suitablity',direction = -1)+
  coord_sf(expand=c(0,0))+
  labs(x='Longitude',y='Latitude',
       title="ANN")+
  theme_void()

ANN


rasdf <- as.data.frame(Ensemble,xy=TRUE)%>%drop_na()
head(rasdf)

Ensemble<-ggplot()+
  geom_raster(aes(x=x,y=y,fill=yellow_EMmeanByTSS_mergedAlgo_RUN1_PA1),data=rasdf)+
  scale_fill_viridis_c('Suitablity',direction = -1)+
  coord_sf(expand=c(0,0))+
  labs(x='Longitude',y='Latitude',
       title="Current")+
  theme_void()

Ensemble

library(ggpubr)

Yrust <- ggarrange(RF, GBM, GLM,ANN,
                        labels = c(""),
                        ncol = 2, nrow = 2)

Yrust

ggsave(file="Yrust.png", dpi=1000)

ggsave(file="Ensemble.png", dpi=1000)


RCP2.6_2050<- raster("C:\\Users\\SONY\\Desktop\\Ayoub Paper\\yellow\\proj_RCP2.6_2050\\proj_RCP2.6_2050_yellow_ensemble.grd")

rasdf <- as.data.frame(RCP2.6_2050,xy=TRUE)%>%drop_na()
head(rasdf)

RCP2.6_2050<-ggplot()+
  geom_raster(aes(x=x,y=y,fill=yellow_EMmeanByTSS_mergedAlgo_RUN1_PA1),data=rasdf)+
  scale_fill_viridis_c('Suitablity',direction = -1)+
  coord_sf(expand=c(0,0))+
  labs(x='Longitude',y='Latitude',
       title="RCP2.6_2050")+
  theme_void()

RCP2.6_2050



ggsave(file="RCP2.6_2050.png", dpi=1000)



RCP8.5_2050<- raster("C:\\Users\\SONY\\Desktop\\Ayoub Paper\\yellow\\proj_RCP8.5_2050\\proj_RCP8.5_2050_yellow_ensemble.grd")

rasdf <- as.data.frame(RCP8.5_2050,xy=TRUE)%>%drop_na()
head(rasdf)

RCP8.5_2050<-ggplot()+
  geom_raster(aes(x=x,y=y,fill=yellow_EMmeanByTSS_mergedAlgo_RUN1_PA1),data=rasdf)+
  scale_fill_viridis_c('Suitablity',direction = -1)+
  coord_sf(expand=c(0,0))+
  labs(x='Longitude',y='Latitude',
       title="RCP8.5_2050")+
  theme_void()

RCP8.5_2050



ggsave(file="RCP8.5_2050.png", dpi=1000)



Future <- ggarrange(RCP2.6_2050,RCP8.5_2050,
                   labels = c(""),
                   ncol = 2, nrow = 1)

Future




Variableimprtance<-ggplot(Importance, aes(x=reorder(Variables, Importance), y=Importance, color=as.factor(Type))) + 
  geom_point() +
  geom_segment(aes(x=Variables,xend=Variables,y=0,yend=Importance)) +
  scale_color_discrete(name="Type") +
  ylab("Relative importance") +
  xlab("Variables") +
  coord_flip()

Variableimprtance
ggsave(file="Variableimprtance.png", dpi=1000)
