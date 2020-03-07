############################
#Load libraries

library("png")
library(maptools)
library(rgeos)
library(tidyverse)
library(rgdal)
library(ggthemes)
library("spatstat")
library(beepr)
library("png")
library("spatstat")
library(parallel)

#https://www.r-bloggers.com/how-to-go-parallel-in-r-basics-tips/
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)
#stopCluster(cl)

#Read PNG, conversion to ppp (point pattern) and calculate Ripley's function (Kest)
generate_ppp <- function(X) {   
  fileNames <- Sys.glob("*.png")
  data <- list()
  for (i in 1:length(fileNames)){
    d <- readPNG(fileNames[i])
    coord <- which(d==1, arr.ind=TRUE)
    colnames(coord) <- c("x","y","int")
    coords <- data.frame(coord)
    x <- coords$x
    y <- coords$y
    d <- ppp(x, y, c(0,1000), c(0,3300))
    basename <-gsub("*.png","", fileNames[i])
    data[[i]] <- d
  }
  return(data)
}  

#Read sample names
root <- "introduce_root_directory"
setwd(root)
#Obtain sample names using folder names
folders <- list.dirs(recursive=FALSE)
samples <- c()
for (i in 1:length(folders)){
  name <- gsub("./","",folders[i])
  samples[i] <- name
}

#PROCESSING LOOP¡
samples_list <- list()
for (i in 1:length(samples)){
  wd <- paste(sep="/", root,samples[i])
  setwd(wd)
  sample <- samples[i]
  data <-generate_ppp()
  samples_list[[i]] <- data 
}

kest_list <- list()
iso <- list()
suma <- list()
porcentajes <- list()
porcentajes_ <- rep(list(list()), length(samples_list))
ratio <- rep(list(list()), length(samples_list))
points <- rep(list(list()), length(samples_list))
total <- list()
data <- NULL
ratio_ <- list()
average <- 0
r <- seq(0,100,by=0.2)
for (i in 1:length(samples_list)){
  kest_list[[i]] <- parLapply(cl, samples_list[[i]], Kest, nlarge=300000, r=r)
  iso[[i]] <- parLapply(cl, kest_list[[i]], function(x) x$iso)
  suma[[i]] <- sum(parSapply(cl, iso[[i]], function(x) {sum(unlist(x))}))
  data <- suma[[i]]
  porcentajes[[i]] <- sapply(iso[[i]], function(x) {x / data})
  
  for(j in 1:length(iso[[i]])){
    porcentajes_[[i]][[j]] <- (porcentajes[[i]])[,j]
  }
  
  points[[i]] <- parLapply(cl, samples_list[[i]], function(x) {x$n})
  total[[i]] <- sum(unlist(points[[i]]))
  
  for(j in 1:length(iso[[i]])){
    ratio[[i]][[j]] <- iso[[i]][[j]] * points[[i]][[j]] / total[[i]]
  }
  data <- 0
  for(j in 1:length(iso[[i]])){
    data <- ratio[[i]][[j]] + data
  }
  ratio_[[i]] <- data
  i
}

#SORT SAMPLES BY GROUP

WT <- grep("WT", samples, value=FALSE)
MUT1 <- grep("MUT1", samples, value=FALSE)
MUT2 <- grep("MUT2", samples, value=FALSE)

#CALCULATE MEAN AND SEM

WT_mean <- NULL
for(i in 1:length(ratio_[WT])){
  for (j in 1:501){
    WT_mean[[j]] <- mean(sapply(ratio_[WT], function(x) (x)[[j]]))
  }} 
MUT1_mean <- NULL
for(i in 1:length(ratio_[MUT1])){
  for (j in 1:501){
    MUT1_mean[[j]] <- mean(sapply(ratio_[MUT1], function(x) (x)[[j]]))
  }} 
MUT2_mean <- NULL
for(i in 1:length(ratio_[MUT2])){
  for (j in 1:501){
    MUT2_mean[[j]] <- mean(sapply(ratio_[MUT2], function(x) (x)[[j]]))
  }} 


sem <- function(x) {sd(x)/sqrt(length(x))}

WT_sem <- NULL
for(i in 1:length(ratio_[WT])){
  for (j in 1:501){
    WT_sem[[j]] <- sem(sapply(ratio_[WT], function(x) (x)[[j]]))
  }}  
MUT1_sem <- NULL
for(i in 1:length(ratio_[MUT1])){
  for (j in 1:501){
    MUT1_sem[[j]] <- sem(sapply(ratio_[MUT1], function(x) (x)[[j]]))
  }} 
MUT2_sem <- NULL
for(i in 1:length(ratio_[MUT2])){
  for (j in 1:501){
    MUT2_sem[[j]] <- sem(sapply(ratio_[MUT2], function(x) (x)[[j]]))
  }} 

#CUMMULATIVE PLOT WITH SEM ENVELOPE

x.index <- c(seq(from=1, to=501, by=1))

plot(WT_mean,type="n")
polygon(c(x.index,rev(x.index)),c((WT_mean + WT_sem),rev((WT_mean - WT_sem))),border=NA,col=blues9[3])
polygon(c(x.index,rev(x.index)),c((MUT1_mean + MUT1_sem),rev((MUT1_mean - MUT1_sem))),border=NA,col=blues9[3])
polygon(c(x.index,rev(x.index)),c((MUT2_mean + MUT2_sem),rev((MUT2_mean - MUT2_sem))),border=NA,col=blues9[3])
lines(WT_mean, col="black")
lines(MUT1_mean, col="red")
lines(MUT2_mean, col="blue")

