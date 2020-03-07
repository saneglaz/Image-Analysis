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


##############
# FUNCTIONS
#############

rescale <- function(x, x0, xm, n) {
  (x - x0)/(xm - x0)*n
}

count_grid <- function(x,y) {
  # ARRAY
  n <-100
  output <- matrix(NA, 100, 100)
  
  xi <- rescale(x, 0, 800, n)
  yi <- rescale(y, 0, 800, n)
  
  plot(xi, yi, xlim=c(0, n), ylim=c(0, n), col="grey", pch=16)
  points(xi, yi)
  abline(h=seq(0, n, 1), col="grey")
  abline(v=seq(0, n, 1), col="grey")
  
  # now we can do something similar, but this time taking the index part as the exact bin
  xii <- ceiling(rescale(x, 0, 679, n))
  yii <- ceiling(rescale(y, 0, 830, n))
  points(xii-0.5, yii-0.5, col="red", pch=16, cex=0.5)
  
  # finally, use table to count the number in occupied cells
  ret <- table(xii, yii)
  N <- 102 - ncol(ret)
  for (i in 1:N){
    ret <- cbind(ret,rep(0, length(nrow(ret))))
  }
  N <- 102 - nrow(ret)
  for (i in 1:N){
    ret <- rbind(ret,rep(0, length(ncol(ret))))
    }
  return(ret)
}

#Export point coordinates from binary mask images
read_coords <- function(X) {
  fileNames <- Sys.glob(X)
  data <- NULL
  for (i in 1:length(fileNames)){
    d <- readPNG(fileNames[i])
    coord <- which(d==1, arr.ind=TRUE)
    colnames(coord) <- c("x","y")
    coords <- data.frame(coord)
    data <- rbind(data,coords)
  }
  return(data)
}

#Quantify signal
count_signal <- function(X,Y) {
  setwd(X)
  fileNames <- Sys.glob(Y)
  sum <- 0
  data <- NULL
  for (i in 1:length(fileNames)){
    input <- paste(sep="/", wd_side, fileNames[i])
    d <- readPNG(fileNames[i])
    sum <- sum + nrow(which(d==1, arr.ind=TRUE))
    data <- sum
  }
  return(data)
}


wd <- paste(sep="", root, samples[i], "/IPSI_clean/")
setwd(wd)
carpetas <- list.dirs(recursive=FALSE)
sections <- c()
for (i in 1:length(carpetas)){
  name <- gsub("./","",carpetas[i])
  sections[i] <- name
}

#################
#PROCESSING LOOP:
################
#browse through folders
data <- NULL
for (i in 1:length(sections)){
  wd2 <- paste(sep="", wd,sections[i])
  sides <- c("green_side", "red_side")
  greenside <- list()
  redside <- list()
  output <- list(greenside, redside)
  #browse through both sides and read point coordinates
  for(j in 1:length(sides)){
    wd_side <- paste(sep="/", wd2,sides[j])
    setwd(wd_side)
    channels <- c("488*.png", "647*.png")
    coords_488 <- count_signal(wd_side, channels[1])
    coords_647 <- count_signal(wd_side, channels[2])
    output[[j]][[1]] <- coords_488
    output[[j]][[2]] <- coords_647
  }
  data[[i]] <- output
}

total <- list(list(0,0),list(0,0))
for (i in 1:length(sections)){
  for (j in 1:2){
    for (k in 1:2){
      total[[j]][[k]] <- data[[i]][[j]][[k]] + total[[j]][[k]]
    }
  }
}

no_cross <- NULL
cross <- NULL
no_cross[[1]] <- sum(as.vector((total[[1]][[1]]+total[[2]][[2]])/2))
cross[[1]] <- sum(as.vector((total[[1]][[2]]+total[[2]][[1]])/2))



#####################

#Read sample names
root <- "introduce_root_directory"
setwd(root)
carpetas <- list.dirs(recursive=FALSE)
samples <- c()
for (i in 1:length(carpetas)){
  name <- gsub("./","",carpetas[i])
  samples[i] <- name
}
sample_list <- NULL

for (i in 1:length(samples)){
#Read section names
wd <- paste(sep="", root, samples[i], "/IPSI_clean/")
setwd(wd)
carpetas <- list.dirs(recursive=FALSE)
sections <- c()
for (i in 1:length(carpetas)){
  name <- gsub("./","",carpetas[i])
  sections[i] <- name
}

#browse through sections
data <- NULL
for (i in 1:length(sections)){
  wd2 <- paste(sep="", wd,sections[i])
  sides <- c("green_side", "red_side")
  greenside <- list()
  redside <- list()
  output <- list(greenside, redside)
  #Navegar por los dos lados y leer puntos
    for(j in 1:length(sides)){
    wd_side <- paste(sep="/", wd2,sides[j])
    setwd(wd_side)
    channels <- c("488*.png", "647*.png")
    coords_488 <- read_coords(channels[1])
    coords_647 <- read_coords(channels[2])
    ret_488 <- count_grid(coords_488$y, coords_488$x)
    ret_647 <- count_grid(coords_647$y, coords_647$x)
    output[[j]][[1]] <- ret_488
    output[[j]][[2]] <- ret_647
  }
  data[[i]] <- output
}

#Aggregate sections
total <- list(list(0,0),list(0,0))
for (i in 1:length(sections)){
  for (j in 1:2){
    for (k in 1:2){
      total[[j]][[k]] <- data[[i]][[j]][[k]] + total[[j]][[k]]
    }
  }
}

no_cross <- NULL
cross <- NULL
no_cross[[1]] <- sum(as.vector((total[[1]][[1]]+total[[2]][[2]])/2))
cross[[1]] <- sum(as.vector((total[[1]][[2]]+total[[2]][[1]])/2))

#RATIO CONVERSION
total_488 <- sum(total[[1]][[1]] + total[[2]][[1]])
total_647 <- sum(total[[1]][[2]] + total[[2]][[2]])  
total[[1]][[1]] <- total[[1]][[1]]/total_488
total[[1]][[2]] <- total[[1]][[2]]/total_647
total[[2]][[1]] <- total[[2]][[1]]/total_488
total[[2]][[2]] <- total[[2]][[2]]/total_647

#ASSIGN SAMPLE NAME
assign(samples[i], total)
sample_list[[i]] <- total
}


#Sacar el nombre de las secciones
wd <- paste(sep="", root, samples[i], "/IPSI_clean/")
setwd(wd)
carpetas <- list.dirs(recursive=FALSE)
sections <- c()
for (i in 1:length(carpetas)){
  name <- gsub("./","",carpetas[i])
  sections[i] <- name
}

#Navegar por las Secciones
data <- NULL
for (i in 1:length(sections)){
  wd2 <- paste(sep="", wd,sections[i])
  sides <- c("green_side", "red_side")
  greenside <- list()
  redside <- list()
  output <- list(greenside, redside)
  #Navegar por los dos lados y leer puntos
  for(j in 1:length(sides)){
    wd_side <- paste(sep="/", wd2,sides[j])
    setwd(wd_side)
    channels <- c("488*.png", "647*.png")
    signal_488 <- count_signal(channels[1])
    signal_647 <- count_signal(channels[2])
    output[[j]][[1]] <- signal_488
    output[[j]][[2]] <- signal_647
  }
  data[[i]] <- output
}




CTRL <- grep("^CTRL", samples)  
MUT <- grep("^MUT", samples)
samples_CTRL <- sample_list[CTRL]
samples_MUT <- sample_list[MUT]

total_CTRL <- list(list(0,0),list(0,0))
  for (i in 1:length(samples_CTRL)){
    for (j in 1:2){
      for (k in 1:2){
        total_CTRL[[j]][[k]] <- samples_CTRL[[i]][[j]][[k]] / length(samples_CTRL) + total_CTRL[[j]][[k]]
      }
    }
  }

total_MUT <- list(list(0,0),list(0,0))
for (i in 1:length(samples_MUT)){
  for (j in 1:2){
    for (k in 1:2){
      total_MUT[[j]][[k]] <- samples_MUT[[i]][[j]][[k]] / length(samples_MUT) + total_MUT[[j]][[k]]
    }
  }
} 

samples_sum <- rep.int(list(list(list(0,0),list(0,0))),length(samples)) 
for (i in 1:length(samples)){
  samples_sum[[i]][[1]] <- (sample_list[[i]][[1]][[2]] + sample_list[[i]][[2]][[1]]) / 2
  samples_sum[[i]][[2]] <- (sample_list[[i]][[1]][[1]] + sample_list[[i]][[2]][[2]]) / 2
}

image(samples_sum[[1]][[1]],  col=fiftyGreys, zlim=range(0, samples_CTRL[[1]][[1]][[1]]))


#PLOT

library(RColorBrewer)
cool = rainbow(25, start=rgb2hsv(col2rgb('purple'))[1], end=rgb2hsv(col2rgb('blue'))[1])
warm = rainbow(25, start=rgb2hsv(col2rgb('red'))[1], end=rgb2hsv(col2rgb('yellow'))[1])
cols = c(rev(cool), rev(warm))
mypalette <- colorRampPalette(cols)(255)


av.colors.ctrl.uni <- (total_CTRL[[1]][[2]] + total_CTRL[[2]][[1]]) / 2
av.colors.ctrl.cross <- (total_CTRL[[1]][[1]] + total_CTRL[[2]][[2]]) / 2
lista <- list(av.colors.ctrl.uni, av.colors.ctrl.cross)
image(lista[[1]] * 0.175, col=mypalette, zlim=range(0,max(c(lista[[1]]), lista[[2]], EH13_list[[1]][[1]][[1]])))

av.colors.mut.uni <- (total_MUT[[1]][[2]] + total_MUT[[2]][[1]]) / 2 
av.colors.mut.cross <- (total_MUT[[1]][[1]] + total_MUT[[2]][[2]]) / 2
lista <- list(av.colors.mut.uni, av.colors.mut.cross)
image(lista[[1]], col=mypalette, zlim=range(0,max(c(lista[[1]]), lista[[2]])))

lista <- list(av.colors.ctrl.uni, av.colors.ctrl.cross, av.colors.mut.uni, av.colors.mut.cross)

par(mfrow=c(1,2))
pie(c(sum(lista[[1]]), sum(lista[[2]])), c("uni","cross"), col=c("white","black"))
pie(c(sum(lista[[3]]), sum(lista[[4]])), c("uni","cross"), col=c("white","black"))

colors <- c("#ED1F24","#F05824","#F8981F","#FED00B","#F9ED32","#C7D92B","#70BE44","#72C16C","#72C385","#6FCCDC","#4A9FD7","#4060AC","#5F52A3","#97509F","#D64296","#EA278F")
mypalette <- rev(colorRampPalette(colors)(255))


plot_ctrl_uncross <- image(lista[[1]], col=mypalette, zlim=range(0,max(c(lista[[1]]), lista[[2]], lista[[3]], lista[[4]])))
plot_ctrl_cross <- image(lista[[2]], col=mypalette, zlim=range(0,max(c(lista[[1]]), lista[[2]], lista[[3]], lista[[4]])))
plot_ctrl_uncross <- image(lista[[3]], col=mypalette, zlim=range(0,max(c(lista[[1]]), lista[[2]], lista[[3]], lista[[4]])))
plot_ctrl_cross <- image(lista[[4]], col=mypalette, zlim=range(0,max(c(lista[[1]]), lista[[2]], lista[[3]], lista[[4]])))

########################### DIAGONAL
diag_CTRL <- rep.int(list(list(list(0,0),list(0,0))),length(CTRL)) 
diag <- NULL
for (i in 1:length(samples_CTRL)){
  for (j in 1:2){
    for (k in 1:2){
      for (m in 1:102){
        diag[m] <- sum(samples_CTRL[[i]][[j]][[k]][(103-m),]) + sum(samples_CTRL[[i]][[j]][[k]][,m]) - samples_CTRL[[i]][[j]][[k]][(103-m),m]
      }
      diag_CTRL[[i]][[j]][[k]] <- diag / sum(diag)
    }}}

diag_MUT <- rep.int(list(list(list(0,0),list(0,0))),length(MUT)) 
diag <- NULL
for (i in 1:length(samples_MUT)){
  for (j in 1:2){
    for (k in 1:2){
      for (m in 1:102){
        diag[m] <- sum(samples_MUT[[i]][[j]][[k]][(103-m),]) + sum(samples_MUT[[i]][[j]][[k]][,m]) - samples_MUT[[i]][[j]][[k]][(103-m),m]
      }
      diag_MUT[[i]][[j]][[k]] <- diag / sum(diag)
      }}}

diag_av_CTRL <- list(list(0,0),list(0,0))
for (i in 1:length(diag_CTRL)){
  for (j in 1:2){
    for (k in 1:2){
      diag_av_CTRL[[j]][[k]] <- diag_CTRL[[i]][[j]][[k]] / length(samples_CTRL) + diag_av_CTRL[[j]][[k]]
    }}}
      
diag_av_MUT <- list(list(0,0),list(0,0))
for (i in 1:length(diag_MUT)){
  for (j in 1:2){
    for (k in 1:2){
      diag_av_MUT[[j]][[k]] <- diag_MUT[[i]][[j]][[k]] / length(samples_MUT) + diag_av_MUT[[j]][[k]]
    }}}


########################### EJE X

X_CTRL <- rep.int(list(list(list(0,0),list(0,0))),length(CTRL)) 
diag <- NULL
for (i in 1:length(CTRL)){
  for (j in 1:2){
    for (k in 1:2){
      for (m in 1:102){
        diag[m] <- sum(samples_CTRL[[i]][[j]][[k]][(103-m),])
      }
      X_CTRL[[i]][[j]][[k]] <- diag / sum(diag)
    }}}

X_av_CTRL <- list(list(0,0),list(0,0))
for (i in 1:length(CTRL)){
  for (j in 1:2){
    for (k in 1:2){
      X_av_CTRL[[j]][[k]] <- X_CTRL[[i]][[j]][[k]] / length(samples_CTRL) + X_av_CTRL[[j]][[k]]
    }}}

X_MUT <- rep.int(list(list(list(0,0),list(0,0))),length(MUT)) 
diag <- NULL
for (i in 1:length(MUT)){
  for (j in 1:2){
    for (k in 1:2){
      for (m in 1:102){
        diag[m] <- sum(samples_MUT[[i]][[j]][[k]][(103-m),])
      }
      X_MUT[[i]][[j]][[k]] <- diag / sum(diag)
    }}}

X_av_MUT <- list(list(0,0),list(0,0))
for (i in 1:length(MUT)){
  for (j in 1:2){
    for (k in 1:2){
      X_av_MUT[[j]][[k]] <- X_MUT[[i]][[j]][[k]] / length(samples_MUT) + X_av_MUT[[j]][[k]]
    }}}


########################### EJE Y

Y_CTRL <- rep.int(list(list(list(0,0),list(0,0))),length(CTRL)) 
diag <- NULL
for (i in 1:length(CTRL)){
  for (j in 1:2){
    for (k in 1:2){
      for (m in 1:102){
        diag[m] <- sum(samples_CTRL[[i]][[j]][[k]][,m])
      }
      Y_CTRL[[i]][[j]][[k]] <- diag / sum(diag)
    }}}

Y_av_CTRL <- list(list(0,0),list(0,0))
for (i in 1:length(CTRL)){
  for (j in 1:2){
    for (k in 1:2){
      Y_av_CTRL[[j]][[k]] <- Y_CTRL[[i]][[j]][[k]] / length(samples_CTRL) + Y_av_CTRL[[j]][[k]]
    }}}

Y_MUT <- rep.int(list(list(list(0,0),list(0,0))),length(MUT)) 
diag <- NULL
for (i in 1:length(MUT)){
  for (j in 1:2){
    for (k in 1:2){
      for (m in 1:102){
        diag[m] <- sum(samples_MUT[[i]][[j]][[k]][,m])
      }
      Y_MUT[[i]][[j]][[k]] <- diag / sum(diag)
    }}}

Y_av_MUT <- list(list(0,0),list(0,0))
for (i in 1:length(MUT)){
  for (j in 1:2){
    for (k in 1:2){
      Y_av_MUT[[j]][[k]] <- Y_MUT[[i]][[j]][[k]] / length(samples_MUT) + Y_av_MUT[[j]][[k]]
    }}}


###### UNIFICACION CANALES DIAGONAL Y EJES

av.ctrl.uni <- (diag_av_CTRL[[1]][[2]] + diag_av_CTRL[[2]][[1]]) / 2
av.ctrl.cross <- (diag_av_CTRL[[1]][[1]] + diag_av_CTRL[[2]][[2]]) / 2
list.ctrl <- list(av.ctrl.uni, av.ctrl.cross)

av.mut.uni <- (diag_av_MUT[[1]][[2]] + diag_av_MUT[[2]][[1]]) / 2
av.mut.cross <- (diag_av_MUT[[1]][[1]] + diag_av_MUT[[2]][[2]]) / 2
list.mut <- list(av.mut.uni, av.mut.cross)

x.av.ctrl.uni <- (X_av_CTRL[[1]][[2]] + X_av_CTRL[[2]][[1]]) / 2
x.av.ctrl.cross <- (X_av_CTRL[[1]][[1]] + X_av_CTRL[[2]][[2]]) / 2
listx.ctrl <- list(avx.ctrl.uni, avx.ctrl.cross)

x.av.mut.uni <- (X_av_MUT[[1]][[2]] + X_av_MUT[[2]][[1]]) / 2
x.av.mut.cross <- (X_av_MUT[[1]][[1]] + X_av_MUT[[2]][[2]]) / 2
listx.mut <- list(avx.mut.uni, avx.mut.cross)

y.av.ctrl.uni <- (Y_av_CTRL[[1]][[2]] + Y_av_CTRL[[2]][[1]]) / 2
y.av.ctrl.cross <- (Y_av_CTRL[[1]][[1]] + Y_av_CTRL[[2]][[2]]) / 2
listy.ctrl <- list(avy.ctrl.uni, avy.ctrl.cross)

y.av.mut.uni <- (Y_av_MUT[[1]][[2]] + Y_av_MUT[[2]][[1]]) / 2
y.av.mut.cross <- (Y_av_MUT[[1]][[1]] + Y_av_MUT[[2]][[2]]) / 2
listy.mut <- list(avy.mut.uni, avy.mut.cross)

plot(NULL, xlim=c(0,100), ylim=c(0,1), ylab="%", xlab="dorsomedial-axis", main="Cummulative distribution", type="l")
lapply(list(WT,EH13), function(x) lines(cumsum(x),col="black"))
lines(cumsum(av.ctrl.uni), col="black")
lines(cumsum(av.mut.uni), col="red")


###### POR REPLICAS

ctrl.uni <- list(NULL)
for (i in 1:length(CTRL)){
      ctrl.uni[[i]] <- (diag_CTRL[[i]][[1]][[2]] + diag_CTRL[[i]][[2]][[1]]) /2
}
mut.uni <- list(NULL)
for (i in 1:length(MUT)){
  mut.uni[[i]] <- (diag_MUT[[i]][[1]][[2]] + diag_MUT[[i]][[2]][[1]]) /2
}

ctrl.cross <- list(NULL)
for (i in 1:length(CTRL)){
  ctrl.cross[[i]] <- (diag_CTRL[[i]][[1]][[1]] + diag_CTRL[[i]][[2]][[2]]) /2
}
mut.cross <- list(NULL)
for (i in 1:length(MUT)){
  mut.cross[[i]] <- (diag_MUT[[i]][[1]][[1]] + diag_MUT[[i]][[2]][[2]]) /2
}


x.ctrl.uni <- list(NULL)
for (i in 1:length(CTRL)){
  x.ctrl.uni[[i]] <- (X_CTRL[[i]][[1]][[2]] + X_CTRL[[i]][[2]][[1]]) /2
}
x.mut.uni <- list(NULL)
for (i in 1:length(MUT)){
  x.mut.uni[[i]] <- (X_MUT[[i]][[1]][[2]] + X_MUT[[i]][[2]][[1]]) /2
}

x.ctrl.cross <- list(NULL)
for (i in 1:length(CTRL)){
  x.ctrl.cross[[i]] <- (X_CTRL[[i]][[1]][[1]] + X_CTRL[[i]][[2]][[2]]) /2
}
x.mut.cross <- list(NULL)
for (i in 1:length(MUT)){
  x.mut.cross[[i]] <- (X_MUT[[i]][[1]][[1]] + X_MUT[[i]][[2]][[2]]) /2
}

y.ctrl.uni <- list(NULL)
for (i in 1:length(CTRL)){
  y.ctrl.uni[[i]] <- (Y_CTRL[[i]][[1]][[2]] + Y_CTRL[[i]][[2]][[1]]) /2
}
y.mut.uni <- list(NULL)
for (i in 1:length(MUT)){
  y.mut.uni[[i]] <- (Y_MUT[[i]][[1]][[2]] + Y_MUT[[i]][[2]][[1]]) /2
}

y.ctrl.cross <- list(NULL)
for (i in 1:length(CTRL)){
  y.ctrl.cross[[i]] <- (CTRL[[i]][[1]][[1]] + Y_CTRL[[i]][[2]][[2]]) /2
}
y.mut.cross <- list(NULL)
for (i in 1:length(MUT)){
  y.mut.cross[[i]] <- (MUT[[i]][[1]][[1]] + Y_MUT[[i]][[2]][[2]]) /2
}

  
data.uni <- data.frame(WT=cumsum(av.ctrl.uni), EH13=cumsum(av.mut.uni), dist=(cumsum(av.mut.uni)-cumsum(av.ctrl.uni)))
data.uni[data.uni$dist==max(data.uni$dist),]
data.cross <- data.frame(WT=cumsum(av.ctrl.cross), EH13=cumsum(av.mut.cross), dist=(cumsum(av.mut.cross)-cumsum(av.ctrl.cross)))
data.cross[data.cross$dist==max(data.cross$dist),]

ctrl.uni.dist <- sapply(ctrl.uni, function(x) (cumsum(x))[56])
mut.uni.dist <- sapply(mut.uni, function(x) (cumsum(x))[56])
ctrl.cross.dist <- sapply(ctrl.cross, function(x) (cumsum(x))[56])
mut.cross.dist <- sapply(mut.cross, function(x) (cumsum(x))[56])


data.uni.X <- data.frame(WT=cumsum(avx.ctrl.uni), EH13=cumsum(avx.mut.uni), dist=(cumsum(avx.mut.uni)-cumsum(avx.ctrl.uni)))
data.uni.X[data.uni.X$dist==max(data.uni.X$dist),]
data.cross.X <- data.frame(WT=cumsum(avx.ctrl.cross), EH13=cumsum(avx.mut.cross), dist=(cumsum(avx.mut.cross)-cumsum(avx.ctrl.cross)))
data.cross.X[data.cross.X$dist==max(data.cross.X$dist),]

ctrl.uni.x <- sapply(x.ctrl.uni, function(x) (cumsum(x))[41])
mut.uni.x <- sapply(x.mut.uni, function(x) (cumsum(x))[41])
ctrl.cross.x <- sapply(x.ctrl.cross, function(x) (cumsum(x))[62])
mut.cross.x <- sapply(x.mut.cross, function(x) (cumsum(x))[62])


data.uni.Y <- data.frame(WT=cumsum(avy.ctrl.uni), EH13=cumsum(avy.mut.uni), dist=(cumsum(avy.mut.uni)-cumsum(avy.ctrl.uni)))
data.uni.Y[data.uni.Y$dist==max(data.uni.Y$dist),]
data.cross.Y <- data.frame(WT=cumsum(avy.ctrl.cross), EH13=cumsum(avy.mut.cross), dist=(cumsum(avy.mut.cross)-cumsum(avy.ctrl.cross)))
data.cross.Y[data.cross.Y$dist==max(data.cross.Y$dist),]

ctrl.uni.y <- sapply(y.ctrl.uni, function(x) (cumsum(x))[66])
mut.uni.y <- sapply(y.mut.uni, function(x) (cumsum(x))[66])
ctrl.cross.y <- sapply(y.ctrl.cross, function(x) (cumsum(x))[56])
mut.cross.y <- sapply(y.mut.cross, function(x) (cumsum(x))[56])

#### DATA TO GENERATE PIE CHART PLOTS: UNI CROSS RATIO

uni.ctrl <- sapply(samples_CTRL, function(x) sum(x[[1]][[2]] + x[[2]][[1]])/2)
cross.ctrl <- sapply(samples_CTRL, function(x) sum(x[[1]][[1]] + x[[2]][[2]])/2)

uni.mut <- sapply(samples_MUT, function(x) sum(x[[1]][[2]] + x[[2]][[1]])/2)
cross.mut <- sapply(samples_MUT, function(x) sum(x[[1]][[2]] + x[[2]][[1]])/2)

source(file="I:/INA/KIR2.1/EH10xAF10+CTX/AF10xEH10xEH17+CTX_P11/SCN_40X/SAMPLES/sample_list.R")

x.ctrl.uni <- x.ctrl.uni[1:4]
x.ctrl.cross <- x.ctrl.cross[1:4]


#############
# SCATTERING ANALYSIS
#############
# Throguh the diagonal

ctrl.uni.cum <- lapply(ctrl.uni, function(x) (cumsum(x)/sum(x)))
ctrl.cross.cum <- lapply(ctrl.cross, function(x) (cumsum(x)/sum(x)))
mut.uni.cum <- lapply(mut.uni, function(x) (cumsum(x)/sum(x)))
mut.cross.cum <- lapply(mut.cross, function(x) (cumsum(x)/sum(x)))



ctrl.uni.sem <- NULL
for (j in 1:102){
  ctrl.uni.sem[[j]] <- sem(sapply(ctrl.uni.cum, function(x) (x)[[j]]))
}
ctrl.cross.sem <- NULL
for (j in 1:102){
  ctrl.cross.sem[[j]] <- sem(sapply(ctrl.cross.cum, function(x) (x)[[j]]))
}

mut.uni.sem <- NULL
for (j in 1:102){
  mut.uni.sem[[j]] <- sem(sapply(mut.uni.cum, function(x) (x)[[j]]))
}
mut.cross.sem <- NULL
for (j in 1:102){
  mut.cross.sem[[j]] <- sem(sapply(mut.cross.cum, function(x) (x)[[j]]))
}


plot(NULL, xlim=c(0,100), ylim=c(0,1), ylab="%", xlab="dorsomedial-axis", main="Cummulative distribution", type="l")
lapply(ctrl.uni.cum, function(x) lines(x,col="black"))
lapply(mut.uni.cum, function(x) lines(x,col="red"))

x.index <- c(seq(from=1, to=102, by=1))
plot(NULL, xlim=c(0,100), ylim=c(0,1), ylab="%", xlab="dorsomedial-axis", main="Cummulative distribution", type="l")
polygon(c(x.index,rev(x.index)),c((cumsum(av.ctrl.uni) + ctrl.uni.sem),rev((cumsum(av.ctrl.uni) - ctrl.uni.sem))),border=NA,col=blues9[3])
polygon(c(x.index,rev(x.index)),c((cumsum(av.mut.uni) + mut.uni.sem),rev((cumsum(av.mut.uni) - mut.uni.sem))),border=NA,col=blues9[3])
lines(cumsum(av.ctrl.uni), col="black")
lines(cumsum(av.mut.uni), col="red")

plot(NULL, xlim=c(0,100), ylim=c(0,1), ylab="%", xlab="dorsomedial-axis", main="Cummulative distribution", type="l")
polygon(c(x.index,rev(x.index)),c((cumsum(av.ctrl.cross) + ctrl.cross.sem),rev((cumsum(av.ctrl.cross) - ctrl.cross.sem))),border=NA,col=blues9[3])
polygon(c(x.index,rev(x.index)),c((cumsum(av.mut.cross) + mut.cross.sem),rev((cumsum(av.mut.cross) - mut.cross.sem))),border=NA,col=blues9[3])
lines(cumsum(av.ctrl.cross), col="black")
lines(cumsum(av.mut.cross), col="red")


uni_p <- NULL
for (j in 1:102){
  ctrl <- sapply(ctrl.uni.cum, function(x) (x[[j]]))
  mut <- sapply(mut.uni.cum, function(x) (x[[j]]))
  if(mean(ctrl)==1){
    uni_p[[j]] <- 1
  }
  else{
    uni_p[[j]] <- (t.test(ctrl, mut))$p.value
  }
}
min(na.omit(uni_p))

cross_p <- NULL
for (j in 1:102){
  ctrl <- sapply(ctrl.uni, function(x) (cumsum(x))[[j]])
  mut <- sapply(mut.uni, function(x) (cumsum(x))[[j]])
  if(mean(ctrl)==1){
    cross_p[[j]] <- 1
  }
  else{
    cross_p[[j]] <- (t.test(ctrl, mut))$p.value
  }
}


uni_sig <- NULL
for (i in 1:102){
  uni_sig[i] <- (uni_p[i] < 0.05)
}
uni_sig

cross_sig <- NULL
for (i in 1:102){
  cross_sig[i] <- (cross_p[i] < 0.05)
}
cross_sig

#############
# SCATTERING ANALYSIS
#############
# Throguh X axis

x.ctrl.uni.cum <- lapply(x.ctrl.uni, function(x) (cumsum(x)/sum(x)))
x.ctrl.cross.cum <- lapply(x.ctrl.cross, function(x) (cumsum(x)/sum(x)))
x.mut.uni.cum <- lapply(x.mut.uni, function(x) (cumsum(x)/sum(x)))
x.mut.cross.cum <- lapply(x.mut.cross, function(x) (cumsum(x)/sum(x)))


x.ctrl.uni.sem <- NULL
for (j in 1:102){
  x.ctrl.uni.sem[[j]] <- sem(sapply(x.ctrl.uni.cum, function(x) (x)[[j]]))
}
x.ctrl.cross.sem <- NULL
for (j in 1:102){
  x.ctrl.cross.sem[[j]] <- sem(sapply(x.ctrl.cross.cum, function(x) (x)[[j]]))
}

x.mut.uni.sem <- NULL
for (j in 1:102){
  x.mut.uni.sem[[j]] <- sem(sapply(x.mut.uni.cum, function(x) (x)[[j]]))
}
x.mut.cross.sem <- NULL
for (j in 1:102){
  x.mut.cross.sem[[j]] <- sem(sapply(x.mut.cross.cum, function(x) (x)[[j]]))
}

x.uni_p <- NULL
for (j in 1:102){
  ctrl <- sapply(x.ctrl.uni.cum, function(x) (x[[j]]))
  mut <- sapply(x.mut.uni.cum, function(x) (x[[j]]))
  if(mean(ctrl)==1){
    x.uni_p[[j]] <- 1
  }
  else{
    x.uni_p[[j]] <- (t.test(ctrl, mut))$p.value
  }
}
min(na.omit(uni_p))

x.cross_p <- NULL
for (j in 1:102){
  ctrl <- sapply(x.ctrl.cross.cum, function(x) (x[[j]]))
  mut <- sapply(x.mut.cross.cum, function(x) (x[[j]]))
  if(mean(ctrl)==1){
    x.cross_p[[j]] <- 1
  }
  else{
    x.cross_p[[j]] <- (t.test(ctrl, mut))$p.value
  }
}
min(na.omit(cross_p))

plot(NULL, xlim=c(0,100), ylim=c(0,1), ylab="%", xlab="dorsomedial-axis", main="Cummulative distribution", type="l")
lapply(x.ctrl.uni.cum, function(x) lines(x,col="black"))
lapply(x.mut.uni.cum, function(x) lines(x,col="red"))

x.index <- c(seq(from=1, to=102, by=1))
plot(NULL, xlim=c(0,100), ylim=c(0,1), ylab="%", xlab="dorsomedial-axis", main="Cummulative distribution", type="l")
polygon(c(x.index,rev(x.index)),c((cumsum(x.av.ctrl.uni) + x.ctrl.uni.sem),rev((cumsum(x.av.ctrl.uni) - x.ctrl.uni.sem))),border=NA,col=blues9[3])
polygon(c(x.index,rev(x.index)),c((cumsum(x.av.mut.uni) + x.mut.uni.sem),rev((cumsum(x.av.mut.uni) - x.mut.uni.sem))),border=NA,col=blues9[3])
lines(cumsum(x.av.ctrl.uni), col="black")
lines(cumsum(x.av.mut.uni), col="red")

plot(NULL, xlim=c(0,100), ylim=c(0,1), ylab="%", xlab="dorsomedial-axis", main="Cummulative distribution", type="l")
polygon(c(x.index,rev(x.index)),c((cumsum(x.av.ctrl.cross) + x.ctrl.cross.sem),rev((cumsum(x.av.ctrl.cross) - x.ctrl.cross.sem))),border=NA,col=blues9[3])
polygon(c(x.index,rev(x.index)),c((cumsum(x.av.mut.cross) + x.mut.cross.sem),rev((cumsum(x.av.mut.cross) - x.mut.cross.sem))),border=NA,col=blues9[3])
lines(cumsum(x.av.ctrl.cross), col="black")
lines(cumsum(x.av.mut.cross), col="red")

#############
# SCATTERING ANALYSIS
#############
# Throguh Y axis

y.ctrl.uni.cum <- lapply(y.ctrl.uni, function(x) (cumsum(x)/sum(x)))
y.ctrl.cross.cum <- lapply(y.ctrl.cross, function(x) (cumsum(x)/sum(x)))
y.mut.uni.cum <- lapply(y.ctrl.uni, function(x) (cumsum(x)/sum(x)))
y.mut.cross.cum <- lapply(y.ctrl.cross, function(x) (cumsum(x)/sum(x)))


y.ctrl.uni.sem <- NULL
for (j in 1:102){
  y.ctrl.uni.sem[[j]] <- sem(sapply(y.ctrl.uni.cum, function(x) (x)[[j]]))
}
y.ctrl.cross.sem <- NULL
for (j in 1:102){
  y.ctrl.cross.sem[[j]] <- sem(sapply(y.ctrl.cross.cum, function(x) (x)[[j]]))
}

y.mut.uni.sem <- NULL
for (j in 1:102){
  y.mut.uni.sem[[j]] <- sem(sapply(y.mut.uni.cum, function(x) (x)[[j]]))
}
y.mut.cross.sem <- NULL
for (j in 1:102){
  y.mut.cross.sem[[j]] <- sem(sapply(y.mut.cross.cum, function(x) (x)[[j]]))
}


y.uni_p <- NULL
for (j in 1:102){
  ctrl <- sapply(y.ctrl.uni.cum, function(x) (x[[j]]))
  mut <- sapply(y.mut.uni.cum, function(x) (x[[j]]))
  if(mean(ctrl)==1){
    y.uni_p[[j]] <- 1
  }
  else{
    y.uni_p[[j]] <- (t.test(ctrl, mut))$p.value
  }
}
min(na.omit(uni_p))

y.cross_p <- NULL
for (j in 1:102){
  ctrl <- sapply(y.ctrl.cross.cum, function(x) (x[[j]]))
  mut <- sapply(y.mut.cross.cum, function(x) (x[[j]]))
  if(mean(ctrl)==1){
    y.cross_p[[j]] <- 1
  }
  else{
    y.cross_p[[j]] <- (t.test(ctrl, mut))$p.value
  }
}
min(na.omit(cross_p))

plot(NULL, xlim=c(0,100), ylim=c(0,1), ylab="%", xlab="dorsomedial-axis", main="Cummulative distribution", type="l")
lapply(y.ctrl.uni.cum, function(x) lines(x,col="black"))
lapply(y.mut.uni.cum, function(x) lines(x,col="red"))

x.index <- c(seq(from=1, to=102, by=1))
plot(NULL, xlim=c(0,100), ylim=c(0,1), ylab="%", xlab="dorsomedial-axis", main="Cummulative distribution", type="l")
polygon(c(x.index,rev(x.index)),c((cumsum(y.av.ctrl.uni) + y.ctrl.uni.sem),rev((cumsum(y.av.ctrl.uni) - y.ctrl.uni.sem))),border=NA,col=blues9[3])
polygon(c(x.index,rev(x.index)),c((cumsum(y.av.mut.uni) + y.mut.uni.sem),rev((cumsum(y.av.mut.uni) - y.mut.uni.sem))),border=NA,col=blues9[3])
lines(cumsum(y.av.ctrl.uni), col="black")
lines(cumsum(y.av.mut.uni), col="red")

plot(NULL, xlim=c(0,100), ylim=c(0,1), ylab="%", xlab="dorsomedial-axis", main="Cummulative distribution", type="l")
polygon(c(x.index,rev(x.index)),c((cumsum(y.av.ctrl.cross) + y.ctrl.uni.sem),rev((cumsum(y.av.ctrl.cross) - y.ctrl.uni.sem))),border=NA,col=blues9[3])
polygon(c(x.index,rev(x.index)),c((cumsum(y.av.mut.cross) + y.ctrl.uni.sem),rev((cumsum(y.av.mut.cross) - y.ctrl.uni.sem))),border=NA,col=blues9[3])
lines(cumsum(y.av.ctrl.cross), col="black")
lines(cumsum(y.av.mut.cross), col="red")
lines((cumsum(y.av.ctrl.cross) + y.ctrl.cross.sem))
lines((cumsum(y.av.ctrl.cross) - y.ctrl.cross.sem))


x.ctrl.CI <- NULL
for (j in 1:102){
  x.ctrl.CI[[j]] <- ci.mean(sapply(x.ctrl.cross, function(x) (cumsum(x))[[j]]))
}

x.mut.CI <- NULL
for (j in 1:102){
  x.mut.CI[[j]] <- ci.mean(sapply(x.mut.cross, function(x) (cumsum(x))[[j]]))
}

sem <- function(x) {sd(x)/sqrt(length(x))}

x.ctrl.sem <- NULL
for (j in 1:102){
  x.ctrl.sem[[j]] <- sem(sapply(x.ctrl.cross, function(x) (cumsum(x))[[j]]))
}

x.mut.sem <- NULL
for (j in 1:102){
  x.mut.sem[[j]] <- sem(sapply(x.mut.cross, function(x) (cumsum(x))[[j]]))
}

x.ctrl_p <- NULL
for (j in 1:102){
  ctrl <- sapply(x.ctrl.uni, function(x) (cumsum(x))[[j]])
  mut <- sapply(x.mut.uni, function(x) (cumsum(x))[[j]])
  if(mean(ctrl)==1){
    x.ctrl_p[[j]] <- 1
    }
  else{
    x.ctrl_p[[j]] <- (t.test(ctrl, mut))$p.value
    }
}
min(na.omit(x.ctrl_p))


ctrl.uni_p <- NULL
for (j in 1:102){
  ctrl <- sapply(ctrl.uni, function(x) (cumsum(x))[[j]])
  mut <- sapply(mut.uni, function(x) (cumsum(x))[[j]])
  if(mean(ctrl)==1){
    ctrl.uni_p[[j]] <- 0
  }
  else{
    ctrl.uni_p[[j]] <- (t.test(ctrl, mut))$p.value
  }
}
min(na.omit(ctrl.uni_p))


x.ctrl.mean <- sapply(x.ctrl.CI, function(x) x$mean)
x.ctrl.upper <- sapply(x.ctrl.CI, function(x) x$upper)
x.ctrl.lower <- sapply(x.ctrl.CI, function(x) x$lower)

x.mut.mean <- sapply(x.mut.CI, function(x) x$mean)
x.mut.upper <- sapply(x.mut.CI, function(x) x$upper)
x.mut.lower <- sapply(x.mut.CI, function(x) x$lower)

x.new <- c(seq(from=1, to=102, by=1))

plot(mean,type="n")
polygon(c(x.new,rev(x.new)),c((x.ctrl.mean + x.ctrl.sem),rev((x.ctrl.mean - x.ctrl.sem))),border=NA,col=blues9[3])
polygon(c(x.new,rev(x.new)),c((x.mut.mean + x.mut.sem),rev((x.mut.mean - x.mut.sem))),border=NA,col="coral")
lines(x.ctrl.mean, col="black")
lines(x.mut.mean, col="red")


