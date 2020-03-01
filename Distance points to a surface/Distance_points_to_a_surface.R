############################
#Insert path to the root directory

root <- c("I:/INA/KIR2.1/SC/CUANTIFICACION/DIST_SUP/P30")
setwd(root)

#Read the name of the samples from the nam of the folders
folders <- list.dirs(recursive=FALSE)
samples <- c()
for (i in 1:length(folders)){
  name <- gsub("./","",folders[i])
  samples[i] <- name
}

#Generate files with the coordinates
for (i in 1:length(samples)){
  wd <- paste(sep="/", root,samples[i])
  setwd(wd)
  get_puntos(wd)
  i
}  

#Read distances
samples_list <- list()
for (i in 1:length(samples)){
  wd <- paste(sep="/", root,samples[i],"DISTANCES/")
  setwd(wd)
  sample <- read_dist()
  #assign(samples[i],sample)
  samples_list[[i]] <- sample
}  

# Sort samples by group
WT <- grep("WT", samples, value=FALSE)
EH13 <- grep("EH13", samples, value=FALSE)
EH10 <- grep("EH10", samples, value=FALSE)


# Average per group
average.WT <- c(unlist(samples_list[WT]))
average.EH13 <- c(unlist(samples_list[EH13]))
average.EH10 <- c(unlist(samples_list[EH10]))


###########
#HISTOGRAM: 
##########

#PLOT DENSITIES
densities_list <- list()
densities_list <- lapply(samples_list, density)
plot(NULL, xlim=c(0,2000), ylim=c(0,0.0045), ylab="%", xlab="distance", main="distancia a superficie", type="l")
polygon(density(average.WT), col=rgb(0, 0, 0,0.25))
polygon(density(average.EH13), col=rgb(1, 0, 0,0.25))
polygon(density(average.EH10), col=rgb(0, 0, 1,0.25))

lapply(densities_list[WT], lines, xlim=c(0,800), col="black")
lapply(densities_list[EH13], lines, xlim=c(0,650), col="red")
lapply(densities_list[EH10], lines, xlim=c(0,650), col="blue")

#Cummulative distribution plot with replicates 
plot(NULL, xlim=c(0,1000), ylim=c(0,1), ylab="%", xlab="distance(px)", main="Cummulative distribution", type="l")
lapply(samples_list[WT], function(x) plot(ecdf(x),add=TRUE,col="black"))
lapply(samples_list[EH13], function(x) plot(ecdf(x),add=TRUE,col="red"))
lapply(samples_list[EH10], function(x) plot(ecdf(x),add=TRUE,col="blue"))

#Rescale distance vectors
hist.list <- lapply(samples_list, function(x) hist(x, breaks=seq(min(x),max(x),l=1001)))
mean.counts <- mean(sapply(hist.list, function(x) sum(x$counts)))
WT.counts <- rowMeans(data.frame(sapply(hist.list[WT], function(x) x$density * mean.counts)))/20
EH10.counts <- rowMeans(data.frame(sapply(hist.list[EH10], function(x) x$density * mean.counts)))/20
EH13.counts <- rowMeans(data.frame(sapply(hist.list[EH13], function(x) x$density * mean.counts)))/20

breaks <- seq(0,1000, by=1)
lines <- list(WT.counts, EH13.counts, EH10.counts)
line_names <- c("WT_cum", "EH13_cum", "EH10_cum")
for (j in 1:(length(lines))){
  total <- NULL
  for (i in 1:(length(breaks)-1)){
    count <- rep.int(breaks[i], (lines[[j]])[i])
    total <- append(total, count, after=length(total))
    assign(line_names[j],total)  
  }
}

samples_resampled <- lapply(hist.list, function(x) (x$density * mean.counts) / 100)
samples_x <- lapply(hist.list, function(x) (x$density))

samples_rework <- list()
for (j in 1:(length(samples_resampled))){
  total <- NULL
  for (i in 1:(length(breaks)-1)){
    count <- rep.int(breaks[i], (samples_resampled[[j]])[i])
    total <- append(total, count, after=length(total))
  }
  samples_rework[[j]] <- total  
}

samples_x <- lapply(samples_resampled, function(x) (cumsum(x)/sum(x)))


#PLOT cummulative distribution by group color
plot(ecdf(WT.counts), col="black")
lines(ecdf(EH13.counts), col="red")  
lines(ecdf(EH10.counts), col="blue")  

plot(samples_freq[[1]],type="n")
lapply(samples_x[WT], lines, col="black")
lapply(samples_x[EH13], lines, col="red")
lapply(samples_x[EH10], lines, col="blue")

plot(WT_mean,type="n")
lines(WT_mean, col="black") 
lines(EH13_mean, col="red") 
lines(EH10_mean, col="blue") 

plot(ecdf(average.WT), col="black")
lines(ecdf(average.EH13), col="red")  
lines(ecdf(average.EH10), col="blue")

# Find the point of the cummulative distribution with the max divergence with the wild-type

cdf.WT <- ecdf(WT_cum)
cdf.EH13 <- ecdf(EH13_cum)
cdf.EH10 <- ecdf(EH10_cum)

genotypes<- c(rep("WT",length(WT)),rep("EH13",length(EH13)),rep("EH10",length(EH10)))
group <- c(rep("WT", length(WT)), rep("EH13", length(EH13)))
dat <- data.frame(KSD = c(WT,EH13), group = group)

minMax <- seq(min(average.WT, average.EH13), max(average.WT, average.EH13), length.out=length(average.WT)) 
x0 <- minMax[which( abs(cdf.WT(minMax) - cdf.EH13(minMax)) == max(abs(cdf.WT(minMax) - cdf.EH13(minMax))) )] 
y0 <- cdf.WT(x0) 
y1 <- cdf.EH13(x0)

minMax.EH10 <- seq(min(average.WT, average.EH10), max(average.WT, average.EH10), length.out=length(average.WT)) 
x0.EH10 <- minMax.EH10[which( abs(cdf.WT(minMax.EH10) - cdf.EH10(minMax.EH10)) == max(abs(cdf.WT(minMax.EH10) - cdf.EH10(minMax.EH10))) )] 
y0.EH10 <- cdf.WT(x0.EH10) 
y1.EH10 <- cdf.EH10(x0.EH10)


sapply(samples_list[WT], function(x) {mean(length(x[x< 586]) / length(x))})

#Calculate the sem (standard error of the mean) for each group

sem <- function(x) {sd(x)/sqrt(length(x))}

WT_sem <- NULL
for(i in length(samples_x[WT]))
  for (j in 1:1000){
    WT_sem[[j]] <- sem(sapply(samples_x[WT], function(x) (x)[[j]]))
  }  
EH13_sem <- NULL
for(i in length(samples_x[EH13]))
  for (j in 1:1000){
    EH13_sem[[j]] <- sem(sapply(samples_x[EH13], function(x) (x)[[j]]))
  } 
EH10_sem <- NULL
for(i in length(samples_x[EH10]))
  for (j in 1:1000){
    EH10_sem[[j]] <- sem(sapply(samples_x[EH10], function(x) (x)[[j]]))
  } 

#Calculate the mean of the groups

WT_mean <- NULL
for(i in length(samples_x[WT]))
  for (j in 1:1000){
    WT_mean[[j]] <- mean(sapply(samples_x[WT], function(x) (x)[[j]]))
  } 
EH13_mean <- NULL
for(i in length(samples_x[EH13]))
  for (j in 1:1000){
    EH13_mean[[j]] <- mean(sapply(samples_x[EH13], function(x) (x)[[j]]))
  } 
EH10_mean <- NULL
for(i in length(samples_x[EH10]))
  for (j in 1:1000){
    EH10_mean[[j]] <- mean(sapply(samples_x[EH10], function(x) (x)[[j]]))
  } 

#CUMMULATIVE PLOT  WITH REPLICATES

x.index <- c(seq(from=1, to=1000, by=1))

plot(WT_mean,type="n")
polygon(c(x.index,rev(x.index)),c((WT_mean + WT_sem),rev((WT_mean - WT_sem))),border=NA,col=blues9[3])
polygon(c(x.index,rev(x.index)),c((EH13_mean + EH13_sem),rev((EH13_mean - EH13_sem))),border=NA,col=blues9[3])
polygon(c(x.index,rev(x.index)),c((EH10_mean + EH10_sem),rev((EH10_mean - EH10_sem))),border=NA,col=blues9[3])
lines(WT_mean, col="black")
lines(EH13_mean, col="black")
lines(EH10_mean, col="black")


# CUMMULATIVE PLOT FOR GROUPS
# The envelope represents the sem of each group

WT_freq <- cumsum(WT_mean)/sum(WT_mean)
EH13_freq <- cumsum(EH13_mean)/sum(EH13_mean)
EH10_freq <- cumsum(EH10_mean)/sum(EH10_mean)

plot(cumsum(WT_mean)/sum(WT_mean),type="n")
polygon(c(x.index,rev(x.index)),c((WT_freq + WT_sem),rev((WT_freq - WT_sem))),border=NA,col=blues9[3])
polygon(c(x.index,rev(x.index)),c((EH13_freq + EH13_sem),rev((EH13_freq - EH13_sem))),border=NA,col=blues9[3])
polygon(c(x.index,rev(x.index)),c((EH10_freq + EH10_sem),rev((EH10_freq - EH10_sem))),border=NA,col=blues9[3])
lines(cumsum(WT_mean)/sum(WT_mean), col="black") 
lines(cumsum(EH13_mean)/sum(EH13_mean), col="red") 
lines(cumsum(EH10_mean)/sum(EH10_mean), col="blue")



#CUMMULATIVE PLOT  FOR GROUPS
plot(cdf.WT, verticals=TRUE, do.points=FALSE, col="black", lwd=4) 
plot(cdf.EH13, verticals=TRUE, do.points=FALSE, col="red", lwd=4, add=TRUE) 
plot(cdf.EH10, verticals=TRUE, do.points=FALSE, col="blue", lwd=4, add=TRUE) 

points(c(x0, x0), c(y0, y1), pch=16, col="red") 
segments(x0, y0, x0, y1, col="black", lty="dotted") 

points(c(x0.EH10, x0.EH10), c(y0.EH10, y1.EH10), pch=16, col="blue") 
segments(x0.EH10, y0.EH10, x0.EH10, y1.EH10, col="black", lty="dotted")


#VIOLIN PLOT

d1 <- data.frame(
  y = c(WT_mean, EH13_mean, EH10_mean),
  group = rep(c("WT", "EH13", "EH10"), c(length(WT_mean), length(EH13_mean), length(EH10_mean)))
)
d1 <- data.frame(
  y = c(log(WT.counts), log(EH13.counts), log(EH10.counts)),
  group = rep(c("WT", "EH13", "EH10"), c(length(WT.counts), length(EH13.counts), length(EH10.counts))))

e <- ggplot(d1, aes(x = group, y = y))
e + geom_violin(aes(fill = group), trim = FALSE) + 
  geom_boxplot(width = 0.2)+
  scale_fill_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"))+
  theme(legend.position = "none")



############################## FUNCTIONS

#Exportar coordenadas puntos para thresholded image
get_puntos <- function(wd) {
  library(png)
  library(jpeg)
  dir.create("PUNTOS", showWarnings = FALSE)  
  fileNames <- Sys.glob("*.png")
  for (i in 1:length(fileNames)){
    d <- readPNG(fileNames[i])
    coord <- which(d==1, arr.ind=TRUE)
    colnames(coord) <- c("x","y","int")
    coords <- data.frame(coord)
    coords <- data.frame(x=coords$x, y=coords$y)
    assign(fileNames[i], coords)
    
    basename <-gsub("*.png","", fileNames[i])
    name <- paste(sep="", getwd(), "/", "PUNTOS/", basename,".png_points.csv")
    
    write.csv(coords, file= name, sep=",", col.names = TRUE, row.names= TRUE)
  }
}

#Leer files y agrupar por muestra ()
read_dist <- function() {
  muestra <- list()
  results <- list()
  data = NULL
  files_results <- Sys.glob("*distances.csv")
  datalist <- list()
  for (i in 1:length(files_results)){
    x = 0
    table <- read.table(files_results[i], header=TRUE, sep=",")
    if (x != 0) {    
      x <- data.frame(table$Mean)
      bad <- is.na(x)
      x <- x[!bad]
      x <- data.frame(x)
    }
    else {  
      x <- table
    }
    data <- c(data, x$X)
  }  
  #data <- rbind(data,x)
  return(data)
}
