############################
#Insert path to the root directory

root <- c("introduce_root_directory")
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
MUT1 <- grep("MUT1", samples, value=FALSE)
MUT2 <- grep("MUT2", samples, value=FALSE)


# Average per group
average.WT <- c(unlist(samples_list[WT]))
average.MUT1 <- c(unlist(samples_list[MUT1]))
average.MUT2 <- c(unlist(samples_list[MUT2]))


###########
#HISTOGRAM: 
##########

#PLOT DENSITIES
densities_list <- list()
densities_list <- lapply(samples_list, density)
plot(NULL, xlim=c(0,2000), ylim=c(0,0.0045), ylab="%", xlab="distance", main="distancia a superficie", type="l")
polygon(density(average.WT), col=rgb(0, 0, 0,0.25))
polygon(density(average.MUT1), col=rgb(1, 0, 0,0.25))
polygon(density(average.MUT2), col=rgb(0, 0, 1,0.25))

lapply(densities_list[WT], lines, xlim=c(0,800), col="black")
lapply(densities_list[MUT1], lines, xlim=c(0,650), col="red")
lapply(densities_list[MUT2], lines, xlim=c(0,650), col="blue")

#Cummulative distribution plot with replicates 
plot(NULL, xlim=c(0,1000), ylim=c(0,1), ylab="%", xlab="distance(px)", main="Cummulative distribution", type="l")
lapply(samples_list[WT], function(x) plot(ecdf(x),add=TRUE,col="black"))
lapply(samples_list[MUT1], function(x) plot(ecdf(x),add=TRUE,col="red"))
lapply(samples_list[MUT2], function(x) plot(ecdf(x),add=TRUE,col="blue"))

#Rescale distance vectors
hist.list <- lapply(samples_list, function(x) hist(x, breaks=seq(min(x),max(x),l=1001)))
mean.counts <- mean(sapply(hist.list, function(x) sum(x$counts)))
WT.counts <- rowMeans(data.frame(sapply(hist.list[WT], function(x) x$density * mean.counts)))/20
MUT2.counts <- rowMeans(data.frame(sapply(hist.list[MUT2], function(x) x$density * mean.counts)))/20
MUT1.counts <- rowMeans(data.frame(sapply(hist.list[MUT1], function(x) x$density * mean.counts)))/20

breaks <- seq(0,1000, by=1)
lines <- list(WT.counts, MUT1.counts, MUT2.counts)
line_names <- c("WT_cum", "MUT1_cum", "MUT2_cum")
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
lines(ecdf(MUT1.counts), col="red")  
lines(ecdf(MUT2.counts), col="blue")  

plot(samples_freq[[1]],type="n")
lapply(samples_x[WT], lines, col="black")
lapply(samples_x[MUT1], lines, col="red")
lapply(samples_x[MUT2], lines, col="blue")

plot(WT_mean,type="n")
lines(WT_mean, col="black") 
lines(MUT1_mean, col="red") 
lines(MUT2_mean, col="blue") 

plot(ecdf(average.WT), col="black")
lines(ecdf(average.MUT1), col="red")  
lines(ecdf(average.MUT2), col="blue")

# Find the point of the cummulative distribution with the max divergence with the wild-type

cdf.WT <- ecdf(WT_cum)
cdf.MUT1 <- ecdf(MUT1_cum)
cdf.MUT2 <- ecdf(MUT2_cum)

genotypes<- c(rep("WT",length(WT)),rep("MUT1",length(MUT1)),rep("MUT2",length(MUT2)))
group <- c(rep("WT", length(WT)), rep("MUT1", length(MUT1)))
dat <- data.frame(KSD = c(WT,MUT1), group = group)

minMax <- seq(min(average.WT, average.MUT1), max(average.WT, average.MUT1), length.out=length(average.WT)) 
x0 <- minMax[which( abs(cdf.WT(minMax) - cdf.MUT1(minMax)) == max(abs(cdf.WT(minMax) - cdf.MUT1(minMax))) )] 
y0 <- cdf.WT(x0) 
y1 <- cdf.MUT1(x0)

minMax.MUT2 <- seq(min(average.WT, average.MUT2), max(average.WT, average.MUT2), length.out=length(average.WT)) 
x0.MUT2 <- minMax.MUT2[which( abs(cdf.WT(minMax.MUT2) - cdf.MUT2(minMax.MUT2)) == max(abs(cdf.WT(minMax.MUT2) - cdf.MUT2(minMax.MUT2))) )] 
y0.MUT2 <- cdf.WT(x0.MUT2) 
y1.MUT2 <- cdf.MUT2(x0.MUT2)


sapply(samples_list[WT], function(x) {mean(length(x[x< 586]) / length(x))})

#Calculate the sem (standard error of the mean) for each group

sem <- function(x) {sd(x)/sqrt(length(x))}

WT_sem <- NULL
for(i in length(samples_x[WT]))
  for (j in 1:1000){
    WT_sem[[j]] <- sem(sapply(samples_x[WT], function(x) (x)[[j]]))
  }  
MUT1_sem <- NULL
for(i in length(samples_x[MUT1]))
  for (j in 1:1000){
    MUT1_sem[[j]] <- sem(sapply(samples_x[MUT1], function(x) (x)[[j]]))
  } 
MUT2_sem <- NULL
for(i in length(samples_x[MUT2]))
  for (j in 1:1000){
    MUT2_sem[[j]] <- sem(sapply(samples_x[MUT2], function(x) (x)[[j]]))
  } 

#Calculate the mean of the groups

WT_mean <- NULL
for(i in length(samples_x[WT]))
  for (j in 1:1000){
    WT_mean[[j]] <- mean(sapply(samples_x[WT], function(x) (x)[[j]]))
  } 
MUT1_mean <- NULL
for(i in length(samples_x[MUT1]))
  for (j in 1:1000){
    MUT1_mean[[j]] <- mean(sapply(samples_x[MUT1], function(x) (x)[[j]]))
  } 
MUT2_mean <- NULL
for(i in length(samples_x[MUT2]))
  for (j in 1:1000){
    MUT2_mean[[j]] <- mean(sapply(samples_x[MUT2], function(x) (x)[[j]]))
  } 

#CUMMULATIVE PLOT  WITH REPLICATES

x.index <- c(seq(from=1, to=1000, by=1))

plot(WT_mean,type="n")
polygon(c(x.index,rev(x.index)),c((WT_mean + WT_sem),rev((WT_mean - WT_sem))),border=NA,col=blues9[3])
polygon(c(x.index,rev(x.index)),c((MUT1_mean + MUT1_sem),rev((MUT1_mean - MUT1_sem))),border=NA,col=blues9[3])
polygon(c(x.index,rev(x.index)),c((MUT2_mean + MUT2_sem),rev((MUT2_mean - MUT2_sem))),border=NA,col=blues9[3])
lines(WT_mean, col="black")
lines(MUT1_mean, col="black")
lines(MUT2_mean, col="black")


# CUMMULATIVE PLOT FOR GROUPS
# The envelope represents the sem of each group

WT_freq <- cumsum(WT_mean)/sum(WT_mean)
MUT1_freq <- cumsum(MUT1_mean)/sum(MUT1_mean)
MUT2_freq <- cumsum(MUT2_mean)/sum(MUT2_mean)

plot(cumsum(WT_mean)/sum(WT_mean),type="n")
polygon(c(x.index,rev(x.index)),c((WT_freq + WT_sem),rev((WT_freq - WT_sem))),border=NA,col=blues9[3])
polygon(c(x.index,rev(x.index)),c((MUT1_freq + MUT1_sem),rev((MUT1_freq - MUT1_sem))),border=NA,col=blues9[3])
polygon(c(x.index,rev(x.index)),c((MUT2_freq + MUT2_sem),rev((MUT2_freq - MUT2_sem))),border=NA,col=blues9[3])
lines(cumsum(WT_mean)/sum(WT_mean), col="black") 
lines(cumsum(MUT1_mean)/sum(MUT1_mean), col="red") 
lines(cumsum(MUT2_mean)/sum(MUT2_mean), col="blue")



#CUMMULATIVE PLOT  FOR GROUPS
plot(cdf.WT, verticals=TRUE, do.points=FALSE, col="black", lwd=4) 
plot(cdf.MUT1, verticals=TRUE, do.points=FALSE, col="red", lwd=4, add=TRUE) 
plot(cdf.MUT2, verticals=TRUE, do.points=FALSE, col="blue", lwd=4, add=TRUE) 

points(c(x0, x0), c(y0, y1), pch=16, col="red") 
segments(x0, y0, x0, y1, col="black", lty="dotted") 

points(c(x0.MUT2, x0.MUT2), c(y0.MUT2, y1.MUT2), pch=16, col="blue") 
segments(x0.MUT2, y0.MUT2, x0.MUT2, y1.MUT2, col="black", lty="dotted")


#VIOLIN PLOT

d1 <- data.frame(
  y = c(WT_mean, MUT1_mean, MUT2_mean),
  group = rep(c("WT", "MUT1", "MUT2"), c(length(WT_mean), length(MUT1_mean), length(MUT2_mean)))
)
d1 <- data.frame(
  y = c(log(WT.counts), log(MUT1.counts), log(MUT2.counts)),
  group = rep(c("WT", "MUT1", "MUT2"), c(length(WT.counts), length(MUT1.counts), length(MUT2.counts))))

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
