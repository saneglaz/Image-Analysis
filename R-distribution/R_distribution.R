############################
#Load libraries
library(prettyR)
library(data.table)
library(Hmisc)
library(gdata)


#Insert path to the root directory
root <- "I:/INA/KIR2.1/dLGN/CUANTIFICACION/AUTO/P11"
setwd(root)
#Obtain sample names using folder names
folders <- list.dirs(recursive=FALSE)
samples <- c()
for (i in 1:length(folders)){
  name <- gsub("./","",folders[i])
  samples[i] <- name
}

#PROCESSING LOOP: GENERATE R DISTRIBUTION
samples_list <- list()
for (i in 1:length(samples)){
  wd <- paste(sep="/", root,samples[i])
  wd_Rdist <- paste(sep="/", wd,"R_DISTRIBUTION")
  setwd(wd_Rdist)
  sample <- samples[i]
  data <-Rdistribution()
  data <- lapply(data, function(x) x[x!=0])
  data <- c(data[[1]], data[[2]], data[[3]])
  assign(sample, data)
  samples_list[[i]] <- data 
}

#############################
#DATA ANALYSIS#
#############################
#GROUP SAMPLES BY CONDITION
rm(samples_WT,samples_EH13,samples_EH10)
objetos <- objects()
samples_WT <- grep("WT", objetos, value=TRUE)
samples_EH13 <- grep("EH13", objetos, value=TRUE)
samples_EH10 <- grep("EH10", objetos, value=TRUE)
samples_list=c(samples_WT,samples_EH13,samples_EH10)

WT <- grep("WT", samples, value=FALSE)
EH13 <- grep("EH13", samples, value=FALSE)
EH10 <- grep("EH10", samples, value=FALSE)


#GENERAR DATA.FRAME CON RESULTADOS
results = data.frame(genotype=genotypes,
                     sample=rep(0,14),
                     overlap=rep(0,14),
                     ipsidom=rep(0,14),
                     contradom=rep(0,14),
                     total_area=rep(0,14),
                     contra_area=rep(0,14))
for(i in 1:length(samples_list)){
  overlap <- lapply(samples_list[[i]], function(x) {length(x[x<0.5 & x>-0.5])})
  total_px <- as.numeric(unlist(lapply(samples_list[i], length)))
  ratio_overlap <- unlist(overlap)/unlist(total_px)
  
  ipsidom <- lapply(samples_list[i], function(x) {length(x[x> 1.75])})
  ratio_ipsidom <- as.numeric(unlist(ipsidom)/unlist(total_px))
  contra_area <- as.numeric((unlist(total_px) - unlist(ipsidom))/unlist(total_px))
  
  contradom<- lapply(samples_list[i], function(x) {(x[abs(x)> 1.75])})
  contradom<- lapply(contradom, function(x) {length(x[x< 0])})
  ratio_contradom <- as.numeric(unlist(contradom)/unlist(total_px))
  
  results[i, ] = c(genotypes[i],samples[i],
                   mean(as.numeric(unlist(ratio_overlap))),
                   mean(as.numeric(unlist(ratio_ipsidom))),
                   mean(as.numeric(unlist(ratio_contradom))),
                   mean(as.numeric(unlist(total_px))),
                   mean(as.numeric(unlist(contra_area))))
}

overlap <- sapply(samples_list, function(x) {length(x[x<0.5 & x>-0.5])})
total_px <- as.numeric(sapply(samples_list, length))
ratio_overlap <- overlap/total_px

ipsidom <- sapply(samples_list, function(x) {length(x[x> 1.75])})
ratio_ipsidom <- ipsidom/total_px
contra_area <- total_px - (ipsidom/total_px)

contradom<- sapply(samples_list, function(x) {(x[abs(x)> 1.75])})
contradom<- sapply(contradom, function(x) {length(x[x< 0])})
ratio_contradom <- contradom/total_px


results$overlap_ <- as.numeric(results$overlap)/rep(mean(as.numeric(results[results=="WT", ]$overlap)), length(as.numeric(results$overlap)))    
results$ipsidom_ <- as.numeric(results$ipsidom)/rep(mean(as.numeric(results[results=="WT", ]$ipsidom)), length(as.numeric(results$ipsidom)))    
results$contra_area_ <- as.numeric(results$contra_area)/rep(mean(as.numeric(results[results=="WT", ]$contra_area)), length(as.numeric(results$contra_area)))    
results$contradom_ <- as.numeric(results$contradom)/rep(mean(as.numeric(results[results=="WT", ]$contradom)), length(as.numeric(results$contradom)))    

overlap_ <- overlap / rep(mean(overlap[WT]), length(overlap))
ipsidom_ <- ipsidom / rep(mean(ipsidom[WT]), length(ipsidom))
contra_area_ <- contra_area / rep(mean(contra_area[WT]), length(contra_area))
contradom_ <- contradom / rep(mean(contradom[WT]), length(contradom))


results$genotype <- reorder(results$genotype, new.order=c("WT","EH13", "EH10"))
results <- results[order(results$genotype), ]

overlap.table <- data.frame(name=c("WT","EH13","EH10"), 
                            value=c(mean(overlap_[WT]), mean(overlap_[EH13]), mean(overlap_[EH10])), 
                            sem=c(se(overlap_[WT]), se(overlap_[EH13]), se(overlap_[EH10])))

ipsidom.table <- data.frame(name=c("WT","EH13","EH10"), 
                            value=c(mean(ipsidom_[WT]), mean(ipsidom_[EH13]), mean(ipsidom_[EH10])), 
                            sem=c(se(ipsidom_[WT]), se(ipsidom_[EH13]), se(ipsidom_[EH10])))


######HISTOGRAM

#Histogram: Average group
WT.list <- grep("WT", samples_list)
EH13.list <- grep("EH13", samples_list)
EH10.list <- grep("EH10", samples_list)
breaks <- seq(from=-2.5, to=2.5, by=0.005)
hist.list <- lapply(samples_list, function(x) hist(get(x), breaks = breaks))

list.index <- list(WT.list, EH13.list, EH10.list)
output <- c("WT.mean", "EH13.mean", "EH10.mean")
for(j in 1:length(list.index)){
  sublist <- list.index[[j]]
  total <- 0
  for(i in 1:length(sublist)){
    total <- hist.list[[sublist[[i]]]]$density + total
  }
  assign(output[j], total / length(sublist))
}

mean.counts <- mean(sapply(hist.list, function(x) sum(x$counts)))
WT.counts <- (WT.mean * mean.counts) / 100
EH13.counts <- (EH13.mean * mean.counts) / 100
EH10.counts <- (EH10.mean * mean.counts) / 100

lines <- list(WT.counts, EH13.counts, EH10.counts)
line_names <- c("WT.group", "EH13.group", "EH10.group")
for (j in 1:(length(lines))){
  total <- NULL
  for (i in 1:(length(breaks)-1)){
    count <- rep.int(breaks[i], (lines[[j]])[i])
    total <- append(total, count, after=length(total))
    assign(line_names[j],total)  
  }
}

plot(density(WT.group, bw=0.25))
lines(density(EH13.group, bw=0.25), col="red")
lines(density(EH10.group, bw=0.25), col="blue")

plot(NULL, xlim=c(-3,3), ylim=c(0,1.2), ylab="%", xlab="distance(px)", main="Cummulative distribution", type="l")
plot(WT)



############  FUNCTIONS

Rdistribution <- function(X) {   
  data <- list()
  filescontra <- Sys.glob("*contra.txt")
  filesipsi <- Sys.glob("*ipsi.txt")
  contra <- NULL
  ipsi <- NULL
  for (i in 1:length(filescontra))
  {
    basename <-gsub("*_contra.txt","", filescontra[i])
    name <- paste(sep="", basename,"_Rdist.txt")
    file_check <- paste(wd_Rdist,name,sep="/")
    
    if(file.exists(file_check)){
      Rdist <- read.delim(file_check, header = FALSE, sep = "\t")
      Rdist <- Rdist[!is.na(Rdist)]
      data[[i]]<- Rdist
    }
    else{
      contra <- read.delim(filescontra[i], header = FALSE, sep = "\t")
      ipsi <- read.delim(filesipsi[i], header = FALSE, sep = "\t")
      
      contra_ <- contra
      ipsi_ <- ipsi
      
      contra_[contra_<5] <- 0
      ipsi_[ipsi_<5] <- 0
      
      suma <- ipsi_ + contra_
      suma[suma==0] <- NA
      
      contra_[is.na(suma)] <- NA
      ipsi_[is.na(suma)] <- NA
      
      contra_[contra_ == 0] <- 0.01
      ipsi_[ipsi_ == 0] <- 0.01
      
      Rdist <- log(ipsi_/contra_, 10)
      Rdist_ <- Rdist
      Rdist_[is.na(Rdist_)] <- 0 
      
      folder <- getwd()
      folder <- paste(sep="", folder,"/")
      output <- paste(sep="", folder,name)
      write.table(Rdist_, output, sep="\t", col.names = FALSE, row.names= FALSE)
      Rdist <- Rdist[!is.na(Rdist)]
      data[[i]]<- Rdist[Rdist!=0]
    }
  }
  return(data)
}

se <- function(x) sqrt(var(x)/length(x))
########################  




