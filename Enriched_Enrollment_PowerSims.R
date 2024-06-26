setwd("~/sandbox/nnlab/clinicaltrial-poweranalyses")
library(dplyr)
library(tidyr)
library(magrittr)
library(ggpubr)

## Full available data
data_full = read.table(file="./data_full.csv",header=T,sep=",",dec=".")

## Filtered data, who would be enrolled according to a certain enrollment model 
## (for e.g using the model that predicts whether an individual would 
## show cognitive decline at the end of the trial, you enroll only those participants)
data_filter = read.table(file="./data_filter.csv",header=T,sep=",",dec=".")

## N fixed as the size of data_filter, given gamma, find power
## For simulations in full sample, we draw a random sample size of N of data_filter
## to compare

powercomp_exp_enroll <- function(NSIMS, gam, data_full, data_filter, boot = F){
  set.seed(12345)
  pb <- txtProgressBar(0, NSIMS,style=3)
    pvals.pred.p<-rep(NA,NSIMS)
    pvals.nopred.p<-rep(NA,NSIMS)
    pvals.lin.p <- rep(NA,NSIMS)
    
    pvals.pred.ti <- rep(NA,NSIMS)
    pvals.nopred.ti<-rep(NA,NSIMS)
    pvals.lin.ti <- rep(NA,NSIMS)
    
    for(i in 1:NSIMS){
      #split ADNI data into 2
      data_matched <- data_full[sample(nrow(data_full),nrow(data_filter),replace=F),]
      tx1 <- data_matched[sample(nrow(data_matched), nrow(data_matched)/2,replace=F),]
      tx2 <- suppressMessages(anti_join(data_matched, tx1))
      n.tx <- nrow(data_matched)/2
      
      #treatment arm 1
      ftx.1 <- tx1 #"bootstrap" a future tx sample of MCI pts, tx arm 1
      ftx.1 %<>% mutate(A = rep(1, n.tx),
                        Y.p = diff_adas_22_2 + A*(gam + rnorm(n.tx, 0, gam*0.1)),
                        Y.t1 = diff_adas_22_2)
      
      #treatment arm 2
      ftx.2 <- tx2 #a future tx sample of MCI pts, tx arm 2
      ftx.2 %<>% mutate(A = rep(0, n.tx),
                        Y.p = diff_adas_22_2 + A*(gam + rnorm(n.tx, 0, gam*0.1)),
                        Y.t1 = diff_adas_22_2)
      
      dat.final = rbind(ftx.1, ftx.2)

      pred.tx1 <- data_filter[sample(nrow(data_filter), nrow(data_filter)/2,replace=F),]
      pred.tx2 <- suppressMessages(anti_join(data_filter, pred.tx1))
      pred.n.tx <- nrow(data_filter)/2
      
      #treatment arm 1
      pred.ftx.1 <- pred.tx1 #"bootstrap" a future tx sample of MCI pts, tx arm 1
      pred.ftx.1 %<>% mutate(A = rep(1, pred.n.tx),
                        Y.p = diff_adas_22_2 + A*(gam + rnorm(pred.n.tx, 0, gam*0.1)),
                        Y.t1 = diff_adas_22_2)
      
      #treatment arm 2
      pred.ftx.2 <- pred.tx2 #a future tx sample of MCI pts, tx arm 2
      pred.ftx.2 %<>% mutate(A = rep(0, pred.n.tx),
                        Y.p = diff_adas_22_2 + A*(gam + rnorm(pred.n.tx, 0, gam*0.1)),
                        Y.t1 = diff_adas_22_2)
      
      pred.dat.final = rbind(pred.ftx.1, pred.ftx.2)
      
      #collecting pvals
      pvals.nopred.p[i] <- summary(lm(Y.p ~ A + age, data = dat.final))$coef[2,4]
      pvals.nopred.ti[i] <- summary(lm(Y.t1 ~ A + age, data = dat.final))$coef[2,4]
      
      pvals.pred.p[i] <- summary(lm(Y.p ~ A + age, data = pred.dat.final))$coef[2,4]
      pvals.pred.ti[i] <- summary(lm(Y.t1 ~ A + age, data = pred.dat.final))$coef[2,4]
    }
    
    pow.pred <- sum(pvals.pred.p<0.05)/NSIMS 
    pow.nopred <- sum(pvals.nopred.p<0.05)/NSIMS 
    
    ti.pred <- sum(pvals.pred.ti<0.05)/NSIMS 
    ti.nopred <- sum(pvals.nopred.ti<0.05)/NSIMS 
    
  close(pb)
  
  print(paste("pow.nopred", pow.nopred))
  print(paste("pow.pred", pow.pred))
  print(paste("ti.nopred", ti.nopred))
  print(paste("ti.pred", ti.pred))
  
  return(list(pow.nopred, pow.pred, ti.nopred, ti.pred))

}

iterate_effectsizes_power <- function(ES_list, NSIMS, data_full, data_filter) {
  # Create an empty list to store the results
  len <- length(ES_list)
  powers.full <- rep(NA,len)
  powers.enroll <- rep(NA,len)
  type1s.full <- rep(NA,len)
  type1s.enroll <- rep(NA,len)
  
  
  # Iterate over the gam_list
  for (i in 1:len) {
    es_value <- ES_list[i]
    # Call the compute function with the current gam_value
    result <- powercomp_exp_enroll(NSIMS, 9*es_value, data_full, data_filter)
    
    # Append the result to the results list
    powers.full[i] <- result[[1]]
    powers.enroll[i] <- result[[2]]
    type1s.full[i] <- result[[3]]
    type1s.enroll[i] <- result[[4]]
  }
  
  pow.p <- c(powers.full, powers.enroll)
  type1 <- c(type1s.full, type1s.enroll)

  method <- c(rep("Without Enriched Enrollment",len), rep("Enriched Enrollment",len))
  
  d.plot<-data.frame( es.plot = rep(ES_list,2), pow.p, method)
  

  
  return(d.plot)
  
}


result_pow <- iterate_effectsizes_power(ES_list = seq(0.3, 0.5, by = 0.025), 1000, data_full, data_filter)

mapping <- c("Without Enriched Enrollment", "Enriched Enrollment")

# Replace factor values with the mapping
#result_pow$method <- factor(result_pow$method, levels = c("Classical Enrollment", "Informed Enrollment"), labels = mapping)


plot.power <- ggplot(aes(x=es.plot, y=pow.p*100, color=method), data=result_pow)+
  geom_line(linewidth=1)+
  geom_point(size = 2) +
  coord_cartesian(ylim=c(0,100))+
  ggtitle("Approach I - During Enrollment")+
  ylab("Power (%)")+
  xlab("Effect Size")+
  scale_x_continuous(breaks = seq(0.2, 0.5, by = 0.05)) +
  theme_classic()+
  theme(text = element_text(size = 14))+
  theme(legend.position = "bottom",
        legend.justification = c(1, 0), # Adjust the x and y justification
        legend.box.just = "left")+
  theme(plot.title = element_text(hjust = 0.5)) +
  theme( axis.title.y = element_text(vjust = 0.5))

plot.power

