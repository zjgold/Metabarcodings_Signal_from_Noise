library(tidyverse)
library(gridExtra)
library(ggplot2)
library(cowplot)
library(mc2d)


length.zero <- function(dat){
  return(length(dat[dat==0]))
}


### START BY JUST LOOKING AT THE NUMBER OF DNA COPIES THAT MAKE IT INTO THE TUBE.
lambda <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,c(seq(0,12,by=0.1) %>% exp()))
V      <- 1

N_samp <- 50000
W <- matrix(0,N_samp,length(lambda))

## poisson version
for(i in 1:length(lambda)){
  W[,i] <- rpois(1000,lambda[i]*V)
}

P.pois <- W / rowSums(W)
p.EV <-  lambda * V / sum(lambda *V)
P.pois.sum <- data.frame(Mean = colMeans(P.pois), Var=  apply(P.pois,2,var))
P.pois.sum$CV <- sqrt(P.pois.sum$Var) / P.pois.sum$Mean


## Neg Binom Version
phi_0_all <- c(seq(0.1,1.0,by=0.1),seq(2,5),seq(10,100,by=10),10000)#(dispersion when expected number of reads is very small)
phi_1_all <- c(0,0.01,0.05,seq(0.1,1,by=0.1),2,5,10) #(positive if dispersion declines with mean)

phi_all <- expand.grid(phi_0=phi_0_all,phi_1=phi_1_all)

#exp(log(phi_0) + phi_1*log(lambda*V) )


W.all <- NULL
for(j in 1:nrow(phi_all)){
  phi_0 <- phi_all$phi_0[j]
  phi_1 <- phi_all$phi_1[j]
  
  for(i in 1:length(lambda)){
    W[,i] <- rnbinom(N_samp,mu = lambda[i]*V,size= exp(log(phi_0) + phi_1*log(lambda[i]*V) ))
  }
  
  W.summary <- t(apply(W,2,quantile,probs=c(0.025,0.05,0.25,0.5,0.75,0.95,0.975)))
  colnames(W.summary) <- c("q.0.025","q.0.05","q.0.25","q.0.50","q.0.75","q.0.95","q.0.975")
  W.summary <- as.data.frame(W.summary)
  W.summary$Mean <- colMeans(W)
  W.summary$Var <- apply(W,2,var)
  W.summary$CV <- sqrt(W.summary$Var) /W.summary$Mean
  W.summary$true.mean <- lambda * V
  W.summary$true.phi <- exp(log(phi_0) + phi_1*log(lambda*V) )
  W.summary$phi_0 <- phi_0
  W.summary$phi_1 <- phi_1
  W.summary$N_zero <- apply(W,2,length.zero)
  W.summary$N_samp <- N_samp

  W.all <- bind_rows(W.all,W.summary)
}

W.all <- W.all %>% mutate(p_zero = N_zero / N_samp)

# subset and plot
  
  p_zero_ylim <- c(0,1)
  BREAKS <- c(0.1,1,10,100,1000,10000,100000)
  
  p_classic <-  ggplot( W.all %>% filter(phi_0 %in% c(1e-01,1,5,10,100,max(phi_0_all)),
                           phi_1 == 0)) +
          geom_line(aes(x=true.mean,y=Mean,color=as.factor(phi_0))) +
          geom_ribbon(aes(x=true.mean,ymax=q.0.95,ymin=q.0.05,fill=as.factor(phi_0),color=as.factor(phi_0)),alpha=0.3) +
          #geom_line(aes(x=true.mean,y=q.0.05,color=as.factor(phi_0))) +
          geom_abline(intercept=0,slope=1,linetype="dashed") +
          scale_color_discrete(expression("Overdispersion ("*phi[0]*")")) +
          scale_fill_discrete(expression("Overdispersion ("*phi[0]*")")) +
          scale_y_continuous("Observation (90% CI)") +
          scale_x_continuous(expression("Mean ("*mu*")")) +
          ggtitle(label=expression("Classic NB ("*mu*", "*phi[0]*","*phi[1]*"=0)")) +
          theme_bw() +
          theme(legend.position = c(0.15,0.8))
  
  p_classic_CV <-  ggplot( W.all %>% filter(phi_0 %in% c(1e-01,0.2,0.5,1,5,10,100,max(phi_0_all)),
                                         phi_1 == 0,
                                         true.mean>=1)) +
    geom_line(aes(x=true.mean,y=CV,color=as.factor(phi_0))) +
    scale_x_continuous(expression("Mean ("*mu*")"),trans="log",
                       breaks=BREAKS) +
    scale_y_continuous("CV") +
    scale_color_discrete(expression("Overdispersion ("*phi[0]*")")) +
    ggtitle(label=expression("CV; Classic NB ("*mu*", "*phi[0]*","*phi[1]*"=0)")) +
    theme_bw() +
    theme(legend.position = c(0.1,0.15))
  
  p_classic_p_zero <-  ggplot( W.all %>% filter(phi_0 %in% c(1e-01,0.2,0.5,1,5,10,100,max(phi_0_all)),
                                            phi_1 == 0,
                                            true.mean>=0.1)) +
    geom_line(aes(x=true.mean,y=p_zero,color=as.factor(phi_0))) +
    scale_x_continuous(expression("Mean ("*mu*")"),trans="log",
                       breaks=BREAKS) +
    scale_y_continuous("P(W=0)",limits = p_zero_ylim) +
    scale_color_discrete(expression("Overdispersion ("*phi[0]*")")) +
    ggtitle(label=expression("P(W=0); Classic NB ("*mu*", "*phi[0]*","*phi[1]*"=0)")) +
    theme_bw() +
    theme(legend.position = c(0.8,0.85))
  
  
  ###
  p_phi_0_is_0.1 <-  ggplot( W.all %>% filter(phi_0 %in% c(0.1),
                                            phi_1 %in% c(0,0.1, 0.5, 1,5,10))) +
    geom_line(aes(x=true.mean,y=Mean,color=as.factor(phi_1))) +
    geom_ribbon(aes(x=true.mean,ymax=q.0.95,ymin=q.0.05,fill=as.factor(phi_1),color=as.factor(phi_1)),alpha=0.3) +
    geom_abline(intercept=0,slope=1,linetype="dashed") +
    scale_color_discrete(expression("Overdispersion slope ("*phi[1]*")")) +
    scale_fill_discrete(expression("Overdispersion slope ("*phi[1]*")")) +
    scale_y_continuous("Observation (90% CI)") +
    scale_x_continuous(expression("Mean ("*mu*")")) +
    ggtitle(label=expression("Mod NB ("*mu*","*phi[0]*"=0.1,"*phi[1]*")")) +
    theme_bw() +
    theme(legend.position = c(0.15,0.8))
  
  p_phi_0_is_0.1_CV <- ggplot( W.all %>% filter(phi_0 %in% c(0.1),
                                              phi_1 %in% c(0,0.1, 0.5, 1,5,10),
                                              true.mean>=1)) +
    geom_line(aes(x=true.mean,y=CV,color=as.factor(phi_1))) +
    scale_x_continuous(expression("Mean ("*mu*")"),trans="log",
                       breaks=BREAKS) +
    scale_y_continuous("CV",trans="sqrt") +
    scale_color_discrete(expression("Overdispersion slope ("*phi[1]*")")) +
    ggtitle(label=expression("CV; Mod NB ("*mu*","*phi[0]*"=0.1,"*phi[1]*")")) +
    theme_bw() +
    theme(legend.position = c(0.1,0.15))
  
  p_phi_0_is_0.1_p_zero <- ggplot( W.all %>% filter(phi_0 %in% c(0.1),
                                                    phi_1 %in% c(0,0.1, 0.5, 1,5,10),
                                                    true.mean>=0.1)) +
    geom_line(aes(x=true.mean,y=p_zero,color=as.factor(phi_1))) +
    scale_x_continuous(expression("Mean ("*mu*")"),trans="log",
                       breaks=BREAKS) +
    scale_y_continuous("P(W=0)",limits = p_zero_ylim) +
    scale_color_discrete(expression("Overdispersion slope ("*phi[1]*")")) +
    ggtitle(label=expression("P(W=0); ("*mu*","*phi[0]*"=0.1,"*phi[1]*")")) +
    theme_bw() +
    theme(legend.position = c(0.8,0.85))
  
  ####
  p_phi_0_is_1 <-  ggplot( W.all %>% filter(phi_0 %in% c(1),
                                            phi_1 %in% c(0,0.1, 0.5, 1,5,10))) +
    geom_line(aes(x=true.mean,y=Mean,color=as.factor(phi_1))) +
    geom_ribbon(aes(x=true.mean,ymax=q.0.95,ymin=q.0.05,fill=as.factor(phi_1),color=as.factor(phi_1)),alpha=0.3) +
    geom_abline(intercept=0,slope=1,linetype="dashed") +
    scale_color_discrete(expression("Overdispersion slope ("*phi[1]*")")) +
    scale_fill_discrete(expression("Overdispersion slope ("*phi[1]*")")) +
    scale_y_continuous("Observation (90% CI)") +
    scale_x_continuous(expression("Mean ("*mu*")")) +
    ggtitle(label=expression("Mod NB ("*mu*","*phi[0]*"=1,"*phi[1]*")")) +
    theme_bw() +
    theme(legend.position = c(0.15,0.8))
  
  p_phi_0_is_1_CV <- ggplot( W.all %>% filter(phi_0 %in% c(1),
                                              phi_1 %in% c(0,0.1, 0.5, 1,5,10),
                                              true.mean>=1)) +
    geom_line(aes(x=true.mean,y=CV,color=as.factor(phi_1))) +
    scale_x_continuous(expression("Mean ("*mu*")"),trans="log",
                       breaks=BREAKS) +
    scale_y_continuous("CV",trans="sqrt") +
    scale_color_discrete(expression("Overdispersion slope ("*phi[1]*")")) +
    ggtitle(label=expression("CV; Mod NB ("*mu*","*phi[0]*"=1,"*phi[1]*")")) +
    theme_bw() +
    theme(legend.position = c(0.1,0.15))
  
  p_phi_0_is_1_p_zero <- ggplot( W.all %>% filter(phi_0 %in% c(1),
                                                    phi_1 %in% c(0,0.1, 0.5, 1,5,10),
                                                    true.mean>=0.1)) +
    geom_line(aes(x=true.mean,y=p_zero,color=as.factor(phi_1))) +
    scale_x_continuous(expression("Mean ("*mu*")"),trans="log",
                       breaks=BREAKS) +
    scale_y_continuous("P(W=0)",limits = p_zero_ylim) +
    scale_color_discrete(expression("Overdispersion slope ("*phi[1]*")")) +
    ggtitle(label=expression("P(W=0); Mod NB ("*mu*","*phi[0]*"=1,"*phi[1]*")")) +
    theme_bw() +
    theme(legend.position = c(0.8,0.85))
  
  
  ##
  
  
  p_phi_0_is_10 <-  ggplot( W.all %>% filter(phi_0 %in% c(10),
                                         phi_1 %in% c(0,0.1, 0.5, 1,5,10))) +
          geom_line(aes(x=true.mean,y=Mean,color=as.factor(phi_1))) +
          geom_ribbon(aes(x=true.mean,ymax=q.0.95,ymin=q.0.05,fill=as.factor(phi_1),color=as.factor(phi_1)),alpha=0.3) +
          geom_abline(intercept=0,slope=1,linetype="dashed") +
          scale_color_discrete(expression("Overdispersion slope ("*phi[1]*")")) +
          scale_fill_discrete(expression("Overdispersion slope ("*phi[1]*")")) +
          scale_y_continuous("Observation (90% CI)") +
          scale_x_continuous(expression("Mean ("*mu*")")) +
          ggtitle(label=expression("Mod NB ("*mu*","*phi[0]*"=10,"*phi[1]*")")) +
          theme_bw() +
          theme(legend.position = c(0.15,0.8))
  
  p_phi_0_is_10_CV <- ggplot( W.all %>% filter(phi_0 %in% c(10),
                                        phi_1 %in% c(0,0.1, 0.5, 1,5,10),
                                        true.mean>=1)) +
            geom_line(aes(x=true.mean,y=CV,color=as.factor(phi_1))) +
            scale_x_continuous(expression("Mean ("*mu*")"),trans="log",
                                breaks=BREAKS) +
            scale_y_continuous("CV",trans="sqrt") +
            scale_color_discrete(expression("Overdispersion slope ("*phi[1]*")")) +
            ggtitle(label=expression("CV; Mod NB ("*mu*","*phi[0]*"=10,"*phi[1]*")")) +
            theme_bw() +
            theme(legend.position = c(0.1,0.15))
  
  
  p_phi_0_is_10_p_zero <- ggplot( W.all %>% filter(phi_0 %in% c(10),
                                                  phi_1 %in% c(0,0.1, 0.5, 1,5,10),
                                                  true.mean>=0.1)) +
    geom_line(aes(x=true.mean,y=p_zero,color=as.factor(phi_1))) +
    scale_x_continuous(expression("Mean ("*mu*")"),trans="log",
                       breaks=BREAKS) +
    scale_y_continuous("P(W=0)",limits = p_zero_ylim) +
    scale_color_discrete(expression("Overdispersion slope ("*phi[1]*")")) +
    ggtitle(label=expression("P(W=0); Mod NB ("*mu*","*phi[0]*"=10,"*phi[1]*")")) +
    theme_bw() +
    theme(legend.position = c(0.8,0.85))
  
  ##
  
  G <- list(p_classic + theme(legend.position = "none"),
            p_classic_CV + theme(legend.position = "none"),
            p_classic_p_zero + theme(legend.position = "none"),
            get_legend(p_classic+theme(legend.position=c(0.5,0.5))),
            p_phi_0_is_0.1 + theme(legend.position = "none"),
            p_phi_0_is_0.1_CV + theme(legend.position = "none"),
            p_phi_0_is_0.1_p_zero + theme(legend.position = "none"),
            get_legend(p_phi_0_is_0.1+theme(legend.position=c(0.5,0.5))),
            p_phi_0_is_1+ theme(legend.position = "none"),
            p_phi_0_is_1_CV+ theme(legend.position = "none"),
            p_phi_0_is_1_p_zero+ theme(legend.position = "none"),
            get_legend(p_phi_0_is_1+theme(legend.position=c(0.5,0.5))),
            p_phi_0_is_10+ theme(legend.position = "none"),
            p_phi_0_is_10_CV+ theme(legend.position = "none"),
            p_phi_0_is_10_p_zero+ theme(legend.position = "none"),
            get_legend(p_phi_0_is_10+theme(legend.position=c(0.5,0.5))))
  
  g.plot <- grid.arrange(grobs=G, nrow=4, ncol=4,widths=c(1,1,1,0.5))
               
  #####################################################################               
  #######################################################################
  #######################################################################
  #######################################################################
  #######################################################################
  #######################################################################
  #######################################################################
  