library(tidyverse)
library(gridExtra)
library(ggplot2)
library(cowplot)
library(mc2d)


length.zero <- function(dat){
  return(length(dat[dat==0]))
}


#######################################################################
## OK. This is the full simulation of sampling molecules, pcr, then sequencing.
#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################
# Double Poisson then multinomial version.

### START BY JUST LOOKING AT THE NUMBER OF DNA COPIES THAT MAKE IT INTO THE TUBE.
N_species <- 20

set.seed(151)
a_amp_vec <- c(1e6,
               1000,1000.001,999.999,
               100,100.001,99.999,
               20, 19.999,20.001,
               10,9.999,10.001,
               5.0001,5,4.9999)

a_amp <- matrix(0,length(a_amp_vec),N_species)
for(i in 1:length(a_amp_vec)){
  a_amp[i,]     <- rbeta(N_species,a_amp_vec[i]*0.7,a_amp_vec[i]*0.3)
}

N_pcr1     <- 35
N_pcr2     <- 10
subsamp <- 0.20 # proportion of PCR 1 product used in PCR 2

# Assume species are equally abundant
lambda <- c(0.5,0.75,0.9,1,2,3,5,10,20,30,50,100,200,300,500,1000,3000,10000)
# lambda here is the total concentration

V      <- 1
N_samp <- 50000
N_multinom_range <- c(60000,140000)
N_multinom_mid <- mean(N_multinom_range)

sim.all <- NULL
sim.prob.dat <- NULL
for(k in 1:nrow(a_amp)){ # loop over amplification variability scenarios.
  
  # First.  Sample Molecules in the tube.
  # W2 columns = species, rows = realizations.
  W2 <- list()
  for(i in 1:length(lambda)){
    W2[[i]] <- matrix(rpois(N_species*N_samp,lambda[i]),N_samp,N_species,byrow = T)
  }
  
  # Then Amplify them stochastically with PCR
  X <- list()
  for(j in 1:length(lambda)){
    X[[j]] <- W2[[j]]
    for(i in 1:N_species){
      X[[j]][,i] <- rpois(N_samp,X[[j]][,i]*(1+a_amp[k,i])^N_pcr1)
      #subsample the reads using a binomial
      X[[j]][,i] <- rbinom(N_samp,prob=subsamp,size=X[[j]][,i])
      #Finish PCR with the second set of PCR cyles.
      X[[j]][,i] <- rpois(N_samp,X[[j]][,i]*(1.9)^N_pcr2)
    }
  }
  
  # Calculate the proportions in each replicate
  X_prop <- list()
  for(j in 1:length(lambda)){
    SUM <- rowSums(X[[j]])
    SUM[SUM ==0] <- 1
    X_prop[[j]] <- X[[j]] / SUM 
    X_prop[[j]] <- X_prop[[j]][which(rowSums(X_prop[[j]])>0),]
  }
  
  # Then sample them with a multinomial, proportionally to their abundance after PCR.
  Y <- list()  
  for(j in 1:length(lambda)){
    Y[[j]] <- rmultinomial(nrow(X_prop[[j]]),
                           size=round(runif(nrow(X_prop[[j]]),N_multinom_range[1],N_multinom_range[2])),
                           prob=X_prop[[j]])
      
      #rep(N_multinom_range,nrow(X_prop[[j]])),prob=X_prop[[j]])
  }
  
  #### HERE IS WHERE YOU START SAVING THINGS THAT CAN BE PLOTTED.
  
  PROBS <- c(0.025,0.05,0.25,0.5,0.75,0.95,0.975)
  
  W.all <- NULL
  Y.all <- NULL
  Y.prob.all <- NULL
  for(i in 1:length(lambda)){
    
    W.summary <- t(apply(W2[[i]],2,quantile,probs=PROBS))
    colnames(W.summary) <- c("q.0.025","q.0.05","q.0.25","q.0.50","q.0.75","q.0.95","q.0.975")
    W.summary <- as.data.frame(W.summary)
    W.summary$Mean <- colMeans(W2[[i]])
    W.summary$Var <- apply(W2[[i]],2,var)
    W.summary$CV <- sqrt(W.summary$Var) /W.summary$Mean
    W.summary$lambda <- lambda[i]
    W.summary$true.mean <- lambda[i] 
    W.summary$sp <- 1:N_species
    W.summary$alpha <- a_amp[k,]
    W.summary$N_zero <- apply(W2[[i]],2,length.zero)
    W.summary$N_samp <- N_samp
    
    W.all <- bind_rows(W.all,W.summary)
    
    # account for sampling depth in the Y observations by standardizing to the mid-point
    # of the observed read sampling depth
    Y_mod <- Y[[i]] * (N_multinom_mid / rowSums(Y[[i]]))
    
    Y.summary <- t(apply(Y_mod,2,quantile,probs=PROBS))
    colnames(Y.summary) <- c("q.0.025","q.0.05","q.0.25","q.0.50","q.0.75","q.0.95","q.0.975")
    Y.summary <- as.data.frame(Y.summary)
    Y.summary$Mean <- colMeans(Y_mod)
    Y.summary$Var <- apply(Y_mod,2,var)
    Y.summary$CV <- sqrt(Y.summary$Var) /Y.summary$Mean
    Y.summary$lambda <- lambda[i] 
    Y.summary$true.mean <- lambda[i] 
    Y.summary$sp <- 1:N_species
    Y.summary$alpha <- a_amp[k,]
    Y.summary$N_zero <- apply(Y[[i]],2,length.zero)
    Y.summary$N_samp <- nrow(Y[[i]])
    
    Y.prob.summary <- t(apply(Y[[i]] / rowSums(Y[[i]]),2,quantile,probs=PROBS))
    colnames(Y.prob.summary) <- c("q.0.025","q.0.05","q.0.25","q.0.50","q.0.75","q.0.95","q.0.975")
    Y.prob.summary <- as.data.frame(Y.prob.summary)
    Y.prob.summary$lambda <- lambda[i] 
    Y.prob.summary$true.mean <- lambda[i] 
    Y.prob.summary$true.mean.prob <- 1/ N_species
    Y.prob.summary$mean.prob <- colMeans(Y[[i]] / rowSums(Y[[i]]))
    Y.prob.summary$sd.prob   <- apply(Y[[i]] / rowSums(Y[[i]]),2,sd)
    Y.prob.summary$sp <- 1:N_species
    Y.prob.summary$alpha <- a_amp[k,]
    
    Y.prob.all <- bind_rows(Y.prob.all,Y.prob.summary)
    Y.all <- bind_rows(Y.all,Y.summary)
  }
  
  Y.all$type <- "Y"
  W.all$type <- "W"
  
  Y.prob.all$a_scenario <- a_amp_vec[k]
  sim.prob.dat <- bind_rows(sim.prob.dat,Y.prob.all)
  
  sim.dat <- bind_rows(W.all,Y.all)  
  sim.dat <- sim.dat %>% mutate(p_zero= N_zero / N_samp)
  sim.dat$a_scenario <- a_amp_vec[k]
  
  sim.all <- bind_rows(sim.all,sim.dat)
  
  print(k)
} # end loop over alpha scenarios


# save these things to file
NOM <- paste0("./signal2noise/code/simulate_overdispersion/Even simulation; N_species= ",N_species," pi = ",subsamp*100," N_pcr=",N_pcr1,".RData")

Output <- list(sim.all = sim.all, # These are the summaries in count space
               sim.prob.dat = sim.prob.dat) # These are the summaries in proportion space
save(Output,file = NOM)


# sim.dat.sum <- sim.dat %>% group_by(sp,alpha,type) %>% 
#                   summarise(mean_p_zero=mean(p_zero),sd_p_zero=sd(p_zero),
#                             grand_mean = mean(Mean),)
# 
##### PLOTS

 BREAKS <- c(0.1,1,10,100,1000,10000)
# # Mean v. Spread for W vs. Y
# p_W<-  ggplot(sim.all %>% filter(sp==1,type=="W")) +
#   geom_line(aes(x=true.mean,y=Mean,color=type)) +
#   geom_ribbon(aes(x=true.mean,ymax=q.0.95,ymin=q.0.05,fill=type,color=type),alpha=0.3) +
#   #geom_line(aes(x=true.mean,y=q.0.05,color=as.factor(phi_0))) +
#   geom_abline(intercept=0,slope=1,linetype="dashed") +
#   scale_color_discrete(expression("Overdispersion ("*phi[0]*")")) +
#   scale_fill_discrete(expression("Overdispersion ("*phi[0]*")")) +
#   scale_y_continuous("Observation (90% CI)") +
#   scale_x_continuous(expression("Mean ("*lambda*")"),trans="log") +
#   ggtitle(label=expression("Classic NB ("*lambda*", "*phi[0]*","*phi[1]*"=0)")) +
#   theme_bw() +
#   theme(legend.position = c(0.15,0.8))

# P_ZERO PLOTS

BREAKS.ZERO <- c(0,0.10,0.20,0.40,0.60,0.80,1.0)
p_Y_zero <- ggplot(sim.all%>% filter(type=="Y")) +
  geom_line(aes(x=true.mean,y=p_zero,color=alpha,group=sp)) +
  scale_x_continuous(expression("Mean("*lambda*")"),trans="log",breaks=BREAKS) +
  scale_y_continuous("P(Y=0)",breaks=BREAKS.ZERO,limits=c(0,1),expand=c(0.01,0.01)) +
  scale_color_viridis_c(option="H") +
  ggtitle("P(0)") +
  facet_wrap(~a_scenario) +
  theme_bw()

BREAKS.ZERO <- c(0,0.01,0.05,0.10,0.20,0.40,0.60,0.80,1.0)
p_W_zero <- ggplot(sim.all%>% filter(type=="W")) +
  geom_line(aes(x=true.mean,y=p_zero,color=alpha,group=sp)) +
  scale_x_continuous(expression("Mean("*lambda*")"),trans="log",breaks=BREAKS) +
  scale_y_continuous("P(W=0)",trans="sqrt",breaks=BREAKS.ZERO,limits=c(0,1)) +
  scale_color_viridis_c(option="H") +
  ggtitle("P(0)") +
  facet_wrap(~a_scenario) +
  theme_bw()


sim.all <- sim.all %>% group_by(true.mean,sp,a_scenario) %>% 
              mutate(max.zero = max(p_zero) ) %>%
              mutate(p_zero_amp = ifelse(type=="W",max.zero-p_zero,NA)) %>%
              mutate(sp.lab = paste0("sp ", sp))

BREAKS.ZERO <- c(0,0.05,0.10,0.20,0.40,0.60,0.80,1.0)
p_zero_both_some <- ggplot(sim.all %>% filter(a_scenario ==5, 
                        sp %in% c(1,2,3,7,11,12,13,14,15),type=="W")) +
  geom_area(aes(x=true.mean,y=max.zero),fill="red",alpha=0.5) +
  geom_area(aes(x=true.mean,y=p_zero_amp),fill="blue",alpha=0.5) +
    facet_wrap(~sp)+
    scale_x_continuous(expression("Mean("*lambda*")"),trans="log",breaks=BREAKS) +
    scale_y_continuous("P(0)",breaks=BREAKS.ZERO,limits=c(0,1)) +
    theme_bw()

BREAKS.ZERO <- c(0,0.10,0.20,0.40,0.60,0.80,1.0)
p_zero_both_all <- ggplot(sim.all %>% filter(a_scenario ==5, 
                                              sp %in% c(1:30),type=="W")) +
  geom_area(aes(x=true.mean,y=max.zero),fill="red",alpha=0.5) +
  geom_area(aes(x=true.mean,y=p_zero_amp),fill="blue",alpha=0.5) +
  facet_wrap(~sp)+
  scale_x_continuous(expression("Mean("*lambda*")"),trans="log",breaks=BREAKS) +
  scale_y_continuous("P(0)",breaks=BREAKS.ZERO,limits=c(0,1)) +
  theme_bw()


# CV PLOTS
BREAKS.CV <- c(0.01,0.03,0.1,0.3,1,3,10,30,100)
p_CV <- ggplot(sim.all%>% filter(type=="Y")) +
  geom_line(aes(x=true.mean,y=CV,color=alpha,group=sp)) +
  scale_x_continuous(expression("Mean("*lambda*")"),trans="log",breaks=BREAKS) +
  scale_y_continuous("CV",trans="log",breaks=BREAKS.CV) +
  scale_color_viridis_c("a",option="H") +
  ggtitle("CV") +
  facet_wrap(~a_scenario) +
  theme_bw()


sim.all <- sim.all %>% mutate(x.loc = (sp - (max(sim.all$sp) / 2) *0.1))

p_var_by_a_scen   <- ggplot(sim.all %>% filter(type=="Y",
                          a_scenario %in% c(5,10,50,100,1000,1e+06),
                          true.mean %in% c(1,10,100,1000,10000))) +
    geom_linerange(aes(ymin=q.0.05,ymax=q.0.95,x=true.mean,group=sp,color=alpha),
                   position = position_dodge(width=1)) +
    scale_y_continuous("Reads") +
    scale_x_continuous(expression("Mean("*lambda*")"),trans="log",breaks=BREAKS) +
    scale_color_viridis_c("a",option="H") +
    facet_wrap(~a_scenario) +  
    ggtitle("Range observed by a_scenario") +
    theme_bw()

p_var_prob_by_a_scen1   <- ggplot(sim.prob.dat %>% filter(
                                               a_scenario %in% c(5,10,50,100,1000,1e+06),
                                               true.mean %in% c(1,10,100,1000,10000))) +
  geom_linerange(aes(ymin=q.0.05,ymax=q.0.95,x=true.mean,group=sp,color=alpha),
                 position = position_dodge(width=1)) +
  scale_y_continuous("Proportion of Read") +
  scale_x_continuous(expression("Mean("*lambda*")"),trans="log",breaks=BREAKS) +
  scale_color_viridis_c("a",option="H") +
  facet_wrap(~a_scenario) +  
  ggtitle("Range observed by a_scenario") +
  theme_bw()

BREAKS.ZERO <- c(0,0.001,0.01,0.05,0.10,0.20,0.40,0.60,0.80,1.0)
p_var_prob_by_a_scen2   <- ggplot(sim.prob.dat %>% 
                            filter(a_scenario %in% c(5,10,50,100,1000,1e+06),
                                  true.mean %in% c(1,10,100,1000,10000))) +
  geom_linerange(aes(ymin=q.0.05,ymax=q.0.95,x=true.mean,group=sp,color=alpha),
                 position = position_dodge(width=1)) +
  scale_y_continuous("Proportion of Reads",trans="sqrt",breaks=BREAKS.ZERO) +
  scale_x_continuous(expression("Mean("*lambda*")"),trans="log",breaks=BREAKS) +
  scale_color_viridis_c("a",option="H") +
  facet_wrap(~a_scenario) +  
  ggtitle("Range observed by a_scenario") +
  theme_bw()

BREAKS.PROBS <- c(0.0001,0.0003,0.001,0.003,0.01,0.03,0.1,0.3,0.5,0.8,1)
p_SD_v_Mean <- ggplot(sim.prob.dat %>% filter( )  ) +
  geom_point(aes(x=mean.prob,y=sd.prob,color=true.mean)) +
  scale_x_continuous(expression("Mean Proportion"),trans="log",breaks=BREAKS.PROBS,limits=c(0.00001,0.5)) +
  scale_y_continuous("SD(Proportion)",trans="log",breaks=BREAKS.PROBS,limits=c(0.00001,NA)) +
  scale_color_viridis_c(expression(lambda),trans="log",option="H",breaks=BREAKS) +
  ggtitle("SD v Mean") +
  # facet_wrap(~true.mean) +
  theme_bw()

p_SD_v_Mean + facet_wrap(~a_scenario)

# Simple marginals 
# 
# lambda.big <- Y[[length(lambda)]][,1:6] %>% as.data.frame() %>% 
#                   pivot_longer(.,cols=1:6,names_to="sp",values_to="counts") %>%
#                   mutate(lambda = lambda[length(lambda)])
# 
# lambda.mid <- Y[[round(length(lambda)/2,0)]][,1:6] %>% as.data.frame() %>% 
#   pivot_longer(.,cols=1:6,names_to="sp",values_to="counts") %>%
#   mutate(lambda = lambda[round(length(lambda)/2,0)])
# 
# lambda.sm <- Y[[1]][,1:6] %>% as.data.frame() %>% 
#   pivot_longer(.,cols=1:6,names_to="sp",values_to="counts") %>%
#   mutate(lambda = lambda[1])
# 
# p_lambda_big <- ggplot(lambda.big) +
#           geom_histogram(aes(counts)) +
#           facet_wrap(~sp) +
#           theme_bw()
# 
# 
# 
# ggplot(Y[[length(lambda)]])

pdf(file=paste0("./signal2noise/code/simulate_overdispersion/Simulate Poisson PCR N_sp=",N_species,"; pi= ",subsamp*100,".pdf"),onefile=T)

    print(p_W_zero)
    print(p_Y_zero)
    print(p_zero_both_all)
    print(p_zero_both_some)
    
    print(p_CV)
    print(p_var_prob_by_a_scen1)
    print(p_var_prob_by_a_scen2)

    X1 <- c(0.0001,0.0003,0.001,0.003,0.01,0.03,0.1,0.3,0.5)
    Y1 <- sqrt(X1)
    xy.dat <- data.frame(x,y)
    print(p_SD_v_Mean) +
      geom_abline(intercept = 0,slope=1,linetype="dashed",color="red")
     #geom_line(data=xy.dat,aes(y=y,x=x),linetype="dashed",color="blue") 
    print(p_SD_v_Mean + facet_wrap(~a_scenario))
    
dev.off()