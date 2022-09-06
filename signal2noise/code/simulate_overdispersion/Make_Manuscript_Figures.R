#### This file is for making pretty figures from the simulations for Zack's overdispersion paper
library(tidyverse)
library(grid)
library(ggplot2)
library(gridExtra)
library(cowplot)
library(here)
library(mc2d)

# Select a simulation read in:
here()
load(here("signal2noise","code","simulate_overdispersion","Even simulation; N_species= 50 pi = 20 N_pcr=35.RData"))
N_PCR <- 35 # needed for labeling below

names(Output)
#sim.all is the summary of the counts.
#sim.prob.dat is the summaries of the simulations in proportions.
sim.all <- Output$sim.all
sim.prob.dat <- Output$sim.prob.dat
BREAKS <- c(0.1,1,10,100,1000,10000,100000)

BREAKS.ZERO <- c(0,0.10,0.20,0.40,0.60,0.80,1.0)
n_scen <- c(5,10,100,1e+06)
p_Y_zero <- list()
p_a_inset <- list()
for(i in 1:length(n_scen)){ 
  dat <- sim.all %>% 
    filter(type=="Y") %>% 
    filter(a_scenario == n_scen[i])
  
  a_scen <- mean(dat$a_scenario)
  
p_Y_zero[[i]] <- 
  ggplot(dat) +
  geom_line(aes(x=true.mean,y=p_zero,color=alpha,group=sp)) +
  scale_color_viridis_c("Amp. Eff. (a)",option="H",breaks=c(0.25,0.4,0.55,0.7,0.85,1),limits=c(0.2,1.0)) +
  #ggtitle("P(0)") +
  theme_bw()

X <- seq(0,1,by=0.001)
Y <- dbeta(X,a_scen*0.7,a_scen*0.3)
BETA <- cbind(X=X,Y=Y) %>% as.data.frame()

p_a_inset[[i]] <-
  ggplot(dat) +
  geom_histogram(aes(x=alpha, y = ..density..),color = 1, fill = "white",binwidth=0.10,boundary=1) +
  geom_line(data=BETA,aes(x=X,y=Y),size=1.2) +
  scale_x_continuous("a",) +
  scale_y_continuous("Probability",expand=c(0,0.1),labels=NULL) +
  theme_classic() +
  theme(axis.title.y = element_text(size=10),
        axis.title.x=element_blank())
}

# Make a grid of plots.

LAY <- rbind(c(1,2,5),c(3,4,5))

quartz(file=here("signal2noise","figures","Figure_1.jpeg"),type="jpeg",dpi=600,height=7,width=9)
grid.arrange(layout_matrix = LAY,
             widths=c(1,1,0.25),
  p_Y_zero[[1]] +
    labs(subtitle="A")+
    scale_y_continuous("p(Y=0)",breaks=BREAKS.ZERO,limits=c(0,1.5),expand=c(0.01,0.01)) +
    scale_x_continuous("",trans="log",breaks=BREAKS) +
    theme(legend.position = "none")+
    annotation_custom(grob = ggplotGrob(p_a_inset[[1]]), 
                      xmin = -0.5, xmax = 5, ymin=1, ymax=1.5),
  p_Y_zero[[2]] +  
    labs(subtitle="B")+
    scale_y_continuous("",breaks=BREAKS.ZERO,limits=c(0,1.5),expand=c(0.01,0.01)) +
    scale_x_continuous("",trans="log",breaks=BREAKS) +
    theme(legend.position = "none")+
    annotation_custom(grob = ggplotGrob(p_a_inset[[2]]), 
                      xmin = -0.5, xmax = 5, ymin=1, ymax=1.5),
  p_Y_zero[[3]] +  
    labs(subtitle="C")+
    scale_y_continuous("p(Y=0)",breaks=BREAKS.ZERO,limits=c(0,1.5),expand=c(0.01,0.01)) +
    scale_x_continuous(expression("DNA concentration ("*lambda*"; copies / "*mu*"L)"),trans="log",breaks=BREAKS) +
        theme(legend.position = "none")+
    annotation_custom(grob = ggplotGrob(p_a_inset[[3]]), 
                      xmin = -0.5, xmax = 5, ymin=1, ymax=1.5),
  p_Y_zero[[4]] +
    labs(subtitle="D")+
    scale_y_continuous("",breaks=BREAKS.ZERO,limits=c(0,1.5),expand=c(0.01,0.01)) +
    scale_x_continuous(expression("DNA concentration ("*lambda*"; copies / "*mu*"L)"),trans="log",breaks=BREAKS) +
    theme(legend.position = "none")+
    annotation_custom(grob = ggplotGrob(p_a_inset[[4]]), 
                      xmin = -0.5, xmax = 5, ymin=1, ymax=1.5),
  get_legend(p_Y_zero[[1]])
)  
  
  dev.off()


  ### SUPPLEMENT FIGURE FOR SKEWED COMMUNITY DISTRIBUTIONS
  # Read in and run Skew plots.

load(here("signal2noise","code","simulate_overdispersion","Skew simulation; N_species= 20 pi = 20 N_pcr=35.RData"))
N_PCR <- 35 # needed for labeling below  

sim.all <- Output$sim.all
sim.prob.dat <- Output$sim.prob.dat

BREAKS.ZERO <- c(0,0.10,0.20,0.40,0.60,0.80,1.0)
n_scen <- c(5,10,100,1e+06)
p_Y_zero <- list()
#p_a_inset <- list()

dat.many <- sim.all %>% 
  filter(type=="Y") %>% 
  filter(a_scenario > 5-1, a_scenario < 5+1) %>%
  mutate(group.sp = paste(sp,a_scenario,sep="_"),
         true_prob = true.mean / lambda) 

a_scen <- mean(dat$a_scenario)

p_single <- ggplot(dat.many %>% filter(a_scenario == "5.0006")) +
  geom_line(aes(x=true.mean,y=p_zero,color=a,group=group.sp)) +
  scale_x_continuous(expression("Mean("*lambda*")"),trans="log",breaks=BREAKS) +
  coord_cartesian(xlim = c(0.1,1.5e+03)) +
  scale_y_continuous("P(Y=0)",breaks=BREAKS.ZERO,limits=c(0,1),expand=c(0.01,0.01)) +
  scale_color_viridis_c("Amp. Eff. (a)",option="H",breaks=c(0.25,0.4,0.55,0.7,0.85,1),limits=c(0.2,1.0)) +
  #ggtitle("P(0)") +
  facet_wrap(~true_prob) +
  theme_bw()

p_overplot <- ggplot(dat.many) +
  geom_line(aes(x=true.mean,y=p_zero,color=a,group=group.sp)) +
  scale_x_continuous(expression("Mean("*lambda*")"),trans="log",breaks=BREAKS) +
  coord_cartesian(xlim = c(0.1,1.5e+03)) +
  scale_y_continuous("P(Y=0)",breaks=BREAKS.ZERO,limits=c(0,1),expand=c(0.01,0.01)) +
  scale_color_viridis_c("Amp. Eff. (a)",option="H",breaks=c(0.25,0.4,0.55,0.7,0.85,1),limits=c(0.2,1.0)) +
  #ggtitle("P(0)") +
  facet_wrap(~true_prob) +
  theme_bw()



# Make a grid of plots.

LAY <- rbind(c(1),c(2))

quartz(file=here("signal2noise","figures","Figure_S1.jpeg"),type="jpeg",dpi=600,height=7,width=9)
grid.arrange(layout_matrix = LAY,
             widths=c(1),
             p_single +
               labs(subtitle="A"),
             p_overplot +  
               labs(subtitle="B")
        
)  

dev.off()



