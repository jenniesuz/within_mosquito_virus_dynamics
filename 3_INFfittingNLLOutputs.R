library(gridExtra)
library(ggplot2)
library(binom)
library(viridis)


#*********Parameter estimates****************


aeg <- readRDS("aegInffits.rds")  # read in parameter value estimates from fits
alb <- readRDS("albInffits.rds")

aegParams <- sapply(aeg,"[[",2)
albParams <- sapply(alb,"[[",2)

aeg <- cbind.data.frame(muV=aegParams[1,]        # convert from list to dataframe
                           ,infRate=aegParams[2,]
                           ,moz="Ae. aegypti"
                           )
alb <- cbind.data.frame(muV=albParams[1,]
                           ,infRate=albParams[2,]
                           ,moz="Ae. albopictus"
                           )

dat <- rbind.data.frame(aeg,alb)


#****************Plot fitted parameter value estimates**********

muVPlot <- ggplot(params,aes(x=escapeRate,fill=species)) + 
  geom_histogram(position="identity", alpha = 0.3) +
  scale_fill_manual(values=c("royalblue4","dodgerblue")) +
  scale_x_log10() +
  xlab("Escape rate") +
  ylab("Number of fits") +
  labs(fill="") +
  ylab("Number of fits") +
  theme_bw() +
  theme(panel.border = element_blank()
        ,axis.line = element_line(color = 'black')
        ,text=element_text(size=10)
        ,plot.margin=unit(c(0.2,0.3,0.1,0.1), "cm")
        ,axis.text=element_text(size=10)
        ,legend.key.size = unit(0.8,"line")
        ,legend.background = element_blank()
        ,legend.text=element_text(size=10)
  )


betaPlot <- ggplot(dat,aes(x=infRate)) + 
  geom_histogram(data=dat,aes(fill = moz), alpha = 0.6) +
  scale_fill_manual(values=c("royalblue4","dodgerblue")) +
  scale_x_log10() +
  labs(fill="") +
  xlab(expression("Fitted rate at which virions infect susceptible cells log"[10]*italic(beta))) +
  ylab("Number of fits") +
  theme_bw() +
  theme(panel.border = element_blank()
        ,axis.line = element_line(color = 'black')
        ,text=element_text(size=10)
        ,plot.margin=unit(c(0.2,0.3,0.1,0.1), "cm")
        ,axis.text=element_text(size=10)
        ,legend.key.size = unit(0.8,"line")
        ,legend.background = element_blank()
        ,legend.text=element_text(size=10)
  )


#***************Simulations**********************

aegInfsims <- readRDS("aegInfSims.rds")
albInfsims  <- readRDS("albInfSims.rds")

aegInfsims <- do.call(rbind.data.frame,aegInfsims)
aegInfsims$Moz <- "Ae. aegypti"
aegInfsims$parmsComb <- paste(aegInfsims$infRate,aegInfsims$muV)
albInfsims <- do.call(rbind.data.frame,albInfsims)
albInfsims$Moz <- "Ae. albopictus"
albInfsims$parmsComb <- paste(albInfsims$infRate,albInfsims$muV)

sims <- rbind.data.frame(aegInfsims,albInfsims)


#******************data to fit to********************
competenceDat <- read.csv("datCiotaOnyango.csv")
competenceDat <- competenceDat[competenceDat$Ref %in% "Ciota 2017",]
bin <- binom.confint(x=competenceDat$NumInf,n=competenceDat$ITotal,method="exact")
competenceDat$meanInf <- bin$mean
competenceDat$lowerInf <- bin$lower
competenceDat$upperInf <- bin$upper
#****************************************************

infFitsPlot <- ggplot(competenceDat) +
  geom_line(data=sims,aes(x=log10(Conc.Min),y=meanInf,col=Moz,group=parmsComb),linetype=3,alpha=0.7) +
  geom_errorbar(aes(x=Conc.Min,ymin=lowerInf,ymax=upperInf),alpha=0.5) +
  geom_point(aes(x=Conc.Min,y=meanInf,fill=Moz),shape=21)+
  # geom_point(data=sims,aes(x=log10(Conc.Min),y=meanInf,col=Moz),size=2)  +
  xlab(expression(paste("Virus concentration (log"[10]," PFU/ml)"))) +
  labs(fill="",title="A") +
  scale_fill_manual(values=c("royalblue4","dodgerblue")) +
  scale_colour_manual(values=c("royalblue4","dodgerblue")) +
  ylab("Proportion of mosquitoes with a midgut infection") +
  theme_set(theme_bw())  +    
  theme(panel.border = element_blank()                   
        ,axis.line = element_line(color = 'black')
        ,text=element_text(size=8)
        ,plot.margin=unit(c(0.2,0.1,0.1,0.1), "cm")
        ,axis.text=element_text(size=6)
        ,legend.key.size = unit(0.5,"line")  
        ,legend.background = element_blank()
        ,legend.text=element_text(size=6)
        ,strip.text=element_text(size=5)
        ,legend.position="none"
  )   
