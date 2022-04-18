library(gridExtra)
library(ggplot2)
library(binom)
library(viridis)


#*********Parameters estimated from midgut infection data****************
aegInf <- readRDS("aegInffits.rds")  # read in parameter value estimates from fits
albInf <- readRDS("albInffits.rds")

aegInfParams <- sapply(aegInf,"[[",2)
albInfParams <- sapply(albInf,"[[",2)

aegInf <- cbind.data.frame(muV=aegInfParams[1,]        # convert from list to dataframe
                           ,probInf=aegInfParams[2,]
                           ,species="Ae. aegypti"
)
albInf <- cbind.data.frame(muV=albInfParams[1,]
                           ,probInf=albInfParams[2,]
                           ,species="Ae. albopictus"
)

datInf <- rbind.data.frame(aegInf,albInf)


#************plot**************
muVPlot <- ggplot(datInf,aes(x=muV,fill=species)) + 
  geom_histogram(position="identity", alpha = 0.3) +
  scale_fill_manual(values=c("royalblue4","dodgerblue")) +
  scale_x_log10() +
  xlab(expression(paste("Virus decay rate (",italic(mu)["V"],")"))) +
  ylab("Number of fits") +
  labs(fill="") +
  ylab("Number of fits") +
  theme_bw() +
  theme(panel.border = element_blank()
        ,axis.line = element_line(color = 'black')
        ,text=element_text(size=7)
        ,plot.margin=unit(c(0.2,0.3,0.1,0.1), "cm")
        ,axis.text=element_text(size=6)
        ,legend.key.size = unit(0.8,"line")
        ,legend.background = element_blank()
        ,legend.text=element_text(size=6)
        ,legend.position=c(0.2,0.9)
  )


betaPlot <- ggplot(datInf,aes(x=probInf,fill=species)) + 
  geom_histogram(position="identity", alpha = 0.3) +
  scale_fill_manual(values=c("royalblue4","dodgerblue")) +
  scale_x_log10() +
  xlab(expression(paste("Rate at which virions infect susceptible cells (log"[10]*italic(beta),")"))) +
  ylab("") +
  labs(fill="") +
  theme_bw() +
  theme(panel.border = element_blank()
        ,axis.line = element_line(color = 'black')
        ,text=element_text(size=7)
        ,plot.margin=unit(c(0.2,0.3,0.1,0.1), "cm")
        ,axis.text=element_text(size=6)
        ,legend.key.size = unit(0.8,"line")
        ,legend.background = element_blank()
        ,legend.text=element_text(size=6)
        ,legend.position="none"
  )


#****************************Dissemination********************************
#*
#*
library(gridExtra)
library(ggplot2)
library(plyr)
library(ggfan)
library(binom)

# read in model fits fits
albF <- readRDS("albDissFits.rds")
aegF <- readRDS("aegDissFits.rds")


# function to change to dataframe and then bind together
listToDat <- function(output=alb, species="alb"){
  params <- lapply(1:length(output),function(x){
    temp <- output[[x]]
    if(is.list(temp)==T){
      temp <- temp[[1]][1]
      temp <- do.call(c,temp)
      names(temp) <- c("prodRate","cellSpread","escapeRate")
      return(exp(temp))
    }else{
      temp <- c(NA,NA,NA)
      names(temp) <- c("prodRate","cellSpread","escapeRate")
      return(temp)
    }
  })
  params <- do.call(rbind.data.frame,params)
  names(params) <- c("prodRate","cellSpread","escapeRate")
  params$species=species
  return(params)
}

albDiss <- listToDat(albF)
aegDiss <- listToDat(output=aegF,species="aeg")

datDiss <- rbind.data.frame(albDiss,aegDiss)

datDiss$species[datDiss$species %in% "aeg"] <- "Ae. aegypti"
datDiss$species[datDiss$species %in% "alb"] <- "Ae. albopictus"


pRPlot <- ggplot(datDiss,aes(x=prodRate,fill=species)) + 
  geom_histogram(position="identity",alpha=0.3) +
  scale_fill_manual(values=c("royalblue4","dodgerblue")) +
  scale_x_log10() +
  xlab(expression(paste("Infected cell virus production rate (",italic(gamma),")"))) +
  ylab("Number of fits") +
  labs(fill="") +
  ylab("Number of fits") +
  theme_bw() +
  theme(panel.border = element_blank()
        ,axis.line = element_line(color = 'black')
        ,text=element_text(size=7)
        ,plot.margin=unit(c(0.2,0.3,0.1,0.1), "cm")
        ,axis.text=element_text(size=6)
        ,legend.key.size = unit(0.8,"line")
        ,legend.background = element_blank()
        ,legend.text=element_text(size=8)
        ,legend.position="none"
  )

csPlot <- ggplot(datDiss,aes(x=cellSpread,fill=species)) + 
  geom_histogram(position="identity", alpha = 0.3) +
  scale_fill_manual(values=c("royalblue4","dodgerblue")) +
  scale_x_log10() +
  xlab(expression(paste("Rate virions spread between midgut cells (",italic(alpha),")"))) +
  ylab("") +
  labs(fill="") +
  theme_bw() +
  theme(panel.border = element_blank()
        ,axis.line = element_line(color = 'black')
        ,text=element_text(size=7)
        ,plot.margin=unit(c(0.2,0.3,0.1,0.1), "cm")
        ,axis.text=element_text(size=6)
        ,legend.key.size = unit(0.8,"line")
        ,legend.background = element_blank()
        ,legend.text=element_text(size=8)
        ,legend.position="none"
  )

erPlot <- ggplot(datDiss,aes(x=escapeRate,fill=species)) + 
  geom_histogram(position="identity", alpha = 0.3) +
  scale_fill_manual(values=c("royalblue4","dodgerblue")) +
  scale_x_log10() +
  xlab(expression(paste("Rate of virion escape into the hemocoel (",italic(rho),")"))) +
  ylab("Number of fits") +
  labs(fill="") +
  theme_bw() +
  theme(panel.border = element_blank()
        ,axis.line = element_line(color = 'black')
        ,text=element_text(size=7)
        ,plot.margin=unit(c(0.2,0.3,0.1,0.1), "cm")
        ,axis.text=element_text(size=6)
        ,legend.key.size = unit(0.8,"line")
        ,legend.background = element_blank()
        ,legend.text=element_text(size=8)
        ,legend.position="none"
  )


pdf(file="fig_modelParameterEstimates.pdf",width=5,height=5)
grid.arrange(muVPlot,betaPlot,pRPlot,csPlot,erPlot,ncol=2)
dev.off()


