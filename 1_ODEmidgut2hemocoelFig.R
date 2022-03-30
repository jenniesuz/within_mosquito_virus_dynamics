source("1_ODEmidgut2hemocoel.R")


sim <- simPop(parms=params())

longFCells <- cbind.data.frame("time"=rep(sim$time,2)
                               ,"cellType"=c(rep("Midgut",length(sim$time))
                                             ,rep("Hemocoel",length(sim$time)))
                               ,"numberCells"=c(sim$Mci,sim$Hci)
)
names(longFCells) <- c("time","cellType","numberCells")


longFVirus <- cbind.data.frame("time"=rep(sim$time,2)
                               ,"cellType"=c(rep("Midgut cells",length(sim$time))
                                             ,rep("Hemocoel cells",length(sim$time)))
                               ,"numberVirus"=c(sim$Mv,sim$Hv)
)
names(longFVirus) <- c("time","cellType","numberVirus")



cellsPlot <- ggplot(longFCells) +
  geom_line(aes(x=time,y=numberCells,col=cellType,group=cellType)) +
  scale_colour_manual(values=c("royalblue4","dodgerblue")) +
  xlab("Time (hours)") +
  ylab("Numbers of infected cells") +
  labs(col="Tissue",title="A") +
  theme_bw() +
  theme(panel.border = element_blank()
        ,axis.line = element_line(color = 'black')
        ,text=element_text(size=8)
        ,plot.margin=unit(c(0.2,0.3,0.1,0.1), "cm")
        ,axis.text=element_text(size=8)
        ,legend.key.size = unit(0.8,"line")
        ,legend.background = element_blank()
        ,legend.text=element_text(size=8)
        ,legend.position="none"
  )

virusPlot <- ggplot(longFVirus) +
  geom_line(aes(x=time,y=numberVirus,col=cellType,group=cellType)) +
  scale_colour_manual(values=c("royalblue4","dodgerblue")) +
  labs(col="Tissue",title="B") +
  xlab("Time (hours)") +
  ylab("Numbers of virus particles") +
  theme_bw() +
  theme(panel.border = element_blank()
        ,axis.line = element_line(color = 'black')
        ,text=element_text(size=8)
        ,plot.margin=unit(c(0.2,0.3,0.1,0.1), "cm")
        ,axis.text=element_text(size=8)
        ,legend.key.size = unit(0.8,"line")
        ,legend.background = element_blank()
        ,legend.text=element_text(size=8)
        ,legend.position=c(0.8,0.2)
  )

#pdf(file="fig_odeModelOutput.pdf",width=5,height=3)
grid.arrange(cellsPlot,virusPlot,ncol=2)
#dev.off()
