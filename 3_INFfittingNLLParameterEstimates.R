
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

aegMuV <- median(dat$muV[dat$moz %in% "Ae. aegypti"])
albMuV <- median(dat$muV[dat$moz %in% "Ae. albopictus"])

aeginfRate <- median(dat$infRate[dat$moz %in% "Ae. aegypti"])
albinfRate <- median(dat$infRate[dat$moz %in% "Ae. albopictus"])

