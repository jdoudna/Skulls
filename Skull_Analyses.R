#####################################################################################
#########################################################
######31 October 2013###########
####seperate analysis of jaw and skull - includes Chariton########
#raw data with metadata
###trial for git

library(geomorph)
library(gplots)
library(ggplot2)
library(gdata)

jawsliders<-matrix(c(1,12,13,12,13,14,13,14,2,7,8,15,8,15,16,15,16,17,16,17,9),7,3,byrow=TRUE)
skullsliders<-matrix(c(1,9,10,9,10,11,10,11,12,11,12,6,4,13,14,13,14,15,14,15,7),7,3,byrow=TRUE)

skulls.meta<-read.csv("Skulls_raw_in.csv")
jaws.meta<-read.csv("Jaws_raw_in.csv")
skullcoords<-skulls.meta[,c(7:38)] #no teeth
jawcoords<-jaws.meta[,c(7:50)]#no teeth
jawcoords.red<-jawcoords[,c(1:2,7:32,39:44)] #no teeth, remove weird curve



#fix atomic error
skulls.numeric<-NULL
for (i in 1:ncol(skullcoords)){
    tmp<-as.numeric(skullcoords[,i])
skulls.numeric<-cbind(skulls.numeric,tmp)
}

jaws.numeric<-NULL
for (i in 1:ncol(jawcoords.red)){
    tmp<-as.numeric(jawcoords.red[,i])
jaws.numeric<-cbind(jaws.numeric,tmp)
}

#make 3D array

jaws.coords<- arrayspecs(jaws.numeric[,1:34], 17, 2, byLand = FALSE)
skulls.coords<- arrayspecs(skulls.numeric[,1:32], 16, 2, byLand = FALSE)

#gpa
jaw.gpa<-gpagen(jaws.coords,curves=jawsliders)
skull.gpa<-gpagen(skulls.coords,curves=skullsliders)

#2D matrix
jaw2d<-two.d.array(jaw.gpa$coords)
skull2d<-two.d.array(skull.gpa$coords)

#analyze

jawlocation<-jaws.meta[,3] 
jawyear<-jaws.meta[,4]


skulllocation<-skulls.meta[,3]
skullyear<-skulls.meta[,4]


procD.lm(jaw2d~jawyear*jawlocation, iter=99)
traj.jaws<-trajectory.analysis(jaw2d~jawlocation*jawyear)

procD.lm(skull2d~skullyear*skulllocation, iter=99)
traj.skulls<-trajectory.analysis(skull2d~skulllocation*skullyear)

skull.allometry<-plotAllometry(skull.gpa$coords[,,-227], skull.gpa$Csize[-227], groups = skulllocation[-227], method = c("CAC", "RegScore", "PredLine"),
warpgrids = TRUE, label = FALSE)
skull.size.shape<-cbind(skullyear, skulllocation, skull.allometry$Csize, skull.allometry$allom.score)
###write.csv(skull.size.shape, file = "Skull CACs.csv")

mod1<-aov(skull.size.shape[,4]~skull.size.shape[,3]*skull.size.shape[,1])
summary(mod1)
mod2<-aov(skull.size.shape[,4]~skull.size.shape[,3]+skull.size.shape[,1])
summary(mod2)
anova(mod1,mod2)

skullreg1<-lm(skull.size.shape[1:150,4]~skull.size.shape[1:150,3])
skullreg2<-lm(skull.size.shape[151:309,4]~skull.size.shape[151:309,3])
plot(skull.size.shape[,4]~skull.size.shape[,3], type = 'n', xlab="Skull Csize", ylab="Allometric Score")
points(skull.size.shape[1:150,3],skull.size.shape[1:150,4], pch=20)
points(skull.size.shape[151:309,3],skull.size.shape[151:309,4], pch=1)
abline(skullreg1, lty=1)
abline(skullreg2, lty=2)
legend("bottomright", c("1900","2012"),lty=c(1,2),pch=c(20,1))

lm(skull.size.shape[1:150,4]~skull.size.shape[1:150,3])
lm(jaw.size.shape[151:310,4]~jaw.size.shape[151:310,3])

jaw.allometry<-plotAllometry(jaw.gpa$coords, jaw.gpa$Csize, groups=jawlocation, method = c( "CAC", "RegScore", "PredLine"),
warpgrids = TRUE, label = FALSE)
jaw.size.shape<-data.frame(jawyear, jawlocation, jaw.allometry$Csize, jaw.allometry$allom.score)
###write.csv(jaw.size.shape, file = "Jaw CACs.csv")

colnames(jaw.size.shape) <- c("Jaw_Year", "Jaw_Location", "Csize", "CAC")
split.jaws<-split(jaw.size.shape, jaw.size.shape$Jaw_Location)

mean.CACs.jaws<-ddply(jaw.size.shape, .(Jaw_Location), summarize, mean.jaw.CAC<-mean(CAC))
mean.Csizes.jaws<-ddply(jaw.size.shape, .(Jaw_Location), summarize, mean.jaw.Csize<-mean(Csize))


mod4<-aov(jaw.size.shape[,4]~jaw.size.shape[,3]*jaw.size.shape[,1])
summary(mod4)
mod5<-aov(jaw.size.shape[,4]~jaw.size.shape[,3]+jaw.size.shape[,1])
summary(mod5)
anova(mod4,mod5)

jawreg1<-lm(jaw.size.shape[1:150,4]~jaw.size.shape[1:150,3])
jawreg2<-lm(jaw.size.shape[151:309,4]~jaw.size.shape[151:309,3])
plot(jaw.size.shape[,4]~jaw.size.shape[,3], type = 'n', xlab="Jaw Csize", ylab="Allometric Score")
points(jaw.size.shape[1:150,3],jaw.size.shape[1:150,4], pch=20)
points(jaw.size.shape[151:309,3],jaw.size.shape[151:309,4], pch=1)
abline(jawreg1, lty=1)
abline(jawreg2, lty=2)
legend("bottomright", c("1900","2012"), lty=c(1,2),pch=c(20,1))




lm(jaw.size.shape[1:150,4]~jaw.size.shape[1:150,3])
lm(jaw.size.shape[151:310,4]~jaw.size.shape[151:310,3])

##add year groups
old<-c(1:150)
new<-c(151:310)
old.allometry.jaw<-plotAllometry(jaw.gpa$coords[,,c(old)], jaw.gpa$Csize[c(old)], method = c("CAC", "RegScore", "PredLine"),
warpgrids = TRUE, label = FALSE)
old.allometry.skull<-plotAllometry(skull.gpa$coords[,,c(old)], skull.gpa$Csize[c(old)], method = c("CAC", "RegScore", "PredLine"),
warpgrids = TRUE, label = FALSE)
new.allometry.jaw<-plotAllometry(jaw.gpa$coords[,,c(new)], jaw.gpa$Csize[c(new)], method = c("CAC", "RegScore", "PredLine"),
warpgrids = TRUE, label = FALSE)
new.allometry.skull<-plotAllometry(skull.gpa$coords[,,c(new)], skull.gpa$Csize[c(new)], method = c("CAC", "RegScore", "PredLine"),
warpgrids = TRUE, label = FALSE)



# Multiple Linear Regression Example 
old.fit.allom.jaw <- lm( old.allometry.jaw$allom.score~old.allometry.jaw$Csize)
summary(old.fit.allom.jaw) # show results
new.fit.allom.jaw <- lm(new.allometry.jaw$allom.score~new.allometry.jaw$Csize)
summary(new.fit.allom.jaw)

old.fit.allom.skull <- lm( old.allometry.skull$allom.score~old.allometry.skull$Csize)
summary(old.fit.allom.skull) # show results
new.fit.allom.skull <- lm(new.allometry.skull$allom.score~new.allometry.skull$Csize)
summary(new.fit.allom.skull)

plotTangentSpace(jaw.gpa$coords, axis1 = 1, axis2 = 2,
    warpgrids = TRUE, label = FALSE)

plotTangentSpace(skull.gpa$coords, axis1 = 1, axis2 = 2,
    warpgrids = TRUE, label = FALSE)

jaw.shape<-prcomp(two.d.array(jaw.gpa$coords))$x[,1:(22*2-4)]

skull.shape<-prcomp(two.d.array(skull.gpa$coords))$x[,1:(16*2-4)]

fit<-aov(jaw.gpa$Csize ~ jawyear*jawlocation)
summary(fit)

fit2<-aov(skull.gpa$Csize ~ skullyear*skulllocation)
summary(fit2)

# Two-way Interaction Plot
interaction.plot(skullyear,skulllocation, skull.gpa$Csize, type="b", col=c(1:3), 
  	 leg.bty="o", leg.bg="beige", lwd=2, pch=c(18,24,22),	
   xlab="Year", 
   ylab="Mean Csize", 
   main="Interaction Plot for Skulls")

# Plot Means with Error Bars
plotmeans(skull.gpa$Csize~skullyear,xlab="Year",
  ylab="Skull Csize", main="Mean Plot\nwith 95% CI", connect=FALSE)

# Two-way Interaction Plot
interaction.plot(jawyear,jawlocation, jaw.gpa$Csize, type="b", col=c(1:3), 
  	 leg.bty="o", leg.bg="beige", lwd=2, pch=c(18,24,22),	
   xlab="Year", 
   ylab="Mean Csize", 
   main="Interaction Plot for Jaws")

# Plot Means with Error Bars
plotmeans(jaw.gpa$Csize~jawyear,xlab="Year",
  ylab="Jaw Csize", main="Mean Plot\nwith 95% CI", connect=FALSE)
###########################################################
##Find shape changes by location##
emnd_old_ref_jaw<-mshape(jaw.gpa$coords[,,c(7:12,29:31,48:51,138:146)])
emnd_new_ref_jaw<-mshape(jaw.gpa$coords[,,280:310])
iowa_old_ref_jaw<-mshape(jaw.gpa$coords[,,c(36:46,51:63,124)])
iowa_new_ref_jaw<-mshape(jaw.gpa$coords[,,151:178])
chariton_ref_jaw<-mshape(jaw.gpa$coords[,,311:338])
nill_old_ref_jaw<-mshape(jaw.gpa$coords[,,c(16:28,135:137)])
nill_new_ref_jaw<-mshape(jaw.gpa$coords[,,179:202])
maka_old_ref_jaw<-mshape(jaw.gpa$coords[,,c(5,6,13:15,64:88,117:118)])
maka_new_ref_jaw<-mshape(jaw.gpa$coords[,,238:253])
otka_old_ref_jaw<-mshape(jaw.gpa$coords[,,c(1,92:109,112:116,150,119:123,125:134)])
otka_new_ref_jaw<-mshape(jaw.gpa$coords[,,203:237])
wamn_old_ref_jaw<-mshape(jaw.gpa$coords[,,c(2:4,32:35,89:91,147:149,110:111)])
wamn_new_ref_jaw<-mshape(jaw.gpa$coords[,,254:279])

emnd_old_ref_skull<-mshape(skull.gpa$coords[,,c(7:12,29:31,48:51,138:146)])
emnd_new_ref_skull<-mshape(skull.gpa$coords[,,280:310])
iowa_old_ref_skull<-mshape(skull.gpa$coords[,,c(36:46,51:63,124)])
iowa_new_ref_skull<-mshape(skull.gpa$coords[,,151:178])
chariton_ref_skull<-mshape(skull.gpa$coords[,,311:338])
nill_old_ref_skull<-mshape(skull.gpa$coords[,,c(16:28,135:137)])
nill_new_ref_skull<-mshape(skull.gpa$coords[,,179:202])
maka_old_ref_skull<-mshape(skull.gpa$coords[,,c(5,6,13:15,64:88,117:118)])
maka_new_ref_skull<-mshape(skull.gpa$coords[,,238:253])
otka_old_ref_skull<-mshape(skull.gpa$coords[,,c(1,92:109,112:116,150,119:123,125:134)])
otka_new_ref_skull<-mshape(skull.gpa$coords[,,203:237])
wamn_old_ref_skull<-mshape(skull.gpa$coords[,,c(2:4,32:35,89:91,147:149,110:111)])
wamn_new_ref_skull<-mshape(skull.gpa$coords[,,254:279])

jawlinks<-matrix(c(1,12,12,13,13,14,14,2,2,3,3,4,4,5,5,6,6,7,7,8,8,15,15,16,16,17,17,9,9,10,10,11,11,1),17,2,byrow=TRUE)
skulllinks<-matrix(c(1,2,2,3,3,7,7,15,15,14,14,13,13,4,4,5,5,6,6,7,8,6,8,3,6,12,12,11,11,10,10,9,9,1),17,2,byrow=TRUE)

old<-matrix(c(0,2,1,0,0,1.75),nrow=3,ncol=2)
new<-matrix(c(0,2,1,0,0,5), nrow=3,ncol=2)
plotRefToTarget(new,old)
##plotreftotarget(old,new) is how shape changed from old to new - what has happened
plotRefToTarget(emnd_old_ref_skull,emnd_new_ref_skull,mag=3, links=skulllinks)
plotRefToTarget(iowa_old_ref_skull,iowa_new_ref_skull,mag=3, links=skulllinks)
plotRefToTarget(nill_old_ref_skull,nill_new_ref_skull,mag=3,links=skulllinks)
plotRefToTarget(maka_old_ref_skull,maka_new_ref_skull,mag=3,links=skulllinks)
plotRefToTarget(otka_old_ref_skull,otka_new_ref_skull,mag=3,links=skulllinks)
plotRefToTarget(wamn_old_ref_skull,wamn_new_ref_skull,mag=3,links=skulllinks)

plotRefToTarget(emnd_old_ref_jaw,emnd_new_ref_jaw,mag=3, links=jawlinks)
plotRefToTarget(iowa_old_ref_jaw,iowa_new_ref_jaw,mag=3, links=jawlinks)
plotRefToTarget(nill_old_ref_jaw,nill_new_ref_jaw,mag=3,links=jawlinks)
plotRefToTarget(maka_old_ref_jaw,maka_new_ref_jaw,mag=3,links=jawlinks)
plotRefToTarget(otka_old_ref_jaw,otka_new_ref_jaw,mag=3,links=jawlinks)
plotRefToTarget(wamn_old_ref_jaw,wamn_new_ref_jaw,mag=3,links=jawlinks)

######################################################
