rm(list=ls())
gc()
graphics.off()

quartz()

library(tidyverse)
library(betapart)
library(BAT)
library(ggplot2)
library(ggrepel)

#######################
##Species communities##
#######################

dat.all <- read.csv("PaciFlora/Species_list_full_2905.csv",sep=";")
dat.all <- dat.all[-which(is.na(dat.all$species)),]
dat.all <- dat.all[-which(is.na(dat.all$island)),]

dat.soc <- dat.all[which(dat.all$islandgroup=="Society"),]
# dat.haw <- dat.all[which(dat.all$islandgroup=="Hawaiian"),]
# dat.sam <- dat.all[which(dat.all$islandgroup=="Samoa"),]
# dat.mar <- dat.all[which(dat.all$islandgroup=="Marquesas"),]
# dat.fij <- dat.all[which(dat.all$islandgroup=="Fiji"),]
# dat.comb <- rbind(dat.soc,dat.haw,dat.sam,dat.mar,dat.fij)

length(unique(dat.soc$island))
# length(unique(dat.haw$island))
# length(unique(dat.sam$island))
# length(unique(dat.mar$island))
# length(unique(dat.fij$island))

dat.soc.red <- dat.soc[,c("species","island")]
# dat.haw.red <- dat.haw[,c("species","island")]
# dat.sam.red <- dat.sam[,c("species","island")]
# dat.mar.red <- dat.mar[,c("species","island")]
# dat.fij.red <- dat.fij[,c("species","island")]
# dat.comb.red <- dat.comb[,c("species","island")]
dat.soc.red$presence <- 1
# dat.haw.red$presence <- 1
# dat.sam.red$presence <- 1
# dat.mar.red$presence <- 1
# dat.fij.red$presence <- 1
# dat.comb.red$presence <- 1

##reshape - pivot matrix
dat.soc.pa <- dat.soc.red %>% 
  pivot_wider(names_from=species,values_from=c(presence))
list0 <- as.list(rep(0,ncol(dat.soc.pa)))
names(list0) <- names(dat.soc.pa)
dat.soc.pa <- as.data.frame(dat.soc.pa %>% replace_na(list0))
row.names(dat.soc.pa) <- dat.soc.pa$island
dat.soc.pa <- dat.soc.pa[,-1]

# dat.haw.pa <- dat.haw.red %>% 
#   pivot_wider(names_from=species,values_from=c(presence))
# list0 <- as.list(rep(0,ncol(dat.haw.pa)))
# names(list0) <- names(dat.haw.pa)
# dat.haw.pa <- as.data.frame(dat.haw.pa %>% replace_na(list0))
# row.names(dat.haw.pa) <- dat.haw.pa$island
# dat.haw.pa <- dat.haw.pa[,-1]
# 
# dat.sam.pa <- dat.sam.red %>% 
#   pivot_wider(names_from=species,values_from=c(presence))
# list0 <- as.list(rep(0,ncol(dat.sam.pa)))
# names(list0) <- names(dat.sam.pa)
# dat.sam.pa <- as.data.frame(dat.sam.pa %>% replace_na(list0))
# row.names(dat.sam.pa) <- dat.sam.pa$island
# dat.sam.pa <- dat.sam.pa[,-1]
# 
# dat.mar.pa <- dat.mar.red %>% 
#   pivot_wider(names_from=species,values_from=c(presence))
# list0 <- as.list(rep(0,ncol(dat.mar.pa)))
# names(list0) <- names(dat.mar.pa)
# dat.mar.pa <- as.data.frame(dat.mar.pa %>% replace_na(list0))
# row.names(dat.mar.pa) <- dat.mar.pa$island
# dat.mar.pa <- dat.mar.pa[,-1]
# 
# dat.fij.pa <- dat.fij.red %>% 
#   pivot_wider(names_from=species,values_from=c(presence))
# list0 <- as.list(rep(0,ncol(dat.fij.pa)))
# names(list0) <- names(dat.fij.pa)
# dat.fij.pa <- as.data.frame(dat.fij.pa %>% replace_na(list0))
# row.names(dat.fij.pa) <- dat.fij.pa$island
# dat.fij.pa <- dat.fij.pa[,-1]
# 
# dat.comb.pa <- dat.comb.red %>% 
#   pivot_wider(names_from=species,values_from=c(presence))
# list0 <- as.list(rep(0,ncol(dat.comb.pa)))
# names(list0) <- names(dat.comb.pa)
# dat.comb.pa <- as.data.frame(dat.comb.pa %>% replace_na(list0))
# row.names(dat.comb.pa) <- dat.comb.pa$island
# dat.comb.pa <- dat.comb.pa[,-1]

dat.soc.pa <- dat.soc.pa[order(row.names(dat.soc.pa)),]

which(rowSums(dat.soc.pa)==0)
which(colSums(dat.soc.pa)==0)

##Gamma
gamma.soc <- ncol(dat.soc.pa)
##Alpha
alpha.soc <- rowSums(dat.soc.pa)
alpha.soc.ratio <- alpha.soc/gamma.soc
alpha.soc.mean <- mean(alpha.soc)
alpha.soc.mean.ratio <- alpha.soc.mean/gamma.soc
##Beta
beta.soc1 <- beta.pair(dat.soc.pa)
beta.soc2 <- beta(dat.soc.pa,func = "sorensen")

plot(beta.soc1$beta.sor,beta.soc2$Btotal,xlim=c(0,1),ylim=c(0,1))
lines(c(0,1),c(0,1),col="red")
plot(beta.soc1$beta.sim,beta.soc2$Brepl,xlim=c(0,1),ylim=c(0,1))
lines(c(0,1),c(0,1),col="red")

##########
##Traits##
##########
library(TR8)

#names.spp.soc <- sub("_"," ",names.spp.soc)

available_traits()

names.spp.soc <- unique(dat.soc$species)

print(Sys.time())
Soc.sp.tr8 <- tr8(species_list = names.spp.soc, download_list=c("seed_mas_cal","seed_mass","seed_wght","SeedMass","Height","h_max","max_height_cal"),allow_persistent=FALSE)
print(Sys.time())

View(Soc.sp.tr8@results)

traits <- Soc.sp.tr8@results

library(stringi)
traits$h_max <- as.numeric(stri_extract_first_regex(traits$h_max, "[0-9]+"))
traits$seed_wght <- as.numeric(stri_extract_first_regex(traits$seed_wght, "[0-9]+"))
traits$Height <- as.numeric(stri_extract_first_regex(traits$Height, "[0-9]+"))
traits$SeedMass <- as.numeric(stri_extract_first_regex(traits$SeedMass, "[0-9]+"))
View(traits)

traits.mean <- data.frame(height=numeric(nrow(traits)),seed.mass=numeric(nrow(traits)))
row.names(traits.mean) <- row.names(traits)
traits.mean$height <- apply(traits[,c(1,4)],1,mean,na.rm=TRUE)
traits.mean$seed.mass <- apply(traits[,c(2,3,5)],1,mean,na.rm=TRUE)

which(!is.na(traits.mean$height) & !is.na(traits.mean$seed.mass))

traits.mean <- traits.mean[-which(is.na(traits.mean$height) | is.na(traits.mean$seed.mass)),]
traits.mean[which(traits.mean==0,arr.ind = TRUE)] <- 0.01

dat.soc.pa.trait <- dat.soc.pa[,which(colnames(dat.soc.pa) %in% row.names(traits.mean))]
colSums(dat.soc.pa.trait)
rowSums(dat.soc.pa.trait)
dat.soc.pa.trait <- dat.soc.pa.trait[-which(rowSums(dat.soc.pa.trait)<3),]

plot(traits.mean)
#traits.mean$height <- log(traits.mean$height)
traits.mean$seed.mass <- log(traits.mean$seed.mass)
plot(traits.mean)

##Convex hull method

Soc.hull <- hull.build(comm = dat.soc.pa.trait, trait = traits.mean)
gamma.trait <- hull.gamma(comm = Soc.hull)
alpha.trait <- hull.alpha(comm = Soc.hull)
alpha.trait.ratio <- alpha.trait/gamma.trait
alpha.trait.mean <- mean(alpha.trait)
alpha.trait.mean.ratio <- alpha.trait.mean/gamma.trait
beta.trait <- hull.beta(comm = Soc.hull)

# plot(alpha.soc[which(names(alpha.soc) %in% names(alpha.trait))],alpha.trait)
# plot(alpha.soc.ratio[which(names(alpha.soc.ratio) %in% names(alpha.trait.ratio))],alpha.trait.ratio,xlim=c(-0.1,1.1),ylim=c(-0.1,1.1))
# lines(c(-10,10),c(-10,10),col="red")
# text(alpha.soc.ratio[which(names(alpha.soc.ratio) %in% names(alpha.trait.ratio))],alpha.trait.ratio,names(alpha.trait.ratio),pos=2)
dat <- data.frame(island=names(alpha.trait.ratio),alpha.tax=alpha.soc.ratio[which(names(alpha.soc.ratio) %in% names(alpha.trait.ratio))],alpha.trait=alpha.trait.ratio)
ggplot(dat,aes(alpha.tax,alpha.trait,label=island))+
  geom_point()+
  geom_text_repel()+
  geom_abline(intercept=0,slope=1,color="red")+
  xlim(0,1)+
  ylim(0,1)

cor(dat$alpha.tax,dat$alpha.trait,method="pearson")
cor(dat$alpha.tax,dat$alpha.trait,method="spearman")

beta.soc2.temp <- lapply(beta.soc2,as.matrix)
beta.soc2.red <- list()
beta.soc2.red$Btotal <- c(as.dist(beta.soc2.temp$Btotal[which(row.names(beta.soc2.temp$Btotal) %in% names(beta.trait$Btotal)),which(colnames(beta.soc2.temp$Btotal) %in% names(beta.trait$Btotal))]))
beta.soc2.red$Brepl <- c(as.dist(beta.soc2.temp$Brepl[which(row.names(beta.soc2.temp$Brepl) %in% names(beta.trait$Brepl)),which(colnames(beta.soc2.temp$Brepl) %in% names(beta.trait$Brepl))]))
beta.soc2.red$Brich <- c(as.dist(beta.soc2.temp$Brich[which(row.names(beta.soc2.temp$Brich) %in% names(beta.trait$Brich)),which(colnames(beta.soc2.temp$Brich) %in% names(beta.trait$Brich))]))

par(mfrow=c(1,2))
plot(beta.soc2.red$Btotal,c(beta.trait$Btotal),xlim=c(0,1),ylim=c(0,1))
lines(c(0,1),c(0,1),col="red")
plot(beta.soc2.red$Brepl,c(beta.trait$Brepl),xlim=c(0,1),ylim=c(0,1))
lines(c(0,1),c(0,1),col="red")

cor(beta.soc2.red$Btotal,c(beta.trait$Btotal),method="pearson")
cor(beta.soc2.red$Btotal,c(beta.trait$Btotal),method="spearman")
cor(beta.soc2.red$Brepl,c(beta.trait$Brepl),method="pearson")
cor(beta.soc2.red$Brepl,c(beta.trait$Brepl),method="spearman")


##Dendrogram method

trait.tree1 <- tree.build(trait = traits.mean,func="nj") ##try different algorithms
trait.tree2 <- tree.build(trait = traits.mean, func="upgma") ##try different algorithms

par(mfrow=c(1,2))
plot(trait.tree1)
plot(trait.tree2)

gamma.dend1 <- gamma(comm=dat.soc.pa.trait,tree=trait.tree1)
alpha.dend1 <- alpha(comm=dat.soc.pa.trait,tree=trait.tree1)
alpha.dend.ratio1 <- alpha.dend1/c(gamma.dend1)
alpha.dend.mean.1 <- mean(alpha.dend1)
alpha.dend.mean.ratio1 <- alpha.dend.mean.1/gamma.dend1
beta.dend1 <- beta(comm=dat.soc.pa.trait,tree=trait.tree1)

gamma.dend2 <- gamma(comm=dat.soc.pa.trait,tree=trait.tree2)
alpha.dend2 <- alpha(comm=dat.soc.pa.trait,tree=trait.tree2)
alpha.dend.ratio2 <- alpha.dend2/c(gamma.dend2)
alpha.dend.mean.2 <- mean(alpha.dend2)
alpha.dend.mean.ratio2 <- alpha.dend.mean.2/gamma.dend2
beta.dend2 <- beta(comm=dat.soc.pa.trait,tree=trait.tree2)


par(mfrow=c(1,2))
plot(c(beta.dend1$Btotal),c(beta.dend2$Btotal),xlim=c(0,1),ylim=c(0,1))
lines(c(0,1),c(0,1),col="red")
plot(c(beta.dend1$Brepl),c(beta.dend2$Brepl),xlim=c(0,1),ylim=c(0,1))
lines(c(0,1),c(0,1),col="red")

dat <- data.frame(island=names(alpha.trait.ratio),alpha.tax=alpha.soc.ratio[which(names(alpha.soc.ratio) %in% names(alpha.trait.ratio))],alpha.trait=c(alpha.dend.ratio1))
ggplot(dat,aes(alpha.tax,alpha.dend1,label=island))+
  geom_point()+
  geom_text_repel()+
  geom_abline(intercept=0,slope=1,color="red")+
  xlim(0,1)+
  ylim(0,1)

cor(dat$alpha.tax,dat$alpha.dend1,method="pearson")
cor(dat$alpha.tax,dat$alpha.dend1,method="spearman")

dat <- data.frame(island=names(alpha.trait.ratio),alpha.tax=alpha.soc.ratio[which(names(alpha.soc.ratio) %in% names(alpha.trait.ratio))],alpha.trait=c(alpha.dend.ratio2))
ggplot(dat,aes(alpha.tax,alpha.dend2,label=island))+
  geom_point()+
  geom_text_repel()+
  geom_abline(intercept=0,slope=1,color="red")+
  xlim(0,1)+
  ylim(0,1)

cor(dat$alpha.tax,dat$alpha.dend2,method="pearson")
cor(dat$alpha.tax,dat$alpha.dend2,method="spearman")


plot(c(beta.dend1$Btotal),c(beta.dend2$Btotal),xlim=c(0,1),ylim=c(0,1))
lines(c(0,1),c(0,1),col="red")
plot(c(beta.dend1$Brepl),c(beta.dend2$Brepl),xlim=c(0,1),ylim=c(0,1))
lines(c(0,1),c(0,1),col="red")


plot(beta.soc2.red$Btotal,c(beta.dend1$Btotal),xlim=c(0,1),ylim=c(0,1))
lines(c(0,1),c(0,1),col="red")
plot(beta.soc2.red$Brepl,c(beta.dend1$Brepl),xlim=c(0,1),ylim=c(0,1))
lines(c(0,1),c(0,1),col="red")

plot(beta.soc2.red$Btotal,c(beta.dend2$Btotal),xlim=c(0,1),ylim=c(0,1))
lines(c(0,1),c(0,1),col="red")
plot(beta.soc2.red$Brepl,c(beta.dend2$Brepl),xlim=c(0,1),ylim=c(0,1))
lines(c(0,1),c(0,1),col="red")


plot(c(beta.trait$Btotal),c(beta.dend1$Btotal),xlim=c(0,1),ylim=c(0,1))
lines(c(0,1),c(0,1),col="red")
plot(c(beta.trait$Brepl),c(beta.dend1$Brepl),xlim=c(0,1),ylim=c(0,1))
lines(c(0,1),c(0,1),col="red")

plot(c(beta.trait$Btotal),c(beta.dend2$Btotal),xlim=c(0,1),ylim=c(0,1))
lines(c(0,1),c(0,1),col="red")
plot(c(beta.trait$Brepl),c(beta.dend2$Brepl),xlim=c(0,1),ylim=c(0,1))
lines(c(0,1),c(0,1),col="red")


cor(c(beta.dend1$Btotal),c(beta.dend2$Btotal))
cor(c(beta.dend1$Brepl),c(beta.dend2$Brepl))

cor(beta.soc2.red$Btotal,c(beta.dend1$Btotal))
cor(beta.soc2.red$Brepl,c(beta.dend1$Brepl))

cor(beta.soc2.red$Btotal,c(beta.dend2$Btotal))
cor(beta.soc2.red$Brepl,c(beta.dend2$Brepl))


cor(c(beta.trait$Btotal),c(beta.dend1$Btotal))
cor(c(beta.trait$Brepl),c(beta.dend1$Brepl))

cor(c(beta.trait$Btotal),c(beta.dend2$Btotal))
cor(c(beta.trait$Brepl),c(beta.dend2$Brepl))


#####################
##Species phylogeny##
#####################

# library(BiocManager)
# BiocManager::install("ggtree")
library(ggtree)
library(ape)


PaciFlora.tree <- read.tree("PaciFlora/Smith_Brown_plus_Wohlwend_tree.tre")

#plot(PaciFlora.tree)
ggtree(PaciFlora.tree,layout="circular")# + geom_tiplab(color='firebrick')


names.spp.soc <- unique(sub(" ", "_", dat.soc$species))
Soc.tree <- keep.tip(PaciFlora.tree,tip=which(PaciFlora.tree$tip.label %in% names.spp.soc))
ggtree(Soc.tree,layout="circular")# + geom_tiplab(color='firebrick') 




# Soc.tree.islands <- list()
# alpha <- numeric(length(unique(dat.soc$island)))
# ii <- 0
# for(i in sort(unique(dat.soc$island))){
#   ii <- ii+1
#   spp.temp <- sub(" ", "_", dat.soc$species[which(dat.soc$island==i)])
#   Soc.tree.islands[[ii]] <- keep.tip(Soc.tree,tip=which(Soc.tree$tip.label %in% spp.temp))
#   alpha[ii] <- sum(Soc.tree.islands[[ii]]$edge.length)
# }



dat.soc.pa.tree <- dat.soc.pa
colnames(dat.soc.pa.tree) <- sub(" ", "_",colnames(dat.soc.pa.tree))
dat.soc.pa.tree <- dat.soc.pa.tree[,which(colnames(dat.soc.pa.tree) %in% Soc.tree$tip.label)]

gamma.phylo <- gamma(comm=dat.soc.pa.tree,tree=Soc.tree)
alpha.phylo <- alpha(comm=dat.soc.pa.tree,tree=Soc.tree)
alpha.phylo.ratio <- alpha.phylo/c(gamma.phylo)
alpha.phylo.mean <- mean(alpha.phylo)
alpha.phylo.mean.ratio <- alpha.phylo.mean/c(gamma.phylo)
beta.phylo <- beta(comm=dat.soc.pa.tree,tree=Soc.tree)


dat <- data.frame(island=names(alpha.soc.ratio),alpha.tax=alpha.soc.ratio,alpha.phylo=c(alpha.phylo.ratio))
ggplot(dat,aes(alpha.tax,alpha.phylo,label=island))+
  geom_point()+
  geom_text_repel()+
  geom_abline(intercept=0,slope=1,color="red")+
  xlim(0,1)+
  ylim(0,1)

cor(alpha.soc.ratio,c(alpha.phylo.ratio))
cor(alpha.soc.ratio,c(alpha.phylo.ratio),method="spearman")


plot(beta.soc2$Btotal,(beta.phylo$Btotal),xlim=c(0,1),ylim=c(0,1))
lines(c(0,1),c(0,1),col="red")
plot(beta.soc2$Brepl,(beta.phylo$Brepl),xlim=c(0,1),ylim=c(0,1))
lines(c(0,1),c(0,1),col="red")


beta.phylo.temp <- lapply(beta.phylo,as.matrix)
beta.phylo.red <- list()
beta.phylo.red$Btotal <- c(as.dist(beta.phylo.temp$Btotal[which(row.names(beta.phylo.temp$Btotal) %in% names(beta.trait$Btotal)),which(colnames(beta.phylo.temp$Btotal) %in% names(beta.trait$Btotal))]))
beta.phylo.red$Brepl <- c(as.dist(beta.phylo.temp$Brepl[which(row.names(beta.phylo.temp$Brepl) %in% names(beta.trait$Brepl)),which(colnames(beta.phylo.temp$Brepl) %in% names(beta.trait$Brepl))]))
beta.phylo.red$Brich <- c(as.dist(beta.phylo.temp$Brich[which(row.names(beta.phylo.temp$Brich) %in% names(beta.trait$Brich)),which(colnames(beta.phylo.temp$Brich) %in% names(beta.trait$Brich))]))


plot(c(beta.phylo.red$Btotal),c(beta.trait$Btotal),xlim=c(0,1),ylim=c(0,1))
lines(c(0,1),c(0,1),col="red")
plot(c(beta.phylo.red$Brepl),c(beta.trait$Brepl),xlim=c(0,1),ylim=c(0,1))
lines(c(0,1),c(0,1),col="red")


##Null models
reps <- 99
gamma.phylo.rand <- alpha.phylo.mean.rand <- alpha.phylo.mean.ratio.rand <- numeric(reps)
alpha.phylo.rand <- alpha.phylo.ratio.rand <- matrix(NA,reps,length(alpha.phylo))
beta.phylo.total.rand <- beta.phylo.repl.rand <- matrix(NA,reps,length(c(beta.phylo$Btotal)))
for(i in 1:reps){
  print(i)
  dat.soc.pa.tree.temp <- dat.soc.pa.tree
  names(dat.soc.pa.tree.temp) <- names(dat.soc.pa.tree.temp)[sample(ncol(dat.soc.pa.tree.temp),ncol(dat.soc.pa.tree.temp))]
  gamma.phylo.rand[i] <- c(gamma(comm=dat.soc.pa.tree.temp,tree=Soc.tree))
  alpha.phylo.rand[i,] <- alpha(comm=dat.soc.pa.tree.temp,tree=Soc.tree)
  alpha.phylo.ratio.rand[i,] <- alpha.phylo.rand[i,]/c(gamma.phylo.rand[i])
  alpha.phylo.mean.rand[i] <- mean(alpha.phylo.rand[i,])
  alpha.phylo.mean.ratio.rand[i] <- alpha.phylo.mean.rand[i]/c(gamma.phylo.rand[i])
  toto <- beta(comm=dat.soc.pa.tree.temp,tree=Soc.tree)
  beta.phylo.total.rand[i,] <- c(toto$Btotal)
  beta.phylo.repl.rand[i,] <- c(toto$Brepl)
}

mean(alpha.phylo)
quantile(c(mean(alpha.phylo),rowMeans(alpha.phylo.rand)),c(0.025,0.05,0.95,0.975))

par(mfrow=c(4,4))
for(i in 1:14){
  hist(alpha.phylo.rand[,i],xlim=c(min(c(alpha.phylo.rand[,i],alpha.phylo[i]))*0.9,max(c(alpha.phylo.rand[,i],alpha.phylo[i]))*1.1),main=row.names(alpha.phylo)[i])
  abline(v=mean(alpha.phylo[i]),col="red")
}

mean(beta.phylo$Btotal)
quantile(c(mean(beta.phylo$Btotal),rowMeans(beta.phylo.total.rand)),c(0.025,0.05,0.95,0.975))

mean(beta.phylo$Brepl)
quantile(c(mean(beta.phylo$Brepl),rowMeans(beta.phylo.repl.rand)),c(0.025,0.05,0.95,0.975))

par(mfrow=c(1,3))

hist(rowMeans(alpha.phylo.rand),xlim=c(5800,7300))
abline(v=mean(alpha.phylo),col="red")

hist(rowMeans(beta.phylo.total.rand),xlim=c(0.65,0.75))
abline(v=mean(beta.phylo$Btotal),col="red")

hist(rowMeans(beta.phylo.repl.rand),xlim=c(0.1,0.2))
abline(v=mean(beta.phylo$Brepl),col="red")



par(mfrow=c(1,2))
dat <- data.frame(island=names(alpha.soc.ratio),alpha.tax=alpha.soc.ratio,alpha.phylo=c(alpha.phylo.ratio))
datt <- data.frame(island=rep(names(alpha.soc.ratio),each=99),alpha.tax=rep(alpha.soc.ratio,each=99),alpha.phylo=c(alpha.phylo.ratio.rand))
ggplot()+
  geom_point(datt,mapping=aes(alpha.tax,alpha.phylo))+
  geom_point(dat,mapping=aes(alpha.tax,alpha.phylo),colour="red")+
  geom_label_repel(dat,mapping=aes(alpha.tax,alpha.phylo,label=island),direction = "x")+
  geom_abline(intercept=0,slope=1,color="red")+
  xlim(0,1)+
  ylim(0,1)


plot(rep(alpha.soc.ratio,each=99),c(alpha.phylo.ratio.rand))
points(alpha.soc.ratio,alpha.phylo.ratio,pch=16,col="red")
abline(a=1,b=0,col="red")
