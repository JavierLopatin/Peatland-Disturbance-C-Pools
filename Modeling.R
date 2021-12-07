
################################################################################
## R-Script: 1_Modeling                                                       ##
## author: Javier Lopatin                                                     ##
## mail: javierlopatin@gmail.com                                              ##
##                                                                            ##
## Manuscript: Disturbance alters relationships between soil carbon pools and ##
##             aboveground vegetation attributes in an anthropogenic peatland ##
##             in Patagonia                                                   ##
##                                                                            ##
##                                                                            ##
## description: This R-code provide the Modeling approach used in the         ##
## manuscript submitted to Ecology and Evolution                              ##
##                                                                            ##
################################################################################

library(isopam)
library(vegan)
library(MASS)
library(ggpubr)
library(rstatix)
library(plspm)
library(plspm.formula)


# funtion to source scripts from GitHub
source_github <- function(u) {
  # load package
  require(RCurl)
  # read script lines from website and evaluate
  script <- getURL(u, ssl.verifypeer = FALSE)
  eval(parse(text = script), envir=.GlobalEnv)
  detach("package:RCurl", unload=TRUE)
} 

setwd('~/Documentos/temp/Peatland_Managements_difference/')

load('management.RData')

# -------------------------------------------
### Load functions
# -------------------------------------------

# Get residuals from PLS-PM object
source_github('https://raw.githubusercontent.com/JavierLopatin/plspmSpatial/master/plspmPredict.R')


# --------------------------------------------------------------
## NMDS CARBON
# --------------------------------------------------------------

nmds.stress <- sapply(1:6, function(x) metaMDS(carbons, k=x)$stress)
plot(1:6, nmds.stress)

pc=metaMDS(carbons, k = 2)

fit = envfit(pc, carbons)
fit2 = envfit(pc, cbind(rich, scores(nmds1)[,1], scores(nmds1)[,2], scores(nmds2)[,1]))

# --------------------------------------------------------------
## NMDS species and PFT level
# --------------------------------------------------------------

# import data
setwd("~/Documentos/PhD/Peatland/")

data <- read.table("data/Peatland1.csv", header=T, sep=",", dec=".")
names(data)

# floristic cover data
sp <- read.table("data/Cover_spp.csv", header=T, sep=",", dec=".")
sp2 <- sp[, -c(1,51,52)]
summary(sp2)

# pft cover data
pft <- read.table("data/PFT1.csv", header=T, sep=",", dec=".")
pft2 <- pft[, -c(1) ]
summary(pft2)

# Applying isopam clustering to check for species aggrupations
ip <- isopam(sp2, c.fix=3)
isotab (ip, 2)
## examine grouping
ip$flat

### To species level

# Selecting the number of k (stress value)
nmds.stress <- sapply(1:6, function(x) metaMDS(sp2, k=x)$stress)
plot(1:6, nmds.stress)

# Selecting a number of dimensions: compromise between as few dimensions and low stress value
nmds1 <- metaMDS(sp2, k=3, trymax=100)
nmds1$stress

## plot
ordiplot(nmds1, choices=c(1,2))
ordiplot(nmds1, choices=c(1,3))
ordiplot(nmds1, choices=c(2,3))

# see order of species on the gradients
ord = data.frame(s=scores(nmds1)[,1], id=ip2$flat)
ord[order(ord$s, decreasing = F),]

### To PFT level

# Applying isopan clustering to check for species agrupations
ip2 <- isopam(pft2, c.fix=3)
isotab(ip2, 2)
ip2$flat

nmds.stress <- sapply(1:6, function(x) metaMDS(pft2, k=x)$stress)
plot(1:6, nmds.stress)
nmds2 <- metaMDS(pft2, k=2, trymax=100)
nmds2$stress
ordiplot(nmds2, choices=c(1,2))

ord = data.frame(s=scores(nmds2)[,1], id=ip2$flat)
ord[order(ord$s, decreasing = F),]

# --------------------------------------------------------------
## Figure 3 and A1
# --------------------------------------------------------------

## fit environmental vectors
carbons = data[, 21:26]
colnames(carbons) = c('Fine', 'Moss', 'R1', 'R2', 'R3', 'Gross')

rich <- data.frame(Shannon_spp=data$shannon,
                   Simpson_spp=data$simpson,
                   Shannon_PFT=data$shannon_pft,
                   Simpson_PFT=data$simpson_pft,
                   H = log(data$Altura_vegetacion_cm),
                   BM_herb=data$Biomasa_herbaceas_kg_m2,
                   BM_shrb=data$Biomasa_arbustivas_kg_m2)

# fit
fit1 = envfit(nmds1, carbons)
fit2 = envfit(nmds1, rich)

conservation = grep("Conservacion", data$Uso)
productive = grep("Productivo", data$Uso)

#### Species-based ordination
d = data.frame(ip=ip$flat, uso=data$Uso)

svg(file = "Managements_difference/Figures2/SP_12_C.svg", width = 5, height = 4.5)
par(mar=c(4,5,1,2))
fig <- plot(nmds1, type = "none", xlab = "NMDS 1", ylab = "NMDS 2", bty = "l", xaxs = "i",
            yaxs = "i",  ylim = c(-1.5, 1), cex.lab = 1.5, las = 1, cex.axis = 1.5)
points(fig$sites[which(d$ip==1 & d$uso=='Conservacion'),], pch = 21, cex = 2, col='black',
       bg=rgb(170,135,222,80, maxColorValue=255))
points(fig$sites[which(d$ip==1 & d$uso=='Productivo'),], pch = 21, cex = 2, col='black',
       bg=rgb(128,128,0,80, maxColorValue=255))
points(fig$sites[which(d$ip==2 & d$uso=='Conservacion'),], pch = 24, cex = 2, col='black',
       bg=rgb(170,135,222,80, maxColorValue=255))
points(fig$sites[which(d$ip==2 & d$uso=='Productivo'),], pch = 24, cex = 2, col='black',
       bg=rgb(128,128,0,80, maxColorValue=255))
points(fig$sites[which(d$ip==3 & d$uso=='Conservacion'),], pch = 22, cex = 2, col='black',
       bg=rgb(170,135,222,80, maxColorValue=255))
points(fig$sites[which(d$ip==3 & d$uso=='Productivo'),], pch = 22, cex = 2, col='black',
       bg=rgb(128,128,0,80, maxColorValue=255))
plot(fit1, col = "red", lty = 2, cex = 1.5)
dev.off()

svg(file = "Managements_difference/Figures2/SP_12_Rich.svg", width = 5, height = 4.5)
fig <- plot(nmds1, type = "none", xlab = "NMDS 1", ylab = "NMDS 2", bty = "l", xaxs = "i",
            yaxs = "i",  ylim = c(-1.5, 1), cex.lab = 1.5, las = 1, cex.axis = 1.5)
points(fig$sites[which(d$ip==1 & d$uso=='Conservacion'),], pch = 21, cex = 2, col='black',
       bg=rgb(170,135,222,80, maxColorValue=255))
points(fig$sites[which(d$ip==1 & d$uso=='Productivo'),], pch = 21, cex = 2, col='black',
       bg=rgb(128,128,0,80, maxColorValue=255))
points(fig$sites[which(d$ip==2 & d$uso=='Conservacion'),], pch = 24, cex = 2, col='black',
       bg=rgb(170,135,222,80, maxColorValue=255))
points(fig$sites[which(d$ip==2 & d$uso=='Productivo'),], pch = 24, cex = 2, col='black',
       bg=rgb(128,128,0,80, maxColorValue=255))
points(fig$sites[which(d$ip==3 & d$uso=='Conservacion'),], pch = 22, cex = 2, col='black',
       bg=rgb(170,135,222,80, maxColorValue=255))
points(fig$sites[which(d$ip==3 & d$uso=='Productivo'),], pch = 22, cex = 2, col='black',
       bg=rgb(128,128,0,80, maxColorValue=255))
plot(fit2, col = "red", lty = 2, cex = 1.5)
dev.off()

#### PFT-based ordination
svg(file = "Managements_difference/Figures2/PFT_12_C.svg", width = 5, height = 4.5)
par(mar=c(4,5,1,2))
fit2 <- envfit(nmds2, rich)
fig <- plot(nmds2, type = "none", xlab = "NMDS 1", ylab = "NMDS 2", bty = "l", xaxs = "i",
            yaxs = "i", cex.lab = 1.5, las = 1, ylim=c(-1.5,1), cex.axis = 1.5)
points(fig$sites[which(d$ip==1 & d$uso=='Conservacion'),], pch = 21, cex = 2, col='black',
       bg=rgb(170,135,222,80, maxColorValue=255))
points(fig$sites[which(d$ip==1 & d$uso=='Productivo'),], pch = 21, cex = 2, col='black',
       bg=rgb(128,128,0,80, maxColorValue=255))
points(fig$sites[which(d$ip==2 & d$uso=='Conservacion'),], pch = 24, cex = 2, col='black',
       bg=rgb(170,135,222,80, maxColorValue=255))
points(fig$sites[which(d$ip==2 & d$uso=='Productivo'),], pch = 24, cex = 2, col='black',
       bg=rgb(128,128,0,80, maxColorValue=255))
points(fig$sites[which(d$ip==3 & d$uso=='Conservacion'),], pch = 22, cex = 2, col='black',
       bg=rgb(170,135,222,80, maxColorValue=255))
points(fig$sites[which(d$ip==3 & d$uso=='Productivo'),], pch = 22, cex = 2, col='black',
       bg=rgb(128,128,0,80, maxColorValue=255))
plot(fit2, col = "red", lty = 2, cex = 1.5)
legend("bottomright", legend = c("class 1", "class 2", "class 3"), pch = c(21, 25,23), cex = 1.3,
       bty = "n")
dev.off()

save.image('ordination.RData')


# --------------------------------------------------------------
## Differences in C pools between managements
# --------------------------------------------------------------

data$Uso = as.factor(data$Uso)
data$MDS_C1

# plot groups
ggboxplot(data, x = "Uso", y = "MDS_C1", add = "jitter")
ggboxplot(data, x = "Uso", y = "MDS_C2", add = "jitter")

# ANOVA y Tukey
# C Comp. 1
aov1 = anova_test(MDS_C1 ~ Uso, data=data); aov1
tuk1 = tukey_hsd(data, MDS_C1 ~ Uso); tuk1
# C Comp. 2
aov2 = anova_test(MDS_C2 ~ Uso, data=data); aov2
tuk2 = tukey_hsd(data, MDS_C2 ~ Uso); tuk2

# agregar coordenadas y posiciones de los valores. Solo necesario para hacer graficos
tuk1 = add_xy_position(tuk1, x = "Uso")
tuk2 = add_xy_position(tuk2, x = "Uso")

# plot
svg('ANOVA_C1.svg',width = 5, height = 4.5)
ggboxplot(data, x = "Uso", y = "MDS_C1", add = "jitter", xlab='Management', ylab='C Comp. 1') +
  stat_pvalue_manual(tuk1, hide.ns = TRUE) +
  labs(
    subtitle = get_test_label(aov1,  detailed = TRUE),
    caption = get_pwc_label(tuk1)
  )
dev.off()

svg('ANOVA_C2.svg', width = 4, height = 3.5)
ggboxplot(data, x = "Uso", y = "MDS_C2", add = "jitter", xlab='Management', ylab='C Comp. 2') +
  stat_pvalue_manual(tuk2, hide.ns = TRUE) +
  labs(
    subtitle = get_test_label(aov2,  detailed = TRUE),
    caption = get_pwc_label(tuk2)
  )

dev.off()

# --------------------------------------------------------------
## Figure 4
# --------------------------------------------------------------

svg(file = "Figures/Supp_MDS_C.svg", width = 5, height = 4.5)
par(mar=c(4,5,1,2))
d = data.frame(ip=ip2$flat, uso=data$Uso)
fig <- plot(pc, type = "none", xlab = c("NMDS 1"), ylab = ("NMDS 2"), 
            bty = "l", cex.lab = 1.5, las = 1, cex.axis = 1.5)
points(fig$sites[which(d$ip==1 & d$uso=='Conservacion'), 1:2], pch = 21, cex = 2, col='black',
       bg=rgb(170,135,222,80, maxColorValue=255))
points(fig$sites[which(d$ip==1 & d$uso=='Productivo'), 1:2], pch = 21, cex = 2, col='black',
       bg=rgb(128,128,0,80, maxColorValue=255))
points(fig$sites[which(d$ip==2 & d$uso=='Conservacion'),], pch = 24, cex = 2, col='black',
       bg=rgb(170,135,222,80, maxColorValue=255))
points(fig$sites[which(d$ip==2 & d$uso=='Productivo'),], pch = 24, cex = 2, col='black',
       bg=rgb(128,128,0,80, maxColorValue=255))
points(fig$sites[which(d$ip==3 & d$uso=='Conservacion'),], pch = 22, cex = 2, col='black',
       bg=rgb(170,135,222,80, maxColorValue=255))
points(fig$sites[which(d$ip==3 & d$uso=='Productivo'),], pch = 22, cex = 2, col='black',
       bg=rgb(128,128,0,80, maxColorValue=255))
plot(fit, col = "red", lty = 2, cex = 1.5)
legend("bottomleft", legend = c("Conservation", "Agriculture"), fill=c(rgb(170,135,222,80, maxColorValue=255),rgb(128,128,0,80, maxColorValue=255)), cex = 1.3,
       bty = "n")
dev.off()

svg(file = "Managements_difference/Figures2/Supp_MDS_plspm.svg", width = 5, height = 4.5)
par(mar=c(4,5,1,2))
fig <- plot(scores(pc)[, 1:2], type = "none", xlab = c("NMDS 1"), ylab = ("NMDS 2"), 
            bty = "l", cex.lab = 1.5, las = 1, cex.axis = 1.5, xlim=c(min(scores(pc)),max(scores(pc))+0.1))
points(scores(pc)[which(d$ip==1 & d$uso=='Conservacion'), 1:2], pch = 21, cex = 2, col='black',
       bg=rgb(170,135,222,80, maxColorValue=255))
points(scores(pc)[which(d$ip==1 & d$uso=='Productivo'), 1:2], pch = 21, cex = 2, col='black',
       bg=rgb(128,128,0,80, maxColorValue=255))
points(scores(pc)[which(d$ip==2 & d$uso=='Conservacion'),], pch = 24, cex = 2, col='black',
       bg=rgb(170,135,222,80, maxColorValue=255))
points(scores(pc)[which(d$ip==2 & d$uso=='Productivo'),], pch = 24, cex = 2, col='black',
       bg=rgb(128,128,0,80, maxColorValue=255))
points(scores(pc)[which(d$ip==3 & d$uso=='Conservacion'),], pch = 22, cex = 2, col='black',
       bg=rgb(170,135,222,80, maxColorValue=255))
points(scores(pc)[which(d$ip==3 & d$uso=='Productivo'),], pch = 22, cex = 2, col='black',
       bg=rgb(128,128,0,80, maxColorValue=255))
plot(fit2, col = "red", lty = 2, cex = 1.5)
legend("bottomright", legend = c("class 1", "class 2", "class 3"), pch = c(1, 2,0), cex = 1.3,
       bty = "n")
dev.off()


# --------------------------------------------------------------
## PLS-PM modeling
# --------------------------------------------------------------

data$MDS_C1 <- scores(pc)[,1]
data$MDS_C2 <- scores(pc)[,2]
data$FC1    <- scores(nmds1)[,1]
data$FC2    <- scores(nmds1)[,2]
data$FC3    <- scores(nmds1)[,3]
data$PFT1   <- scores(nmds2)[,1]
data$PFT2   <- scores(nmds2)[,2]

# ------------------------------------------
## Select floristic and PFT components
# -----------------------------------------
fit1 <- lm(MDS_C1 ~ FC1+FC2+FC3, data=data)
summary(fit1)
fit1$coefficients[2]

fit2 <- lm(MDS_C2 ~ FC1+FC2+FC3, data=data)
summary(fit2)

fit3 <- lm(MDS_C1 ~ PFT1+PFT2, data=data)
summary(fit3)

fit4 <- lm(MDS_C2 ~ PFT1+PFT2, data=data)
summary(fit4)


### Inner formula to predict C Comp. 1
formulaMDS1 <-"
  FC1 =~ FC1
  FC2 =~ FC2
  PFT =~ PFT1
  H =~ Altura_vegetacion_cm
  BM =~ Biomasa_herbaceas_kg_m2 + Biomasa_arbustivas_kg_m2
  DV =~ shannon + simpson + shannon_pft + simpson_pft
  C =~ MDS_C1

  H  ~~ FC1 + FC2 + PFT
  BM ~~ FC1 + FC2 + PFT + H 
  DV ~~ FC1 + FC2 + PFT + H
  SD ~~ FC1 + FC2 + PFT + BM + DV
  C ~~ FC1 + FC2 + PFT + H + BM + DV 
  "

# predict overall model
set.seed(123)
plsPCA1 = plspm.formula(formulaMDS1, data, modes=rep('A',7), maxiter= 500, boot.val =T, br = 500,
                        scheme = "factor", scaled = T)
plsPCA1$outer_model

# predict undisturbed model
set.seed(123)
plsPCA1_cons = plspm.formula(formulaMDS1, data[conservation,], modes=rep('A',7), maxiter= 500, boot.val =F, br = 500,
                             scheme = "factor", scaled = T)

# predict disturbed model 
set.seed(123)
plsPCA1_agri = plspm.formula(formulaMDS1, data[productive,], modes=rep('A',7), maxiter= 500, boot.val = T, br = 500,
                             scheme = "factor", scaled = T)

plsPCA1$gof; plsPCA1$boot$rsq 
plsPCA1_cons$gof; plsPCA1_cons$boot$rsq 
plsPCA1_agri$gof; plsPCA1_agri$boot$rsq 

### Inner formula to predict C Comp. 2
formulaMDS2 <-"
  FC1 =~ FC1
  FC2 =~ FC2
  PFT =~ PFT1
  H =~ Altura_vegetacion_cm
  BM =~ Biomasa_herbaceas_kg_m2 + Biomasa_arbustivas_kg_m2
  DV =~ shannon + simpson + shannon_pft + simpson_pft
  C =~ MDS_C2

  H  ~~ FC1 + FC2 + PFT
  BM ~~ FC1 + FC2 + PFT + H 
  DV ~~ FC1 + FC2 + PFT + H
  SD ~~ FC1 + FC2 + PFT + BM + DV
  C ~~ FC1 + FC2 + PFT + H + BM + DV
  "

# predict overall model
set.seed(123)
plsPCA2 = plspm.formula(formulaMDS2, data, modes=rep('A',7), maxiter= 500, boot.val = T, br = 500,
                        scheme = "factor", scaled = T)
plsPCA2$outer_model

# predict undisturbed model
set.seed(123)
plsPCA2_cons = plspm.formula(formulaMDS2, data[conservation,], modes=rep('A',7), maxiter= 500, boot.val = F, br = 500,
                             scheme = "factor", scaled = T)

# predict distuerbed model
set.seed(123)
plsPCA2_agri = plspm.formula(formulaMDS2, data[productive,], modes=rep('A',7), maxiter= 500, boot.val = T, br = 500,
                             scheme = "factor", scaled = T)

plsPCA2$gof; plsPCA2$boot$rsq 
plsPCA2_cons$gof; plsPCA2_cons$boot$rsq 
plsPCA2_agri$gof; plsPCA2_agri$boot$rsq 


# function to get significances from plspm directly
sig.effects <- function(PLS){ 
  path = PLS$boot$paths
  path.out = matrix(nrow=nrow(path), ncol=2)
  colnames(path.out) <- c('Mean.Boot', 'sig. 0.05')
  rownames(path.out) <- rownames(path)
  path.out[,1] = path$Mean.Boot
  for (i in 1:nrow(path)){
    if (sign(path$perc.025[i]) == sign(path$perc.975[i])) path.out[i,2] = 'yes'
    if (sign(path$perc.025[i]) != sign(path$perc.975[i])) path.out[i,2] = 'no'
  }
  
  path = PLS$boot$total.efs
  path.out2 = matrix(nrow=nrow(path), ncol=2)
  colnames(path.out2) <- c('Mean.Boot', 'sig. 0.05')
  rownames(path.out2) <- rownames(path)
  path.out2[,1] = path$Mean.Boot
  for (i in 1:nrow(path)){
    if (sign(path$perc.025[i]) == sign(path$perc.975[i])) path.out2[i,2] = 'yes'
    if (sign(path$perc.025[i]) != sign(path$perc.975[i])) path.out2[i,2] = 'no'
  }
  
  out <- list(path.out, path.out2)
  names(out) <- c('path.coefficients', 'total.effects')
  return(out)
}

sig1_all  <- sig.effects(plsPCA1) 
sig1_cons <- sig.effects(plsPCA1_cons) 
sig1_agri <- sig.effects(plsPCA1_agri) 

sig2_all  <- sig.effects(plsPCA2) 
sig2_cons <- sig.effects(plsPCA2_cons) 
sig2_agri <- sig.effects(plsPCA2_agri) 

# check for the significant direct path coefficients
sig1_all$path.coefficients
sig2_all$path.coefficients

sig1_cons$path.coefficients
sig2_cons$path.coefficients

sig1_agri$path.coefficients
sig2_agri$path.coefficients

# check for the total effects on C components
sig1_all$total.effects
sig2_all$total.effects

sig1_cons$total.effects
sig2_cons$total.effects

sig1_agri$total.effects
sig2_agri$total.effects

##############
plsPCA1_cons$effects
plsPCA2_cons$effects

# ---------------------------------------------------------
# Chack for Spatial Autocorrelation in the models
# ---------------------------------------------------------

# check for spatial autocorrelation
library(ape) # Moran´s I

# function to get distance matrix of coordinates
dist_mat <- function(X, Y){
  xy.dists <- as.matrix(dist(cbind(X, Y)))
  xy.dists.inv <- 1/xy.dists # invers
  diag(xy.dists.inv) <- 0
  return (xy.dists.inv)
}

# obtain distance matrix of XY coordinates
dist_all = dist_mat(xy$x, xy$y)
dist_cons = dist_mat(xy$x[data$Uso == 'Conservacion'], xy$y[data$Uso == 'Conservacion'])
dist_agri = dist_mat(xy$x[data$Uso == 'Productivo'], xy$y[data$Uso == 'Productivo'])

# get model residuals
res_PLCA1 = plspmResiduals(plsPCA1)
res_PLCA1_cons = plspmResiduals(plsPCA1_cons)
res_PLCA1_agri = plspmResiduals(plsPCA1_agri)

res_PLCA2 = plspmResiduals(plsPCA2)
res_PLCA2_cons = plspmResiduals(plsPCA2_cons)
res_PLCA2_agri = plspmResiduals(plsPCA2_agri)

#### Moran's I

# Vegetation height
moran.H      <- Moran.I(res_PLCA1$inner_residuals[,1], dist_all)
moran.H_cons <- Moran.I(res_PLCA1_cons$inner_residuals[,1], dist_cons)
moran.H_agri <- Moran.I(res_PLCA1_agri$inner_residuals[,1], dist_agri)

# Aboveground biomass
moran.BM      <- Moran.I(res_PLCA1$inner_residuals[,2], dist_all)
moran.BM_cons <- Moran.I(res_PLCA1_cons$inner_residuals[,2], dist_cons)
moran.BM_agri <- Moran.I(res_PLCA1_agri$inner_residuals[,2], dist_agri)

# Species richness
moran.DV      <- Moran.I(res_PLCA1$inner_residuals[,3], dist_all)
moran.DV_cons <- Moran.I(res_PLCA1_cons$inner_residuals[,3], dist_cons)
moran.DV_agri <- Moran.I(res_PLCA1_agri$inner_residuals[,3], dist_agri)

# Belowground C stock
moran.C      <- Moran.I(res_PLCA1$inner_residuals[,4], dist_all)
moran.C_cons <- Moran.I(res_PLCA1_cons$inner_residuals[,4], dist_cons)
moran.C_agri <- Moran.I(res_PLCA1_agri$inner_residuals[,4], dist_agri)

moran.C2      <- Moran.I(res_PLCA2$inner_residuals[,4], dist_all)
moran.C2_cons <- Moran.I(res_PLCA2_cons$inner_residuals[,4], dist_cons)
moran.C2_agri <- Moran.I(res_PLCA2_agri$inner_residuals[,4], dist_agri)


# --------------------------------
# Plot variable interactions
# --------------------------------

setwd('~/Documentos/temp/Peatland_Managements_difference/')

Spointtype <- factor(data$Uso)

plot_manage <- function(x, y, xlab, ylab='First C component', managements=1, ...){ 
  library(mgcv) 
  if(managements==1){
    bg="mediumpurple4"
    id = which(Spointtype=='Conservacion')
  } else if(managements==2){
    bg="lavender"
    id = which(Spointtype=='Productivo')
  }
  dat = data.frame(x=x,y=y)
  plot(100, bty="l", xlim=c(min(x),max(x)), ylim=c(min(y),max(y)), las=1, xlab=xlab, ylab=ylab)
  points(dat[ip$flat[id]==1, 1:2], pch=21, col="black", bg=bg, cex = 1.5)
  points(dat[ip$flat[id]==2, 1:2], pch=24, col="black", bg=bg, cex = 1.5)
  points(dat[ip$flat[id]==3, 1:2], pch=22, col="black", bg=bg, cex = 1.5)
  lm <- loess(y ~ x)
  curve.dat = data.frame(x=x, y=predict(lm))
  curve.dat = curve.dat[order(curve.dat$x),]
  lines(curve.dat, lwd=2, lty=2)  
  r1 = round((cor(x,y,method = 'spearman')), 2)
  nRSE = round( (lm$s/(max(y)-min(y)))*100, 2)
  mtext(bquote(paste(r == .(r1),';',' ', nRSE == .(nRSE),'%')), side = 3, line = 0.5, adj = 0, cex = 1.3, font = 2)
}

#### First C component - conservation
y = plsPCA1_cons$scores[,7]
# C1 - FC1
svg(file = "Figures2/C_FC1.svg", width=4, height=4)
plot_manage(x=plsPCA1_cons$scores[,1], y, xlab = "Spp. Comp. 1")
dev.off()
# C1 - FC2
svg(file = "Figures2/C_FC2.svg", width=4, height=4)
plot_manage(plsPCA1_cons$scores[,2], y, xlab = "Spp. Comp. 2")
dev.off()
# C1 - PFT1
svg(file = "Figures2/C_PFT1.svg", width=4, height=4)
plot_manage(plsPCA1_cons$scores[,3], y, xlab = "PFT Comp. 1")
dev.off()
# C1 - Diversity
svg(file = "Figures2/C_Div.svg", width=4, height=4)
plot_manage(plsPCA1_cons$scores[,6], y, xlab = "Diversity")
dev.off()
# C1 - PFT1
svg(file = "Figures2/C_H.svg", width=4, height=4)
plot_manage(plsPCA1_cons$scores[,4], y, xlab = "Vegetation height")
dev.off()
# C1 - PFT1
svg(file = "Figures2/C_BM.svg", width=4, height=4)
plot_manage(plsPCA1_cons$scores[,5], y, xlab = "Biomass")
dev.off()

#### Second C component - conservation
y = plsPCA2_cons$scores[,7]
# C1 - FC1
svg(file = "Figures2/C2_FC1.svg", width=4, height=4)
plot_manage(plsPCA2_cons$scores[,1], y, xlab = "Spp. Comp. 1", ylab='Second C component')
dev.off()
# C1 - FC2
svg(file = "Figures2/C2_FC2.svg", width=4, height=4)
plot_manage(plsPCA2_cons$scores[,2], y, xlab = "Spp. Comp. 2", ylab='Second C component')
dev.off()
# C1 - PFT1
svg(file = "Figures2/C2_PFT1.svg", width=4, height=4)
plot_manage(plsPCA2_cons$scores[,3], y, xlab = "PFT Comp. 1", ylab='Second C component')
dev.off()
# C1 - Diversity
svg(file = "Figures2/C2_Div.svg", width=4, height=4)
plot_manage(plsPCA2_cons$scores[,6], y, xlab = "Diversity", ylab='Second C component')
dev.off()
# C1 - H
svg(file = "Figures2/C2_H.svg", width=4, height=4)
plot_manage(plsPCA2_cons$scores[,4], y, xlab = "Vegetation height", ylab='Second C component')
dev.off()
# C1 - BM
svg(file = "Figures2/C2_BM.svg", width=4, height=4)
plot_manage(plsPCA2_cons$scores[,5], y, xlab = "Biomass", ylab='Second C component')
dev.off()


#### First C component - Agriculture
y = plsPCA1_agri$scores[,7]
# C1 - FC1
svg(file = "Figures2/C_FC1_agri.svg", width=4, height=4)
plot_manage(plsPCA1_agri$scores[,1], y, xlab = "Spp. Comp. 1", managements = 2)
dev.off()
# C1 - FC2
svg(file = "Figures2/C_FC2_agri.svg", width=4, height=4)
plot_manage(x=plsPCA1_agri$scores[,2], y, xlab = "Spp. Comp. 2", managements = 2)
dev.off()
# C1 - PFT1
svg(file = "Figures2/C_PFT1_agri.svg", width=4, height=4)
plot_manage(plsPCA1_agri$scores[,3], y, xlab = "PFT Comp. 1", managements = 2)
dev.off()
# C1 - Diversity
svg(file = "Figures2/C_Div_agri.svg", width=4, height=4)
plot_manage(plsPCA1_agri$scores[,6], y, xlab = "Diversity", managements = 2)
dev.off()
# C1 - PFT1
svg(file = "Figures2/C_H_agri.svg", width=4, height=4)
plot_manage(plsPCA1_agri$scores[,4], y, xlab = "Vegetation height", managements = 2)
dev.off()
# C1 - PFT1
svg(file = "Figures2/C_BM_agri.svg", width=4, height=4)
plot_manage(plsPCA2_agri$scores[,5], y, xlab = "Biomass", managements = 2)
dev.off()

#### Second C component - Agriculture
y = plsPCA1_agri$scores[,7]
# C1 - FC1
svg(file = "Figures2/C2_FC1_agri.svg", width=4, height=4)
plot_manage(plsPCA2_agri$scores[,1], y, xlab = "Spp. Comp. 1" , ylab='Second C component',managements = 2)
dev.off()
# C1 - FC2
svg(file = "Figures2/C2_FC2_agri.svg", width=4, height=4)
plot_manage(plsPCA2_agri$scores[,2], y, xlab = "Spp. Comp. 2", ylab='Second C component', managements = 2)
dev.off()
# C1 - PFT1
svg(file = "Figures2/C2_PFT1_agri.svg", width=4, height=4)
plot_manage(plsPCA2_agri$scores[,3], y, xlab = "PFT Comp. 1", ylab='Second C component', managements = 2)
dev.off()
# C1 - Diversity
svg(file = "Figures2/C2_Div_agri.svg", width=4, height=4)
plot_manage(plsPCA2_agri$scores[,6], y, xlab = "Diversity", ylab='Second C component' ,managements = 2)
dev.off()
# C1 - PFT1
svg(file = "Figures2/C2_H_agri.svg", width=4, height=4)
plot_manage(plsPCA2_agri$scores[,4], y, xlab = "Vegetation height", ylab='Second C component', managements = 2)
dev.off()
# C1 - PFT1
svg(file = "Figures2/C2_BM_agri.svg", width=4, height=4)
plot_manage(plsPCA2_agri$scores[,5], y, xlab = "Biomass", ylab='Second C component', managements = 2)
dev.off()


# ------------------------------------------------------------------------
# Significant diff, between models with floristic  and PFT compositions
# ------------------------------------------------------------------------

SifTest <- function(data1=plsPCA2$scores){
  fit1 <- lm(C ~., data=as.data.frame(data1))
  fit2 <- lm(C ~ H+BM+DV+PFT, data=as.data.frame(data1) )
  fit3 <- lm(C ~ H+BM+DV+FC1+FC2, data=as.data.frame(data1) )
  
  a1 = anova(fit1, fit2)
  a2 = anova(fit1, fit3)
  aic1 = AIC(fit1, fit2)
  aic2 = AIC(fit1, fit3)
  
  out <- list(a1$F[2], a2$F[2], a1$`Pr(>F)`[2], a2$`Pr(>F)`[2], aic1$AIC, aic2$AIC)
  names(out) <- c('Fvalue_spp', 'Fvalue_pft', 'Pvalue_spp', 'Pvalue_pft', 'AIC_spp', 'AIC_pft')
  out
}

# C1-overall
SifTest(plsPCA1$scores)
# C1-Cons
SifTest(plsPCA1_cons$scores)
# C1-Agri
SifTest(plsPCA1_agri$scores)

# C2-overall
SifTest(plsPCA2$scores)
# C2-Cons
SifTest(plsPCA2_cons$scores)
# C2-Agri
SifTest(plsPCA2_agri$scores)

save.image('management.RData')


