## Consensus RDA ##
# *************** #

# Useful packages:
# ****************

library(adespatial)
library(usdm)
library(ordiconsensus)

# Data input:
# ***********
# Total abundance below which a species is not considered: 
# We consider rare species too (interesting for the consensus-RDA)
ab <- 2

C <- read.table("spa.txt", h = T, sep = "\t", row.names = 1)

# Some quadrats must be removed because they are empty (when ab <- 15):
# C <- C[-c(60, 123, 166, 175, 198, 220, 242, 245, 274, 303, 334, 354, 412, 546, 548, 640), ]
# If ab != 15, check that C and Y have the same number and identity of rows.

env <- read.table("env.txt", h = T, sep = "\t", row.names = 1)

vif(env)
vifcor(env, th = 0.7)
env <- env[, -c(4, 6, 7, 9, 14, 20, 24)]   # Pour cor < 70 %

spe <- read.table("spe-vivant.txt", h = T, sep = "\t", row.names = 1)
sum <- c() ; for(i in 1:ncol(spe)) sum <- c(sum, sum(spe[,i]))
spe <- spe[, which(sum >= ab)]

Y <- spe

# coeffCompare() construit le MST des différentes "association coefficients" --> permet de 
# décider quels coefficients on garde
# consensusRDA() à utiliser après avoir utilisé coeffCompare() pour effectuer la consRDA
# sur les rda des coefficients sélectionnés

### Construct results object for the rda on abundances
ndis_abun <- 10
ordiRes_abun <- vector("list", length = ndis_abun)

#---------------------------------------------
### Perform the various constrained ordination
#---------------------------------------------
### RDA species profile
sp <- Y / apply(Y, 1, sum)
ordiRes_abun[[1]] <- rda(sp ~., data = env)

### RDA chord
chord <- Y/sqrt(apply(Y^2, 1, sum))
ordiRes_abun[[2]] <- rda(chord ~., data = env)

### RDA Hellinger
hell <- decostand(Y, method = "hellinger")
ordiRes_abun[[3]] <- rda(hell ~., data = env)

### RDA chi2
chisq <- decostand(Y, method = "chi.square")
ordiRes_abun[[4]] <- rda(chisq ~., data = env)

### db-RDA Bray-Curtis
ordiRes_abun[[5]] <- capscale(sqrt(vegdist(Y, method = "bray")) ~., data = env, comm = Y)

### db-RDA square-root Bray-Curtis
ordiRes_abun[[6]] <- capscale(sqrt(vegdist(Y^0.5, method = "bray")) ~., data = env, comm = Y)

### db-RDA fourth-root Bray-Curtis
ordiRes_abun[[7]] <- capscale(sqrt(vegdist(Y^0.25, method = "bray")) ~., data = env, comm = Y)

### db-RDA modified Gower log 2
ordiRes_abun[[8]] <- capscale(vegdist(decostand(Y, "log", logbase = 2), "altGower") ~., 
                              data = env, comm = Y) ### Warning message stem from 
# log transformation of 0

### db-RDA modified Gower log 5
ordiRes_abun[[9]] <- capscale(vegdist(decostand(Y, "log", logbase = 5), "altGower") ~., 
                              data = env, comm = Y) ### Warning message stem from 
# log transformation of 0

### db-RDA modified Gower log 10
ordiRes_abun[[10]] <- capscale(vegdist(decostand(Y, "log", logbase = 10), "altGower") ~., 
                               data = env, comm = Y) ### Warning message stem from 
# log transformation of 0

### Construct results object for the rda on the presence-absence data
ndis_bin <- 6
ordiRes_bin <- vector("list", length = ndis_bin)
# Transformation of abundances in presence-absence data:
Ybin <- Y
for (i in 1:nrow(Y)) {
  for (j in 1:ncol(Y)) if (Y[i, j] != 0) Ybin[i, j] = 1
}

Y[1:3, 1:10]
Ybin[1:3, 1:10]

### RDA species profile - binary
spbin <- Ybin / apply(Ybin, 1, sum)
ordiRes_bin[[1]] <- rda(spbin ~., data = env)

### db-RDA Ochiai
library(ade4)
och <- dist.binary(Ybin, method = 7)
ordiRes_bin[[2]] <- capscale(och ~., data = env, comm = Ybin)

### db-RDA Raup-Crick
ordiRes_bin[[3]] <- capscale(sqrt(vegdist(Y, method = "raup", binary = T)) ~., data = env, 
                          comm = Ybin)

### RDA chi2
chisq <- decostand(Ybin, method = "chi.square")
ordiRes_bin[[4]] <- rda(chisq ~., data = env)

### db-RDA Jaccard
ordiRes_bin[[5]] <- capscale(sqrt(vegdist(Y, method = "jaccard", binary = T)) ~., data = env, 
                          comm = Ybin)

### db-RDA Soerensen (Bray-Curtis)
ordiRes_bin[[6]] <- capscale(sqrt(vegdist(Y, method = "bray", binary = T)) ~., data = env, 
                          comm = Ybin)

# Significant axes of each RDA:
rda_test_abun <- vector("list", ndis_abun)
for (i in 1:ndis_abun) rda_test_abun[[i]] <- anova.cca(ordiRes_abun[[i]], by = "axis")
nb_sign_ax_abun <- sapply(rda_test_abun, function (x) length(which(x$Pr <= 0.05)))

rda_test_bin <- vector("list", ndis_bin)
for (i in 1:ndis_bin) rda_test_bin[[i]] <- anova.cca(ordiRes_bin[[i]], by = "axis")
nb_sign_ax_bin <- sapply(rda_test_bin, function (x) length(which(x$Pr <= 0.05)))

### Compare association coefficients
AssoComp_abun <- coeffCompare(ordiRes_abun, nb_sign_ax_abun)
AssoComp_bin  <- coeffCompare(ordiRes_bin, nb_sign_ax_bin)

#---------------------------------------------
### Draw a graphic to visualize the comparison
#---------------------------------------------
### Name of association coefficient compared
name_abun <- c("Species profiles", "Chord", "Hellinger", "Chi2", "Bray-Curtis", 
               "(Bray-Curtis)^0.5", "(Bray-Curtis)^0.25", "mGowerlog2", "mGowerlog5", 
               "mGowerlog10")
name_bin <- c("Species profiles", "Ochiai", "Raup", "chi2", "Jaccard", "Soerensen")

plot(AssoComp_abun$mst, type = "t", labels = name, xlab = "", ylab = "", 
     main = "MST Sites scores")
plot(AssoComp_bin$mst, type = "t", labels = name, xlab = "", ylab = "", 
     main = "MST Sites scores")

# We remove Chi2 (abund et bin), Chord and Species profile (abund):
ordiRes_abun_final <- ordiRes_abun[-c(1, 2, 4)]
rda_test_abun_final <- rda_test_abun[-c(1, 2, 4)]
(ndis_abun <- length(ordiRes_abun_final))

ordiRes_bin_final <- ordiRes_bin[-4]
rda_test_bin_final <- rda_test_bin[-4]
(ndis_bin <- length(ordiRes_bin_final))

source("consensusRDA_DB.R")
consRDA_abun <- consensusRDA_DB(ordiRes_abun_final, rda_test_abun_final, Y, env, scaling = 2)
summary(consRDA_abun)
# How many axis are worth being considered:
sum_eig_abun <- sum(consRDA_abun$values)
ax_sel_abun <- which(as.numeric(consRDA_abun$values) >= mean(as.numeric(consRDA_abun$values)))
round(consRDA_abun$values[ax_sel_abun] / sum_eig_abun, 4)

axisLabels <- c(paste("Axis 1 - ", 
                      round(consRDA_abun$values[1]/sum(consRDA_abun$values), 4) * 100, 
                      sep = ""), 
                paste("Axis 2 - ", 
                      round(consRDA_abun$values[2]/sum(consRDA_abun$values), 4) * 100, 
                      sep = ""))

plot(consRDA_abun$siteConsensus[, 1:2], pch = 19, xlab = axisLabels[1], ylab = axisLabels[2], 
     las = 1)
abline(h = 0, v = 0, lty = 3)

arrows(0, 0, consRDA_abun$spConsensus[, 1] * 0.1, consRDA_abun$spConsensus[, 2] * 0.1, 
       col = "red", length = 0.1, angle = 12)
text(consRDA_abun$spConsensus[, 1:2] * 0.1, labels = rownames(consRDA_abun$spConsensus), 
     col = "red")

arrows(0, 0, consRDA_abun$descConsensus[, 1] * 0.6, consRDA_abun$descConsensus[, 2] * 0.6, 
       col = "blue", length = 0.1, angle = 12)
text(consRDA_abun$descConsensus[, 1:2] * 0.6, labels = rownames(consRDA_abun$descConsensus), 
     col = "blue")

consRDA_bin <- consensusRDA_DB(ordiRes_bin_final, rda_test_bin_final, Ybin, env, scaling = 1)
summary(consRDA_bin)
# How many axis are worth being considered:
sum_eig_bin <- sum(consRDA_bin$values)
ax_sel_bin <- which(as.numeric(consRDA_bin$values) >= mean(as.numeric(consRDA_bin$values)))
round(consRDA_bin$values[ax_sel_bin] / sum_eig_bin, 4)

axisLabels <- c(paste("Axis 1 - ", 
                      round(consRDA_bin$values[1]/sum(consRDA_bin$values), 4) * 100, sep = ""), 
                paste("Axis 2 - ", 
                      round(consRDA_bin$values[2]/sum(consRDA_bin$values), 4) * 100, sep = ""))

plot(consRDA_bin$siteConsensus[, 1:2], pch = 19, xlab = axisLabels[1], ylab = axisLabels[2], 
     las = 1)
abline(h = 0, v = 0, lty = 3)

arrows(0, 0, consRDA_bin$spConsensus[, 1] * 0.6, consRDA_bin$spConsensus[, 2] * 0.6, 
       col = "red", length = 0.1, angle = 12)
text(consRDA_bin$spConsensus[, 1:2] * 0.6, labels = rownames(consRDA_bin$spConsensus), 
     col = "red")

arrows(0, 0, consRDA_bin$descConsensus[, 1] * 0.6, consRDA_bin$descConsensus[, 2] * 0.6, 
       col = "blue", length = 0.1, angle = 12)
text(consRDA_bin$descConsensus[, 1:2] * 0.6, labels = rownames(consRDA_bin$descConsensus), 
     col = "blue")

