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

# consensusRDA() à utiliser après avoir utilisé coeffCompare()

### Construct results object
ndis <- 10
ordiRes <- vector("list", length = ndis)

#---------------------------------------------
### Perform the various constrained ordination
#---------------------------------------------
### RDA species profile
sp <- Y / apply(Y, 1, sum)
ordiRes[[1]] <- rda(sp~., data = env)

### RDA chord
chord <- Y/sqrt(apply(Y^2,1,sum))
ordiRes[[2]] <- rda(chord~ ., data = env)

### RDA Hellinger
hell <- decostand(Y, method = "hellinger")
ordiRes[[3]] <- rda(hell~ ., data = env)

### RDA chi2
chisq <- decostand(Y, method = "chi.square")
ordiRes[[4]] <- rda(chisq~ ., data = env)

### db-RDA Bray-Curtis
ordiRes[[5]] <- capscale(sqrt(vegdist(Y, method = "bray"))~ ., data = env, comm = Y)

### db-RDA square-root Bray-Curtis
ordiRes[[6]] <- capscale(sqrt(vegdist(Y^0.5, method = "bray"))~ ., data = env, comm = Y)

### db-RDA fourth-root Bray-Curtis
ordiRes[[7]] <- capscale(sqrt(vegdist(Y^0.25, method = "bray"))~ ., data = env, comm = Y)

### db-RDA modified Gower log 2
ordiRes[[8]] <- capscale(vegdist(decostand(Y, "log", logbase = 2), "altGower")~ ., data = env, 
                        comm = Y) ### Warning message stem from log transformation of 0

### db-RDA modified Gower log 5
ordiRes[[9]] <- capscale(vegdist(decostand(Y, "log", logbase = 5), "altGower")~ ., data = env,
                        comm = Y) ### Warning message stem from log transformation of 0

### db-RDA modified Gower log 10
ordiRes[[10]] <- capscale(vegdist(decostand(Y, "log", logbase = 10), "altGower")~ ., data = env,
                          comm = Y) ### Warning message stem from log transformation of 0

rda_test <- vector("list", ndis)
for (i in 1:ndis) rda_test[[i]] <- anova.cca(ordiRes[[i]], by = "axis")

### Compare association coefficients
AssoComp <- coeffCompare(ordiRes, rep(7, ndis))
