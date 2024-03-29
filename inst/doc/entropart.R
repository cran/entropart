## ----global_options, include=FALSE--------------------------------------------
set.seed(97310)

## ----LoadParacou18------------------------------------------------------------
library("entropart")
data("Paracou618")
N18 <- Paracou618.MC$Nsi[, "P018"]

## ----PlotN18------------------------------------------------------------------
Abd18 <- as.AbdVector(N18)
autoplot(Abd18, Distribution="lnorm")

## ----PN18---------------------------------------------------------------------
P18 <- as.ProbaVector(N18)

## ----rCommunity---------------------------------------------------------------
rc <- rCommunity(1, size=10000, Distribution = "lseries", alpha = 30)
autoplot(rc, Distribution="lseries")

## ----IndicesP-----------------------------------------------------------------
Richness(P18)
Shannon(P18)
Simpson(P18)

## ----IndicesAbd---------------------------------------------------------------
Richness(Abd18)
Shannon(Abd18)
Simpson(Abd18)

## ----Tsallis------------------------------------------------------------------
Tsallis(Abd18, q=1)

## ----Diversity----------------------------------------------------------------
Diversity(Abd18, q=1)

## ----lnq----------------------------------------------------------------------
(d2 <- Diversity(Abd18,q=2))
lnq(d2, q=2)
(e2 <-Tsallis(Abd18,q=2))
expq(e2, q=2)

## ----DiversityProfile---------------------------------------------------------
DP <- CommunityProfile(Diversity, Abd18)
autoplot(DP)

## ----PhyloDiversity-----------------------------------------------------------
summary(PhyloDiversity(Abd18,q=1,Tree=Paracou618.Taxonomy))

## ----SBDiversity--------------------------------------------------------------
# Prepare the similarity matrix
DistanceMatrix <- as.matrix(Paracou618.dist)
# Similarity can be 1 minus normalized distances between species
Z <- 1 - DistanceMatrix/max(DistanceMatrix)
# Calculate diversity of order 2
Dqz(Abd18, q=2, Z)

## ----PDiversityProfile--------------------------------------------------------
sbDP <- CommunityProfile(Dqz, Abd18, Z=Z)
pDP <- CommunityProfile(function(X, ...) PhyloDiversity(X, ...)$Total, Abd18, Tree=Paracou618.Taxonomy)
autoplot(pDP)

## ----MetaCommunitydf----------------------------------------------------------
library("entropart")
(df <- data.frame(C1 = c(10, 10, 10, 10), C2 = c(0, 20, 35, 5), C3 = c(25, 15, 0, 2), row.names = c("sp1", "sp2", "sp3", "sp4")))
w <- c(1, 2, 1)

## ----MetaCommunityMC----------------------------------------------------------
MC <- MetaCommunity(Abundances = df, Weights = w)
plot(MC)

## ----DivPart------------------------------------------------------------------
p <- DivPart(q = 1, MC = Paracou618.MC)
summary(p)

## ----DivEst, fig.width=6, fig.height=6----------------------------------------
de <- DivEst(q = 1, Paracou618.MC, Simulations = 50)
# Margin adjustment required for this vignette html output
par(mar=c(1,1,2.2,1))
autoplot(de)

## ----DivProfile, fig.width=6, fig.height=6------------------------------------
dp <- DivProfile(q.seq = seq(0, 2, 0.1), Paracou618.MC)
autoplot(dp)

