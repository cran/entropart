### R code from vignette source 'entropart.Rnw'
### Encoding: ISO8859-1

###################################################
### code chunk number 1: Declarations
###################################################
options(prompt = "R> ", continue = "+  ", width = 70, useFancyQuotes = FALSE)


###################################################
### code chunk number 2: MetaCommunity
###################################################
library("entropart")
(df <- data.frame(C1 = c(10, 10, 10, 10), C2 = c(0, 20, 35, 5), C3 = c(25, 15, 0, 2), row.names = c("sp1", "sp2", "sp3", "sp4")))
w <- c(1, 2, 1)
MC <- MetaCommunity(Abundances = df, Weights = w)


###################################################
### code chunk number 3: MCPlot
###################################################
par(mar=c(2.5, 4, 1, 0))
plot(MC)


###################################################
### code chunk number 4: Ps
###################################################
MC$Ps


###################################################
### code chunk number 5: data
###################################################
data("Paracou618")
summary(Paracou618.MC)


###################################################
### code chunk number 6: ShannonP6.1
###################################################
Tsallis(Ps = Paracou618.MC$Ps, q = 1)


###################################################
### code chunk number 7: ShannonP6.2
###################################################
Shannon(Ps = Paracou618.MC$Ps)


###################################################
### code chunk number 8: Coverage
###################################################
Coverage(Ns = Paracou618.MC$Ns)


###################################################
### code chunk number 9: bcShannonP6.1
###################################################
bcTsallis(Ns = Paracou618.MC$Ns, q = 1)


###################################################
### code chunk number 10: bcShannonP6.2
###################################################
bcShannon(Ns = Paracou618.MC$Ns)


###################################################
### code chunk number 11: Diversity
###################################################
expq(Simpson(Ps = Paracou618.MC$Ps), q = 2)
Diversity(Ps = Paracou618.MC$Ps, q = 2)
expq(bcTsallis(Ns = Paracou618.MC$Ns, q = 2), q = 2)
bcDiversity(Ns = Paracou618.MC$Ns, q = 2)


###################################################
### code chunk number 12: AlphaEntropy
###################################################
e <- AlphaEntropy(Paracou618.MC, q = 1)
summary(e)


###################################################
### code chunk number 13: DivPart
###################################################
p <- DivPart(q = 1, MC = Paracou618.MC, Biased = FALSE)
summary(p)
p$CommunityAlphaEntropies


###################################################
### code chunk number 14: DivPartPlot
###################################################
par(mar=c(4, 4, 2, 1))
plot(p)


###################################################
### code chunk number 15: DivEst
###################################################
de <- DivEst(q = 1, Paracou618.MC, Biased = FALSE, Correction = "Best", Simulations = 1000)
summary(de)


###################################################
### code chunk number 16: DivEstPlot
###################################################
par(mar=c(4, 4, 2, 1))
plot(de)


###################################################
### code chunk number 17: DivProfile
###################################################
dp <- DivProfile(seq(0, 2, 0.2), Paracou618.MC, Biased = FALSE)
summary(dp)


###################################################
### code chunk number 18: DivProfilePlot
###################################################
plot(dp)


###################################################
### code chunk number 19: ShannonBeta
###################################################
ShannonBeta(Paracou618.MC$Psi[, 1], Paracou618.MC$Ps)


###################################################
### code chunk number 20: PhyloDiversity
###################################################
phd <- bcPhyloDiversity(Paracou618.MC$Ns, q = 1, Tree = Paracou618.Taxonomy, Normalize = TRUE)
summary(phd)


###################################################
### code chunk number 21: PhyloDiversityPlot
###################################################
par(mar=c(4, 4, 2, 1))
plot(phd, main = "")


###################################################
### code chunk number 22: PhylodivPart
###################################################
dp <- DivPart(q = 1, Paracou618.MC, Biased = FALSE, Correction = "Best", Tree = Paracou618.Taxonomy)
summary(dp)


###################################################
### code chunk number 23: PhyloBetaEntropy
###################################################
summary(BetaEntropy(Paracou618.MC, q = 2, Tree = Paracou618.Taxonomy, Correction = "None", Normalize = FALSE))


###################################################
### code chunk number 24: divc
###################################################
library("ade4")
divc(as.data.frame(Paracou618.MC$Wi), disc(as.data.frame(Paracou618.MC$Nsi), Paracou618.Taxonomy$Wdist))


###################################################
### code chunk number 25: Dqz
###################################################
DistanceMatrix <- as.matrix(Paracou618.dist)
Z <- 1 - DistanceMatrix/max(DistanceMatrix)
bcDqz(Paracou618.MC$Ns, q = 2, Z)


###################################################
### code chunk number 26: DqzHCDT
###################################################
Dqz(Paracou618.MC$Ps, q = 2, Z = diag(length(Paracou618.MC$Ps)))
Diversity(Paracou618.MC$Ps, q = 2)


###################################################
### code chunk number 27: Hqz
###################################################
Hqz(Paracou618.MC$Ps, q = 2, Z)
lnq(Dqz(Paracou618.MC$Ps, q = 2, Z), q = 2)


###################################################
### code chunk number 28: AlphaEntropy
###################################################
e <- AlphaEntropy(Paracou618.MC, q = 1, Z = Z)
summary(e)


###################################################
### code chunk number 29: EntropyCI
###################################################
SimulatedDiversity <- expq(EntropyCI(FUN = bcTsallis, Simulations = 1000, Ns = Paracou618.MC$Ns, q = 1), q = 1)
bcDiversity(Paracou618.MC$Ns, q = 1)
quantile(SimulatedDiversity, probs = c(0.025, 0.975))


###################################################
### code chunk number 30: CommunityProfile
###################################################
bcProfile <- CommunityProfile(bcDiversity, Paracou618.MC$Ns)
Profile <- CommunityProfile(Diversity, Paracou618.MC$Ps)


###################################################
### code chunk number 31: CommunityProfileFigCode (eval = FALSE)
###################################################
## plot(bcProfile, type="l", main="", xlab="q", ylab="Diversity")
## lines(Profile$y~Profile$x, lty=3)
## legend("topright", c("Bias Corrected", "Biased"), lty=c(1,3), inset=0.02)


###################################################
### code chunk number 32: CommunityProfileFig
###################################################
par(mar=c(4, 4, 2, 1))
plot(bcProfile, type="l", main="", xlab="q", ylab="Diversity")
lines(Profile$y~Profile$x, lty=3)
legend("topright", c("Bias Corrected", "Biased"), lty=c(1,3), inset=0.02)


###################################################
### code chunk number 33: PhyloApply
###################################################
pa <- PhyloApply(Tree=Paracou618.Taxonomy, FUN=bcTsallis, NorP=Paracou618.MC$Ns)
summary(pa)
exp(pa$Cuts)
exp(pa$Total)


###################################################
### code chunk number 34: MC1
###################################################
(df <- data.frame(C1 = c(10, 10, 10, 10), C2 = c(0, 20, 35, 5),
  C3 = c(25, 15, 0, 2), row.names = c("sp1", "sp2", "sp3", "sp4")))
w <- colSums(df)
MC1 <- MetaCommunity(Abundances = df, Weights = w)


###################################################
### code chunk number 35: MC1
###################################################
(df <- data.frame(C1 = c(10, 4), C2 = c(3, 4), row.names = c("sp1", "sp5")))
w <- colSums(df)
MC2 <- MetaCommunity(Abundances = df, Weights = w)


###################################################
### code chunk number 36: MCMergeC
###################################################
mergedMC1 <- MergeC(list(MC1, MC2))
mergedMC1$Nsi


###################################################
### code chunk number 37: MCMergeMC
###################################################
mergedMC2 <- MergeMC(list(MC1, MC2), Weights = sapply(list(MC1, MC2), function(x) (x$N)))
mergedMC2$Nsi


###################################################
### code chunk number 38: MCMergeHi
###################################################
dpAll <- DivPart(q=1, MC=mergedMC2)
summary(dpAll)


###################################################
### code chunk number 39: MCMergeLo1
###################################################
dpMC1 <- DivPart(q=1, MC=MC1)
summary(dpMC1)


