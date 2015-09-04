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
### code chunk number 6: AbdVector
###################################################
data("Paracou618")
PAbd <- as.AbdVector(Paracou618.MC$Ns)


###################################################
### code chunk number 7: ProbaVector
###################################################
PProba <- as.ProbaVector(Paracou618.MC$Ns)


###################################################
### code chunk number 8: ProbaVectorPlot
###################################################
plot(PAbd)


###################################################
### code chunk number 9: lnq
###################################################
curve(log(x), 0, 1, lty=1, ylab = expression(ln[q](x)))
curve(lnq(x, 0), 0, 1, lty = 2, add = TRUE)
curve(lnq(x, 2), 0, 1, lty = 3, add = TRUE)
curve(lnq(x, 3), 0, 1, lty = 4, add = TRUE)  
legend("bottomright", legend = c(expression(ln[0](x)), "ln(x)", expression(ln[2](x)), expression(ln[3](x))), lty = c(2, 1, 3, 4), inset=  0.02)


###################################################
### code chunk number 10: ShannonP6.1
###################################################
Tsallis(Ps = Paracou618.MC$Ps, q = 1)


###################################################
### code chunk number 11: ShannonP6.2
###################################################
Shannon(Ps = Paracou618.MC$Ps)


###################################################
### code chunk number 12: Coverage
###################################################
Coverage(Ns = Paracou618.MC$Ns)


###################################################
### code chunk number 13: bcShannonP6.1
###################################################
bcTsallis(Ns = Paracou618.MC$Ns, q = 1)


###################################################
### code chunk number 14: bcShannonP6.2
###################################################
bcTsallis(Paracou618.MC$Ns, q = 1)


###################################################
### code chunk number 15: bcShannonP6.3
###################################################
Tsallis(PAbd, q = 1)


###################################################
### code chunk number 16: bcShannonP6.4
###################################################
Tsallis(PProba, q = 1)


###################################################
### code chunk number 17: Diversity
###################################################
expq(Simpson(Ps = Paracou618.MC$Ps), q = 2)
Diversity(Ps = Paracou618.MC$Ps, q = 2)
expq(bcTsallis(Ns = Paracou618.MC$Ns, q = 2), q = 2)
bcDiversity(Ns = Paracou618.MC$Ns, q = 2)


###################################################
### code chunk number 18: Hurlbert
###################################################
Hurlbert(Ps = Paracou618.MC$Ps, k = 2)
bcHurlbert(Ns = Paracou618.MC$Ns, k = 2)
bcHurlbertD(Ns = Paracou618.MC$Ns, k = 2)


###################################################
### code chunk number 19: AlphaEntropy
###################################################
e <- AlphaEntropy(Paracou618.MC, q = 1)
summary(e)


###################################################
### code chunk number 20: DivPart
###################################################
p <- DivPart(q = 1, MC = Paracou618.MC, Biased = FALSE)
summary(p)
p$CommunityAlphaEntropies


###################################################
### code chunk number 21: DivPartPlot
###################################################
par(mar=c(4, 4, 2, 1))
plot(p)


###################################################
### code chunk number 22: DivEst
###################################################
de <- DivEst(q = 1, Paracou618.MC, Biased = FALSE, Correction = "Best", Simulations = 100)
summary(de)


###################################################
### code chunk number 23: DivEstPlot
###################################################
par(mar=c(4, 4, 2, 1))
plot(de)


###################################################
### code chunk number 24: DivProfile
###################################################
dp <- DivProfile(seq(0, 2, 0.2), Paracou618.MC, Biased = FALSE)
summary(dp)


###################################################
### code chunk number 25: DivProfilePlot
###################################################
plot(dp)


###################################################
### code chunk number 26: ShannonBeta
###################################################
ShannonBeta(Paracou618.MC$Psi[, 1], Paracou618.MC$Ps)


###################################################
### code chunk number 27: PhyloDiversity
###################################################
phd <- bcPhyloDiversity(Paracou618.MC$Ns, q = 1, Tree = Paracou618.Taxonomy, Normalize = TRUE)
summary(phd)


###################################################
### code chunk number 28: PhyloDiversityPlot
###################################################
par(mar=c(4, 4, 2, 1))
plot(phd, main = "")


###################################################
### code chunk number 29: PhylodivPart
###################################################
dp <- DivPart(q = 1, Paracou618.MC, Biased = FALSE, Correction = "Best", Tree = Paracou618.Taxonomy)
summary(dp)


###################################################
### code chunk number 30: PhyloBetaEntropy
###################################################
summary(BetaEntropy(Paracou618.MC, q = 2, Tree = Paracou618.Taxonomy, Correction = "None", Normalize = FALSE))


###################################################
### code chunk number 31: divc
###################################################
library("ade4")
divc(as.data.frame(Paracou618.MC$Wi), disc(as.data.frame(Paracou618.MC$Nsi), Paracou618.Taxonomy$Wdist))


###################################################
### code chunk number 32: Dqz
###################################################
DistanceMatrix <- as.matrix(Paracou618.dist)
Z <- 1 - DistanceMatrix/max(DistanceMatrix)
bcDqz(Paracou618.MC$Ns, q = 2, Z)


###################################################
### code chunk number 33: DqzHCDT
###################################################
Dqz(Paracou618.MC$Ps, q = 2, Z = diag(length(Paracou618.MC$Ps)))
Diversity(Paracou618.MC$Ps, q = 2)


###################################################
### code chunk number 34: Hqz
###################################################
Hqz(Paracou618.MC$Ps, q = 2, Z)
lnq(Dqz(Paracou618.MC$Ps, q = 2, Z), q = 2)


###################################################
### code chunk number 35: AlphaEntropy
###################################################
e <- AlphaEntropy(Paracou618.MC, q = 1, Z = Z)
summary(e)


###################################################
### code chunk number 36: rCommunity
###################################################
rCommunity(n = 1, size = 1000, S=300, Distribution = "lnorm", sd=1) -> NsRef
Richness(as.ProbaVector(NsRef))
bcRichness(NsRef)
plot(NsRef, Distribution="lnorm")


###################################################
### code chunk number 37: rCommunityFig
###################################################
plot(NsRef, Distribution="lnorm") -> fit


###################################################
### code chunk number 38: EntropyCI
###################################################
SimulatedDiversity <- expq(EntropyCI(FUN = Tsallis, Simulations = 100, Ns = PAbd, q = 1), q = 1)
Diversity(PAbd, q = 1)
quantile(SimulatedDiversity, probs = c(0.025, 0.975))


###################################################
### code chunk number 39: CommunityProfile
###################################################
bcProfile <- CommunityProfile(Diversity, PAbd, NumberOfSimulations = 10)
Profile <- CommunityProfile(Diversity, PProba)


###################################################
### code chunk number 40: CommunityProfileFigCode (eval = FALSE)
###################################################
## plot(bcProfile)
## CEnvelope(Profile, lty=3)
## legend("topright", c("Bias Corrected", "Biased"), lty=c(1,3), inset=0.02)


###################################################
### code chunk number 41: CommunityProfileFig
###################################################
par(mar=c(4, 4, 2, 1))
plot(bcProfile)
CEnvelope(Profile, lty=3)
legend("topright", c("Bias Corrected", "Biased"), lty=c(1,3), inset=0.02)


###################################################
### code chunk number 42: PhyloApply
###################################################
pa <- PhyloApply(Tree=Paracou618.Taxonomy, FUN=bcTsallis, NorP=Paracou618.MC$Ns)
summary(pa)
exp(pa$Cuts)
exp(pa$Total)


###################################################
### code chunk number 43: MC1
###################################################
(df <- data.frame(C1 = c(10, 10, 10, 10), C2 = c(0, 20, 35, 5),
  C3 = c(25, 15, 0, 2), row.names = c("sp1", "sp2", "sp3", "sp4")))
w <- colSums(df)
MC1 <- MetaCommunity(Abundances = df, Weights = w)


###################################################
### code chunk number 44: MC1
###################################################
(df <- data.frame(C1 = c(10, 4), C2 = c(3, 4), row.names = c("sp1", "sp5")))
w <- colSums(df)
MC2 <- MetaCommunity(Abundances = df, Weights = w)


###################################################
### code chunk number 45: MCMergeC
###################################################
mergedMC1 <- MergeC(list(MC1, MC2))
mergedMC1$Nsi


###################################################
### code chunk number 46: MCMergeMC
###################################################
mergedMC2 <- MergeMC(list(MC1, MC2), Weights = sapply(list(MC1, MC2), function(x) (x$N)))
mergedMC2$Nsi


###################################################
### code chunk number 47: MCMergeHi
###################################################
dpAll <- DivPart(q=1, MC=mergedMC2)
summary(dpAll)


###################################################
### code chunk number 48: MCMergeLo1
###################################################
dpMC1 <- DivPart(q=1, MC=MC1)
summary(dpMC1)


