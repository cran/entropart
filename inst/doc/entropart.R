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
### code chunk number 6: lnq
###################################################
curve(log(x), 0, 1, lty=1, ylab = expression(ln[q](x)))
curve(lnq(x, 0), 0, 1, lty = 2, add = TRUE)
curve(lnq(x, 2), 0, 1, lty = 3, add = TRUE)
curve(lnq(x, 3), 0, 1, lty = 4, add = TRUE)  
legend("bottomright", legend = c(expression(ln[0](x)), "ln(x)", expression(ln[2](x)), expression(ln[3](x))), lty = c(2, 1, 3, 4), inset=  0.02)


###################################################
### code chunk number 7: ShannonP6.1
###################################################
Tsallis(Ps = Paracou618.MC$Ps, q = 1)


###################################################
### code chunk number 8: ShannonP6.2
###################################################
Shannon(Ps = Paracou618.MC$Ps)


###################################################
### code chunk number 9: Coverage
###################################################
Coverage(Ns = Paracou618.MC$Ns)


###################################################
### code chunk number 10: bcShannonP6.1
###################################################
bcTsallis(Ns = Paracou618.MC$Ns, q = 1)


###################################################
### code chunk number 11: bcShannonP6.2
###################################################
bcShannon(Ns = Paracou618.MC$Ns)


###################################################
### code chunk number 12: Diversity
###################################################
expq(Simpson(Ps = Paracou618.MC$Ps), q = 2)
Diversity(Ps = Paracou618.MC$Ps, q = 2)
expq(bcTsallis(Ns = Paracou618.MC$Ns, q = 2), q = 2)
bcDiversity(Ns = Paracou618.MC$Ns, q = 2)


###################################################
### code chunk number 13: AlphaEntropy
###################################################
e <- AlphaEntropy(Paracou618.MC, q = 1)
summary(e)


###################################################
### code chunk number 14: DivPart
###################################################
p <- DivPart(q = 1, MC = Paracou618.MC, Biased = FALSE)
summary(p)
p$CommunityAlphaEntropies


###################################################
### code chunk number 15: DivPartPlot
###################################################
par(mar=c(4, 4, 2, 1))
plot(p)


###################################################
### code chunk number 16: DivEst
###################################################
de <- DivEst(q = 1, Paracou618.MC, Biased = FALSE, Correction = "Best", Simulations = 1000)
summary(de)


###################################################
### code chunk number 17: DivEstPlot
###################################################
par(mar=c(4, 4, 2, 1))
plot(de)


###################################################
### code chunk number 18: DivProfile
###################################################
dp <- DivProfile(seq(0, 2, 0.2), Paracou618.MC, Biased = FALSE)
summary(dp)


###################################################
### code chunk number 19: DivProfilePlot
###################################################
plot(dp)


###################################################
### code chunk number 20: ShannonBeta
###################################################
ShannonBeta(Paracou618.MC$Psi[, 1], Paracou618.MC$Ps)


###################################################
### code chunk number 21: PhyloDiversity
###################################################
phd <- bcPhyloDiversity(Paracou618.MC$Ns, q = 1, Tree = Paracou618.Taxonomy, Normalize = TRUE)
summary(phd)


###################################################
### code chunk number 22: PhyloDiversityPlot
###################################################
par(mar=c(4, 4, 2, 1))
plot(phd, main = "")


###################################################
### code chunk number 23: PhylodivPart
###################################################
dp <- DivPart(q = 1, Paracou618.MC, Biased = FALSE, Correction = "Best", Tree = Paracou618.Taxonomy)
summary(dp)


###################################################
### code chunk number 24: PhyloBetaEntropy
###################################################
summary(BetaEntropy(Paracou618.MC, q = 2, Tree = Paracou618.Taxonomy, Correction = "None", Normalize = FALSE))


###################################################
### code chunk number 25: divc
###################################################
library("ade4")
divc(as.data.frame(Paracou618.MC$Wi), disc(as.data.frame(Paracou618.MC$Nsi), Paracou618.Taxonomy$Wdist))


###################################################
### code chunk number 26: Dqz
###################################################
DistanceMatrix <- as.matrix(Paracou618.dist)
Z <- 1 - DistanceMatrix/max(DistanceMatrix)
bcDqz(Paracou618.MC$Ns, q = 2, Z)


###################################################
### code chunk number 27: DqzHCDT
###################################################
Dqz(Paracou618.MC$Ps, q = 2, Z = diag(length(Paracou618.MC$Ps)))
Diversity(Paracou618.MC$Ps, q = 2)


###################################################
### code chunk number 28: Hqz
###################################################
Hqz(Paracou618.MC$Ps, q = 2, Z)
lnq(Dqz(Paracou618.MC$Ps, q = 2, Z), q = 2)


###################################################
### code chunk number 29: AlphaEntropy
###################################################
e <- AlphaEntropy(Paracou618.MC, q = 1, Z = Z)
summary(e)


###################################################
### code chunk number 30: EntropyCI
###################################################
SimulatedDiversity <- expq(EntropyCI(FUN = bcTsallis, Simulations = 1000, Ns = Paracou618.MC$Ns, q = 2), q = 2)
bcDiversity(Paracou618.MC$Ns, q = 2)
quantile(SimulatedDiversity, probs = c(0.025, 0.975))


