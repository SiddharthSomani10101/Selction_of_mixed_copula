# Copula package
######################## Plotting the single copula models in ggplot
library(copula)
library(ggplot2)
library(grid)

install.packages("fCopulae")
library(fCopulae)
library(QRM)
set.seed(235)
# Build and plotting a Frank, a Gumbel and a Clayton single copula
frank <- frankCopula(dim = 2, param = 8)
gumbel <- gumbelCopula(dim = 3, param = 5.6)
clayton <- claytonCopula(dim = 4, param = 19)
# Select the copula
cp <- claytonCopula(param = c(3.4), dim = 2)
normal <- normalCopula(param = 0.7, dim = 2)
fr <- rCopula(2000, frank)
gu <- rCopula(2000, gumbel)
cl <- rCopula(2000, clayton)


# Plot the samples
p1 <- qplot(fr[,1], fr[,2], colour = fr[,1], main="Frank copula", xlab = "u", ylab = "v")
p2 <- qplot(gu[,1], gu[,2], colour = gu[,1], main="Gumbel copula", xlab = "u", ylab = "v") 
p3 <- qplot(cl[,1], cl[,2], colour = cl[,1], main="Clayton copula", xlab = "u", ylab = "v")

# Define grid layout to locate plots and print each graph^(1)
pushViewport(viewport(layout = grid.layout(1, 3)))
print(p1, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(p2, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
print(p3, vp = viewport(layout.pos.row = 1, layout.pos.col = 3))


###########This is mixed copula code - Clayton-Gumbel
cc <- claytonCopula(iTau(claytonCopula(), tau = 0.50)) # the first component
gc <- gumbelCopula(iTau(gumbelCopula(),   tau = 0.50)) # the second component
wts <- c(1/3, 2/3) # the corresponding weights
(mcg <- mixCopula(list(cc, gc), w = wts)) # the mixture copula

stopifnot(
  all.equal(   rho(mcg), wts[1] *    rho(cc) + wts[2] *    rho(gc)),
  all.equal(lambda(mcg), wts[1] * lambda(cc) + wts[2] * lambda(gc)))
lambda(mcg)
set.seed(127)
U <- rCopula(1000, copula = mcg) # sample from the mixture
wireframe2(mcg, FUN = dCopula, delta = 0.050) # density
contourplot2(mcg, FUN = pCopula) # copula
contourplot2(mcg, FUN = dCopula, cuts = 32, # density
             n.grid = 50, pretty = FALSE,
             col = adjustcolor(1, 1/3), alpha.regions = 3/4)
plot(U, xlab = quote(U[1]), ylab = quote(U[2])) # scatter plot


###########This is mixed copula code - normal-Gumbel
cc <- normalCopula(iTau(normalCopula(), tau = 0.50)) # the first component
gc <- gumbelCopula(iTau(gumbelCopula(),   tau = 0.50)) # the second component
wts <- c(1/3, 2/3) # the corresponding weights
(mcg <- mixCopula(list(cc, gc), w = wts)) # the mixture copula

stopifnot(
  all.equal(   rho(mcg), wts[1] *    rho(cc) + wts[2] *    rho(gc)),
  all.equal(lambda(mcg), wts[1] * lambda(cc) + wts[2] * lambda(gc)))
lambda(mcg)
set.seed(127)
U <- rCopula(1000, copula = mcg) # sample from the mixture
wireframe2(mcg, FUN = dCopula, delta = 0.050) # density
contourplot2(mcg, FUN = pCopula) # copula
contourplot2(mcg, FUN = dCopula, cuts = 32, # density
             n.grid = 50, pretty = FALSE,
             col = adjustcolor(1, 1/3), alpha.regions = 3/4)
plot(U, xlab = quote(U[1]), ylab = quote(U[2])) # scatter plot
############################# Real data example

library(RColorBrewer)
cols <- brewer.pal(3, "BuGn")
pal <- colorRampPalette(cols)
exch<- read.csv(file="D:\\SPY.csv") 
usd<-exch$Adj_Close_SPY 
uk<-exch$Adj_Close_UK 
jp<-exch$Adj_Close_JP 
############ plotting closing prices of markets
plot(usd, type = "l", col = "black") 
############ calculating returns of jpy
usd_ret <- diff(usd) / usd[- length(usd)]  # Calculate returns
uk_ret <- diff(uk) / uk[- length(uk)]  # Calculate returns
jp_ret <- diff(jp) / jp[- length(jp)]  # Calculate returns
########################Plotting the returns of markets
plot(usd_ret, type = "l", col = "black") 
plot(jp_ret, type = "l", col = "blue")  
############## Finding out correlation
cor(usd,uk)
cor(usd,jp)
cor(uk,jp)
################## Finding kendall's tau correlation
res<-cor.test(usd,uk, method="kendall")
res
res<-cor.test(uk,jp, method="kendall")
res
res<-cor.test(usd,jp, method="kendall")
res
##################################### estimating parameter and log likelihood of copulas of international market data
####################USD-JP copula
val.ln <- cbind(exch[2],exch[3])
val.ln<- as.matrix(val.ln) 
n<-nrow(val.ln)
summary(val.ln)
ro<-cor(val.ln) 
ro<-ro[1,2]
ro
Udata <- pobs(val.ln)
rotau<-Kendall(val.ln)
rotau<-rotau[1,2]
rotau
ParGum<-1/(1-rotau)
ParClay<-(2*rotau)/(1-rotau)
norm.cop <- normalCopula(ro, dim = 2, dispstr = "un")
###
EstimatedNormCop<-fitCopula(norm.cop,Udata, method="mpl")
EstimatedNormCop
logLik(EstimatedNormCop) #estimating log likelihood
AIC(EstimatedNormCop) ## estimatingthe AIC
BIC(EstimatedNormCop) #estimating BIC
##gaussian copula estimation
norm.cop <- normalCopula(0.50688, dim = 2, dispstr = "un")
TailDep<-tailIndex(norm.cop)
TailDep

###### Cópula Gumbel #

gumb.cop0 <- gumbelCopula(ParGum, dim =2)
gumbCopEst<-fitCopula(gumb.cop0,Udata, method="mpl")
gumbCopEst
logLik(gumbCopEst)
AIC(gumbCopEst)
BIC(gumbCopEst)

# Dependencies on Gumbel copula
gumb.cop <- gumbelCopula(1.44915, dim =2)
TailDep<-tailIndex(gumb.cop)
TailDep

##############################
# Ajust clayton copula #

clay.cop0<- claytonCopula(param =ParClay, dim = 2)    
ClayCopEst<-fitCopula(clay.cop0,Udata, method="mpl")
ClayCopEst
logLik(ClayCopEst)
AIC(ClayCopEst)
BIC(ClayCopEst)
#Dependencies on clayton copula
clay.cop<- claytonCopula(0.8450, dim = 2)    
TailDep<-tailIndex(clay.cop)
TailDep

############################
# Ajuste frank copula #

frank.cop0<-frankCopula(param = NA_real_, dim = 2)
frank.cop0
FrankCopEst<-fitCopula(frank.cop0,Udata,method="mpl")
FrankCopEst
logLik(FrankCopEst)
AIC(FrankCopEst)
BIC(FrankCopEst)
#Dependencies on Frank copula
frank.cop<-frankCopula(3.3007, dim = 2)
TailDep<-tailIndex(frank.cop)
TailDep

########################## Estimate parameters for US-UK copula models
val.ln <- cbind(exch[2],exch[4])
val.ln<- as.matrix(val.ln) 
n<-nrow(val.ln)
summary(val.ln)
ro<-cor(val.ln) 
ro<-ro[1,2]
Udata <- pobs(val.ln)
rotau<-Kendall(val.ln)
rotau<-rotau[1,2]
ParGum<-1/(1-rotau)
ParClay<-(2*rotau)/(1-rotau)
############# Gaussian copula
norm.cop <- normalCopula(ro, dim = 2, dispstr = "un")
NormCopEst<-fitCopula(norm.cop,Udata, method="mpl")
NormCopEst
logLik(NormCopEst) 
AIC(NormCopEst) 
BIC(NormCopEst) 


# Ajust - Cópula Gumbel #
gumb.cop0 <- gumbelCopula(ParGum, dim =2)
gumbCopEst<-fitCopula(gumb.cop0,Udata, method="mpl")
gumbCopEst
logLik(gumbCopEst)
AIC(gumbCopEst)
BIC(gumbCopEst)


##############################
# Ajust - Cópula Clayton #
clay.cop0<- claytonCopula(param =ParClay, dim = 2)    
ClayCopEst<-fitCopula(clay.cop0,Udata, method="mpl")
ClayCopEst
logLik(ClayCopEst)
AIC(ClayCopEst)
BIC(ClayCopEst)


############################
# Ajust - Cópula  Frank #
frank.cop0<-frankCopula(param = NA_real_, dim = 2)
FrankCopEst<-fitCopula(frank.cop0,Udata,method="mpl")
FrankCopEst
logLik(FrankCopEst)
AIC(FrankCopEst)
BIC(FrankCopEst)


####################################### Estimate parameters for UK-JP copulas
val.ln <- cbind(exch[3],exch[4])
val.ln<- as.matrix(val.ln)
n<-nrow(val.ln)
summary(val.ln)
ro<-cor(val.ln) 
ro<-ro[1,2]
Udata <- pobs(val.ln)
rotau<-Kendall(val.ln)
rotau<-rotau[1,2]
rotau
ParGum<-1/(1-rotau)
ParClay<-(2*rotau)/(1-rotau)
############## Normal copula
norm.cop <- normalCopula(ro, dim = 2, dispstr = "un")
NormCopEst<-fitCopula(norm.cop,Udata, method="mpl")
NormCopEst
logLik(NormCopEst)
AIC(NormCopEst) 
BIC(NormCopEst) 







# Ajust - Cópula Gumbel #
gumb.cop0 <- gumbelCopula(ParGum, dim =2)
gumbCopEst<-fitCopula(gumb.cop0,Udata, method="mpl")
gumbCopEst
logLik(gumbCopEst)
AIC(gumbCopEst)
BIC(gumbCopEst)



##############################
# Ajust Cópula  Clayton #

clay.cop0<- claytonCopula(param =ParClay, dim = 2)    
ClayCopEst<-fitCopula(clay.cop0,Udata, method="mpl")
ClayCopEst
logLik(ClayCopEst)
AIC(ClayCopEst)
BIC(ClayCopEst)



############################
# Ajust Frank Copula
frank.cop0<-frankCopula(param = NA_real_, dim = 2)
FrankCopEst<-fitCopula(frank.cop0,Udata,method="mpl")
FrankCopEst
logLik(FrankCopEst)
AIC(FrankCopEst)
BIC(FrankCopEst)


