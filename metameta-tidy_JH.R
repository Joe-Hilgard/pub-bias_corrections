library(pwr)
library(truncnorm)
library(psych)
library(grid)
library(ggplot2)
library(metafor)




# build a function to run a meta-analysis and store some values. 
#You can tell it whether or not you want publication bias (dropping all nonsignificant studies), 
#and whether or not you want heterogeneity in effects

meta <- function(k, pubbias = c("absolute", "small", "none"), het, 
                 klarge = 0, nlarge = 280) { ## this is an example for publication bias and homogenous ES
  nobs <- numeric(k)  # set up some empty bins to store outputs
  tobs <- numeric(k)
  pobs <- numeric(k)
  dfobs <- numeric(k)
  d <- numeric(k)
  dobs <- numeric(k)
  sedobs <- numeric(k)
  vardobs <- numeric(k)
  nNeed <- numeric(k)
  mx <- numeric(k)
  my <- numeric(k)
  sdx <- numeric(k)
  sdy <- numeric(k)
  nx <-nobs
 
  #simulate data for each experiment
  
  for(i in 1:k){ #for each simulated experiment
    dtrue <-   if (het==TRUE) {
      rtruncnorm(1, a=.2, b=.4, mean = .2, sd = .12) #under heterogeneity it picks realistic values between .2 and .4
    } else {
      runif(1,.27,.28) # under homogeneity, it sticks close to the mean of the het. distribution
    }  
    nreq <-  round(power.t.test(d=dtrue, n=NULL, power=.8)$n, 0) #pick n for .8 power
    nused <- round(rtruncnorm(1, a=20, b=200, mean = 30, sd = 50),0) # actual N, again picked at random. I'm using a truncated normal so that n peaks around 30 and tapers from there. max 200, min 20
    if (i <= klarge) nused = nlarge # klarge studies have sample size nlarge
    x <- rnorm(n = nused, mean = 0, sd = 1) # create data
    y <- rnorm(n = nused, mean = dtrue, sd = 1)
    nobs[i] <- nused  # record n
    nNeed[i] <- nreq
    d[i] <- dtrue #record true d
    mx[i] <- mean(x) # store the stuff needed to calculate d and se(d)
    my[i] <- mean(y)
    sdx[i] <- sd(x)
    sdy[i] <- sd(x)
    test <- t.test(y,x) # run a t test
    tobs[i] <- test$statistic # grab t
    pobs[i] <- test$p.value # grab the p-value
    dfobs[i] <- test$parameter
  }
  
  
  
  meandiff = my-mx # again, just use standard formulae to calculate things
  sdpool = sqrt(((sdx^2+sdy^2)/2))
  dobs <- meandiff/sdpool # here's cohen's d
  vardobs <- ((2*nobs/nobs^2)+(dobs^2/(2*(2*nobs-2))))*(2*nobs/(nobs+nobs-2)) # calculate the variance from standard formula
  sedobs <- sqrt(vardobs) #ditto
  
  

  ### publication bias. Toss out nonsignificant studies. 
  # This section is ommitted in the conditions with no publication bias
  

  simdat <- data.frame(cbind(dobs, sedobs, vardobs, nNeed, nobs, pobs, d)) #gather required inputs
  
  if (pubbias=="absolute") {
	  simdat.selected <- subset(simdat, pobs <=.05) # if T, toss all nonsignificant data. Harsh.
  } else if (pubbias == "small") {
  	simdat.selected <- subset(simdat, pobs <= .05 | nobs > 200) # if F, use the full set
  } else if (pubbias == "none") {
    simdat.selected <- simdat
  }
  
  
  ### run PET-PEESE
  
  pet <- lm(dobs ~ sedobs, data=simdat.selected, weights= 1/vardobs) #run pet, store values
  pet.d <- summary(pet)$coefficients[1]
  pet.ld <- confint(pet)[1,1]
  pet.ud <- confint(pet)[1,2]
  

  peese <- lm(dobs ~ vardobs, data=simdat.selected, weights = 1/vardobs)  #run peese, store values
  peese.d <- summary(peese)$coefficients[1]
  peese.ld <- confint(peese)[1,1]
  peese.ud <- confint(peese)[1,2]
  
  # observed mean ES store this
  
  mean.obs = mean(simdat.selected$dobs)
  
  # real mean ES from ALL studies run (regardless of pub. bias). store this
 
  mean.real = mean(simdat$d)
  
  # give full pet-peese estimates. basically tell it to use pet if can't reject nil, use peese if can
  
  pp.d <- ifelse(pet.ld < 0 & pet.ud > 0, pet.d, peese.d)
  pp.ld <- ifelse(pet.ld < 0 & pet.ud > 0, pet.ld, peese.ld)
  pp.ud <- ifelse(pet.ld < 0 & pet.ud > 0, pet.ud, peese.ud)
  
  
  # return how many studies actually made the cut
  
  kper <- nrow(simdat.selected)
  
  # return total n and median n per observed study
  
  ntot <- sum(simdat.selected$nobs)
  nmed <- median(simdat.selected$nobs)
  nmax <- max(simdat.selected$nobs)
  
  # calculate correlation between actual N and required N among selected studies
  
  nCorr <- cor(simdat.selected$nobs, simdat.selected$nNeed)
  
  # calculate bias in pet-peese estimate
  
  pp.bias <- pp.d-mean.real
  
  return(data.frame(pet.d, pet.ld, pet.ud, peese.d, peese.ld, peese.ud, 
                    mean.obs, mean.real, pp.d, pp.ld, pp.ud, 
                    kper, ntot, nmed, nmax, nCorr, pp.bias)) # this tells it what to spit out for each meta-analysis
}


### so, that was creating a function for generating 100 studies 
# and using pet-peese to meta-analyze them. 

# to test it, you can just run the function once
meta(k = 100, 
     pubbias = "small", 
     het = F,
     klarge = 1,
     nlarge = 280) 






#### now do the meta-metas. simulate a bunch of meta-analyses under different conditions ####

set.seed(9999) # this should ensure similar results

nMeta <- 1000 # how many meta-analyses do you want?

phet <- data.frame() #this is the meta-meta simulation for publication bias + heterogeneity
for (j in 1:nMeta) {
	print(paste0("Running ", j, "/", nMeta))
	phet <- rbind(phet, meta(k=100, pubbias="small", het=T, klarge = 1))
}
summary(phet)


phom <- data.frame() #this is the meta-meta simulation for publication bias + homogeneity. it also makes me want to eat pho
for (j in 1:nMeta) {
  print(paste0("Running ", j, "/", nMeta))
  phom <- rbind(phom, meta(k=100, pubbias = "small", het=F, 
                           klarge = 1, nlarge = 280))
}
summary(phom)


nhet <- data.frame() #this is the meta-meta simulation for NO publication bias + heterogeneity
for (j in 1:nMeta) {
  print(paste0("Running ", j, "/", nMeta))
  nhet <- rbind(nhet, meta(k=100, pubbias=F, het=T))
}
summary(nhet)


nhom <- data.frame() #this is the meta-meta simulation for publication bias + homogeneity
for (j in 1:nMeta) {
  print(paste0("Running ", j, "/", nMeta))
  nhom <- rbind(nhom, meta(k=100, pubbias=F, het=F))
}
summary(nhom)



#### plot densities of bias in each condition ####

ggplot(phet, aes(x=pp.bias)) + 
  geom_vline(xintercept = 0, lty= 3, color = "darkgrey")+
  geom_density(fill="purple", colour="black", alpha=.2) +
  theme_bw() +
  geom_vline(xintercept = median(phet$pp.bias), lty= 1, color = "purple") +
  coord_cartesian(xlim = c(-.5, .2), ylim = c(0,11)) +
  labs(x="\nbias in PET-PEESE estimate (d)", y = "Density\n", title = "Heterogenous ES\n+ Publication Bias\n") +
  geom_segment(aes(x = 0, y = 1, xend = median(phet$pp.bias), yend = 1), arrow = arrow(length = unit(0.4, "cm")), color="red") +
  annotate("text", x = median(phet$pp.bias)/2, y = 1.4, label = round(median(phet$pp.bias), 3))


ggplot(phom, aes(x=pp.bias)) + 
  geom_vline(xintercept = 0, lty= 3, color = "darkgrey")+
  geom_density(fill="purple", colour="black", alpha=.2) +
  theme_bw() +
  geom_vline(xintercept = median(phom$pp.bias), lty= 1, color = "purple") +
  #coord_cartesian(xlim = c(-.5, .2), ylim = c(0,11)) +
  labs(x="\nbias in PET-PEESE estimate (d)", y = "Density\n", title = "Homogenous ES\n+ Publication Bias\n") +
  geom_segment(aes(x = 0, y = 1, xend = median(phom$pp.bias), yend = 1), arrow = arrow(length = unit(0.4, "cm")), color="red") +
  annotate("text", x = median(phom$pp.bias)/2, y = 1.4, label = round(median(phom$pp.bias), 3))


ggplot(nhet, aes(x=pp.bias)) + 
  geom_vline(xintercept = 0, lty= 3, color = "darkgrey")+
  geom_density(fill="purple", colour="black", alpha=.2) +
  theme_bw() +
  geom_vline(xintercept = median(nhet$pp.bias), lty= 1, color = "purple") +
  coord_cartesian(xlim = c(-.5, .2), ylim = c(0,11)) +
  labs(x="\nbias in PET-PEESE estimate (d)", y = "Density\n", title = "Heterogenous ES\nNo Publication Bias\n") +
  geom_segment(aes(x = 0, y = 1, xend = median(nhet$pp.bias), yend = 1), arrow = arrow(length = unit(0.4, "cm")), color="red") +
  annotate("text", x = median(nhet$pp.bias)/2, y = 1.4, label = round(median(nhet$pp.bias), 3))

ggplot(nhom, aes(x=pp.bias)) + 
  geom_vline(xintercept = 0, lty= 3, color = "darkgrey")+
  geom_density(fill="purple", colour="black", alpha=.2) +
  theme_bw() +
  geom_vline(xintercept = median(nhom$pp.bias), lty= 1, color = "purple") +
  coord_cartesian(xlim = c(-.5, .2), ylim = c(0,11)) +
  labs(x="\nbias in PET-PEESE estimate (d)", y = "Density\n", title = "Homogenous ES\nNo Publication Bias\n") +
  geom_segment(aes(x = 0, y = 1, xend = median(nhom$pp.bias), yend = 1), arrow = arrow(length = unit(0.4, "cm")), color="red") +
  annotate("text", x = median(nhom$pp.bias)/2, y = 1.4, label = round(median(nhom$pp.bias), 3))



#### check out the relationships between sample size accuracy and pet-peese bias ####

cor(phet$nCorr, phet$pp.bias)
cor(phom$nCorr, phom$pp.bias)
cor(nhet$nCorr, nhet$pp.bias)
cor(nhom$nCorr, nhom$pp.bias)

# what does this look like with no bias?
nhet$cor.1 <- nhet$nCorr*10  #rescale the correlation to units of .1

acc.lm <- lm(pp.bias ~ cor.1, data=nhet)
summary(acc.lm)

# looks bad for heterogenous effect sizes. Plot this


ggplot(phet, aes( x = nCorr, y = pp.bias)) +
  theme_bw() +
  geom_point(alpha = .15, color = "black") +
  geom_smooth(method=lm, se = F, color = "purple", size = 2) +
  labs( x = "\nCorrelation between required N and actual N", y = "Bias in PET-PEESE Estimate\n", title="Heterogenous ES\n+ Publication Bias")

ggplot(nhet, aes( x = nCorr, y = pp.bias)) +
  theme_bw() +
  geom_point(alpha = .15, color = "black") +
  geom_smooth(method=lm, se = F, color = "purple", size = 2) +
  labs( x = "\nCorrelation between required N and actual N", y = "Bias in PET-PEESE Estimate\n", title="Heterogenous ES\nNo Publication Bias")



#### put them all together, and run some descriptives ####

phet$cond <- factor("PBHet")
phom$cond <- factor("PBHom")
nhet$cond <- factor("NBHet")
nhom$cond <- factor("NBHom")

full <- rbind(phet,phom,nhet,nhom)

full$miss.real <- ifelse(full$pp.ld <= full$mean.real & full$pp.ud >= full$mean.real, "Hit", "Miss") # does pet-peese fail to include the actual effect size?
full$nofx <- ifelse(full$pp.ld <= 0 & full$pp.ud >= 0, "No Effect", "Effect") # does pet-peese mistakenly claim that no effect exists?
full$neg.bias <- ifelse(full$pp.bias < 0, "Negative", "Nonnegative") # is it negatively biased?


prop.table(table(full$cond, full$miss.real),1)

prop.table(table(full$cond, full$nofx),1)

prop.table(table(full$cond, full$neg.bias),1)

summs <- full[c("mean.real", "pp.d", "pp.ld", "pp.ud", "pp.bias", "cond")]
describeBy(x=summs, group=summs$cond) # grab sume descriptives



## Finally some plots for the sample size and effect size distributions the sim samples from

ggplot(data.frame(x=c(.2, .4)), aes(x)) +
  stat_function(fun = dnorm, args = list(mean = .2, sd = .12), geom="density", fill = "purple", alpha=.2) +
  theme_bw() +
  labs( x= "\nEffect size (d)", y = "density\n")

ggplot(data.frame(x=c(20, 200)), aes(x)) +
  stat_function(fun = dnorm, args = list(mean = 30, sd = 50), geom="density", fill = "purple", alpha=.2) +
  theme_bw() +
  labs( x= "\nSample size", y = "density\n")

# Hilgard is going ham ----
phom_k1_pow80 <- data.frame() #Publication bias, one guaranteed 80% pow study
for (j in 1:nMeta) {
  print(paste0("Running ", j, "/", nMeta))
  phom_k1_pow80 <- rbind(phom_k1_pow80, meta(k=100, pubbias = "small", het=F, 
                           klarge = 1, nlarge = 209))
}
summary(phom_k1_pow80)

phom_k1_pow90 <- data.frame() #Publication bias, one guaranteed 90% pow study
for (j in 1:nMeta) {
  print(paste0("Running ", j, "/", nMeta))
  phom_k1_pow90 <- rbind(phom_k1_pow90, meta(k=100, pubbias = "small", het=F, 
                           klarge = 1, nlarge = 280))
}
summary(phom_k1_pow90)

phom_k3_pow80 <- data.frame() #Publication bias, three guaranteed 80% pow studies
for (j in 1:nMeta) {
  print(paste0("Running ", j, "/", nMeta))
  phom_k3_pow80 <- rbind(phom_k3_pow80, meta(k=100, pubbias = "small", het=F, 
                           klarge = 3, nlarge = 209))
}
summary(phom_k3_pow80)

phom_k10_pow80 <- data.frame() #Publication bias, three guaranteed 80% pow studies
for (j in 1:nMeta) {
  print(paste0("Running ", j, "/", nMeta))
  phom_k10_pow80 <- rbind(phom_k10_pow80, meta(k=100, pubbias = "small", het=F, 
                                             klarge = 10, nlarge = 209))
}
summary(phom_k10_pow80)

# Hilgard plots
ggplot(phom_k1_pow80, aes(x=pp.bias)) + 
  geom_vline(xintercept = 0, lty= 3, color = "darkgrey")+
  geom_density(fill="purple", colour="black", alpha=.2) +
  theme_bw() +
  geom_vline(xintercept = median(phom$pp.bias), lty= 1, color = "purple") +
  #coord_cartesian(xlim = c(-.5, .2), ylim = c(0,11)) +
  labs(x="\nbias in PET-PEESE estimate (d)", y = "Density\n", title = "Homogenous ES\n+ Publication Bias\n+ Guaranteed k=1 @ 80% power") +
  geom_segment(aes(x = 0, y = 1, xend = median(phom$pp.bias), yend = 1), arrow = arrow(length = unit(0.4, "cm")), color="red") +
  annotate("text", x = median(phom$pp.bias)/2, y = 1.4, label = round(median(phom$pp.bias), 3))

ggplot(phom_k1_pow90, aes(x=pp.bias)) + 
  geom_vline(xintercept = 0, lty= 3, color = "darkgrey")+
  geom_density(fill="purple", colour="black", alpha=.2) +
  theme_bw() +
  geom_vline(xintercept = median(phom$pp.bias), lty= 1, color = "purple") +
  #coord_cartesian(xlim = c(-.5, .2), ylim = c(0,11)) +
  labs(x="\nbias in PET-PEESE estimate (d)", y = "Density\n", title = "Homogenous ES\n+ Publication Bias\n+ Guaranteed k=1 @ 90% power") +
  geom_segment(aes(x = 0, y = 1, xend = median(phom$pp.bias), yend = 1), arrow = arrow(length = unit(0.4, "cm")), color="red") +
  annotate("text", x = median(phom$pp.bias)/2, y = 1.4, label = round(median(phom$pp.bias), 3))

ggplot(phom_k3_pow80, aes(x=pp.bias)) + 
  geom_vline(xintercept = 0, lty= 3, color = "darkgrey")+
  geom_density(fill="purple", colour="black", alpha=.2) +
  theme_bw() +
  geom_vline(xintercept = median(phom$pp.bias), lty= 1, color = "purple") +
  #coord_cartesian(xlim = c(-.5, .2), ylim = c(0,11)) +
  labs(x="\nbias in PET-PEESE estimate (d)", y = "Density\n", title = "Homogenous ES\n+ Publication Bias\n+ Guaranteed k=3 @ 80% power") +
  geom_segment(aes(x = 0, y = 1, xend = median(phom$pp.bias), yend = 1), arrow = arrow(length = unit(0.4, "cm")), color="red") +
  annotate("text", x = median(phom$pp.bias)/2, y = 1.4, label = round(median(phom$pp.bias), 3))

ggplot(phom_k10_pow80, aes(x=pp.bias)) + 
  geom_vline(xintercept = 0, lty= 3, color = "darkgrey")+
  geom_density(fill="purple", colour="black", alpha=.2) +
  theme_bw() +
  geom_vline(xintercept = median(phom$pp.bias), lty= 1, color = "purple") +
  #coord_cartesian(xlim = c(-.5, .2), ylim = c(0,11)) +
  labs(x="\nbias in PET-PEESE estimate (d)", y = "Density\n", title = "Homogenous ES\n+ Publication Bias\n+ Guaranteed k=10 @ 80% power") +
  geom_segment(aes(x = 0, y = 1, xend = median(phom$pp.bias), yend = 1), arrow = arrow(length = unit(0.4, "cm")), color="red") +
  annotate("text", x = median(phom$pp.bias)/2, y = 1.4, label = round(median(phom$pp.bias), 3))

### Credits ####
# initial sloppy work by Will Gervais. 
# Suggestions for more realistic effect and sample size functions by Evan Carter. 
# Fixing all Will's mistakes and streamlining by Felix Schönbrodt
