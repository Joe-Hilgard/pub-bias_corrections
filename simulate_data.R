# NOTE: To-do --------------------
#  Current funnel plot looks preposterous. How do I fix it?
# One possibility: distribution is too restricted in range on n.
# 2nd possibility: publication should favor z>0 over z<0
# 3rd possibility: perhaps I'm going about this backwards by simulating
#   k studies to get some random number of published results.
#   Should I instead start with k published results?
# END: To-do ---------------------

library(metafor)
library(dplyr)
# Loss function for Simonsohn's p-curve
# borrowed in turn from Lakens blog post
# does this mean I have to do everything in d?
loss=function(t_obs,df_obs,d_est) {  
  ncp_est=sqrt((df_obs+2)/4)*d_est
  tc=qt(.975,df_obs)              
  power_est=1-pt(tc,df_obs,ncp_est)
  p_larger=pt(t_obs,df=df_obs,ncp=ncp_est)  
  ppr=(p_larger-(1-power_est))/power_est 	
  KSD=ks.test(ppr,punif)$statistic       	
  return(KSD)     
}

# Select true effect size (Fisher's z)
z = 0
# Select degree of heterogeneity
sd = 0
# Select number of tests
k = 300
# Select average sample size and negbinom dispersion parameter
  # Note meanlog param gives median sample size
plot(1:600, dlnorm(1:600, 5, .6), typ='l')
n_median = 80
n_sd = .6
# Select probability of null or negative result being censored
# This isn't so simple either:
  # I think sig results get published with small n, but
  # nonsig results only get published with larger n
  # Note that in lieu of some p-hacking mechanism, 
    # p_censor will have to be quite large to account for base rate!
  # Simplest to say that p_censor includes self-censorship of null analyses
    # e.g. p-hacking away analyses that didn't turn out.
p_censor = .99

# Begin data frame by generating study sizes
temp = data.frame("size" = rlnorm(k, meanlog = log(n_median), 
                                   sdlog = n_sd))
# (may wish to impose some absolute minimum on sample size)
# say no study with less than n = 18. not even Psych Science would allow.
temp$size[temp$size < 18] = 18
# Next, determine true latent effect sizes (random effects)
# Note that if sd == 0, all z_true == z (fixed effects)
temp$z_true = z + rnorm(k, 0, sd)
# Sampling precision is determined by study size
# Recall that this may be an approximation!
temp$se = 1/sqrt(temp$size-3)
# Get observed effect sizes as z_true + sampling error
temp$z_hat = rnorm(k, temp$z_true, temp$se)
# Make t-value
temp$t = temp$z_hat / temp$se
# Make p-value
temp$p = 2*pt(temp$t, 
              df = temp$size - 2,
              lower.tail = F)
# Censor or publish null results
temp$pubbed = rbinom(k, 1, 1-p_censor)
# Publish all positive results
temp$pubbed[temp$p < .05] = 1
# Transform to logical for rma()
temp$pubbed = as.logical(temp$pubbed)

# Conduct fixed-effects naive estimate on only pubbed tests
fe_naive = rma(yi = z_hat, 
               sei = se, 
               method = "FE", 
               data=temp, 
               subset = pubbed)
# Conduct random-effects naive estimate on only pubbed tests
re_naive = rma(yi = z_hat, 
               sei = se, 
               method = "REML",
               data = temp, 
               subset = pubbed)
# Apply trim&fill to each
fe_taf = trimfill(fe_naive)
re_taf = trimfill(re_naive)
# Apply PET-PEESE
PET = rma(yi = z_hat, 
          sei = se, 
          mods = ~se,
          method = "FE",
          data = temp,
          subset = pubbed)
PEESE = rma(yi = z_hat, 
            sei = se, 
            mods = ~I(se^2),
            method = "FE",
            data = temp,
            subset = pubbed)
# p-curve estimate

# Rosenthal's Fail-safe N, just for laughs
failsafe_n = fsn(yi = z_hat,
                 sei = se,
                 type = "Rosenthal",
                 data = temp,
                 subset = pubbed)

# R-Index? Would be interesting to put Uli to the test.
# Percentage of significant tests - median(observed power)
# Get observed powers via inspection of naive meta-analysis
# then calculate R-index
percent_sig = sum(temp$p < .05 & temp$pubbed == T) / sum(temp$pubbed == T) 
# R_index = percent_sig - median(temp$observed_power)

# Retrieve meta-analysis results
results = data.frame("z_true" = z,
                     "count_pub" = sum(temp$pubbed),
                     "percent_sig" = percent_sig,
                     "fe_naive" = fe_naive$b, 
                     "re_naive" = re_naive$b,
                     "fe_taf" = fe_taf$b,
                     "re_taf" = re_taf$b,
                     "PET" = PET$b[1,1],
                     "PEESE" = PEESE$b[1,1],
                     "PETPEESE" = ifelse(PET$pval[1] < .05, 
                                         PEESE$b[1,1],
                                         PET$b[1,1]),
                     "failsafe_n" = failsafe_n$fsnum)

# Checks to see if pub-bias model working as intended
rma(yi = z_hat, 
    sei = se, 
    method = "FE", 
    data=temp) %>% 
  funnel(main = "All tests", 
         pch = ifelse(temp$pubbed==T, 19, 1))
rma(yi = z_hat, 
    sei = se, 
    method = "FE", 
    data=temp, 
    subset = !pubbed) %>% 
  funnel(main = "Censored tests")
funnel(fe_naive, main = "Published tests")
table(temp$p<.05, temp$pubbed, dnn=c("significant", "pubbed"))
