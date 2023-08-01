devtools::install_github("Cole-Monnahan-NOAA/wham", ref='add_recruit_reporting')
library(tmbstan)
library(shinystan)
library(wham)
library(TMB)
library(tidyverse)

## modified from example 1 of WHAM to test bias adjustment for
## recruitment. Had to add ADREPORT for recruitment to the
## template, hence the branch install above.

# create directory for analysis, E.g.,
#write.dir <- "/path/to/save/output"
if(!exists("write.dir")) write.dir = getwd()
if(!dir.exists(write.dir)) dir.create(write.dir)
setwd(write.dir)
# copy asap3 data file to working directory
wham.dir <- find.package("wham")
file.copy(from=file.path(wham.dir,"extdata","ex1_SNEMAYT.dat"), to=write.dir, overwrite=TRUE)
# read asap3 data file and convert to input list for wham
asap3 <- read_asap3_dat("ex1_SNEMAYT.dat")

# ---------------------------------------------------------------
# model 1
#   recruitment expectation (recruit_model): random about mean (no S-R function)
#   recruitment deviations (NAA_re): independent random effects
#   selectivity: age-specific (fix sel=1 for ages 4-5 in fishery, age 4 in index1, and ages 2-4 in index2)
input1 <- prepare_wham_input(asap3, recruit_model=2, model_name="Ex 1: SNEMA Yellowtail Flounder",
	                            selectivity=list(model=rep("age-specific",3),
                                	re=rep("none",3),
                                	initial_pars=list(c(0.5,0.5,0.5,1,1,0.5),c(0.5,0.5,0.5,1,0.5,0.5),c(0.5,1,1,1,0.5,0.5)),
                                	fix_pars=list(4:5,4,2:4)),
	                            NAA_re = list(sigma="rec", cor="iid"))
#input1$map$logit_selpars[14] <- NA
m1 <- fit_wham(input1, do.osa = FALSE, do.retro=FALSE, )

input2 <- input1
input2$par <- m1$parList
input2$map$logit_selpars <- as.factor(NA*input2$map$logit_selpars)
effN <- floor(seq(5,100, len=44))
effN[10:20] <- 1
input2$data$catch_Neff[,1] <- effN
input2$data$index_Neff <- cbind(effN,effN)*0+1

## simulate from previous fit but different effN
m2 <- fit_wham(input2, do.osa = FALSE, do.retro=FALSE, do.fit=FALSE)
input3 <- input2
input3$data <- m2$simulate(complete=TRUE)
input3$data$index_paa
input3$data$bias_correct_pe <- 1

m3 <- fit_wham(input3, do.osa=FALSE, do.retro=FALSE)
sdrep.corrected <- sdreport(m3, bias.correct=TRUE)
est <- summary(sdrep.corrected) |> as.data.frame() %>% rownames_to_column(var='par') %>%
  janitor::clean_names() %>% filter(grepl('recruits', par)) %>%
  mutate(CV=std_error/estimate)
est

mcmc <- tmbstan(obj=m3, iter=1000, init='last.par.best', chains=4, cores=4, open_progress=FALSE)
## launch_shinystan(mcmc)

## get posterior means
post <- as.data.frame(mcmc)
recruits <- matrix(NA, nrow=nrow(post), ncol=44)
for(ii in 1:nrow(post)){
  recruits[ii,] <- m3$report(par=post[ii,-ncol(post)])$recruits
}
recruits.mcmc <- colMeans(recruits)

sigmaR <- exp(m2$parList$log_NAA_sigma)
recruits.true <- exp(m2$parList$log_NAA[,1])
## compare them all
plot(est$estimate)
lines(est$est_bias_correct)
lines(colMeans(recruits), col=2, lty=3)

plot(est$estimate-recruits.true, ylab='Difference from true (recruits)')
lines(est$est_bias_correct-recruits.true)
lines(recruits.mcmc-recruits.true, col=2, lty=3)

cbind(uncorrected=est$estimate, corrected=est$est_bias_correct,
      MCMC=recruits.mcmc) |> round(1) |> head()
