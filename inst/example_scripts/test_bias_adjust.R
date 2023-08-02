## devtools::install_github("Cole-Monnahan-NOAA/wham", ref='add_recruit_reporting')
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
input1 <- prepare_wham_input(asap3, recruit_model=2,
                             model_name="Ex 1: SNEMA Yellowtail Flounder",
                             selectivity=list(model=rep("age-specific",3),
                                              re=rep("none",3),
                                              initial_pars=list(c(0.5,0.5,0.5,1,1,0.5),c(0.5,0.5,0.5,1,0.5,0.5),c(0.5,1,1,1,0.5,0.5)),
                                              fix_pars=list(4:5,4,2:4)),
                             NAA_re = list(sigma="rec", cor="iid"),
                             )
#input1$map$logit_selpars[14] <- NA
m1 <- fit_wham(input1, do.osa = FALSE, do.retro=FALSE)

input2 <- input1
input2$par <- m1$parList
input2$map$logit_selpars <- as.factor(NA*input2$map$logit_selpars)
input2$map$log_N1_pars<- as.factor(NA*input2$par$log_N1_pars)
## input2$map$F_devs<- as.factor(NA*input2$par$F_devs)
effN <- c(100, floor(seq(1,20, len=43)))
#input2$data$use_index_paa[2:30,] <- 0
#input2$data$index_Neff <- cbind(effN,effN)*0+1
input2$data$use_index_paa <- 0*input2$data$use_index_paa
input2$data$use_catch_paa[2:30] <- 0
input2$data$use_indices[5:15,] <- 0
input2$data$catch_Neff[,1] <- effN
input2$data$agg_index_sigma <- input2$data$agg_index_sigma

## simulate from previous fit but different effN
m2 <- fit_wham(input2, do.osa = FALSE, do.retro=FALSE, do.fit=FALSE)
set.seed(2314)
input3 <- input2
input3$data <- m2$simulate(complete=TRUE)
input3$data$bias_correct_pe <- 1

m3 <- fit_wham(input3, do.osa=FALSE, do.retro=FALSE)
sdrep.corrected <- sdreport(m3, bias.correct=TRUE)
est <- summary(sdrep.corrected) |> as.data.frame() %>% rownames_to_column(var='par') %>%
  janitor::clean_names() %>% filter(grepl('recruits', par)) %>%
  mutate(CV=std_error/estimate)
recdevs <- summary(sdrep.corrected) |> as.data.frame() %>% rownames_to_column(var='par') %>%
  janitor::clean_names() %>% filter(grepl('log_NAA', par)) %>% filter(!grepl('rep|sigma',par))%>%
  mutate(CV=std_error/estimate)
est
ggplot(recdevs, aes(x=1:43, y=estimate, ymin=estimate-2*std_error,
                ymax=estimate+2*std_error)) + geom_pointrange() +
  ylab('log recruits')

mcmc <- tmbstan(obj=m3, iter=1000, init='last.par.best',
                chains=4, cores=4, open_progress=FALSE,
                control=list(max_treedepth=8))
## launch_shinystan(mcmc)
## get posterior means
post <- as.data.frame(mcmc)
recruits <- matrix(NA, nrow=nrow(post), ncol=44)
for(ii in 1:nrow(post)){
  recruits[ii,] <- m3$report(par=post[ii,-ncol(post)])$recruits
}
recruits.mcmc <- apply(recruits, 2,median)

sigmaR <- exp(m2$parList$log_NAA_sigma)
recruits.true <- exp(m2$parList$log_NAA[,1])
## compare them all
plot(est$estimate)
lines(est$est_bias_correct)
lines(recruits.mcmc, col=2, lty=3)



xx <- data.frame(uncorrected=est$estimate,
      corrected=est$est_bias_correct, mle_se=est$std_error,
      MCMC=recruits.mcmc) %>%
  mutate(diff_uncorr=uncorrected-MCMC, diff_corr=corrected-MCMC, corr_better=abs(diff_corr)<abs(diff_uncorr))
xx |> round(0) |> select(1:3)
