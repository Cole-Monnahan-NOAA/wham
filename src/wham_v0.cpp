#define TMB_LIB_INIT R_init_wham
#include <TMB.hpp>
#include <iostream>
#include "helper_functions.hpp"

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_INTEGER(n_years_catch);
  DATA_INTEGER(n_years_indices);
  DATA_INTEGER(n_years_model);
  DATA_INTEGER(n_ages);
  DATA_INTEGER(n_fleets);
  DATA_INTEGER(n_indices);
  DATA_INTEGER(n_selblocks);
  DATA_IVECTOR(selblock_models);
  DATA_IMATRIX(selblock_pointer_fleets);
  DATA_IMATRIX(selblock_pointer_indices);
  DATA_IVECTOR(age_comp_model_fleets);
  DATA_IVECTOR(n_age_comp_pars_fleets);
  DATA_IVECTOR(age_comp_model_indices);
  DATA_IVECTOR(n_age_comp_pars_indices);
  DATA_VECTOR(fracyr_SSB);
  DATA_MATRIX(mature);
  DATA_IVECTOR(waa_pointer_fleets);
  DATA_INTEGER(waa_pointer_totcatch);
  DATA_IVECTOR(waa_pointer_indices);
  DATA_INTEGER(waa_pointer_ssb);
  DATA_INTEGER(waa_pointer_jan1);
  DATA_ARRAY(waa);
  DATA_MATRIX(agg_catch);
  DATA_MATRIX(agg_catch_sigma);
  DATA_ARRAY(catch_paa); //n_fleets x n_years x n_ages
  DATA_IMATRIX(use_catch_paa);
  DATA_MATRIX(catch_Neff);
  DATA_IMATRIX(catch_aref);
  DATA_IVECTOR(units_indices);
  DATA_MATRIX(fracyr_indices);
  DATA_MATRIX(agg_indices);
  DATA_IMATRIX(use_indices);
  DATA_MATRIX(agg_index_sigma);
  DATA_IVECTOR(units_index_paa);
  DATA_ARRAY(index_paa); //n_indices x n_years x n_ages
  DATA_IMATRIX(use_index_paa);
  DATA_MATRIX(index_Neff);
  DATA_IMATRIX(index_aref);
  DATA_VECTOR(q_lower);
  DATA_VECTOR(q_upper);
  DATA_MATRIX(selpars_lower);
  DATA_MATRIX(selpars_upper);
  DATA_INTEGER(n_NAA_sigma);
  DATA_IVECTOR(NAA_sigma_pointers);
  DATA_INTEGER(recruit_model);
  DATA_INTEGER(n_M_re);
  DATA_IVECTOR(MAA_pointer); //n_ages
  DATA_IVECTOR(M_sigma_par_pointer); //n_M_re
  DATA_INTEGER(M_model); //0: just age-specific M, 1: Lorenzen M decline with age/weight
  DATA_INTEGER(N1_model); //0: just age-specific numbers at age, 1: 2 pars: log_N_{1,1}, log_F0, age-structure defined by equilibrium NAA calculations
  DATA_INTEGER(use_M_re);
  DATA_INTEGER(use_NAA_re);
  DATA_INTEGER(use_b_prior);
  DATA_INTEGER(random_recruitment);
  DATA_INTEGER(which_F_age); //which age of F to use for full total F for msy/ypr calculations
  DATA_INTEGER(use_steepness); // which parameterization to use for BH/Ricker S-R, if needed.
  DATA_INTEGER(bias_correct_pe); //bias correct lognormal process error?
  DATA_INTEGER(bias_correct_oe); //bias correct lognormal observation error?
  DATA_IVECTOR(Fbar_ages);
  DATA_INTEGER(simulate_state); //if 1 then state parameters will be simulated
  DATA_SCALAR(percentSPR); //percentage to use for SPR-based reference points

  // data for one-step-ahead (OSA) residuals
  DATA_VECTOR(obsvec); // vector of all observations for OSA residuals
  DATA_VECTOR_INDICATOR(keep, obsvec); // for OSA residuals
  DATA_IMATRIX(keep_C); // indices for catch obs, can loop years/fleets with keep(keep_C(y,f))
  DATA_IMATRIX(keep_I);
  DATA_IARRAY(keep_Cpaa);
  // DATA_IARRAY(keep_Ipaa);

  // data for environmental covariate(s), Ecov
  DATA_INTEGER(n_Ecov); // also = 1 if no Ecov
  DATA_INTEGER(n_years_Ecov); // num years in Ecov  process model
  DATA_IMATRIX(Ecov_use_obs); // all 0 if no Ecov
  DATA_MATRIX(Ecov_obs);
  DATA_MATRIX(Ecov_obs_sigma);
  DATA_IVECTOR(Ecov_lag);
  DATA_IVECTOR(Ecov_how); // 1 = controlling, 2 = limiting, 3 = lethal, 4 = masking, 5 = directive
  DATA_IVECTOR(Ecov_where); // 1 = recruit, 2 = growth, 3 = mortality
  DATA_IVECTOR(Ecov_model); // 0 = no Ecov, 1 = RW, 2 = AR1
  DATA_IVECTOR(n_Ecov_pars); // RW: 1 par (sig), AR1: 2 par (phi, sig)
  DATA_INTEGER(Ecov_recruit); // Ecov index to use for recruitment
  DATA_INTEGER(Ecov_growth); // Ecov index to use for growth
  DATA_INTEGER(Ecov_mortality); // Ecov index to use for mortality
  DATA_INTEGER(year1_Ecov); // first year Ecov
  DATA_INTEGER(year1_model); // first year model
  DATA_IVECTOR(ind_Ecov_out_start); // index of Ecov_x to use for Ecov_out (operates on pop model, lagged)
  DATA_IVECTOR(ind_Ecov_out_end); // index of Ecov_x to use for Ecov_out (operates on pop model, lagged)

  // parameters - general
  PARAMETER_VECTOR(mean_rec_pars);
  PARAMETER_VECTOR(logit_q);
  PARAMETER_VECTOR(log_F1);
  PARAMETER_MATRIX(F_devs);
  PARAMETER_VECTOR(log_N1_pars); //length = n_ages or 2
  PARAMETER_VECTOR(log_NAA_sigma);
  PARAMETER_MATRIX(logit_selpars); //n_selblocks x n_ages + 6 (n_ages for by age, 2 for logistic, 4 for double-logistic)
  PARAMETER_VECTOR(catch_paa_pars);
  PARAMETER_VECTOR(index_paa_pars);
  PARAMETER_MATRIX(log_NAA);
  PARAMETER_VECTOR(M_pars1);
  PARAMETER_VECTOR(M_sigma_pars); //up to n_M_re
  PARAMETER_MATRIX(M_re); //n_years-1 x n_M_re
  PARAMETER(log_b);
  PARAMETER_VECTOR(log_R); //n_years-1, if used.
  PARAMETER(log_R_sigma);
  PARAMETER_VECTOR(log_catch_sig_scale) //n_fleets
  PARAMETER_VECTOR(log_index_sig_scale) //n_indices

  // parameters - environmental covariate ("Ecov")
  PARAMETER_MATRIX(Ecov_re); // nrows = n_years_Ecov, ncol = N_Ecov
  PARAMETER_VECTOR(Ecov_beta); // one for each ecov, beta_R in eqns 4-5, Miller et al. (2016)
  PARAMETER_MATRIX(Ecov_process_pars); // nrows = RW: 2 par (sig, Ecov1), AR1: 3 par (mu, phi, sig); ncol = N_ecov

  Type nll= 0.0; //negative log-likelihood
  vector<int> any_index_age_comp(n_indices);
  vector<int> any_fleet_age_comp(n_fleets);
  vector<Type> SSB(n_years_model);
  matrix<Type> F(n_years_model,n_fleets);
  matrix<Type> log_F(n_years_model,n_fleets);
  array<Type> pred_CAA(n_years_model,n_fleets,n_ages);
  array<Type> pred_catch_paa(n_years_model,n_fleets,n_ages);
  matrix<Type> pred_catch(n_years_model,n_fleets);
  matrix<Type> log_pred_catch(n_years_model,n_fleets);
  array<Type> pred_IAA(n_years_model,n_indices,n_ages);
  array<Type> pred_index_paa(n_years_model,n_indices,n_ages);
  matrix<Type> pred_indices(n_years_model,n_indices);
  matrix<Type> NAA(n_years_model,n_ages);
  matrix<Type> pred_NAA(n_years_model,n_ages);
  array<Type> FAA(n_years_model,n_fleets,n_ages);
  array<Type> log_FAA(n_years_model,n_fleets,n_ages);
  matrix<Type> FAA_tot(n_years_model,n_ages);
  matrix<Type> ZAA(n_years_model,n_ages);
  array<Type> QAA(n_years_model,n_indices,n_ages);
  matrix<Type> selblocks(n_selblocks,n_ages);
  vector<Type> q(n_indices);
  vector<Type> t_paa(n_ages);
  vector<Type> t_pred_paa(n_ages);
  matrix<Type> selpars(n_selblocks,n_ages+6);
  //Type SR_a, SR_b, SR_R0, SR_h;
  for(int i = 0; i < n_selblocks; i++) for(int j = 0; j < n_ages + 6; j++)
  {
    selpars(i,j) = selpars_lower(i,j) +  (selpars_upper(i,j)-selpars_lower(i,j))/(1.0+exp(-logit_selpars(i,j)));
  }
  REPORT(selpars);
  selblocks = get_selblocks(n_ages, n_selblocks, selpars, selblock_models);

  for(int i = 0; i < n_indices; i++)
  {
    any_index_age_comp(i) = 0;
    for(int y = 0; y < n_years_indices; y++) if(use_index_paa(y,i) == 1) any_index_age_comp(i) = 1;
  }
  for(int i = 0; i < n_fleets; i++)
  {
    any_fleet_age_comp(i) = 0;
    for(int y = 0; y < n_years_catch; y++) if(use_catch_paa(y,i) == 1) any_fleet_age_comp(i) = 1;
  }
  vector<Type> sigma2_log_NAA = exp(log_NAA_sigma*2.0);

  // Environmental covariate process model(s)
  matrix<Type> Ecov_x(n_years_Ecov, n_Ecov); // 'true' estimated Ecov (x_t in Miller et al. 2016 CJFAS)
  matrix<Type> Ecov_out(n_years_model, n_Ecov); // Pop model uses Ecov_out(t) for processes in year t (Ecov_x shifted by lag and padded)
  Type nll_Ecov = 0.0;

  if(Ecov_model(0) == 0){ // no Ecov
    for(int y = 0; y < n_years_model; y++) Ecov_out(y,0) = Type(0); // set Ecov_out = 0
  } else { // yes Ecov
    for(int i = 0; i < n_Ecov; i++){ // loop over Ecovs
      // Ecov model option 1: RW
      if(Ecov_model(i) == 1){
        Type Ecov_sig; // sd (sig_x in Eq1, pg 1262, Miller et al. 2016)
        Ecov_sig = exp(Ecov_process_pars(0,i));
        Type Ecov1; // Ecov_x in year 1 (fixed effect)
        Ecov1 = Ecov_process_pars(1,i);

        Ecov_x(0,i) = Ecov1;
        nll_Ecov -= dnorm(Ecov_re(0,i), Ecov1, Ecov_sig, 1);
        // Ecov_x(0,i) = Ecov_re(0,i); // initial year value (x_1, pg 1262, Miller et al. 2016)
        // nll_Ecov -= dnorm(Ecov_x(0,i), Type(0), Type(1000), 1);
        for(int y = 1; y < n_years_Ecov; ++y){
          nll_Ecov -= dnorm(Ecov_re(y,i), Ecov_re(y-1,i), Ecov_sig, 1);
          Ecov_x(y,i) = Ecov_re(y,i);
        }

        // Type Ecov1; // initial year value ("given a fixed effect x_1", pg 1262, Miller et al. 2016)
        // Type Ecov_sig; // sd (sig_x in Eq1, pg 1262, Miller et al. 2016)
        // Ecov1 = Ecov_process_pars(0);
        // Ecov_sig = exp(Ecov_process_pars(1));

        // Ecov_y(0) = Ecov1;
        // nll_Ecov -= dnorm(Ecov_y(0), Ecov1, Ecov_sig, 1);
        // // nll_Ecov -= -(half * (log(two*M_PI) + square(Ecov(0,i) - Ecov1(i))/sigma2_Ecov(i)) + log_Ecov_sigma(i));
        // for(int y = 1; y < n_years_Ecov - 1; y++)
        // {
        //   Ecov_y(y) = Ecov_re(y-1); // Ecov_re is one shorter than Ecov_y bc first value is Ecov1, fixed effect
        //   nll_Ecov -= dnorm(Ecov_re(y), Ecov_re(y-1), Ecov_sig, 1);
        //   // nll_Ecov -= -(half * (log(two*M_PI) + square(Ecov(y,i) - Ecov(y-1,i))/sigma2_Ecov(i)) + log_Ecov_sigma(i));
        // }
        // Ecov_y(n_years_Ecov-1) = Ecov_re(n_years_Ecov-2); // why does the last year not contribute to nll?
      }

      // Ecov model option 2: AR1
      if(Ecov_model(i) == 2){
        Type Ecov_mu; // mean
        Type Ecov_phi; // autocorrelation
        Type Ecov_sig; // conditional sd
        Ecov_mu = Ecov_process_pars(0,i);
        Ecov_phi = -Type(1) + Type(2)/(Type(1) + exp(-Ecov_process_pars(1,i)));
        Ecov_sig = exp(Ecov_process_pars(2,i));
        for(int y = 0; y < n_years_Ecov; y++) Ecov_x(y,i) = Ecov_mu + Ecov_re(y,i);

        nll_Ecov -= dnorm(Ecov_re(0,i), Type(0), Ecov_sig*exp(-Type(0.5) * log(Type(1) - pow(Ecov_phi,Type(2)))), 1);
        for(int y = 1; y < n_years_Ecov; y++) nll_Ecov -= dnorm(Ecov_re(y,i), Ecov_phi * Ecov_re(y-1,i), Ecov_sig, 1);
      }
    } // end loop over Ecovs
  }
  nll += nll_Ecov;

  // Environmental covariate observation model
  Type nll_Ecov_obs = Type(0);
  for(int i = 0; i < n_Ecov; i++){
    for(int y = 0; y < n_years_model; y++){
      if(Ecov_use_obs(y,i) == 1){
        nll_Ecov_obs -= dnorm(Ecov_obs(y,i), Ecov_x(y,i), Ecov_obs_sigma(y,i), 1);
      }
    }
  }
  nll += nll_Ecov_obs;

  // Lag environmental covariates
  // Then use Ecov_out(t) for processes in year t, instead of Ecov_x
  // matrix<Type> Ecov_out(n_years_model, n_Ecov);
  for(int i=0; i < n_Ecov; i++){
    int ct = 0;
    for(int y=ind_Ecov_out_start(i); y < ind_Ecov_out_end(i)+1; y++){
      Ecov_out(ct,i) = Ecov_x(y,i);
      ct++;
    }
  }

  // Natural mortality estimation
  if(use_M_re == 1)
  {
    vector<Type> sigma_M(n_M_re);
    for(int i = 0; i < n_M_re; i++) sigma_M(i) = exp(M_sigma_pars(M_sigma_par_pointer(i)-1));
    matrix<Type> nll_M(n_years_model-1, n_M_re);
    nll_M.setZero();

    if(M_model == 0) for(int i = 0; i < n_M_re; i++)
    {
      Type mu = M_pars1(i);
      if(bias_correct_pe == 1) mu -= 0.5*exp(2*log(sigma_M(i)));
      nll_M(0,i) -= dnorm(M_re(0,i), mu, sigma_M(i), 1);
      SIMULATE if(simulate_state == 1) M_re(0,i) = rnorm(mu, sigma_M(i));
      for(int y = 1; y < n_years_model - 1; y++)
      {
        mu = M_re(y-1,i);
        if(bias_correct_pe == 1) mu -= 0.5*exp(2*log(sigma_M(i)));
        nll_M(y,i) -= dnorm(M_re(y,i), mu, sigma_M(i), 1);
        SIMULATE if(simulate_state == 1) M_re(y,i) = rnorm(mu, sigma_M(i));
      }
    }
    else
    {
      for(int i = 0; i < n_M_re; i++)
      {
        Type mu = M_pars1(i);
        if(bias_correct_pe == 1) mu -= 0.5*exp(2*log(sigma_M(i)));
        nll_M(0,i) -= dnorm(M_re(0,i), mu, sigma_M(i), 1);
        SIMULATE if(simulate_state == 1) M_re(0,i) = rnorm(mu, sigma_M(i));
      }
      for(int i = 0; i < n_M_re; i++) for(int y = 1; y < n_years_model - 1; y++)
      {
        Type mu = M_re(y-1,i);
        if(bias_correct_pe == 1) mu -= 0.5*exp(2*log(sigma_M(i)));
        nll_M(y,i) -= dnorm(M_re(y,i), mu, sigma_M(i), 1);
        SIMULATE if(simulate_state == 1) M_re(y,i) = rnorm(mu, sigma_M(i));
      }
    }
    SIMULATE REPORT(M_re);
    REPORT(nll_M);
    nll += nll_M.sum();
  }
  //see(nll);

  matrix<Type> MAA(n_years_model,n_ages);
  for(int i = 0; i < n_ages; i++)
  {
    if(M_model == 0) //M by age
    {
      MAA(0,i) = exp(M_pars1(MAA_pointer(i)-1));
      for(int y = 1; y < n_years_model; y++) MAA(y,i) = exp(M_re(y-1,MAA_pointer(i)-1));
    }
    else //M_model == 1, allometric function of weight
    {
      MAA(0,i) = exp(M_pars1(0) - exp(log_b) * log(waa(waa_pointer_jan1-1,0,i)));
      if(use_M_re == 1) for(int y = 1; y < n_years_model; y++) MAA(y,i) = exp(M_re(y-1,0) - exp(log_b) * log(waa(waa_pointer_jan1-1,y,i)));
      else for(int y = 1; y < n_years_model; y++) MAA(y,i) = exp(M_pars1(0) - exp(log_b) * log(waa(waa_pointer_jan1-1,y,i)));
    }
  }

  for(int i = 0; i < n_indices; i++)
  {
    q(i) = q_lower(i) + (q_upper(i) - q_lower(i))/(1 + exp(-logit_q(i)));
    for(int y = 0; y < n_years_model; y++)
    {
      for(int a = 0; a < n_ages; a++) QAA(y,i,a) = q(i) * selblocks(selblock_pointer_indices(y,i)-1,a);
    }
  }
  FAA_tot.setZero();
  for(int f = 0; f < n_fleets; f++)
  {
    log_F(0,f) = log_F1(f);
    F(0,f) = exp(log_F(0,f));
    for(int a = 0; a < n_ages; a++)
    {
      FAA(0,f,a) = F(0,f) * selblocks(selblock_pointer_fleets(0,f)-1,a);
      log_FAA(0,f,a) = log(FAA(0,f,a));
      FAA_tot(0,a) = FAA_tot(0,a) + FAA(0,f,a);
    }
    for(int y = 1; y < n_years_model; y++)
    {
      log_F(y,f) = log_F(y-1,f) + F_devs(y-1,f);
      F(y,f) = exp(log_F(y,f));
      for(int a = 0; a < n_ages; a++)
      {
        FAA(y,f,a) = F(y,f) * selblocks(selblock_pointer_fleets(y,f)-1,a);
        log_FAA(y,f,a) = log(FAA(y,f,a));
        FAA_tot(y,a) = FAA_tot(y,a) + FAA(y,f,a);
      }
    }
  }

  ZAA = FAA_tot + MAA;

  SSB.setZero();
  //year 1
  for(int a = 0; a < n_ages; a++)
  {
    if(N1_model == 0) NAA(0,a) = exp(log_N1_pars(a));
    else
    {
      if(a==0) NAA(0,0) = exp(log_N1_pars(0));
      else
      {
        if(a == n_ages-1) NAA(0,a) = NAA(0,a-1)/(1.0 + exp(-MAA(0,a) - exp(log_N1_pars(1)) * FAA_tot(0,a)/FAA_tot(0,n_ages-1)));
        else NAA(0,a) = NAA(0,a-1)* exp(-MAA(0,a) -  exp(log_N1_pars(1)) * FAA_tot(0,a)/FAA_tot(0,n_ages-1));
      }
    }
    SSB(0) += NAA(0,a) * waa(waa_pointer_ssb-1,0,a) * mature(0,a) * exp(-ZAA(0,a)*fracyr_SSB(0));
    pred_NAA(0,a) = NAA(0,a);
  }

  //after year 1
  //get predicted numbers at age
  vector<Type> M(n_ages), mat(n_ages), waassb(n_ages), log_SPR0(n_years_model);
  int nh = 1, na = n_years_model;
  if(use_steepness != 0)
  {
    nh = n_years_model;
    na = 1;
  }
  vector<Type> log_SR_a(na), log_SR_b(na), SR_h(nh), SR_R0(nh);

  for(int y = 0; y < n_years_model; y++)
  {
    for(int a = 0; a < n_ages; a++)
    {
      M(a) = MAA(y,a);
      waassb(a) = waa(waa_pointer_ssb-1,y,a);
      mat(a) = mature(y,a);
    }
    log_SPR0(y) = log(get_SPR_0(M, mat, waassb, fracyr_SSB(y)));
  }
  REPORT(log_SPR0);
  ADREPORT(log_SPR0);
  if(recruit_model > 2) //BH or Ricker SR
  {
    if(recruit_model == 3) //BH stock recruit
    {
      if(use_steepness == 1)
      {
        SR_h.fill(0.2 + 0.8/(1+exp(-mean_rec_pars(0)))); //SR_a * SPR0/(4.0 + SR_a*SPR0);
        SR_R0.fill(exp(mean_rec_pars(1))); //(SR_a - 1/SPR0) / SR_b;
        log_SR_a = log(4 * SR_h/(exp(log_SPR0)*(1 - SR_h)));
        log_SR_b = log((5*SR_h - 1)/((1-SR_h)*SR_R0*exp(log_SPR0)));
      }
      else
      {
        log_SR_a.fill(mean_rec_pars(0));
        log_SR_b.fill(mean_rec_pars(1));
        SR_h = exp(log_SR_a) * exp(log_SPR0)/(4.0 + exp(log_SR_a + log_SPR0));
        SR_R0 = (exp(log_SR_a) - 1/exp(log_SPR0)) / exp(log_SR_b);
      }
    }
    if(recruit_model>3) //Ricker stock recruit
    {
      if(use_steepness == 1)
      {
        SR_h.fill(0.2 + exp(mean_rec_pars(0)));
        SR_R0.fill(exp(mean_rec_pars(1)));
        log_SR_a = 1.25*log(5*SR_h) - log_SPR0;
        log_SR_b = log(1.25*log(5*SR_h)/(SR_R0*exp(log_SPR0)));
      }
      else
      {
        log_SR_a.fill(mean_rec_pars(0));
        log_SR_b.fill(mean_rec_pars(1));
        SR_h = 0.2 * exp(0.8*log(exp(log_SR_a) * exp(log_SPR0)));
        SR_R0 = log(exp(log_SR_a + log_SPR0))/(exp(log_SR_b + log_SPR0));
      }
    }
    ADREPORT(log_SR_a);
    ADREPORT(log_SR_b);
    vector<Type> logit_SR_h = log(SR_h - 0.2) - log(1 - SR_h);
    vector<Type> log_SR_R0 = log(SR_R0);
    ADREPORT(logit_SR_h);
    ADREPORT(log_SR_R0);
    REPORT(log_SR_a);
    REPORT(log_SR_b);
    REPORT(logit_SR_h);
    REPORT(log_SR_R0);
  }
  matrix<Type> nll_NAA(n_years_model-1,n_ages);
  nll_NAA.setZero();
  vector<Type> nll_recruit(n_years_model-1);
  nll_recruit.setZero();

  // Population model (get NAA, numbers-at-age, for all years)
  for(int y = 1; y < n_years_model; y++)
  {
    // Expected recruitment
    if(recruit_model == 1) pred_NAA(y,0) = NAA(y-1,0) * exp(Ecov_beta(Ecov_recruit-1) * Ecov_out(y,Ecov_recruit-1)); // random walk
    else
    {
      if(recruit_model == 2) pred_NAA(y,0) = exp(mean_rec_pars(0) + Ecov_beta(Ecov_recruit-1) * Ecov_out(y,Ecov_recruit-1)); // random about mean
      else
      {
        if(recruit_model == 3) // BH stock recruit
        {
          if(Ecov_how(Ecov_recruit-1) == 1) // "controlling" = dens-indep mortality
          {
            if(use_steepness == 1) pred_NAA(y,0) = exp(log_SR_a(y-1) + Ecov_beta(Ecov_recruit-1) * Ecov_out(y,Ecov_recruit-1)) * SSB(y-1)/(1 + exp(log_SR_b(y-1))*SSB(y-1));
            else pred_NAA(y,0) = exp(log_SR_a(0) + Ecov_beta(Ecov_recruit-1) * Ecov_out(y,Ecov_recruit-1)) * SSB(y-1)/(1 + exp(log_SR_b(0))*SSB(y-1));
          }
          if(Ecov_how(Ecov_recruit-1) == 2) // "limiting" = carrying capacity
          {
            if(use_steepness == 1) pred_NAA(y,0) = exp(log_SR_a(y-1)) * SSB(y-1)/(1 + exp(log_SR_b(y-1) + Ecov_beta(Ecov_recruit-1) * Ecov_out(y,Ecov_recruit-1))*SSB(y-1));
            else pred_NAA(y,0) = exp(log_SR_a(0)) * SSB(y-1)/(1 + exp(log_SR_b(0) + Ecov_beta(Ecov_recruit-1) * Ecov_out(y,Ecov_recruit-1))*SSB(y-1));
          }
          if(Ecov_how(Ecov_recruit-1) == 4) // "masking" = metabolic/growth (decreases dR/dS)
          {
            if(use_steepness == 1) pred_NAA(y,0) = exp(log_SR_a(y-1)) * SSB(y-1)/(exp(Ecov_beta(Ecov_recruit-1) * Ecov_out(y,Ecov_recruit-1)) + exp(log_SR_b(y-1))*SSB(y-1));
            else pred_NAA(y,0) = exp(log_SR_a(0)) * SSB(y-1)/(exp(Ecov_beta(Ecov_recruit-1) * Ecov_out(y,Ecov_recruit-1)) + exp(log_SR_b(0))*SSB(y-1));
          }
        }
        else // recruit_model = 4, Ricker stock recruit
        {
          if(Ecov_how(Ecov_recruit-1) == 1) // "controlling" = dens-indep mortality
          {
            if(use_steepness == 1) pred_NAA(y,0) = exp(log_SR_a(y-1)) * SSB(y-1) * exp(-exp(log_SR_b(y-1)) * SSB(y-1) + Ecov_beta(Ecov_recruit-1) * Ecov_out(y,Ecov_recruit-1));
            else pred_NAA(y,0) = exp(log_SR_a(0)) * SSB(y-1) * exp(-exp(log_SR_b(0)) * SSB(y-1) + Ecov_beta(Ecov_recruit-1) * Ecov_out(y,Ecov_recruit-1));
          }
          if(Ecov_how(Ecov_recruit-1) == 4) // "masking" = metabolic/growth (decreases dR/dS)
          {
            if(use_steepness == 1) pred_NAA(y,0) = exp(log_SR_a(y-1)) * SSB(y-1) * exp(-exp(log_SR_b(y-1)) * SSB(y-1) * (1 + Ecov_beta(Ecov_recruit-1) * Ecov_out(y,Ecov_recruit-1)));
            else pred_NAA(y,0) = exp(log_SR_a(0)) * SSB(y-1) * exp(-exp(log_SR_b(0)) * SSB(y-1) * (1 + Ecov_beta(Ecov_recruit-1) * Ecov_out(y,Ecov_recruit-1)));
          }
        }
      }
    }

    //expected numbers at age after recruitment
    for(int a = 1; a < n_ages-1; a++) pred_NAA(y,a) = NAA(y-1,a-1) * exp(-ZAA(y-1,a-1));
    pred_NAA(y,n_ages-1) = NAA(y-1,n_ages-2) * exp(-ZAA(y-1,n_ages-2)) + NAA(y-1,n_ages-1) * exp(-ZAA(y-1,n_ages-1));

    if(use_NAA_re == 1) //random effects NAA, state-space model for all numbers at age
    {
      for(int a = 0; a < n_ages; a++)
      {
        Type mu = log(pred_NAA(y,a));
        if(bias_correct_pe == 1) mu -= 0.5*exp(2*log_NAA_sigma(NAA_sigma_pointers(a)-1));
        nll_NAA(y-1,a) -= dnorm(log_NAA(y-1,a), mu, exp(log_NAA_sigma(NAA_sigma_pointers(a)-1)), 1);
        SIMULATE
        {
          if(simulate_state == 1) log_NAA(y-1,a) = rnorm(mu, exp(log_NAA_sigma(NAA_sigma_pointers(a)-1)));
        }
        NAA(y,a) = exp(log_NAA(y-1,a));
      }
    }
    else //recruitments are still estimated parameters, fixed or random effects
    {
      if(random_recruitment == 1) //estimate recruitment as random effects, otherwise as fixed effects. pred_NAA(y,0) should be properly specified above in any case.
      {
        Type mu = log(pred_NAA(y,0));
        if(bias_correct_pe == 1) mu -= 0.5*exp(2*log_R_sigma);
        nll_recruit(y-1) -= dnorm(log_R(y-1), mu, exp(log_R_sigma), 1);
        SIMULATE if(simulate_state == 1) log_R(y-1) = rnorm(mu, exp(log_R_sigma));
      }
      NAA(y,0) = exp(log_R(y-1));
      //when random effects not used for all numbers at age, survival is deterministic.
      for(int a = 1; a < n_ages; a++) NAA(y,a) = pred_NAA(y,a);
    }
    for(int a = 0; a < n_ages; a++) SSB(y) += NAA(y,a) * waa(waa_pointer_ssb-1,y,a) * mature(y,a) * exp(-ZAA(y,a)*fracyr_SSB(y));
  }
  if(use_NAA_re == 1)
  {
    REPORT(nll_NAA);
    nll += nll_NAA.sum();
    SIMULATE REPORT(log_NAA);
  }

  if(random_recruitment == 1)
  {
    REPORT(nll_recruit);
    nll += nll_recruit.sum();
    SIMULATE REPORT(log_R);
  }

  if(use_b_prior == 1)
  {
    Type mu = log(0.305);
    if(bias_correct_pe == 1) mu -= 0.5*exp(2*log(0.08));
    Type lprior_b = dnorm(log_b, mu, Type(0.08), 1);
    SIMULATE
    {
      if(simulate_state == 1) log_b = rnorm(mu, Type(0.08));
      REPORT(log_b);
    }
    //see(lprior_b);
    REPORT(lprior_b);
    nll -= lprior_b;
  }
  //see(nll);

  matrix<Type> nll_agg_catch(n_years_catch,n_fleets), nll_catch_acomp(n_years_catch,n_fleets);
  nll_agg_catch.setZero();
  nll_catch_acomp.setZero();
  for(int y = 0; y < n_years_catch; y++)
  {
    int acomp_par_count = 0;
    for(int f = 0; f < n_fleets; f++)
    {
      pred_catch(y,f) = 0.0;
      Type tsum = 0.0;
      for(int a = 0; a < n_ages; a++)
      {
        pred_CAA(y,f,a) =  NAA(y,a) * FAA(y,f,a) * (1 - exp(-ZAA(y,a)))/ZAA(y,a);
        pred_catch(y,f) += waa(waa_pointer_fleets(f)-1,y,a) * pred_CAA(y,f,a);
        tsum += pred_CAA(y,f,a);
      }
      Type mu = log(pred_catch(y,f));
      Type sig = agg_catch_sigma(y,f)*exp(log_catch_sig_scale(f));
      if(bias_correct_oe == 1) mu -= 0.5*exp(2*log(sig));
      nll_agg_catch(y,f) -= keep(keep_C(y,f)) * dnorm(obsvec(keep_C(y,f)), mu, sig,1);
      // nll_agg_catch(y,f) -= keep(keep_C(y,f)) * dnorm(log(agg_catch(y,f)), mu, sig,1);
      SIMULATE agg_catch(y,f) = exp(rnorm(mu, sig));
      log_pred_catch(y,f) = log(pred_catch(y,f));
      if(any_fleet_age_comp(f) == 1)
      {
        vector<Type> acomp_pars(n_age_comp_pars_fleets(f));
        for(int j = 0; j < n_age_comp_pars_fleets(f); j++)
        {
          acomp_pars(j) = catch_paa_pars(acomp_par_count);
          acomp_par_count++;
        }
        if(use_catch_paa(y,f) == 1)
        {
          // vector<Type> t_keep(n_ages);
          for(int a = 0; a < n_ages; a++)
          {
            pred_catch_paa(y,f,a) = pred_CAA(y,f,a)/tsum;
            t_pred_paa(a) = pred_catch_paa(y,f,a);
            t_paa(a) = catch_paa(f,y,a);
            // t_paa(a) = obsvec(keep_Cpaa(f,y,a));
            // t_keep(a) = keep(keep_Cpaa(f,y,a));
          }
          nll_catch_acomp(y,f) -= get_acomp_ll(y, n_ages, catch_Neff(y,f), age_comp_model_fleets(f), t_paa, t_pred_paa, acomp_pars, catch_aref(y,f));
          // nll_catch_acomp(y,f) -= get_acomp_ll_osa(y, n_ages, catch_Neff(y,f), age_comp_model_fleets(f), t_paa, t_pred_paa, acomp_pars, catch_aref(y,f), t_keep);
          SIMULATE
          {
            t_paa = sim_acomp(y, n_ages, catch_Neff(y,f), age_comp_model_fleets(f), t_paa, t_pred_paa, acomp_pars, catch_aref(y,f));
            for(int a = 0; a < n_ages; a++) catch_paa(f,y,a) = t_paa(a);
          }
        }
      }
    }
  }
  SIMULATE REPORT(agg_catch);
  SIMULATE REPORT(catch_paa);
  REPORT(nll_agg_catch);
  nll += nll_agg_catch.sum();
  //see(nll);
  REPORT(nll_catch_acomp);
  nll += nll_catch_acomp.sum();
  //see(nll);

  matrix<Type> nll_agg_indices(n_years_catch,n_indices), nll_index_acomp(n_years_catch,n_indices);
  nll_agg_indices.setZero();
  nll_index_acomp.setZero();
  pred_indices.setZero();
  for(int y = 0; y < n_years_indices; y++)
  {
    int acomp_par_count = 0;
    for(int i = 0; i < n_indices; i++)
    {
      Type tsum = 0.0;
      for(int a = 0; a < n_ages; a++)
      {
        pred_IAA(y,i,a) =  NAA(y,a) * QAA(y,i,a) * exp(-ZAA(y,a) * fracyr_indices(y,i));
        if(units_indices(i) == 1) pred_indices(y,i) += waa(waa_pointer_indices(i)-1,y,a) * pred_IAA(y,i,a);
        else pred_indices(y,i) += pred_IAA(y,i,a);
      }
      for(int a = 0; a < n_ages; a++)
      {
        if(units_index_paa(i) == 1) pred_IAA(y,i,a) = waa(waa_pointer_indices(i)-1,y,a) * pred_IAA(y,i,a);
        tsum += pred_IAA(y,i,a);
      }

      if(use_indices(y,i) == 1)
      {
        Type mu = log(pred_indices(y,i));
        Type sig = agg_index_sigma(y,i)*exp(log_index_sig_scale(i));
        if(bias_correct_oe == 1) mu -= 0.5*exp(2*log(sig));
        // nll_agg_indices(y,i) -= dnorm(log(agg_indices(y,i)), mu, sig, 1);
        // nll_agg_indices(y,i) -= keep(keep_I(y,i)) * dnorm(log(agg_indices(y,i)), mu, sig, 1);
        nll_agg_indices(y,i) -= keep(keep_I(y,i)) * dnorm(obsvec(keep_I(y,i)), mu, sig, 1);
        SIMULATE agg_indices(y,i) = exp(rnorm(mu, sig));
      }
      if(any_index_age_comp(i) == 1)
      {
        vector<Type> acomp_pars(n_age_comp_pars_indices(i));
        for(int j = 0; j < n_age_comp_pars_indices(i); j++)
        {
          acomp_pars(j) = index_paa_pars(acomp_par_count);
          acomp_par_count++;
        }
        if(use_index_paa(y,i) > 0)
        {
          // vector<Type> t_keep(n_ages);
          for(int a = 0; a < n_ages; a++)
          {
            pred_index_paa(y,i,a) = pred_IAA(y,i,a)/tsum;
            t_pred_paa(a) = pred_index_paa(y,i,a);
            t_paa(a) = index_paa(i, y, a);
            // t_keep(a) = keep(keep_Ipaa(i,y,a));
          }
          nll_index_acomp(y,i) -= get_acomp_ll(y, n_ages, index_Neff(y,i), age_comp_model_indices(i), t_paa, t_pred_paa, acomp_pars, index_aref(y,i));
          // nll_index_acomp(y,i) -= get_acomp_ll_osa(y, n_ages, index_Neff(y,i), age_comp_model_indices(i), t_paa, t_pred_paa, acomp_pars, index_aref(y,i), t_keep);
          SIMULATE
          {
            t_paa = sim_acomp(y, n_ages, index_Neff(y,i), age_comp_model_indices(i), t_paa, t_pred_paa, acomp_pars, index_aref(y,i));
            for(int a = 0; a < n_ages; a++) index_paa(i,y,a) = t_paa(a);
          }
        }
      }
    }
  }
  SIMULATE REPORT(agg_indices);
  SIMULATE REPORT(index_paa);
  REPORT(nll_agg_indices);
  nll += nll_agg_indices.sum();
  //see(nll);
  REPORT(nll_index_acomp);
  nll += nll_index_acomp.sum();
  //see(nll);


  //////////////////////////////////////////
  //Still need to add in yearly vectors of biological inputs, make sure to calculate SR_a,SR_b vector or otherwise.
  //////////////////////////////////////////
  //calculate BRPs
  //First SPR-based proxies
  //Type percentSPR = 40;
  vector<Type> predR = pred_NAA.col(0);
  matrix<Type> SPR_res = get_SPR_res(MAA, FAA_tot, which_F_age, waa, waa_pointer_ssb, waa_pointer_totcatch, mature, percentSPR, predR, fracyr_SSB, log_SPR0);
  vector<Type> log_FXSPR = SPR_res.col(0);
  vector<Type> log_SSB_FXSPR = SPR_res.col(1);
  vector<Type> log_Y_FXSPR = SPR_res.col(2);
  vector<Type> log_SPR_FXSPR = SPR_res.col(3);
  vector<Type> log_YPR_FXSPR = SPR_res.col(4);
  matrix<Type> log_FXSPR_iter = SPR_res.block(0,5,n_years_model,10);
  REPORT(log_FXSPR_iter);
  REPORT(log_FXSPR);
  REPORT(log_SSB_FXSPR);
  REPORT(log_Y_FXSPR);
  REPORT(log_SPR_FXSPR);
  REPORT(log_YPR_FXSPR);
  ADREPORT(log_FXSPR);
  ADREPORT(log_SSB_FXSPR);
  ADREPORT(log_Y_FXSPR);

  //If stock-recruit models
  if(recruit_model > 2) //Beverton-Holt or Ricker selected
  {
    int n = 10;
    vector<Type> log_FMSY(n_years_model), log_FMSY_i(1), waacatch(n_ages), sel(n_ages);
    matrix<Type> log_FMSY_iter(n_years_model,n);
    log_FMSY_iter.col(0).fill(log(0.2)); //starting value
    vector<Type> log_YPR_MSY(n_years_model), log_SPR_MSY(n_years_model), log_R_MSY(n_years_model);
    Type SR_a, SR_b;
    for(int y = 0; y < n_years_model; y++)
    {
      for(int a = 0; a < n_ages; a++)
      {
        M(a) = MAA(y,a);
        sel(a) = FAA_tot(y,a)/FAA_tot(y,which_F_age-1); //have to look at FAA_tot to see where max F is.
        waassb(a) = waa(waa_pointer_ssb-1,y,a);
        waacatch(a) = waa(waa_pointer_totcatch-1, y, a);
        mat(a) = mature(y,a);
      }
      if(use_steepness == 1)
      {
        SR_a = exp(log_SR_a(y));
        SR_b = exp(log_SR_b(y));
      }
      else
      {
        SR_a = exp(log_SR_a(0));
        SR_b = exp(log_SR_b(0));
      }
      if(recruit_model == 3) //Beverton-Holt selected
      {
        sr_yield<Type> sryield(SR_a, SR_b, M, sel, mat, waassb, waacatch,fracyr_SSB(y),0);
        for (int i=0; i<n-1; i++)
        {
          log_FMSY_i(0) = log_FMSY_iter(y,i);
          vector<Type> grad_sr_yield = autodiff::gradient(sryield,log_FMSY_i);
          matrix<Type> hess_sr_yield = autodiff::hessian(sryield,log_FMSY_i);
          log_FMSY_iter(y,i+1) = log_FMSY_iter(y,i) - grad_sr_yield(0)/hess_sr_yield(0,0);
        }
      }
      else //Ricker selected
      {
        sr_yield<Type> sryield(SR_a, SR_b, M, sel, mat, waassb, waacatch,fracyr_SSB(y),1);
        for (int i=0; i<n-1; i++)
        {
          log_FMSY_i(0) = log_FMSY_iter(y,i);
          vector<Type> grad_sr_yield = autodiff::gradient(sryield,log_FMSY_i);
          matrix<Type> hess_sr_yield = autodiff::hessian(sryield,log_FMSY_i);
          log_FMSY_iter(y,i+1) = log_FMSY_iter(y,i) - grad_sr_yield(0)/hess_sr_yield(0,0);
        }
      }
      log_FMSY(y) = log_FMSY_iter(y,n-1);
      log_SPR_MSY(y) = log(get_SPR(log_FMSY(y), M, sel, mat, waassb, fracyr_SSB(y)));
      log_YPR_MSY(y) = log(get_YPR(log_FMSY(y), M, sel, waacatch));
      if(recruit_model == 3) log_R_MSY(y) = log((SR_a - 1/exp(log_SPR_MSY(y))) / SR_b); //bh
      else log_R_MSY(y) = log(log(SR_a) + log_SPR_MSY(y)) - log(SR_b) - log_SPR_MSY(y); //ricker
    }
    vector<Type> log_SSB_MSY = log_R_MSY + log_SPR_MSY;
    vector<Type> log_MSY = log_R_MSY + log_YPR_MSY;

    ADREPORT(log_FMSY);
    ADREPORT(log_SSB_MSY);
    ADREPORT(log_R_MSY);
    ADREPORT(log_MSY);
    ADREPORT(log_SPR_MSY);
    ADREPORT(log_YPR_MSY);
    REPORT(log_FMSY);
    REPORT(log_FMSY_iter);
    REPORT(log_SSB_MSY);
    REPORT(log_R_MSY);
    REPORT(log_MSY);
    REPORT(log_SPR_MSY);
    REPORT(log_YPR_MSY);
  }

  matrix<Type> log_FAA_tot = log(FAA_tot.array());
  matrix<Type> log_index_resid = log(agg_indices.block(0,0,n_years_model,n_indices).array()) - log(pred_indices.array());
  matrix<Type> log_catch_resid = log(agg_catch.block(0,0,n_years_model,n_fleets).array()) - log(pred_catch.array());
  vector<Type> log_SSB =  log(SSB);
  vector<Type> Fbar(n_years_model);
  Fbar.setZero();
  int n_Fbar_ages = Fbar_ages.size();
  for(int y = 0; y < n_years_model; y++) for(int a = 0; a < n_Fbar_ages; a++) Fbar(y) += FAA_tot(y,Fbar_ages(a)-1)/n_Fbar_ages;

  vector<Type> log_Fbar = log(Fbar);
  matrix<Type> log_NAA_rep = log(NAA.array());

  REPORT(NAA);
  REPORT(pred_NAA);
  REPORT(SSB);
  REPORT(selblocks);
  REPORT(MAA);
  REPORT(q);
  REPORT(QAA);
  REPORT(F);
  REPORT(FAA);
  REPORT(FAA_tot);
  REPORT(Fbar);
  REPORT(pred_catch);
  REPORT(pred_catch_paa);
  REPORT(pred_CAA);
  REPORT(pred_indices);
  REPORT(pred_index_paa);
  REPORT(pred_IAA);
  REPORT(Ecov_x);
  REPORT(Ecov_process_pars);
  REPORT(Ecov_re);
  REPORT(Ecov_beta);
  REPORT(mean_rec_pars);

  ADREPORT(log_F);
  ADREPORT(log_FAA);
  ADREPORT(log_FAA_tot);
  ADREPORT(log_Fbar);
  ADREPORT(log_NAA_rep);
  ADREPORT(log_SSB);
  ADREPORT(log_pred_catch);
  ADREPORT(pred_IAA);
  ADREPORT(log_index_resid);
  ADREPORT(log_catch_resid);
  ADREPORT(Ecov_x);
  ADREPORT(Ecov_process_pars);
  ADREPORT(Ecov_re);
  ADREPORT(Ecov_beta);
  ADREPORT(mean_rec_pars);

  REPORT(nll);
  REPORT(nll_Ecov);
  REPORT(nll_Ecov_obs);

  return nll;
}

