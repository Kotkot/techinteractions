
//******************************************************************************************
// This is a two-sexes, multi-population, multi-fishery model fitted to 
// length frequency data (sex-aggregated and separated by population), 
// catch-rate indices (for both populations), 
// spawning abundance indices (for both populations) and recruitment index
// using Bayesian methods (MCMC). Catchabilities and sampling variances are set to their MLEs.
//
// cab.dat has the data 
// cab.ctl has the controls
// cab.pin has the initial values for the parameters
//
//**************************************************************************************
//
DATA_SECTION

  //*********COUNTERS*******************************************************************
  // All counters are declared globally (here) in this version
  int z // counters for size
  int l
  int g // counter for sex
  int a // counter for ages
  int p // counter for populations
  int t // counter for time
  int i
  int j
  int f // counter for fleet
  int pop // count for population

 // int Diags; //CRM commented this out- Diags does not seem to do anything anymore.
 // !! Diags = 1; //CRM commented this out.
 
 // *=========================================================================*
 
  
  //*****************MODEL DIMENSIONS ************************************************
  init_int styr                        // start year of the model
  !! log_input(styr);
  init_int endyr                       // end year of the model
  !!log_input(endyr);
  init_int npop                        // pops
  !!log_input(npop);
  init_int fleet                       // number of fleets (each one with a different selectivity function)
  !!log_input(fleet);
  init_int survey                      // number of surveys (each one with a different selectivity function)
  !!log_input(survey);
  int TotFleet;
  !! TotFleet = fleet + survey;        // total number of fleets (survey and fishery)
  init_ivector fleetpop(1,TotFleet)    // pointers for fleets to populations
  !!log_input(fleetpop);
  init_int sex                      // number of sex
  !!log_input(sex);
  init_int nages                       // plus group
  !!log_input(nages);
  init_int nlength                     // number of length classes
  !!log_input(nlength);
  init_number len_start                // first length bin
  !!log_input(len_start);
  init_number len_step                 // size of the length bins
  !!log_input(len_step);
  init_int age_specific_sel            // selectivity is age-specifc (1=Yes;0=No)?
  !!log_input(age_specific_sel);
  init_int wt_sel_bias                 // weight-at-age accounts for selectivity bias
  !!log_input(wt_sel_bias);

  int NvalMax;
  !! if (age_specific_sel == 0) NvalMax = nlength;
  !! if (age_specific_sel == 1) NvalMax = nages;
  
  // Vectors of ages, years and sizes
  ivector age_vector(0,nages) 
  ivector years(styr,endyr)            // vector of the years of the model
  vector size_vector(1,nlength)  //vector with all length bins
  
 LOCAL_CALCS
  years(styr) = styr;
  for (i=styr+1;i<=endyr;i++) years(i)=years(i-1)+1;
  for (a=0;a<=nages;a++) age_vector(a) = a;
  size_vector(1) = len_start;
  for (z=2;z<=nlength;z++) 
   size_vector(z) = size_vector(z-1) + len_step;
 END_CALCS
  
  //************BIOLOGICAL INFORMATION *************************
  // growth curve parameters, first females sex=1, then males sex=2
  // Von Bertalanffy growth function reparametrized as Lmin, Lmax, K intead of Linf, t0 and K
  init_vector Lmin(1,sex)
  !!log_input(Lmin);  
  init_vector Lmax(1,sex)
  !!log_input(Lmax);
  init_vector K(1,sex)
  !!log_input(K);
  init_vector CVLmin(1,sex)
  !!log_input(CVLmin);
  init_vector CVLmax(1,sex)
  !!log_input(CVLmax);
  
   // Age-length key
  3darray age_length_pop(1,sex,0,nages,1,nlength) 

   // length-weight parameters, same for both sex and populations
  init_number WL_intercept
  !!log_input(WL_intercept);
  init_number WL_slope
  !!log_input(WL_slope);
  init_3darray wt_input(1,2,1,sex,0,nages)               // 3d array: beginning of the year weight, middle of the year weight, sex, ages (females then males) 
  !!log_input(wt_input);
  vector wt_at_length(1,nlength)                            // weight as a function of length

  //CareyFix 28March2013: This used to be fecundity and now it is maturity

  //Relative fecundity at age
  init_vector Maturity(0,nages)
  !!log_input(Maturity);
     
  //Define fecundity as weight-at-age*Maturity-at-age:
  vector Fec(0,nages)
 LOCAL_CALCS
  for (a=0;a<=nages;a++) Fec(a) = Maturity(a)*wt_input(1,1,a); //Maturity * weight-at-age for females at the beginning of the year  
 END_CALCS
 //End CareyFix: Now we have Fec, which is Fecundity = Maturity-at-age x Weight-at-age
  !!cout << "Fec = " <<Fec << endl;
  !!cout << "Maturity = " <<Maturity <<endl;
  !!cout << "wt_input = " << wt_input(1,1) << endl;
  //*******CATCH biomass
  init_matrix catch_bio(styr,endyr,0,fleet)
  !!log_input(catch_bio);
  !! CheckFile << "Catches" << endl << catch_bio << endl;
  //!!cout<<"catch_bio = "<<catch_bio<<endl;
  // *************** Abundance data ***********************************************
  // *** Indx_type determines the nature of the index: 1+ CPUE; -1 SSB; -2 Recruitment
  init_int ncpue;
  !!log_input(ncpue);
  init_ivector indx_fleet(1,ncpue);
  !!log_input(indx_fleet);
  init_ivector indx_type(1,ncpue);
  !!log_input(indx_type);
  init_ivector nyrs_cpue(1,ncpue);
  !!log_input(nyrs_cpue);
  
  init_3darray TheAbund(1,ncpue,1,nyrs_cpue,1,4);
  !!log_input(TheAbund);

  // ************* Length-frequency data *****************************************
  // *** Lenf_type determines the nature of the length data
  init_int nlenf;
  !!log_input(nlenf);
  init_ivector lenf_fleet(1,nlenf);
  !!log_input(lenf_fleet);
  init_ivector lenf_type(1,nlenf);
  !!log_input(lenf_type);
  init_ivector nyrs_lenf(1,nlenf);
  !!log_input(nyrs_lenf);
  init_3darray TheLengths(1,nlenf,1,nyrs_lenf,-4,nlength);
  !!log_input(TheLengths);
  matrix nlensamp(1,nlenf,1,nyrs_lenf)               // Effective sample size for multinomial, ragged array
  vector offleng(1,nlenf)                            // Compute OFFSET for multinomial (i.e, value for the multinonial function for a perfect fit, or observed length frequency equal expected length frequency  
                                                     
  // ************* Age-structure data ********************************************
  // *** Agef_type determines the nature of the age data
  init_int nagef;
  !!log_input(nagef);
  init_ivector agef_fleet(1,nagef);
  !!log_input(agef_fleet);
  init_ivector agef_type(1,nagef);
  !!log_input(agef_type);
  init_ivector nyrs_agef(1,nagef);
  !!log_input(nyrs_agef);
  init_3darray TheAges(1,nagef,1,nyrs_agef,-5,nages);
  !!log_input(TheAges);
  matrix nagesamp(1,nagef,1,nyrs_agef);              // Effective sample size for multinomial, ragged array
  vector offage(1,nagef)                             // Compute OFFSET for multinomial (i.e, value for the multinonial function for a perfect fit, or observed length frequency equal expected length frequency
  
  //End of file indicator
  init_int endoffile
  !!log_input(endoffile);
  !!mark_input("END OF DATA FILE");
  !!cout<<"If you see 999, we got to the end of the data input sucessfully! "<<endoffile<<endl;
  !!cout<<"Otherwise, CAB will exit abnormally with code 1 so that you can correct the DAT file"<<endl;
  !! if (endoffile != 999) exit(1);
  
  //********CONTROLS********************************************************************
  //  Lets change the main datafile to a control data file for 
  //  specifying controls over the estimation (range of parameters, weights of data sets,
  //  switch on and off data sets,etc.
  
  !! ad_comm::change_datafile_name("cab.ctl");
  !!mark_input("BEGIN CTL FILE");
  init_int dummy;    //dummy==1 turns off all parameters, makes determinitic projections
  !!log_input(dummy);
  int phase_dummy;
  
  //estimate the extra the variability around the observed abundance indices (yes==1; no==0)
  init_int do_var; 
  !!log_input(do_var);  
  // Proportion of recruitment by area
  init_vector Prop_rec(1,npop);
  !!log_input(Prop_rec);  
  // 3 controls: Lower_limit, Upper_limit, phase_of_estimation
  init_number low_M;       init_number upp_M;       init_int phase_M; 
  init_number low_ln_S0;   init_number upp_ln_S0;   init_int phase_ln_S0;
  init_number low_h;       init_number upp_h;       init_int phase_h;
  init_number low_s2age;   init_number upp_s2age;   init_int phase_s2age;
  
  !!log_input(low_M);      log_input(upp_M);        log_input(phase_M);
  !!log_input(low_ln_S0);  log_input(upp_ln_S0);    log_input(phase_ln_S0);
  !!log_input(low_h);      log_input(upp_h);        log_input(phase_h);
  !!log_input(low_s2age);  log_input(upp_s2age);    log_input(phase_s2age); 

  //Recruitment residuals:
  // 5 controls: start_year, end year,Lower_limit, Upper_limit, phase_of_estimation           
  init_int start_rec; init_int end_rec;   init_number low_rec;   init_number upp_rec;  init_int phase_rec;
  !!log_input(start_rec); log_input(end_rec); log_input(low_rec); log_input(upp_rec); log_input(phase_rec);
  // Selectivity
  init_ivector SelexType(1,TotFleet); //Need to specify what these numbers mean - right now set at 1 for fishery and survey for am.
  !!log_input(SelexType);
  int NselexPar; //Number of selectivity parameters (calculated, not input)
  int NselexVec; //+1 each time SelexType(f) >0: number of selex types that are NOT mirrored
  imatrix SelexPnt(1,TotFleet,0,3); //Something about mirrored selectivity
  imatrix SelexSpec(1,TotFleet,1,3); //Same as SelexType, but just for non-mirrored fleets
 LOCAL_CALCS
   NselexPar = 0;
   NselexVec = 0;
   for (f=1;f<=TotFleet;f++)
    {
     
     // Mirror selectivity
     if (SelexType(f) < 0)
      {
       SelexPnt(f,0) = -1;
       SelexPnt(f,1) = SelexPnt(-1*SelexType(f),1); //What is this about? CRM
       SelexPnt(f,2) = SelexPnt(-1*SelexType(f),2);
       SelexPnt(f,3) = SelexPnt(-1*SelexType(f),3);
      } 
     
     if (SelexType(f) > 0)
      {
       NselexVec += 1;
       SelexSpec(NselexVec,1) = SelexType(f);
       SelexSpec(NselexVec,2) = NselexPar+1;
       SelexPnt(f,0) = SelexType(f);
       SelexPnt(f,1) = SelexType(f);
       SelexPnt(f,2) = NselexVec;
       SelexPnt(f,3) = NselexPar+1;
       if (SelexType(f) == 1) NselexPar += 2; //CRM: Is this logistic selex?
       if (SelexType(f) == 2) NselexPar += 7;
       if (SelexType(f) == 3) NselexPar += 4;
       if (SelexType(f) == 4) NselexPar += 1;
       if (SelexType(f) > 4) { cout << "Selex Type " << SelexType(f) << " out of range" << endl; exit(1); }
      }
      
    }
   CheckFile << endl << "Selectivity overview" << endl; 
   CheckFile << SelexType << endl; 
   CheckFile << SelexPnt << endl; 
   CheckFile << NselexPar << " " << NselexVec << endl;
   
   
 END_CALCS
  
  // Check the number of salectivity parameters
  init_int NselexCheck
  !!log_input(NselexCheck);
  !! if (NselexCheck != NselexPar) { cout << "Number of selectivity parameters does not match " << NselexPar << endl; exit(1); }
 
  // Read in 
  vector SelParLow(1,NselexPar);
  vector SelParHi(1,NselexPar);
  vector SelParInit(1,NselexPar);
  ivector SelParPh(1,NselexPar);
  !! for (i=1;i<=NselexPar;i++)
  !!  *(ad_comm::global_datafile) >> SelParLow(i) >> SelParHi(i) >> SelParInit(i) >> SelParPh(i);
  !!log_input(SelParLow); log_input(SelParHi); log_input(SelParInit); log_input(SelParPh);

  init_vector surv_lambda(1,ncpue)
  !!log_input(surv_lambda);
  init_vector length_lambda(1,nlenf)
  !!log_input(length_lambda);
  init_vector age_lambda(1,nagef)
  !!log_input(age_lambda);
  init_number lambda_rec
  !!log_input(lambda_rec);
  init_ivector effectleng(1,nlenf) // switch to use 0 - number of trips or number of anglers or  1  to use number of samples,
                                   // or another value that will then be used as the effective sample size for all the years
  !!log_input(effectleng);
  init_ivector effectage(1,nagef)  // switch to use 0 - number of trips or number of anglers or  1  to use number of samples,
                                   // or another value that will then be used as the effective sample size for all the years
  !!log_input(effectage);
  init_int Do_rec_Bias;
  !!log_input(Do_rec_Bias);

 //effective sample size for multinomial, maximum index for each row of the ragged array "nlensamp"
 LOCAL_CALCS
      
  //Turn off all the parameters if debbuging is on. 
  if (dummy==1) 
  {
   phase_dummy = 1;
   for (i=1;i<=NselexPar;i++) SelParPh(i) = -1;
   phase_M            =-6; phase_ln_S0        =-1; 
   phase_h            =-5; 
   phase_s2age        =-6; phase_rec          =-3; 
  }
  
  //Turn off the dummy_par is the debugging is off, o/w the hessian will have a 0
  if (dummy==0)  phase_dummy=-1;
    
 END_CALCS 
  
  init_int fim // end of file indicator
  !!mark_input("If you see fim = 999 then reached end of ctl file successfully!");
  !!log_input(fim);
  
  !!cout<<"If you see 999, we got to the end of the control file sucessfully! "<<fim<<endl;
  !!cout<<"Otherwise, Cab will exit abnormally with code 1 so you can correct the CTL file"<<endl;
  !! if (fim != 999) exit(1);

//===================================================================================================
//===================================================================================================
 
PARAMETER_SECTION
  // Dummy parameter for debugging (if dummy== 1, then turn all the parameters off)
  init_number dummy_par(phase_dummy)
  
   //recruitment and initial conditions
  init_bounded_number M(low_M,upp_M,phase_M)
  init_bounded_number ln_S0(low_ln_S0,upp_ln_S0,phase_ln_S0)
  init_bounded_number h_steep(low_h,upp_h,phase_h)
  init_bounded_number sigmasq_rec(low_s2age,upp_s2age,phase_s2age)
  init_bounded_dev_vector rec_dev1(start_rec,end_rec,low_rec,upp_rec,phase_rec) // recruitment residuals population 1
  init_bounded_number_vector SelexPar(1,NselexPar,SelParLow,SelParHi,SelParPh)  // Selectivity parameters
  
  number S0                                                       // Virgin spawning biomass (all stocks) - kind of - SSB0 = S0/2 for a 2 sex model.
  number R0_pop                                                   // Virgin recruitment
  vector log_R0_pop(1,sex)
  vector Spbio(styr,endyr+1)                                      // Spawning biomass
  vector exp_rec(start_rec,end_rec)                               // expected value for recruitment (deterministic)
  vector pred_rec(start_rec,end_rec)                              // predicted value for recruitment (stochastic)
 
  4darray natage(1,npop,styr,endyr+1,1,sex,0,nages)            // Numbers-at-age
  4darray catage(styr,endyr,1,fleet,1,sex,0,nages)             // Catch-at-age by fleet
  4darray catage_tot(1,npop,styr,endyr,1,sex,0,nages)          // Catch-at-age by population
  matrix Hrate(1,fleet,styr,endyr)                                // Harvest by fleet
  
  // Selectivity-related parameters
  matrix Sel_len(1,NselexVec,1,nlength)                           // Selectivity (by length and fleet)
  3darray Sel_age(1,NselexVec,1,sex,0,nages)                   // Selectivity (by age and fleet)
  3darray Sel_wt_age(1,NselexVec,1,sex,0,nages)                // Weight times selectivity
  
  // weight at the middle of the year
  3darray wt_age_middle(1,TotFleet,1,sex,0,nages)      

  // Parameters related to indices
  vector log_cpue(1,ncpue);
  vector sigma_cpue(1,ncpue);
  matrix exp_cpue(1,ncpue,1,nyrs_cpue); 
  matrix sd_cpue(1,ncpue,1,nyrs_cpue);
  
  // Expected age- and length-frequency
  3darray exp_lenf(1,nlenf,1,nyrs_lenf,1,nlength)
  matrix neffleng(1,nlenf,1,nyrs_lenf)           // Estimated effective sample size for multinomial
  3darray exp_agef(1,nagef,1,nyrs_agef,0,nages)
  matrix neffage(1,nagef,1,nyrs_agef)            // Estimated effective sample size for multinomial
    
  vector surv_like(1,ncpue)                      // likelihood of the indices
  vector length_like(1,nlenf)                    // likelihood of the length-frequency data
  vector age_like(1,nagef)                       // likelihood of the age-frequency data
  number prior_rec
  number CrashPen;
  
  objective_function_value obj_fun 
  !!cout<<"end of parameter section"<<endl;
 
  //likelihood profile numbers
  likeprof_number S0_lprof
  sdreport_number Depl;
  sdreport_matrix LogTheNatAge(1,sex,0,nages);  //Carey added this
  //sdreport_vector LogTheNatAgeMale(0,nages);  //Carey added this
  
// *=========================================================================*

PRELIMINARY_CALCS_SECTION
  int Iseries;
  float total,length;
  
  //This will guarantee that the vectors are set = 0  at the beginning of the run
  age_length_pop.initialize();
  wt_age_middle.initialize();

  //GROWTH 
  for (g=1;g<=sex;g++)
   age_length_pop(g) = SizeTrans(Lmin(g),Lmax(g),K(g),CVLmin(g),CVLmax(g),1);
 
  // wt_input(1) is beginning of the year weight, wt_input(2) is middle of the year weight, from age 0 to nages
  for (f=1;f<=TotFleet;f++)
   for (g=1;g<=sex;g++) 
    for (a=0;a<=nages;a++)
     wt_age_middle(f,g,a)=wt_input(2,g,a);

  // Weight-at-length
  length = len_start+ len_step/2.0;
  for (l=1;l<=nlength;l++)
   {
    wt_at_length(l) = WL_intercept * pow(length,WL_slope);
    length += len_step;
   }

  // Specify the effective sample sizes (length)
  for (Iseries=1;Iseries<=nlenf;Iseries++)
   for (i=1;i<=nyrs_lenf(Iseries);i++)
    {
     if (effectleng(Iseries) == 0)
      nlensamp(Iseries,i) = TheLengths(Iseries,i,-1);
     else if (effectleng(Iseries) == -1)
      nlensamp(Iseries,i) = TheLengths(Iseries,i,0);
     else
      nlensamp(Iseries,i) = effectleng(Iseries); 
    }
    
  // Specify the effective sample sizes (age)
  for (Iseries=1;Iseries<=nagef;Iseries++)
   for (i=1;i<=nyrs_agef(Iseries);i++)
    {
     if (effectage(Iseries) == 0)
      nagesamp(Iseries,i) = TheAges(Iseries,i,-2);
     else if (effectage(Iseries) == -1)
      nagesamp(Iseries,i) = TheAges(Iseries,i,-1);
     else
      nagesamp(Iseries,i) = effectage(Iseries); 
    }
    
  // Calculate the OFFSET and renormalize the data
  offleng.initialize();
  for (Iseries=1;Iseries<=nlenf;Iseries++)
   for (i=1;i<=nyrs_lenf(Iseries);i++)
    {
     total = 0;
     for (l=1;l<=nlength;l++) total += TheLengths(Iseries,i,l);
     for (l=1;l<=nlength;l++)
      {
       TheLengths(Iseries,i,l) = TheLengths(Iseries,i,l)/total;
       offleng(Iseries) -= nlensamp(Iseries,i)*TheLengths(Iseries,i,l)*log(0.0001+TheLengths(Iseries,i,l));
      } 
    }

  // Calculate the OFFSET and renormalize the data
  offage.initialize();
  for (Iseries=1;Iseries<=nagef;Iseries++)
   for (i=1;i<=nyrs_agef(Iseries);i++)
    {
     total = 0;
     for (l=0;l<=nages;l++) total += TheAges(Iseries,i,l);
     for (l=0;l<=nages;l++)
      {
       TheAges(Iseries,i,l) = TheAges(Iseries,i,l)/total;
       offage(Iseries) -= nagesamp(Iseries,i)*TheAges(Iseries,i,l)*log(0.0001+TheAges(Iseries,i,l));
      } 
    }

  // Initial values for bias-correction-2
  sigmasq_rec = (low_s2age+upp_s2age)/2.0;
  //cout<<"CRM: what is rec_dev1 before and after Do_rec_Bias==2? Before: "<<endl<<rec_dev1<<endl;
  if (Do_rec_Bias == 2) rec_dev1 = -sigmasq_rec/2.0;
  //cout<<"CRM: After: "<<rec_dev1<<endl;

  if (dummy==1) 
     {cout<<endl<<endl<<endl<<"Debugging is on, no parameters are being estimated"<<endl;
      cout<<"model being projected using pin file values"<<endl<<endl<<endl;}
  else
     {cout<<endl<<endl<<endl<<"Estimating...please wait..."<<endl; 
      if (do_var==0) cout<<endl<<"Please note: extra variability around observations is not been estimated"<<endl;} 
      
  cout << "Completed Preliminary Calculations" << endl;
     
// *=========================================================================*
  
PROCEDURE_SECTION

  // Reset the crash penalty
  CrashPen.initialize();
 
  // selectivity does not change over time, it can be compute only once in each iteration
  get_selectivity(); 
  //cout<<"end of get selectivity"<<endl;

  get_initial_conditions(); 
  //cout<<"end of get initial conditions"<<endl;
  
  get_numbers_at_age();  
  //cout<<"end of get numbers at age"<<endl;

  get_recruitment_prior();   
  //cout<<"end of get recruitment prior"<<endl;

  // compute the expected values for the indices and the length frequency
  get_predictions(); 
  //cout << "end of Get_predictions" << endl;
  
  // compute likelihood functions
  evaluate_the_objective_function(); 
  //cout << "end of calculate obj fn" << endl;
  
  Depl = Spbio(endyr+1)/Spbio(styr)*100;
  LogTheNatAge = log(natage(1,endyr+1));  //Carey added this
  //LogTheNatAgeMale = log(natage(1,endyr+1,2));  //Carey added this

  if (mceval_phase())
   {
    cout << obj_fun << " " << Depl << " " << ln_S0 << " ";
    cout << SelexPar << endl;
    cout << rec_dev1 << " ";
    cout << Spbio << endl;
   }


//=========================================================================================
 
FUNCTION get_recruitment_prior
  dvariable chi,tmp;

  //The recruitment prior is assumed to be a lognormal pdf with expected
  // value equal to the deterministic stock-recruitment curve  
  chi = 0.;
  //if (phase_rec > 0) 
   for (i=start_rec;i<=end_rec;i++) //CRM changed this from i = start_rec to i = start_rec+1 and then back
    {
     if (Do_rec_Bias==2)
      tmp = log( pred_rec(i) / exp_rec(i) ) + sigmasq_rec/2;
     if (Do_rec_Bias==1)
      tmp = log( pred_rec(i) / exp_rec(i) );
     if (Do_rec_Bias==0)
      tmp = log( pred_rec(i) / exp_rec(i) );
     chi += square(tmp)/(2*sigmasq_rec) + log(sqrt(sigmasq_rec));
    } 
  prior_rec = chi;
  
  // Recruitment variability: EARLY PHASES ONLY
  if (!last_phase() ) 
    prior_rec += 1. * norm2(rec_dev1);
     
  // Adjust to weight
  prior_rec *= lambda_rec ;
   
//====================================================================================

FUNCTION get_numbers_at_age
  dvariable vul_bio;                                  // Vulnerable biomass
  dvariable harvest_rate;                             // Harvest rate
  dvariable Spaw_bio;                                 // Spawning biomass
  dvariable Recruits;                                 // Age0 Recruits 
  int SelPnt;                                         // Pointer
  
  // Reset variables
  catage.initialize();  catage_tot.initialize(); 

  //loop over time
  for (t=styr;t<=endyr;t++)
   {

    //loop over fleets
    for (f=1;f<=fleet;f++) 
     if (catch_bio(t,f) > 0)
      {         
       pop = fleetpop(f);
       SelPnt = SelexPnt(f,2);
       
       // vul_bio for each fleet, need this to calculate the harvest rate
       vul_bio = 0.0; 
       for (g=1;g<=sex;g++)
        vul_bio += (natage(pop,t,g)*mfexp(-M/2))*Sel_wt_age(SelPnt,g);
       
       // Compute the harvest rate and store it
       harvest_rate = (1.- posfun(1.-catch_bio(t,f)/vul_bio,.0001,CrashPen));
       Hrate(f,t)   = harvest_rate;
          
       // Compute the predicted catch at age 
       for (g=1;g<=sex;g++)
        {
         catage(t,f,g) = harvest_rate* elem_prod( (natage(pop,t,g)*mfexp(-M/2) ) ,Sel_age(SelPnt,g));
         // for (a=1;a<=nages;a++)
          // if (catage(t,f,g,a)<=0.0) {catage(t,f,g,a) = 0.0;}       // avoid negative catches
         catage_tot(pop,t,g) += catage(t,f,g);                     //catch at age for all fleets, this is summing the two fleets
        }
       cout << "Find glitch in catch at age " << endl;
       cout << "Year = " << t << ", fleet = "<< f <<", Hrate(f,t) = "<< Hrate(f,t) << ", vul_bio = " << vul_bio << endl;

      }    
     else
      Hrate(f,t) = 0.0;
    // Update the dynamics------------------------------------------------------------------------
    for (pop=1;pop<=npop;pop++) 
     for (g=1;g<=sex;g++)
      {               
       // Recruitment (by sex)
       natage(pop,t+1,g,0) = 0;
                   
       // the rest of the ages
       for (a=1;a<nages;a++)
        natage(pop,t+1,g,a) = natage(pop,t,g,a-1)*mfexp(-M)-catage_tot(pop,t,g,a-1)*mfexp(-M/2);
         
       // plus group
       natage(pop,t+1,g,nages) = natage(pop,t,g,nages-1)*mfexp(-M) - catage_tot(pop,t,g,nages-1)*mfexp(-M/2);
       natage(pop,t+1,g,nages) += natage(pop,t,g,nages)*mfexp(-M) - catage_tot(pop,t,g,nages)*mfexp(-M/2);
                     
       // now make sure all numbers at age are above 0
       for (a=0;a<=nages;a++)
        if (natage(pop,t+1,g,a)<=0.0001) {natage(pop,t+1,g,a) = 0.0001;} 
      }
       
    // Compute the spawning biomass (males and females)-----------------------------------------
    Spaw_bio = 0.0; 
    for (pop=1;pop<=npop;pop++) Spaw_bio += Fec*natage(pop,t+1,1); 
    if (Spaw_bio < 0.0001) {Spaw_bio= 0.0001;}
    Spbio(t+1) = Spaw_bio;                                  //store it for the report
       
    // Compute recruitment------------------------------------------------------------------------ 
    Recruits =  (4*h_steep*R0_pop*Spaw_bio) / ((S0/sex)*(1-h_steep)+(5*h_steep-1)*Spaw_bio) ; //deterministic 
    // add stochastics bits and store the quantities we need for the recruitment prior 
    cout << "h_steep = " << h_steep << ", t = " << t << ", Spaw_bio = " << Spaw_bio << ", Fec = " << Fec << ", R0_pop = " << R0_pop << ", S0/sex = " << S0/sex << endl;
    cout << "Spbio(t+1) = " << Spbio(t+1) << endl;
    //CareyFix: changed below from t+1 to t throughout if statement to match with get_recruitment_prior section
    if (t >= start_rec & t <=end_rec) //CRM: Recruits is defined by the population dynamics at the end of the first year and applies as year t+1; likewise, Recruits in the last year (t+1 = end_rec+1) does not correspond to any data, so it is not included in the likelihood component for recruitment later on.
     {
      exp_rec(t)  = Recruits;                         //store deterministic
      Recruits *= mfexp(rec_dev1(t));
      if (Do_rec_Bias == 1) Recruits *= mfexp(-sigmasq_rec*0.5);
      pred_rec(t) = Recruits;
     }  
    //End CareyFix   
    
    // Recruitment (by sex)
    for (pop=1;pop<=npop;pop++)
     for (g=1;g<=sex;g++) natage(pop,t+1,g,0) =  Prop_rec(pop)*Recruits/sex;      


   } //close time loop
     
//==================================================================================

FUNCTION get_selectivity
  dvariable total;
  int Ipnt,SelType;
  dvariable par1,par2,par3,par4,par5,par6,par7;

  //Reset variables
  Sel_len.initialize(); Sel_age.initialize();Sel_wt_age.initialize();
  
  //----------------------------Selectivity at length------------------------------
  if (age_specific_sel == 0)
   {
    for (f=1;f<=NselexVec;f++)
     {
      Ipnt =  SelexSpec(f,2); 
      SelType =  SelexSpec(f,1);

      // Logistic
      if (SelType == 1)
       {
        par1 = SelexPar(Ipnt); par2 = SelexPar(Ipnt+1);
        Sel_len(f) = 1/(1+mfexp(-log(19)*(size_vector-par1)/par2));
       } 

      // Double logistic
      if (SelType == 2)
       {
        par1 = SelexPar(Ipnt); par2 = SelexPar(Ipnt+1); par3 = SelexPar(Ipnt+2); par4 = SelexPar(Ipnt+3);
        par5 = SelexPar(Ipnt+4); par6 = SelexPar(Ipnt+5); par7 = SelexPar(Ipnt+6);
        Sel_len(f) = DoubLogistic(par1,par2,par3,par4,par5,par6,par7);
       } 

      // Double normal        
      if (SelType == 3)
       {
        par1 = SelexPar(Ipnt); par2 = SelexPar(Ipnt+1); par3 = SelexPar(Ipnt+2); par4 = SelexPar(Ipnt+3);
        Sel_len(f) = DoubNormal(par1,par2,par3,par4);
       }

      // Knife-edged
      if (SelType == 4)
       {
        par1 = SelexPar(Ipnt); 
        for (l=1;l<=nlength;l++)
         if (size_vector(l) < par1) Sel_len(f,l) = 0; else Sel_len(f,l) = 1;
       } 

     }
   
    //--------------Selectivity at age, is the selectivity at length times the Age-Length Key
    //  sel_wt_age is the selectivity at age times the weight at age at the middle of the year 
    if (wt_sel_bias == 0)
     {
      for (f=1;f<=NselexVec;f++)
       for(g=1;g<=sex;g++)
        {
         for (a=0;a<=nages;a++)
          Sel_age(f,g,a)= Sel_len(f)*age_length_pop(g,a);
          Sel_wt_age(f,g) = elem_prod(Sel_age(f,g), wt_age_middle(f,g));    
        }
     }   
    if (wt_sel_bias == 1)
     {
      for (f=1;f<=NselexVec;f++)
       for(g=1;g<=sex;g++)
        for (a=0;a<=nages;a++)
         {
          for (l=1;l<=nlength;l++)
           {
            Sel_age(f,g,a) += Sel_len(f,l)*age_length_pop(g,a,l);
            Sel_wt_age(f,g,a) += Sel_len(f,l)*age_length_pop(g,a,l)*wt_at_length(l);
           } 
         } 
     } 
   }
   
 if (age_specific_sel == 1)
  {
   cout <<" starting age_specific_sel = 1 "<< endl;
   for (g=1;g<=sex;g++)
    {
     for (f=1;f<=NselexVec;f++)
      {
       Ipnt =  SelexSpec(f,2); 
       SelType =  SelexSpec(f,1);
       //cout << "Ipnt = " << Ipnt <<endl;
       //cout << "SelType (1 is logistic) = " << SelType << endl;
       // Logistic
       if (SelType == 1)
        {
         par1 = SelexPar(Ipnt); par2 = SelexPar(Ipnt+1);
         Sel_age(f,g) = 1/(1+mfexp(-log(19)*(age_vector-par1)/par2)); 
        } 
        
        //cout << "Doing logistic selex " << endl;
        //cout << "a50 = " << par1 << endl;
        //cout << "slope = " << par2 << endl;
       // Double logistic
       if (SelType == 2)
        {
         cout << "Selex Type 2 is not yet available for age-based selex" << endl;
         exit(1);
        } 

       // Double normal        
       if (SelType == 3)
        {
         par1 = SelexPar(Ipnt); par2 = SelexPar(Ipnt+1); par3 = SelexPar(Ipnt+2); par4 = SelexPar(Ipnt+3);
         Sel_age(f,g) = DoubNormalAge(par1,par2,par3,par4);
        }
     
      // Knife-edged
      if (SelType == 4)
       {
        par1 = SelexPar(Ipnt);
        for (a=0;a<=nages;a++)
         {
         if (age_vector(a) < par1) Sel_age(f,g,a) = 0; else Sel_age(f,g,a) = 1;
         }
       } 

     }
      
     //  sel_wt_age 
     for (f=1;f<=NselexVec;f++)
      Sel_wt_age(f,g) = elem_prod(Sel_age(f,g), wt_age_middle(f,g));    
    }  
  }  
 
//=====================================================================

FUNCTION get_initial_conditions
  dvariable sum_fec=0.0;
  
  // reset
  natage.initialize();
  
  //Virgin_Recruitment (by population)
  S0 = mfexp(ln_S0);
  S0_lprof = S0; //for likelihood profile
   
  //Calculate R0 from S0 and fecundity at age
  sum_fec = Fec(nages)*mfexp(-M*double(nages))/(1-mfexp(-M));
  for (a=1;a<nages;a++)
   sum_fec  += Fec(a)*mfexp(-M*double(age_vector(a)));
  R0_pop = S0/sum_fec; 
  Spbio(styr) = S0/sex;
  cout << "S0 = " << S0 << ", sum_fec = " << sum_fec <<", R0_pop = "<< R0_pop;
  //Allocate half of the recruitment for each sex
  for (g=1;g<=sex;g++)
   if (R0_pop > 0) 
    log_R0_pop(g)= log (R0_pop/sex); 
   else
    log_R0_pop(g) = 0.0;  
   
  //Initial age structure
  for (pop=1;pop<=npop;pop++)
   for (g=1;g<=sex;g++)
    {
     natage(pop,styr,g,0) = Prop_rec(pop)*mfexp(log_R0_pop(g));
     for (j=1;j<=nages-1;j++) 
      natage(pop,styr,g,j) = natage(pop,styr,g,j-1) * mfexp(-M); 
     natage(pop,styr,g,nages) = natage(pop,styr,g,nages-1) * mfexp(-M)/(1-mfexp(-M));
    }

// ==============================================================================

FUNCTION get_predictions 
  int iyr,Ip,Iseries,Yr,Sex,SelPnt;
  dvariable log_tot,total,tot1,tot2,term;

  // Catch-rate data
  for (Iseries=1;Iseries<=ncpue;Iseries++)
   {
    log_cpue(Iseries) = 0;
    log_tot = 0;
    f = indx_fleet(Iseries);
    pop = fleetpop(f);
    SelPnt = SelexPnt(f,2);
    for (i=1;i<=nyrs_cpue(Iseries);i++)
     {
      Yr = int(TheAbund(Iseries,i,2));
      exp_cpue(Iseries,i) = 0;
      
      // CPUE [or survey] index {in biomass}
      if (indx_type(Iseries) == 1)
       for (g=1;g<=sex;g++)
        exp_cpue(Iseries,i) += Sel_wt_age(SelPnt,g)*(mfexp(-M/2)*natage(pop,Yr,g)-catage_tot(pop,Yr,g) );
        
      // CPUE [or survey] index {in numbers}
      else if (indx_type(Iseries) == 2)
       for (g=1;g<=sex;g++)
        exp_cpue(Iseries,i) += Sel_age(SelPnt,g)*(mfexp(-M/2)*natage(pop,Yr,g)-catage_tot(pop,Yr,g) );

      // SSB index  
      else if (indx_type(Iseries) == -1)
       for (pop=1;pop<=npop;pop++)
        for (g=1;g<=sex;g++)
         exp_cpue(Iseries,i) += Fec*natage(pop,Yr-1,g);
        
      // Recruitment index  
      else if (indx_type(Iseries) == -2)  
       for (pop=1;pop<=npop;pop++)
        for (g=1;g<=sex;g++)
         exp_cpue(Iseries,i) += natage(pop,Yr,g,0) + 0.5*natage(pop,Yr,g,1);
      
      // Ugh problem
      else
       {cout << "Undefined index type " << indx_type(Iseries) << endl; exit(1); }
        
      if (exp_cpue(Iseries,i) <= 0) exp_cpue(Iseries,i) = 0.001;
      log_cpue(Iseries) += log(TheAbund(Iseries,i,3)/exp_cpue(Iseries,i))/square(TheAbund(Iseries,i,4)); 
      log_tot += 1.0/square(TheAbund(Iseries,i,4));
     } 
    log_cpue(Iseries) /= log_tot;  
    
    // ML estimate for sigma
    if (do_var == 1)
     {
      sigma_cpue(Iseries) = 0;
      for (i=1;i<=nyrs_cpue(Iseries);i++)
       // GML 17Apr13: Bug: Need to close brackets after first log.
       // AEP: Yup, the brackets were wrong
       //sigma_cpue(Iseries) += square(log(TheAbund(Iseries,i,3)-log(exp_cpue(Iseries,i))-log_cpue(Iseries))/TheAbund(Iseries,i,4));
       sigma_cpue(Iseries) += square((log(TheAbund(Iseries,i,3))-log(exp_cpue(Iseries,i))-log_cpue(Iseries))/TheAbund(Iseries,i,4));
       sigma_cpue(Iseries) = sqrt(sigma_cpue(Iseries)/nyrs_cpue(Iseries));
     }
    for (i=1;i<=nyrs_cpue(Iseries);i++)
     if (do_var == 1)
      sd_cpue(Iseries,i) = sigma_cpue(Iseries)*TheAbund(Iseries,i,4);
     else
      sd_cpue(Iseries,i) = TheAbund(Iseries,i,4);
   }
  //if (Diags >= 2) cout << "Predictions for the indices done" << endl; 
  // Length-frequency data
  neffleng.initialize();
  exp_lenf.initialize();
  for (Iseries=1;Iseries<=nlenf;Iseries++)
   {
    f = lenf_fleet(Iseries);
    pop = fleetpop(f);
    SelPnt = SelexPnt(f,2);

    for (i=1;i<=nyrs_lenf(Iseries);i++)
     {
      iyr = int(TheLengths(Iseries,i,-3));
      Sex = int(TheLengths(Iseries,i,-2));
     
      for (g=1;g<=sex;g++) 
       if (Sex == 0 || Sex ==g)
        for (l=1;l<=nlength;l++)
         {
          for (j=0;j<=nages;j++)
           {
            if (age_specific_sel == 0) term = natage(pop,iyr,g,j)*age_length_pop(g,j,l)*Sel_len(SelPnt,l);
            if (age_specific_sel == 1) term = natage(pop,iyr,g,j)*age_length_pop(g,j,l)*Sel_age(SelPnt,g,j);
            exp_lenf(Iseries,i,l) += term;
           }
         }    
      for (l=1;l<=nlength;l++)
       if (exp_lenf(Iseries,i,l) <= 0) exp_lenf(Iseries,i,l)= 0.00001;
      total = sum(exp_lenf(Iseries,i)); 
      tot1 = 0; tot2 = 0;
      for (l=1;l<=nlength;l++)
       {
        exp_lenf(Iseries,i,l) /= total;
        tot1 += exp_lenf(Iseries,i,l)*(1.0-exp_lenf(Iseries,i,l));
        tot2 += square(TheLengths(Iseries,i,l)-exp_lenf(Iseries,i,l));
        }
      neffleng(Iseries,i) = tot1/tot2;
     }
   }  
  //if (Diags >= 2) cout << "Predictions for the lengths done" << endl; 
  
  // Age-frequency data
  neffage.initialize();
  exp_agef.initialize();
  for (Iseries=1;Iseries<=nagef;Iseries++)
   {
    f = agef_fleet(Iseries);
    pop = fleetpop(f);
    SelPnt = SelexPnt(f,2);

    for (i=1;i<=nyrs_agef(Iseries);i++)
     {
      iyr = int(TheAges(Iseries,i,-4));
      Sex = int(TheAges(Iseries,i,-3));
      for (g=1;g<=sex;g++)
       if (Sex == 0 || g == Sex)
        for (j=0;j<=nages;j++)
         {
          if (age_specific_sel == 0)
           for (l=1;l<=nlength;l++)
            exp_agef(Iseries,i,j) += natage(pop,iyr,g,j)*age_length_pop(g,j,l)*Sel_len(SelPnt,l);
          if (age_specific_sel == 1)  
           exp_agef(Iseries,i,j) += natage(pop,iyr,g,j)*Sel_age(SelPnt,g,j);
         }
      for (j=0;j<=nages;j++)   
       if (exp_agef(Iseries,i,j) <= 0) exp_agef(Iseries,i,j) = 0.00001;
      total = sum(exp_agef(Iseries,i)); 
      tot1 = 0; tot2 = 0;
      for (j=0;j<=nages;j++)
       {
        exp_agef(Iseries,i,j) /= total;
        tot1 += exp_agef(Iseries,i,j)*(1.0-exp_agef(Iseries,i,j));
        tot2 += square(TheAges(Iseries,i,j)-exp_agef(Iseries,i,j));
       }
      neffage(Iseries,i) = tot1/tot2;
     }
   }  
  //if (Diags >= 2) cout << "Predictions for the ages done" << endl; 
  
//=========================================================================================

FUNCTION evaluate_the_objective_function
  int Iseries;

  // Fit to indices (lognormal) 
  surv_like.initialize();   
  for (Iseries=1;Iseries<=ncpue;Iseries++)
   if (do_var == 1)
    {
     surv_like(Iseries) = nyrs_cpue(Iseries)*log(sigma_cpue(Iseries)) + nyrs_cpue(Iseries)/2.0;
    }
   else
    {
     surv_like(Iseries) = 0;
     for (i=1;i<=nyrs_cpue(Iseries);i++)
      surv_like(Iseries) += square((log(TheAbund(Iseries,i,3))-log(exp_cpue(Iseries,i))-log_cpue(Iseries))/TheAbund(Iseries,i,4))/2.0;
    } 
  
  // Fit to length data (multinomial)
  length_like.initialize();
  for (Iseries=1;Iseries<=nlenf;Iseries++)
   {
    for (i=1;i<=nyrs_lenf(Iseries);i++)
      length_like(Iseries) -= nlensamp(Iseries,i)*TheLengths(Iseries,i)(1,nlength)*log(exp_lenf(Iseries,i)+0.0001);
     // for (l=1;l<=nlength;l++) length_like(Iseries) -= nlensamp(Iseries,i)*TheLengths(Iseries,i,l)*log(exp_lenf(Iseries,i,l)+0.0001);
    length_like(Iseries) -= offleng(Iseries);
   }
  
  // Fit to age data (multinomial)
  age_like.initialize();
  for (Iseries=1;Iseries<=nagef;Iseries++)
   {
    for (i=1;i<=nyrs_agef(Iseries);i++)
      age_like(Iseries) -= nagesamp(Iseries,i)*TheAges(Iseries,i)(0,nages)*log(exp_agef(Iseries,i)+0.0001);
     // for (l=0;l<=nages;l++) age_like(Iseries) -= nagesamp(Iseries,i)*TheAges(Iseries,i,l)*log(exp_agef(Iseries,i,l)+0.0001);
    age_like(Iseries) -= offage(Iseries);
   }

  // if (dummy == 1) //debugging mode, turn off all parameters obj_fun= dummy_par*dummy_par;
  // else
  {
   obj_fun += surv_lambda*surv_like;       
   obj_fun += length_lambda*length_like;
   obj_fun += age_lambda*age_like;          
   obj_fun += prior_rec + CrashPen;
   // if (!mceval_phase()) cout << obj_fun << endl;
  }

// *=========================================================================*

REPORT_SECTION
  int Iseries;
  dvariable BioOut,Recruit;
  

  // Likelihood summary
  report << "obj_fun "<<endl<<obj_fun<<endl; //CRM added
  report << "Likelihood components " <<endl;
  report <<" indices "<<endl<<surv_like<<endl << elem_prod(surv_like,surv_lambda)<<endl;
  report <<" length frequency "<< endl<<length_like <<endl << elem_prod(length_like,length_lambda)<<endl;
  report <<" age frequency "<< endl << age_like<<endl<< elem_prod(age_like,age_lambda)<<endl;
  report <<" penalties: prior_rec, CrashPen "<<prior_rec<< " " << CrashPen <<endl;
  report << endl;
  //cout<<"Report: before summary stats"<<endl;
  //cout<<"Spbio(endyr+1) = "<<Spbio(endyr+1)<<endl;
  //cout<<"Spbio(styr) = "<<Spbio(styr)<<endl;
  //cout<<"Spbio = "<<Spbio<<endl;
  // Summary stats
  //CRM added the following descriptive line:
  report <<"Year, SpBio(Year), SpBio(Year)/SpBio(styr), BioOut, Recruit, Hrate "<<endl;
  for (i=styr;i<=endyr+1;i++)
   {
    report << i << " ";
    report << Spbio(i) << " " << Spbio(i)/Spbio(styr) << " ";
    //cout<<"npop = "<<npop<<endl;
    //cout<<"nages = "<<nages<<endl;
    //cout<<"sex = "<<sex<<endl;
    for (pop=1;pop<=npop;pop++)
     {
      BioOut = 0;
      for (j=0;j<=nages;j++)
       for (g=1;g<=sex;g++)
       {//cout<<"pop, i, g, j = "<<pop<<", "<<i<<", "<<g<<", "<<j<<endl;
        //cout<<"natage(pop,i,g,j) = "<<natage(pop,i,g,j)<<endl;
        BioOut += natage(pop,i,g,j)*wt_age_middle(1,g,j);
       }
      Recruit = 0;
      for (g = 1;g<=sex;g++)
      Recruit += natage(pop,i,g,0); //Carey's fix for next line problem
      //Recruit = natage(pop,i,1,0) + natage(pop,i,g,0) //This only works if you have 2 sexes in the model!
      report << BioOut << " " << Recruit << " "; 
      //cout<<"i, pop, j, g = "<<i<<", "<<pop<<", "<<j<<", "<<g<<endl;
      } 
    //cout<<"fleet = "<<fleet<<endl;
    for (f=1;f<=fleet;f++)
     if (i != endyr+1) report << Hrate(f,i) << " "; else report << "-1" << " ";
    report << endl;
   } 
   //cout<<"Report: After looping over years"<<endl;
  report << endl;
  report << setw(5) << h_steep << " " 
         << setw(5) << sigmasq_rec << " ";
  report << setw(12) << S0/2 << endl;
  report << " " << endl;
  
  // Abundance index information
  //cout<<"Report: just before Fits to indices"<<endl;
 report << "Fits to indices" << endl; 
 for (Iseries=1;Iseries<=ncpue;Iseries++)
  for (t=1;t<=nyrs_cpue(Iseries);t++)
   report << "I"<< Iseries << " " << indx_fleet(Iseries) << " " << indx_type(Iseries) << " " << TheAbund(Iseries,t,2) << " " << TheAbund(Iseries,t,3) << " " << mfexp(log_cpue(Iseries))*exp_cpue(Iseries,t) << " " << sd_cpue(Iseries,t) << endl;
  if (do_var==1)
   for (Iseries=1;Iseries<=ncpue;Iseries++)
    report << "Cpue-sigma " << Iseries << " " << sigma_cpue(Iseries) << endl;
  for (Iseries=1;Iseries<=ncpue;Iseries++)
   report << "log_q_cpue " << Iseries << " " << log_cpue(Iseries) << endl;

  //cout<<"Report: just before Length-frequency information"<<endl;
  // Length-frequency information
 report << endl << "length-frequency information" << endl;
 for (Iseries=1;Iseries<=nlenf;Iseries++)
  for (t=1;t<=nyrs_lenf(Iseries);t++)
   for (l=1;l<=nlength;l++)
    {
     report<<"lf "<<Iseries<<" "<<lenf_fleet(Iseries) << " " << lenf_type(Iseries) << " " << TheLengths(Iseries,t,-3)<<" "<<nlensamp(Iseries,t)<<" "<<neffleng(Iseries,t)<<" ";
     report<<size_vector(l)<<" "<<TheLengths(Iseries,t,l)<<" "<<exp_lenf(Iseries,t,l)<<endl;
    }
   //cout<<"Report: just before Age-frequency information"<<endl;
  // Age-frequency information
 report << endl << "Age-frequency information" << endl;
 for (Iseries=1;Iseries<=nagef;Iseries++)
  for (t=1;t<=nyrs_agef(Iseries);t++)
   for (l=0;l<=nages;l++)
    {
     report<<"Ag "<<Iseries<<" "<<agef_fleet(Iseries) << " " << agef_type(Iseries) <<" " << TheAges(Iseries,t,-4)<<" "<<nagesamp(Iseries,t)<<" "<<neffage(Iseries,t)<<" ";
     report<<age_vector(l)<<" "<<TheAges(Iseries,t,l)<<" "<<exp_agef(Iseries,t,l)<<endl;
    }
  
  //Selectivity at length
  report << endl << "Selectivity-at-age" << endl;
  report<<"fleet"<<" "<<"length"<<" "<<"selL"<<endl;
  for(f=1;f<=NselexVec;f++)
   for(l=1;l<=nlength;l++)
    report<<f<<" "<<size_vector(l)<<" "<<Sel_len(f,l)<<endl;
  
  //Selectivity at age
  report << endl << "Selectivity-at-age" << endl;
  report<<"fleet"<<" "<<"sex"<<" "<<"age"<<" "<<"selA"<<endl;
  for(f=1;f<=NselexVec;f++)
   for(g=1;g<=sex;g++)
    for(a=0;a<=nages;a++)
     report<<f<<" "<<g<<" "<<a<<" "<<Sel_age(f,g,a)<<endl;

   //Number and catch at age for quick look
  report << endl << "Estimated numbers of fish " << endl;
  for (pop=1;pop<=npop;pop++)
   for (g=1;g<=sex;g++)
    {report << "Population "<<pop<<" sex "<<g<<endl;
     report << "Age "<<"0 "<<age_vector <<endl;
     for (i=styr;i<=endyr;i++) 
      report <<i<<" "<< natage(pop,i,g) << endl;}

  report <<endl<< "Estimated catch at age " << endl;
  for (g=1;g<=sex;g++)
   for (f=1;f<=fleet;f++)
    {
     pop = fleetpop(f);
     report << "Population "<<pop<<" fleet "<<f<<" sex "<<g<<endl;
     report << "Age "<<"0 "<<age_vector <<endl;
     for (i=styr;i<=endyr;i++) 
      report <<i<<" "<<catage(i,f,g)  << endl;
    }

  report <<endl<< "Estimated total catch at age " << endl;
  for (pop=1;pop<=npop;pop++)
   for (g=1;g<=sex;g++)
    {report << "Population "<<pop <<" sex "<<g<<endl;
     report << "Age "<<"0 "<<age_vector <<endl;
     for (i=styr;i<=endyr;i++) 
      report <<i<<" "<<catage_tot(pop,i,g)  << endl;}
  ////cout << "Here" << endl;    

  //cout << "Done report Section" << endl;


// *=====================================================================================*
// *=====================================================================================*

RUNTIME_SECTION
  maximum_function_evaluations 1000 1000 1000 2000;
  convergence_criteria 0.01,0.01,1e-7;

// *======================================================================================*
// *======================================================================================*

TOP_OF_MAIN_SECTION
  arrmblsize = 5000000;
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(560000);
  gradient_structure::set_CMPDIF_BUFFER_SIZE(1500000);
  gradient_structure::set_MAX_NVAR_OFFSET(500);
  gradient_structure::set_NUM_DEPENDENT_VARIABLES(500);
  time(&start); //this is to see how long it takes to run
  cout<<endl<<"Start time : "<<ctime(&start)<<endl; 
  CheckFile.open("Check.Out");


// *======================================================================================*

GLOBALS_SECTION
  #include <admodel.h>
  #include <time.h>
  time_t start,finish;
  long hour,minute,second;
  double elapsed_time;
  ofstream CheckFile;

   #undef log_input
   #define log_input(object) write_input_log << "# " #object "\n" << object << endl;
   #undef log_param

   #undef mark_input
   #define mark_input(object) write_input_log <<"\n" "# " #object "\n\n\n" <<endl;
   
  ofstream write_input_log("input.log");
  
//---------------------------------------------------------------------------------------

FUNCTION dmatrix SizeTrans(const double& Lbeg,const double& Lmax,const double& K,const double& CVLmin,const double& CVLmax,const int& m)
  { 
   RETURN_ARRAYS_INCREMENT(); //Need this statement because the function
   // m is a switch, if m==0, the function will calculate the length transition for the beginning of the year,
   //		     if m==1, the function will calculate it for the middle of the year; 
   dmatrix Size_Trans(0,nages,1,nlength);
   dvector Average_Size(0,nages);
   dvector Sd(0,nages);
   double age; 
   for (i=0; i<=nages;i++) 
    {
     //first calculate average values...
     if(m==0) age=double(i);//beginning of the year
     else age=double(i)+0.5;//middle of the year
     Average_Size(i) = Lmax + (Lbeg - Lmax) * mfexp (- K * (age-1));
     Sd(i)= (CVLmin+(age-1)*(CVLmax-CVLmin)/(nages-1))*Average_Size(i);

     //...then calculate the distribution arround those values
     // first bin, note: need to standarize before using cumd_norm;
     Size_Trans(i,1)=((size_vector(1)+(size_vector(2)-size_vector(1))/2)-Average_Size(i))/Sd(i);
     Size_Trans(i,1)=cumd_norm(Size_Trans(i,1));

     //other bins but the last;
     for (j=2;j<=nlength-1;j++)
      {
       Size_Trans(i,j)= 0;
       Size_Trans(i,j)= cumd_norm(((size_vector(j)+(size_vector(j+1)-size_vector(j))/2)-Average_Size(i))/Sd(i));
       Size_Trans(i,j)-= cumd_norm(((size_vector(j)-(size_vector(j)-size_vector(j-1))/2)-Average_Size(i))/Sd(i));
      }
     //last bin;
     Size_Trans(i,nlength)= 1 - cumd_norm(((size_vector(nlength)-(size_vector(nlength)-size_vector(nlength-1))/2)-Average_Size(i))/Sd(i));
    }
    
   RETURN_ARRAYS_DECREMENT(); // Need this to decrement the stack increment
                       // caused by RETURN_ARRAYS_INCREMENT();
   return(Size_Trans);
   
  }
  
//-----------------------------------------------------------------------------------
// SEE if we have any problem because I didn't declare the arguments to be const
// According to Jim this may cause the arguments to be changed within the function
// But there is a limit to the number of strings I can pass, so if I include the const will exceed the limit

FUNCTION dvar_vector DoubLogistic(dvariable& pk,dvariable& in,dvariable& infl,dvariable& sl,dvariable& fin,dvariable& infl2,dvariable& sl2)
  {
   //This code is based on POP model (from Ian J. Stewart, SAFS-UW) 
   RETURN_ARRAYS_INCREMENT(); //Need this statement because the function
   dvar_vector sel_at_length(1,nlength);
                                   
   for (j=1; j<=nlength; j++)  //calculate the value over length bins
	{
	  if (double(j) < pk) // ascending limb
	  {
	     sel_at_length(j)= in + 
		    (1 - in)/((1 / (1 + (mfexp(-1 * sl * (pk - infl))))) - 
	            (1 / (1 + (mfexp(-1 * sl * (1 - infl)))))) *
                    ((1/(1+mfexp(-1*sl*(size_vector(j)-infl)))) - 
		    (1 / (1 + (mfexp(-1 * sl * (1 - infl))))));
	  }
	  else
	  {
	      if (double(j) > (pk + 1)) // descending limb
	      {
	          sel_at_length(j) = 1 + 
		        (fin - 1)/((1 / (1 + (mfexp(-1 * sl2 * (size_vector(nlength) - infl2))))) - 
	                (1 / (1 + (mfexp(-1 * sl2 * ((pk+1) - infl2)))))) *
	                ((1/(1+mfexp(-1*sl2*(size_vector(j)-infl2)))) - 
		        (1 / (1 + (mfexp(-1 * sl2 * ((pk+1) - infl2))))));
	      }
	      else // between the peaks
	      {
	   	  sel_at_length(j) = 1.0;
	      };
           };
        };

   RETURN_ARRAYS_DECREMENT(); // Need this to decrement the stack increment
                       // caused by RETURN_ARRAYS_INCREMENT();
   return(sel_at_length);
  }   
   
FUNCTION dvar_vector DoubNormal(dvariable& mean1,dvariable& mean2,dvariable& width1,dvariable& width2)
  {
   RETURN_ARRAYS_INCREMENT(); 
   
   dvar_vector sel_at_length(1,nlength);
   dvariable MaxSel;
   
   MaxSel = 0;
   for (j=1; j<=nlength; j++)  //calculate the value over length bins
	{
     sel_at_length(j) = 1.0/(1.0+mfexp(-1*(size_vector(j)-mean1)/width1))
                           /(1.0+mfexp(-1*(mean2-size_vector(j))/width2));
     if (sel_at_length(j) > MaxSel) MaxSel = sel_at_length(j);
    }
   for (j=1; j<=nlength; j++) sel_at_length(j) = sel_at_length(j) / MaxSel;

   RETURN_ARRAYS_DECREMENT();
   return(sel_at_length);
  }   
   
FUNCTION dvar_vector DoubNormalAge(dvariable& mean1,dvariable& mean2,dvariable& width1,dvariable& width2)
  {
   RETURN_ARRAYS_INCREMENT(); 
   
   dvar_vector sel_at_age(0,nages);
   dvariable MaxSel;
                                   
   MaxSel = 0;                                
   for (j=0; j<=nages; j++)  //calculate the value over length bins
	{
     sel_at_age(j) = 1.0/(1.0+mfexp(-1*(age_vector(j)-mean1)/width1))
                        /(1.0+mfexp(-1*(mean2-age_vector(j))/width2));
     if (sel_at_age(j) > MaxSel) MaxSel = sel_at_age(j);
    }
   for (j=0; j<=nages; j++) sel_at_age(j) = sel_at_age(j) / MaxSel;

   RETURN_ARRAYS_DECREMENT(); 
   
   return(sel_at_age);
  }   
   
//*======================================================================================*=/

FINAL_SECTION
  //Calculates how long is taking to run
  // this code is based on the Widow Rockfish model (from Erik H. Williams, NMFS-Santa Cruz)  
  time(&finish);
  elapsed_time = difftime(finish,start);
  hour = long(elapsed_time)/3600;
  minute = long(elapsed_time)%3600/60;
  second = (long(elapsed_time)%3600)%60;
  cout<<endl<<endl<<"starting time: "<<ctime(&start);
  cout<<"finishing time: "<<ctime(&finish);
  cout<<"This run took: ";
  cout<<hour<<" hours, "<<minute<<" minutes, "<<second<<" seconds."<<endl<<endl<<endl;


