C========================================================================
C This is the multispecies technical interaction MSE model adapted from Carey McGilliard's code
C 
C General model description:
C This model sets the OM, EM, CR for the multispecies tech int model and uses Jim Ianelli's (modified) tech int admb code to implement the multispecies catch tradeoffs
C The estimation model is based on CAB (very simple age structured assessment model)     
C
C Authors:
C Kotaro Ono
C Carey McGilliard      
      
C========================================================================

	MODULE ShareDataRatioSubs
	IMPLICIT NONE
	SAVE


      INTEGER MaxSpec,MinPYr,MaxPYr,MaxAge,Nsex,MaxTime
      INTEGER MaxNSBlk,MaxEWBlk,MaxSim,MaxLengthBins
      INTEGER MaxAssArea,MaxOpenArea,Max3StrataArea,MaxConfigs,MaxMVN
      PARAMETER (MaxSpec=3,MinPYr=1,MaxPYr=200,MaxTime=1)
      PARAMETER (MaxAge=40,Nsex=2,MaxNSBlk=30,MaxEWBlk=5)
      PARAMETER (MaxAssArea=10,MaxOpenArea = 12, Max3StrataArea = 20)
	PARAMETER (MaxSim = 102)
      PARAMETER (MaxLengthBins = 100)
      PARAMETER (MaxConfigs = 5)
      PARAMETER (MaxMVN = 1001)
C
C     Major parameters
      INTEGER :: Isim
      INTEGER Nspec,Amax(MaxSpec)
      INTEGER Ntime,NNSBlock,NEWBlock,UsePrevVals
      INTEGER BaseYr,FinDepYr,FirstFishYr
	REAL*8  InitDep(MaxSpec),FinalDep(MaxSpec)
      REAL*8  Qvalseed(MaxSpec),Rtime,SigmaR1(MaxSpec)
      REAL*8  SigmaR2(MaxSpec),SigmaSpat(MaxSpec)


C      Nspec        - Number of species
C      Amax(s)      - Maximum age for species s
C      InitDep(s)   - Initial depletion at year 1 for species s (Carey added this)
C      FinalDep(s)  - Final depletion at year FinDepYr
C      Qvalseed(s)  - Seed for Qval (catchability coefficient) for species s (Carey added)
C      Ntime        - Number of time-steps during the year
C      Rtime        - Width of each time step
C      NNSBlock     - Number of blocks (North / South)
C      NEWBlock     - Number of blocks (East / West)
C      SigmaR1(s)   - SigmaR (global) by species
C      SigmaR2(s)   - SigmaR (local) by species
C      UsePrevVals  - Read the q-values from SAVE.ALL  (if 1)
C      SigmaSpat(s) - Spatial variation in catchability
C      BaseYr       - Base year for simulation
C      FinDepYr     - Year in which 2nd depletion is pre-specified
C      
C     Specifications for closed areas
      INTEGER:: NClosedAreas
      INTEGER MinCEW(MaxAssArea),MaxCEW(MaxAssArea)
      INTEGER MinCNS(MaxAssArea),MaxCNS(MaxAssArea)
      INTEGER MinOstrataNS(MaxOpenArea),MaxOstrataNS(MaxOpenArea)
      INTEGER MinOstrataEW(MaxOpenArea),MaxOstrataEW(MaxOpenArea)

      INTEGER NoSurveyClosed
      REAL*8 Closed(MaxNSBlk,MaxEWBlk)


C      NclosedAreas   - Number of closed areas
C      MinCEW(a)      - boundary of closed areas      
C      MaxCEW(a)      - boundary of closed areas      
C      MinCNS(a)      - boundary of closed areas      
C      MaxCNS(a)      - boundary of closed areas    
C	 MinOStrataNS(a)		- boundary of open areas
C	 MaxOStrataNS(a)		- boundary of open areas
C	 MinOStrataEW(a)		- boundary of open areas
C	 MaxOStrataEW(a)		- boundary of open areas
C      Closed(n,e)    - is the current cell (North/South, East/West) closed?
C      NoSurveyClosed - Set to 1 for there to be no surveys in closed areas
C
C     Specifications for effort allocation blocks
      INTEGER NEffortAreas
      INTEGER MinEEW(0:MaxAssArea),MaxEEW(0:MaxAssArea)
      INTEGER MinENS(0:MaxAssArea),MaxENS(0:MaxAssArea)


C      NeffortAreas - Number of effort areas
C      MinEEW(a)    - boundary of effort areas      
C      MaxEEW(a)    - boundary of effort areas      
C      MinENS(a)    - boundary of effort areas      
C      MaxENS(a)    - boundary of effort areas      
C
C     Global N variables
      REAL*8 N(MaxSpec,Nsex,MaxNSBlk,MaxEWBlk,MinPYr:MaxPYr,0:MaxAge)
      REAL*8 Neq(MaxSpec,Nsex,MaxNSBlk,MaxEWBlk,0:MaxAge)
      REAL*8 SpBio(MaxSpec,MaxNSBlk,MaxEWBlk,MinPYr:MaxPYr)  
      REAL*8 SpBioTot(MaxSpec,MinPYr:MaxPYr)
      REAL*8 SpBioArea(MaxSpec,MaxAssArea,MinPyr:MaxPYr)
      REAL*8 RecSum(MaxSpec,MinPYr:MaxPYr,MaxNSBlk,MaxEWBlk)
      REAL*8 RecVirg(MaxSpec,MaxNSBlk,MaxEWBlk)
      REAL*8 SSBEq(MaxSpec)
      REAL*8 Nequ(MaxSpec,Nsex,0:MaxAge)
      REAL*8 Steep(MaxSpec)
      REAL*8 SRAlpha(MaxSpec,MaxNSBlk,MaxEWBlk)
      REAL*8 SRBeta(MaxSpec,MaxNSBlk,MaxEWBlk)
      REAL*8 SRR0(MaxSpec,MaxNSBlk,MaxEWBlk)
      INTEGER SrTYpe(MaxSpec)
	INTEGER RecruitOpt(MaxSpec)


C      N(s,)              - Numbers-at-age
C      SpBio(s,y)         - Spawning biomass
C      Nequ(s,g,a)        - Equilibrium numbers (by age)
C      Neq(s,g,ans,aew,a) - Equilibrium numbers (by age and space)
C      SSBEq(s)           - Spawning biomass-per-recruit 
C      SpBioTot(s,y)      - Spawning biomass (total)
C      Steep(s)           - Steepness of the stock-recruitment relationship
C      SRAlpha(s,ans,aew) - SR parameter
C      SRBeta(s,ans,aew)  - SR parameter
C      SrType(s)          - SR type
C      SrR0(s,ans,aew)    - SR virgin recruitment
C      RecSum(s,y,n,s)    - Recruitment (age 0)
C      RecVirg(s,n,s)     - Distribution of recruits
C      SpBioArea(s,a,y)   - Spawning biomass (by assessment area)
C
C     Biological parameters
      REAL*8 M(MaxSpec),SBPR(MaxSpec),Fec(MaxSpec,0:MaxAge)
      REAL*8 Wght(MaxSpec,Nsex,0:MaxAge,0:MaxTime)
      REAL*8 MixN(MaxSpec),MixE(MaxSpec)
      REAL*8 X1N(MaxSpec),X2N(MaxSpec),V1N(MaxSpec),V2n(MaxSpec)
      REAL*8 X1E(MaxSpec),X2E(MaxSpec),V1E(MaxSpec),V2E(MaxSpec)
      REAL*8 MixMat(MaxSpec,MaxNSBlk,MaxEWBlk,MaxNSBlk,MaxEWBlk)
      REAL*8 RecRes(MaxSim,MaxSpec,MinPYr:MaxPYr,MaxNSBlk,MaxEWBlk)
      REAL*8 RecQ(MaxSim,MaxSpec,MaxNSBlk,MaxEWBlk)
      REAL*8 RecX(MaxSpec),RecY(MaxSpec)
      REAL*8 VRecXL(MaxSpec),VRecYL(MaxSpec)
      REAL*8 VRecXH(MaxSpec),VRecYH(MaxSpec)
      REAL*8 APar(MaxSpec,Nsex),Bpar(MaxSpec,Nsex),kpar(MaxSpec,Nsex)


C      M(s)              - Natural mortality for species s     
C      SBPR(s)           - Spawner biomass-per-recruit
C      Fec(s,a)          - Fecundity
C      Wght(s,g,a,i)     - Weight-at-age (by time-step)
C      MixN(s)           - North-south mixing rate for species s
C      MixE(s)           - East-west mixing rate for species s
C      MaxMat(s,n,s,n,s) - Mixing matrix
C      RecRes(s,y,n,s)   - Spatially-varying local recruitment residual
C      RecQ(s,n,s)       - Spatially-varying catchability
C      RecX/Y(s)         - Centre of recruitment for species S
C      VRecX/Y&LH         - How recruitment of species S is distributed
C      X1/2E/N           - Attractiveness of location (X-values)
C      V1/2E/N           - Attractiveness of location (Y-values)
C      Apar, Bpar        - Length-weight parameters
C      kpar              - growth parameter
C
C     Catch, Effort, etc.
      REAL*8 SCL50(MaxSpec),SCL95(MaxSpec),SSL50(MaxSpec),SSL95(MaxSpec)
      REAL*8 TotEffort(MinPYr:MaxPyr),Qval(MaxSpec)
	REAL*8 TotEffortSeed
      REAL*8 ExpBio(MaxSpec,MaxNSBlk,MaxEWBlk,MinPYr:MaxPYr,MaxTime)
      REAL*8 Sel(MaxSpec,0:MaxAge)
      REAL*8 FullFish(MaxSpec,MaxNSBlk,MaxEWBlk,MinPYr:MaxPYr,MaxTime)
      REAL*8 CatBio(MaxSpec,MaxAssArea,MinPYr:MaxPYr)
      REAL*8 CatNum(MaxSpec,MaxAssArea,MinPYr:MaxPYr,Nsex,0:MaxAge)
      REAL*8 CatArea(MaxSpec,MaxNSBlk,MaxEWBlk,MinPYr:MaxPYr)
      REAL*8 TruExp(MaxSpec,MaxAssArea,MinPYr:MaxPYr)
      REAL*8 TruSuv(MaxSpec,MaxAssArea,MinPYr:MaxPYr)
      REAL*8 SuvNum(MaxSpec,MaxAssArea,MinPYr:MaxPYr,Nsex,0:MaxAge)
      REAL*8 Effort(MaxAssArea,MinPYr:MaxPYr)
      REAL*8 NetEffort(MaxNSBlk,MaxEWBlk,MinPYr:MaxPYr)



C      TotEffort(y)        - Effort by year (summed over all cells in the model)                 
C	 TotEffortSeed	     - Seed for the Total Effort from year - 19 to 0.
C      Qval(y)             - q value
C      Sel(s,a)            - Selectivity-at-age
C      ExpBio(s,g,s,b,y)   - Exploitable biomass (mid-season)
C      FullFish(s,n,w,y,t) - Fully-specified F
C      CatBio(s,a,y)       - Catch biomass
C      CatNum(s,a,y,g,a)   - Catch numbers
C      CatArea(s,i,j,y)    - Catch biomass (by cell)
C      TruExp(s,a,y)       - Exploitable biomass
C      TruSuv(s,a,y)       - Survey biomass
C      SuvNum(s,a,y,g,a)   - Survey numbers
C      Effort(a,y)         - Effort by year
C      NetEffort(n,s,y)    - Effort by year (and spatial cell)
C      SCL50,SCL95         - Parameters to determine selectivity (fishery)
C      SSL50,SSL95         - Parameters to determine selectivity (survey)

C
C     Data for assessment purposes
      REAL*8 ObsCatch(MaxSpec,0:MaxAssArea,MinPyr:MaxPyr)
      REAL*8 Cpue(MaxSpec,0:MaxAssArea,MinPYr:MaxPyr)
      REAL*8 CVCpue(MaxSpec,0:MaxAssArea,MinPYr:MaxPyr)
      REAL*8 SurveyEst(MaxSpec,0:MaxAssArea,MinPYr:MaxPyr)
      REAL*8 SurveyCV(MaxSpec,0:MaxAssArea,MinPYr:MaxPyr)
      INTEGER ComSS(MaxSpec,0:MaxAssArea,MinPYr:MaxPyr)
      INTEGER SuvSS(MaxSpec,0:MaxAssArea,MinPYr:MaxPyr)
      INTEGER EffectSampCatch1,EffectSampSurvey1
      INTEGER EffectSampCatch2,EffectSampSurvey2
      REAL*8 CNObs(MaxSpec,0:MaxAssArea,MinPyr:MaxPyr,MaxSpec,0:MaxAge)
      REAL*8 SNObs(MaxSpec,0:MaxAssArea,MinPyr:MaxPyr,MaxSpec,0:MaxAge)
      REAL*8 Ssel(MaxSpec,0:MaxAge)
      REAL*8 CpueCV,SurvCV
      INTEGER IsDeterm,FirstSurvYr,FreqSurv


C      ObsCatch(s,a,y)     - Observed catch by species, area and year
C      CPUE(s,a,y)         - Cpue data by species, area and year
C      CVCPUE(s,a,y)       - CV of the CPUE data
C      SurveyEst(s,a,y)    - Survey estimates by species, area and year
C      SurveyCV(s,a,y)     - Survey CVs by species, area and year
C      ComSS(s,a,y)        - Effective sample size (Commercial data)
C      SuvSS(s,a,y)        - Effective sample size (Survey data)
C      CNObs(s,a,y,g,a)    - Commercial catch-at-age data
C      SNObs(s,a,y,g,a)    - Survey catch-at-age data
C      Ssel(s,a)           - Survey selectivity
C      IsDeterm            - Are the data deterministic?
C      FirstSurvYr         - First year with a survey
C      FreqSurv            - Frequency of surveys
C      CpueCV              - CV of CPUE data
C      SurveyCV            - CV of Survey data
C      EffectSampCatch1/2  - Effective Sample Size (Catch) [Pre/with MPAs]
C      EffectSampSurvey1/2 - Effective Sample Size (Survey) [Pre/with MPAs]
C      
C     Pre-specified parameters
      INTEGER AmaxEst(MaxSpec),IbiasType,IestIgR(MaxSpec)
      INTEGER Narea,AssOpt
      INTEGER MinEW(MaxAssArea),MaxEW(MaxAssArea)
      INTEGER MinNS(MaxAssArea),MaxNS(MaxAssArea)
      INTEGER PhaseM,PhaseS0,PhaseH,PhaseSel
      INTEGER NhistAss,HistAssYr(100)
      INTEGER NoCPUEMPA
      REAL*8 EstM(MaxSpec,0:MaxConfigs),IestSteep(MaxSpec,0:MaxConfigs)
      REAL*8 Iestln_S0(MaxSpec)
      REAL*8 IestSigmaR(MaxSpec),WghtEst(2,0:MaxAge)
      REAL*8 BoundM,BoundS0,BoundH,BoundSel,BoundSela50,BoundSelSlope
      REAL*8 EstCpueCV,EstEffCatchSamp,EstSurveyCV,EstEffSurveySamp


C      AmaxEst(s)          - Maximum age for assessment purposes
C      EstM(s,Iconfig)   - Assumed value for M in the assessment
C      IestSteep(s,Iconfigs) - Steepness
C      Iestln_S0(s)        - Log of something related to SSB0: SSB0 if single sex model; SSB0*2 if 2-sex model
C      IestSigmaR(s)       - SigmaR
C      IbiasType           - Recruitment bias-correction
C      IestIgR             - Recruitments to Ignore
C      Min/MaxEW(a)        - EW bounds for area a
C      Min/MaxNS(a)        - NS bounds for area a
C      Narea               - Number of assessment areas
C      WghtEst(s,a)        - Estimated weight-at-age
C      Phase/BoundM        - Phase and bound for M
C      Phase/BoundS0       - Phase and bound for S0
C      Phase/Boundh        - Phase and bound for H
C      Phase/BoundSel      - Phase and bound for Selectivity
C      BoundSela50,BoundSelSlope    - Bound for a50 and slope for SS3 runs (selectivity params)
C      EstCpueCV           - Cpue CV used in CAB
C      EstSurveyCV         - Survey CV used in CAB
C      EstEffCatchSamp     - Effective sample size (commercial catch-at-age)
C      EstEffSurveySamp    - Effective sample size (survey catch-at-age)
C      NHistAss            - Number of historical assessment years
C      HistAssYr(i)        - Years in which assessments are conducted
C      AssOpt              - Assessment option to run (0=all; 1=single area; 2=N areas; 3=Nfleets)
C      NoCPUEMPA           - Don't produce CPUE data after year 1
C
C     Estimated parameters
      REAL*8 EstSB(MaxSpec,MaxAssArea,MinPyr:MaxPyr)
      REAL*8 EstRec(MaxSpec,MaxAssArea,0:MaxConfigs,MinPyr:MaxPyr)
      REAL*8 EstN(MaxSpec,MaxAssArea,Nsex,MinPyr:MaxPyr,0:MaxAge)
      REAL*8 EstS(MaxSpec,MaxAssArea,Nsex,0:MaxAge,0:MaxConfigs)
      REAL*8 SteepCAEst(MaxSpec,MaxAssArea)
      REAL*8 SigmaRCAest(MaxSpec,MaxAssArea)
      REAL*8 R0CAEst(MaxSpec,MaxAssArea)
      REAL*8 S0CAEst(MaxSpec,MaxAssArea,0:MaxConfigs)
      REAL*8 PassEstVar(MaxSpec,MinPyr:MaxPyr,1:8)


C      EstSB(s,a,y)        - Estimated SSB by species, area and year      
C      EstRec(s,a,c,y)       - Estimated Recruitment by species, area and year and assessment configuration     
C      EstN(s,a,g,y,a)     - Estimated numbers-at-age
C      EstS(s,a,g,a,c)       - Estimated selectivity-at-age
C      EstSB(s,a,y)        - Estimated SSB by species, area and year      
C      SteepCAEst          - Estimated steepness
C      SigmaREst           - Estimated Sigma(R)
C      R0CAEst             - Estimated R0
C      S0CAEst             - Estimated S0
C      PassEstVar          - Material at the stock level
C
C     Misc quantities
      INTEGER TheProjYear,MinOpt,Nsim,NprojYr,ISEEDO,ISEEDP
      INTEGER SolFinDepl,ImplOpt
      REAL*8 Quota(MaxSpec),QuotaAct(MaxSpec)

C      TheProjYear         - Year of projection
C      MinOpt              - Option for minimization
C      Quota(s)            - Species quota      
C      QuotaAct(s)         - Quota for model run (sum over areas)
C      Nsim                - Number of simulations
C      NprojYr             - Number of projection years for MSE only (Carey thinks)
C      ISEEDO              - Random number seed
C      ISEEDP              - Process error
C      SolFinDepl          - Solve for the effort that gets the final depletion
C      ImplOpt             - Option for Implementation uncertainty


C==========================================================================================
C	Carey's variables in module form only: Variables that Carey added to the program
C	Shared data declarations
	 REAL(KIND=8), DIMENSION(Maxspec,Nsex,MaxNSBlk,MaxEWBlk,MinPYr:MaxPYr,
     +	 0:Maxage) :: SpatialSuvNum !Leaving this variable when taking out density ratio project - might be generally useful.
	 REAL(KIND=8), DIMENSION(Maxspec,MinPYr:MaxPYr,MaxNSBlk,MaxEWBlk) :: 
     +	 SpatialSuvNumSum !Leaving this variable when taking out density ratio project
     	REAL(8),DIMENSION(MinPYr:MaxPYr) :: GlobalEffort !Leaving this variable when taking out density ratio project
	REAL(8) :: CareyGRD  !Leaving this variable when taking out density ratio project
	REAL(8) :: ConstantEff,MSYEffort !Leaving this variable when taking out density ratio project
	INTEGER :: MPAYear !Leaving this variable when taking out density ratio project - how did this one get all the way down to this area??
    	REAL(8), DIMENSION(MinPYr:MaxPYr) :: NSum,CatNumSum !Leaving this variable when taking out density ratio project
	REAL(8) :: MSYEffort1 !Leaving this variable when taking out density ratio project
	REAL(8), DIMENSION(MinPYr:MaxPYr) :: TotSpatialSuvNum !Leaving this variable when taking out density ratio project
	REAL(8),DIMENSION(MaxSpec,2,MinPYr:MaxPYr) :: SpBioMPA !Leaving this variable when taking out density ratio project
	INTEGER :: QValOne !Leaving this variable when taking out density ratio project
!      INTEGER :: SuperSim !Leaving this variable when taking out density ratio project


C	ConstantEff		- Constant Effort level
C	MSYEffort			- The effort level that produces the maximum sum of the long term discounted catch with an MPA 
C     MPAYear				- Year that MPA is implemented
C	MSYEffort1			- Effort at MSY under a particular value for steepness (buried within steepness optimization
C	TotSpatialSuvNum	- Numbers summed over age (vulnerable to survey and currently same as fishery vulnerability)
C	SpBioMPA			- Spawning biomass split by species and reserve status
C	QValOne				- 1 if letting Q = 1
C     #End original section that Carey added (all DR variables removed)
C #=================================================================================


      REAL*8 Linf(MaxSpec,Nsex),t0(MaxSpec,Nsex)
      REAL*8 L50(MaxSpec),L95(MaxSpec)
      REAL(8),DIMENSION(MaxSpec,0:Maxage) :: Mature
	INTEGER :: TotYrs, OptCRrecord

C	Linf				- Maximum length
C	t0					- Von Bertalanffy t0
C	L50, L95			- Length at 50% maturity; length at 95% maturity.
C     Mature              - Maturity-at-age vector based on L50 and L95 (which are ages, not lengths)

      REAL*8 MixNLarval(MaxSpec),MixELarval(MaxSpec)
      REAL*8 X1NLarval(MaxSpec),X2NLarval(MaxSpec),V1NLarval(MaxSpec)
      REAL*8 V2NLarval(MaxSpec)
      REAL*8 X1ELarval(MaxSpec),X2ELarval(MaxSpec)
	REAL*8 V1ELarval(MaxSpec),V2ELarval(MaxSpec)
      REAL*8 MixMatLarval(MaxSpec,MaxNSBlk,MaxEWBlk,MaxNSBlk,MaxEWBlk)

C      MixN(s)           - North-south mixing rate for species s
C      MixE(s)           - East-west mixing rate for species s
C      MaxMat(s,n,s,n,s) - Mixing matrix


C	

	REAL(8),DIMENSION(MaxSpec,MaxNSBlk,MaxEWBlk,MinPYr:MaxPYr)::
     +	TotRecruits
	REAL(8),DIMENSION(MaxSpec,MaxNSBlk,MaxEWBlk,MinPYr:MaxPYr)::Larvae
	REAL(8) :: FinDeplEffort,FinDepRecord

C	TotRecruits(s,i,j,y)	- Total number of recruits per cell post movement and recruitment
C	Larvae(s,i,j,y)		- Total number of larvae by cell prior to movement or recruitment
C	FinDeplEffort			- Specified constant effort value when not solving for Fin Depl.
C	FinDepRecord			- The Depletion found when trying to force the depletion level of FinDep.


	REAL(8) :: InitDeplRecord, BestInitEffort,BestFinalEffort
	REAL(8) :: 	mysystemtime

C	InitDeplRecord		- A record of the InitDep that results from the best initial effort level applied from year -19 to BaseYr
C	BestInitEffort		- Effort level that produces InitDepl (specified depletion level in BaseYr)
C	BestFinEffort			- Effort level that produces FinDep (specified depletion level in FinDepYr)
C	mysystemtime			- The amount of time it took to run the whole program

	REAL(8) :: AvgDepletion,AvgSSB,AvgCatch
      REAL(8) :: PretendVar

C	AvgDepletion			- Depletion at the end of NProjYrs, averaged over simulations in the function for finding MSY
C	AvgSSB					- SSB "" "" ""
C	AvgCatch				- Catch "" "" ""

	REAL(8) :: EffortMSY
	INTEGER :: DoMSYOnly
	INTEGER :: WriteOpt

C	EffortMSY				- Effort that produces MSY, averaged over simulations
C	DoMSYOnly				- Equal to 1 if choosing to find MSY for the input file setup and then exit out of the program entirely.
C	WriteOpt				- Equal to 1 if info from an MSY run should be written out: this is set correctly within the function for finding MSY and doesn't need to be manually set; isn't in an input file.

	REAL(8), DIMENSION(MaxSim,MaxSpec,MinPYr:MaxPYr) :: ObsResOpen,
     +	ObsResClosed

C	ObsResOpen, ObsResClosed	- Arrays to specify residual values for surveys in open and closed areas, respectively, under the density ratio control rule

C     Debugging variables:
      REAL(8), DIMENSION(MaxSpec,Nsex,0:Maxage) :: TestNeq
      INTEGER :: StartYrRecDevs
      REAL(8),DIMENSION(MaxSpec) :: Kvec
      INTEGER :: DebugAssessment
      INTEGER :: LastSuperSim
      INTEGER :: StartYr,EndYr
      REAL(8),DIMENSION(100) :: ConstantEffVec
      INTEGER :: NTargetDepl
      REAL(8),DIMENSION(10,MaxSpec) :: TargetDeplVec
      INTEGER, DIMENSION(10) ::DeplYr
    !  REAL(8),DIMENSION(MaxSpec) :: Depletion
      INTEGER :: PlusGroup
      INTEGER :: DoProcessError

      ! Kvec                - Carrying capacity in the model (refers to SSB for females only, I think)
      !                     - This replaces just specifying Kvals = 10000.0 inside of FUNK, which is impossible to find
      ! StartYrRecDevs      - This specifies the year that recruitment devs should start being estimated in the assessment model
      ! DebugAssessment     - This is a switch: 1 = run the assessment deterministically (don't estimate anything, just project); 0 = estimate parameters in the assessment
	! StartYr,EndYr       - Changes for each run of GetEffLev such that can find a constant effort level used from STartYr to EndYr such that the depletion level at EndYr is some specified value.
      ! Depletion           - Reports the depletion calcuated in FUNK back to GetEffLev and the main program.
      ! PlusGroup           - 1 = yes, do a plus group; 2 = no, don't do a plus group (e.g. if there are local recruitment dynamics AND adult movement, as for density ratio project)
      ! DoProcessError      - 1 = yes, include the process error; 0 = no, run the model with no process error for the purpose of debugging and comparing to deterministic version of estimation model (just projecting, not estimating)

! Stuff below was added because of SS3:
      REAL(8),DIMENSION(MaxSpec,Nsex) :: Lamin
      REAL(8),DIMENSION(Maxspec) :: MaturitySlope
      REAL(8),DIMENSION(Maxspec) :: R0forSS3

      !Lamin                - Size at min. age: used instead of t0 in the Von-Bert in SS3. conversion is: Lamin = Linf*(1-exp(-kpar*(amin-t0))); Final eqn is : L(age) = Linf+(Lamin-Linf)*exp(-kpar*(a-amin)))
      INTEGER :: PhaseRecDevs
      REAL(8),DIMENSION(MaxSpec) :: SlopeFishery,SlopeSurvey
      INTEGER,DIMENSION(MaxSpec) :: LengthBinIntervals
      INTEGER,DIMENSION(MaxSpec) :: NLengthBins
      INTEGER,DIMENSION(MaxSpec,0:MaxLengthBins) :: LengthBins
      INTEGER,DIMENSION(MaxSpec,0:Maxage) :: BinnedLength
      INTEGER,DIMENSION(MaxSpec,0:Maxage) :: BinnedLengthBeg
      REAL(8),DIMENSION(MaxSpec,NSex,MinPYr:MaxPYr,0:Maxage,0:MaxAss
     +Area) :: NatAgeCompact
      REAL(8),DIMENSION(MaxSpec,Nsex,0:Maxage,0:MaxTime) :: WghtBeg
      REAL(8),DIMENSION(MaxSpec,MaxAssArea,MaxAssArea) :: AreaMove
      REAL(8),DIMENSION(MinPYr:MaxPYr) :: DeplSS
      REAL(8),DIMENSION(MinPYr-3:MaxPYr) ::TotBioSS,ssbSS,RecruitsSS
      INTEGER :: StartFishing
      REAL(8),DIMENSION(MaxSpec,MinPYr:MaxPYr) :: TotalBiomass
      REAL(8),DIMENSION(MaxSpec,MinPYr:MaxPYr) :: TotalRecruits
      INTEGER :: UseAreaMove
      INTEGER,DIMENSION(18) :: PhaseMove
      REAL(8),DIMENSION(18) :: InitMove,BoundMove
      INTEGER,DIMENSION(9) :: MinMoveAge,MaxMoveAge
      INTEGER,DIMENSION(MaxAssArea) :: PhaseRAllocate
      REAL(8),DIMENSION(MaxAssArea) :: BoundRAllocate,InitRAllocate

      !If switching to MSE this would be MinPYr:MaxPYr instead of 1:3 (NHistAss)
      REAL(8),DIMENSION(1:3,1:MaxSim)::Bias1,Bias2,Bias3,Bias4,
     + Bias5
       REAL(8),DIMENSION(1:3,1:MaxAssArea)::BBias1,BBias2,BBias3,BBias4,
     + BBias5
      INTEGER :: ReadBiasPars
      REAL(8),DIMENSION(0:MaxAssArea)::TheBias1,TheBias2,TheBias3,
     + TheBias4,TheBias5
      INTEGER :: OMonly
      REAL(8),DIMENSION(MaxSpec) :: TheTarg
      REAL(8),DIMENSION(MaxAssArea,MaxAssArea,0:Maxage) :: PropN
 
      INTEGER :: DevsInRecAlloc,RecAllocMinYr,
     +    RecAllocMaxYr
      REAL(8) :: RecAllocStddev
      INTEGER :: PhaseRecAllocDevs
      REAL(8),DIMENSION(MaxSpec,MaxAssArea) :: R0forSS3B
      REAL(8) :: loBoundSigmaR,hiBoundSigmaR
      INTEGER :: PhaseSigmaR
      REAL(8) :: DeplFound
      INTEGER :: Hyperstability

      !For MPA and assessment paper review responses:
      INTEGER :: UseFisheryCPUE
      
      !PhaseRecDevs                 - The estimation phase for recruitment deviations within SS3
      !SlopeFishery,SlopeSurvey     - The slope in the logistic selectivity function (parameterized by slope and a50)
      !LengthBinIntervals           - The width of the length bins (2 means length bins are 0 2 4 6 8, etc.)
      !NLengthBins                  - Number of length bins to input to SS3 (and for binning lengths before calculating mid-year weight)
      !LengthBins                   - The length bins themselves (input to SS3)
      !BinnedLength                 - The binned version of the lengths at age (this is what's used to calculate weight at age)
      !NatAgeCompact                - Numbers at age, not by spatial area
      !WghtBeg                      - Weight at age at the beginning of the year (Wght is the middle of the year)
      !AreaMove                     - Probability of movement from one assessment area to another AreaMove(Iarea,Jarea) = pr(move from Iarea to Jarea)
      !DeplSS                       - Est. Depletion in each year from SS3 results
      !TotBioSS                     - Estimated total biomass in each year from SS3 results
      !ssbSS                        - Est. total spawning bimass each year from SS3 results
      !RecruitsSS                   - Est. total recruits each year from SS3 results
      !StartFishing                 - The first year that catches > 0
      !TotalBiomass                 - Total biomass in the OM in each year
      !TotalRecruits                - Total recruits in teh OM in each year
      !UseAreaMove                  - 1 = use AreaMove matrix to assign seeds for movement parameters for SS3
      !PhaseMove                    - Phase for each movement parameter (see ESTMODEL.dAT for order) for SS3
      !InitMove,BoundMove           - Initial values and bound (+- this amount) for movement parameters for SS3
      !MinMoveAge,MaxMoveAge        - Parameters of the SS movement function: minimum movement age and max. movement age.
      !PhaseRAllocate               - Phase for recruitment allocation parameters in SS3
      !BoundRAllocate,InitRAllocate - Bounds (+- this amount) and initial values for recruitment allocation parameters for SS3
      !Bias1 to Bias5               - Correct bias adjustment parameters for the DoAssess.ctl file (from running r4ss code after starter simulations)
      !OMonly                       - 1 = don't do any assessments, just run the OM portion of the program (useful for generating Qvals and Constant Effort levels for each sim that can be used for multiple management strategies on the same operating model - speeds things up)
      !DevsInRecAlloc               - Do yearly deviations in recruitment allocations?
      !RecAllocMinYr                - Minimum year for the devs. in rec. allocations
      !RecAllocMaxYr                - Maximum year for the devs. in rec. allocations
      !RecAllocStddev               - Standard deviation for annual devs. in recruitment allocations
      !PhaseRecAllocDevs            - Phase for deviations in recruitment allocations (if there are any)
      !R0forSS3B                    - R0 by assessment area (used to calc. lnR0 for Method B (three sep. assessments)
      !DeplFound                    - The depletion level found and read in from SaveConstantEff.All
      !Hyperstability               - 1 if effort occurs in areas of highest biomass and 0 if effort is evenly distributed among open areas
      !UseFisheryCPUE               - 1 if using the fishery cpue in the "single" configuration assessments and 0 to leave out fishery cpue

      
      !Added when updating to new cab
      REAL(8),DIMENSION(MaxConfigs) :: cab_obj_fun  
      
      !cab_obj_fun                - The objective function from a cab assessment run
      
      
      !Harvest control rules for technical interactions model
      INTEGER :: UsePstar
      INTEGER :: Nmvn
      REAL(8) :: ABC
      REAL(8),DIMENSION(1:MaxSpec,1:Nsex,1:MaxMVN,0:MaxAge,
     +  0:MaxConfigs) :: Nvec
      REAL(8),DIMENSION(1:MaxSpec,1:MaxMVN,0:MaxConfigs) :: 
     +SSB_MVNvec
      REAL(8) :: pCatBio
      REAL, DIMENSION(1:MaxSpec,1:Nsex,0:Maxage) :: MyNs
      REAL(8), DIMENSION(1:MaxSpec) :: TrialCatch
!      REAL(8),DIMENSION(1:MaxSpec,1:MaxConfigs) :: EstSPR
!      REAL(8),DIMENSION(1:MaxSpec,1:MaxConfigs) :: EstB_SPR
      REAL(8) :: EstSPR, EstB_SPR !can't have Ispec,Iconfig dimensions b.c. referred to within FUNK.
      INTEGER :: GlobalProjYr
      REAL(8),DIMENSION(1:MaxSpec,1:MaxConfigs) :: EstF35,EstB35
      REAL(8),DIMENSION(1:MaxSpec,1:MaxConfigs) :: EstF40,EstB40
      REAL(8) :: AFSC_alpha
      INTEGER :: AFSC_MinusYrs
      REAL(8),DIMENSION(1:MaxSpec,1:MaxConfigs,1:MaxMVN) 
     +:: EstFOFL,SSEstFOFL
      REAL(8),DIMENSION(1:MaxSpec,MinPYr:MaxPYr) :: TrueFOFL
      REAL(8),DIMENSION(1:MaxSpec,MinPYr:MaxPYr) :: TrueB35,TrueB40
      REAL(8),DIMENSION(1:MaxSpec) :: TrueF35,TrueF40
!      REAL(8),DIMENSION(1:MaxSpec) :: SeedQuota
      REAL(8) :: SeedQuota
      REAL(8) :: Pstar
      REAL(8) :: EstimatedP
      INTEGER :: Nconfigs
      REAL(8),DIMENSION(1:MaxSpec,1:MaxConfigs) :: ConfigWts !,hConfigs, MConfigs
      REAL(8) :: TheM,TheX
      REAL(8), DIMENSION(1:Nsex,0:Maxage) :: ThoseNs
      REAL(8), DIMENSION(1:Nsex,0:Maxage) :: TheSelex
      REAL(8), DIMENSION(1:Nsex,0:Maxage) :: MyWeight
!      REAL(8) :: LoBoundM,HiBoundM,LoBoundh,HiBoundh
      
       REAL(8),DIMENSION(1:MaxSpec,MinPYr:MaxPYr,0:MaxConfigs,
     +1:MaxMVN) :: FDistn,SSFDistn
       REAL(8),DIMENSION(1:MaxSpec,MinPYr:MaxPYr) :: TrueF,SSTrueF
       REAL(8),DIMENSION(0:Maxage) :: TheFecVec
       REAL(8),DIMENSION(MinPYr:MaxPYr) :: TheRecVec
       INTEGER :: ThePlusAge
       INTEGER :: GlobalIspec
       INTEGER :: Funk2_MinOpt
       REAL(8) :: TrialQuota
       REAL(8),DIMENSION(1:MaxSpec,1:MaxConfigs) :: EstSPR35,SSEstspr35
       REAL(8),DIMENSION(1:MaxSpec,1:MaxConfigs) :: EstSPR40,SSEstspr40
       REAL(8),DIMENSION(1:MaxSpec) :: TrueSPR35,SSTrueSPR35,
     +TrueSPR40,SSTrueSPR40
      LOGICAL :: cor_exists
      CHARACTER(11) :: TestMVN
      
      !Pstar                      - (1 = use the pstar approach, 0 = use AFSC Tier 3 HCR)
      !Nmvn                       - Number of multivariate normal draws of log(N-at-age) used for pstar control rule 
      !SSB_MVNvec                 - The current SSB associated with each vector of random Numbers-at-age generated using the MVN distribution
      !ABC                        - The ABC for the following year
      !MyNs                       - This is needed to define true N for most recent year within FUNK (which isn't passed info on year)
      !TrialCatch                 - The trial ABC or OFL or MSY for the P* control rule; needed as a starting point for figuring out the F that would produce the TrialCatch and more generally the F that would result in the catch associated with a P* value for each random vector of Ns and the true N.
      !GlobalProjYr               - Globally defines what year of the MSE loop is the current year so as to use the year within the FUNK function when doing optimizations for the p* HCR.
      !EstF35,EstB35,EstF40,EstB40   - estimated reference pts needed to calculate the Estimated FOFL (EstFOFL) for the AFSC OFL control rule. B's are SSB and will change by year within the MSE due to the years included in the calc for avg recruitment
      !AFSC_alpha                 - The alpha used in the AFSC OFL control rule (usually 0.05)
      !AFSC_MinusYrs              - The number of years to leave out of calculating Avg Recruitment for the AFSC OFL control rule
      !EstFOFL                    - The estimated FOFL for each configuration and random Natage vector for use in the P* control rule.
      !TrueFOFL                   - The true AFSC FOFL based on the OM
      !TrueF35,TrueF40,TrueB35,TrueB40 - true ref pts needed to calculate TrueFOFL
      !SeedQuota                  - The seed value (a trial ABC value to start with for doing the calculations necessary to come up with a quota with a pr(overfishing) = p* as specified.
      !Pstar                      - The p* specified in CR.dat that represents a (manager-determined) acceptable p* (probability of overfishing).
      !Nconfigs                   - Number of alternative assessment configurations (combos of Ms and hs) to use in assessments when using the p* approach (an assessment is done assuming each configuration)
      !ConfigWts                  - Weights (probabilities) assigned to each config combo for a species (used in p* control rule)
 !canceled idea     !LoBoundM,HiBoundM          - Lower bound and higher bound on estimated M (for p* approach needed to have bounds that would work for all configs, rather than the bound being subtracted and added to the value, which might not make sense).
 !canceled idea     !LoBoundh,HiBoundh          - As for LoBoundM, etc., but for steepness.
      !TheM,TheX                  - Iconfig and Ivec (for TheX) specific M and F values that are input into funk for finding the F(t) that matches a trial ABC value as part of the P* control rule.
      !ThoseNs                      - Iconfig and Ivec specific values for Numbers-at-age to input to funk (MinOpt = 10) to find the F(t) that matches a trial aBC.
      !TheSelex                   - as for TheNs, but for estimated (or true selectivity) MinOpt = 10 (or 11 for true F(t), rather than estimated F(t)).
      !MyWeight                   - Weight at age for a particular species for entry to CalcF...for a particular species (without needing Ispec)
      !FDistn,SSFDistn            - F resulting in the trial ABC for each iteration of the random numbers-at-age vector generated as part of the p* approach nad associated SS value (shoudl be near 0 if the optimization to find the F resulting in the trial abc worked.
      !TrueF,SSTrueF              - True fishing mort resulting in a particular trial abc and associated SS value (should be near 0)
      !TheFecVec                  - Fec at age defined for a particular species for application within FUNK, where Ispec, etc. are not directly recognized (to go into AFSC_SPR subroutine)
      !TheRecVec                  - Recruitment over time for a particular species and configuration for application within FUNK, where Ispec, etc. are not directly recognized (to go into AFSC_SPR subroutine)
      !ThePlusAge                 - The Amax(Ispec) without having to refer directly to Ispec within FUNK.
      !GlobalIspec                - globally defines progress through the big Ispec loop so that Ispec is defined within subroutines used in the P* optimization
      !Funk2_MinOpt               - For use with FUNK2 and FIT2, which is needed to avoid conflicts when doing the outer layer of the nested optimization required (p = p* by changing SeedQuota) so as to avoid conflicts with nested optimizations (to determine the joint distn of (F(t),FOFL)
      !TrialQuota                 - Quota that changes with the step of the P* optimization; SeedQuota is the initial value for TrialQuota.
      !EstSpr35,EstSPR40          - Should be equal to 0.35 and 0.40, respectively
      !SSEstspr35,SSEstspr40      - Should be close to 0; SS values from optimizing to get SPR = 0.35 and 0.40
      !TrueSPR35,TrueSPR40        -  Should be equal to 0.35 and 0.40, respectively
      !SSTrueSPR35,SSTrueSPR40    - Should be close to 0; SS values from optimizing to get SPR = 0.35 and 0.40
      !cor_exists                 - True = cab.cor exists (used if skipping over a run that didn't converge).
      !TestMVN                    - Will equal "Correlation" if MVN ran incorrectly.
      END MODULE ShareDataRatioSubs
C======================================================================




C Output files
C =============
C SUMM1.OUT   Specifications for the simulation (also pasted at the bottom of other output files)
C SUMM.OUT    Biological parameters
C SUMM2.OUT   Total: SSB, Depl, Catch (in biomass), Quota (Unit 14)
C SUMM3.OUT   Depl (by region) (Unit 15)     
C SUMM4.OUT   Recruits (scaled to the maximum recruits over blocks) (Unit 16)
C SUMM5.OUT   Depl (scaled to maximum SSB over blocks) (Unit 17)     
C SUMM6.OUT   Effort (scaled to maximum effort over blocks) (Unit 18)        
C SUMM7.OUT   Mixing matrix (Unit 19)
C SUMM8.OUT   Numbers-at-age by block (Unit 26)
C SUMM9.OUT   Year, Catch, Effort, Quota (Unit 27)
C SUMM10.OUT  Stuff for quota setting (Unit 29)
C SUMM11.OUT  Effort (unscaled) (Unit=32)
C SUMM12.OUT  Distribution of recruitment spatially (Unit=33)
C SUMM13.OUT  Spawning biomass by block (Unit 34)     
C SUMM14.OUT  Catch by block (Unit 35)
C SUMM17.OUT  GlobalEffort not divided by blocks, but overall (Unit 38)
C SUMM18.OUT  Optimization results: Qval found by solving for the Q that produces a particular depletion level given the effort series (Unit 39)
C SUMM19.OUT  Total numbers, total catch in numbers (2 columns)
C SUMM20.OUT  Spawning biomass by MPA status (Unit 41)
C SUMM27.OUT  Ispec, IBlkNS,IBlkEW,Year,Larvae,Age 0s (used for testing larval dispersal)
C SUMM29.OUT  Constant effort for years 1 to 30 (so you can do subsequent runs setting this and not having to find it again) (Unit 51)
C SUMM30.OUT  MSY information from true MSY runs
C CatchCPUE.OUT  Sim,Yr,Ispec,Iarea,CPUE; Iarea = 0 is aggregated across areas
C SurveyCPUE.OUT  Sim,Yr,Ispec,Iarea,Survey CPUE; Iarea = 0 is aggregated across areas
C RECRUITS.OUT Recruits in each cell (males and females) made to be smaller than SUMM8.OUT, but the same info for age 0s. (Unit 260)


      PROGRAM MPAMODEL

	USE ShareDataRatioSubs
C	USE mpamodelsubs
      IMPLICIT NONE
C
      !------------------------------------------------------------
C     Local variables
      INTEGER Yr,ISpec,Iarea,IBlkNS,IBlkEW,SEED1,SEED2,counter,counter1
      INTEGER :: itarg,mystyr,myendyr,IAssess
      REAL(8),DIMENSION(MaxSpec) :: mytarg

C	Time counter variables
      real(8) :: etime          ! Declare the type of etime()
      real(8),DIMENSION(2) :: elapsed     ! For receiving user and system time
      real(8) :: total          ! For receiving total time

C	INTEGER Isim
      CHARACTER*1 MYINDC
!     CHARACTER*2 HistChar1,HistChar2,HistChar3
      INTEGER Iass,U,VV,i,iu,TheSim
      REAL*8 X(100),SS
      REAL*8 FUNK !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      REAL*8 FUNKSCALAR
C     Carey: EXTERNAL means to use the user-defined definition of this function (FUNK), rather than an intrinsically defined function, like what's written into the FORTRAN programming language
      EXTERNAL FUNK
      EXTERNAL FUNKSCALAR
      !End local variables and definition of externals
C     !---------------------------------------------------------------------

      

! What time did we start?
      OPEN(UNIT=2001,FILE='results/TimeDateStamp.OUT')
      CALL mydate('Start time and date      ')

C     Clear some files
	CALL SYSTEM("Del SUMM1.OUT")
      CALL SYSTEM("Del SUMM1.EST") 
      CALL SYSTEM("Del SUMM2.EST") 
      CALL SYSTEM("Del REC.EST")
	CALL SYSTEM("Del SUMM29.OUT")
C
C     Read in the specifications for the run
      CALL ReadSpec()
      CALL Summary() !reads out fecundity and weight at age
      CALL ReadCntl()
      CALL ReadCR()


!	OPEN(UNIT=100, FILE = "SUMM1.OUT") !Acts like an echoinput file
!	CALL WriteSpecs(100) !Acts like an echoinput file, but would need updating to list variables by name. Kind of useless/outdated, but could be made useful, therefore kept.


C     Set up specifications that are independent of recruitment residuals
      CALL Equn()             ! this is unfished equilibrium number at age per recruit
      CALL InitProj()         ! this is the unfished spawning stock number and the number at age per recruit
C

C	Input Files

!     OPEN(UNIT=2,FILE="RANDOM6.NUM")
      OPEN(UNIT=2,FILE="RANDOMTRY.NUM") !This is used - seeds for each sim come from this file
	OPEN(UNIT=3,FILE="RANDOM4.NUM")

C	Output Files
	OPEN(UNIT=39,FILE="results/SUMM18.OUT")

!     LastSuperSim = 1
!	DO 10001 SuperSim = 1,LastSuperSim

C     Open output files 
      OPEN(UNIT=14,FILE="results/SUMM2.OUT")
      OPEN(UNIT=15,FILE="results/SUMM3.OUT")
      OPEN(UNIT=16,FILE="results/SUMM4.OUT")
      OPEN(UNIT=17,FILE="results/SUMM5.OUT")
      OPEN(UNIT=18,FILE="results/SUMM6.OUT")
C     Carey: Add a file for analyzing movement characteristics
      OPEN(UNIT=19,FILE="results/SUMM7.OUT")
      OPEN(UNIT=26,FILE="results/SUMM8.OUT")
      OPEN(UNIT=27,FILE="results/SUMM9.OUT")
      OPEN(UNIT=28,FILE="results/SelAfor.OUT")
      OPEN(UNIT=29,FILE="results/SUMM10.OUT")
      OPEN(UNIT=32,FILE="results/SUMM11.OUT")
      OPEN(UNIT=33,FILE="results/SUMM12.OUT")
      OPEN(UNIT=34,FILE="results/SUMM13.OUT")
      OPEN(UNIT=35,FILE="results/SUMM14.OUT")
	OPEN(UNIT=38,FILE="results/SUMM17.OUT")
      OPEN(UNIT=22,FILE="results/SAVE.ALL")
	OPEN(UNIT=40,FILE="results/SUMM19.OUT")
	OPEN(UNIT=41,FILE="results/SUMM20.OUT")
	OPEN(UNIT=49,FILE="results/SUMM27.OUT")
      OPEN(UNIT = 1500,FILE = "results/CatchCPUE.OUT")
      OPEN(UNIT = 1600,FILE = "results/SurveyCPUE.OUT")
      OPEN(UNIT = 78,FILE = 'results/SaveConstantEff.All')
      OPEN(UNIT = 79,FILE = 'results/Move4SS.OUT')
      OPEN(UNIT = 260,FILE="results/RECRUITS.OUT")
    
      !Pstar-related files
      OPEN(UNIT = 400,FILE = "results/Pstar1.OUT") !FDistn,FOFL information by vector, configuration, year,simulation
      OPEN(UNIT = 401,FILE = "results/Pstar2.OUT") !EstF35,EstF40 and associated optimization info to check if FOFL-related optimizations worked. QAQC, mostly.
      OPEN(UNIT = 402,FILE = "results/Pstar3.OUT") !Flags simulation, year, and configuration that didn't converge (defined as no .cor file)
      OPEN(UNIT = 4000,FILE = "results/Pstar1_Old.OUT")
      OPEN(UNIT = 4001,FILE = "results/Pstar2_Old.OUT")
      OPEN(UNIT = 4002,FILE = "results/Pstar3_Old.OUT")

      
      WRITE(402,*) "The_following_assessments_didnt_produce_cor_or_MVNn"
      WRITE(402,*) "Isim, GlobalProjYr, Ispec, Iconfig, TestMVN"
      !Make character variables of historical assessment years to then string together file names containing info on assessment year
!      Write(HistChar1,'(I2)') HistAssYr(1)
!      Write(HistChar2,'(I2)') HistAssYr(2)
!      Write(HistChar3,'(I2)') HistAssYr(3)

!      OPEN(UNIT = 4989,File ='MainResults.Yr'//HistChar1//'Method1')
!      OPEN(UNIT = 4990,File ='MainResults.Yr'//HistChar1//'Method2')
!      OPEN(UNIT = 4991,File ='MainResults.Yr'//HistChar1//'Method3')
!      OPEN(UNIT = 4992,File ='MainResults.Yr'//HistChar1//'Method4')
!
!      OPEN(UNIT = 4993,File ='MainResults.Yr'//HistChar2//'Method1')
!      OPEN(UNIT = 4994,File ='MainResults.Yr'//HistChar2//'Method2')
!      OPEN(UNIT = 4995,File ='MainResults.Yr'//HistChar2//'Method3')
!      OPEN(UNIT = 4996,File ='MainResults.Yr'//HistChar2//'Method4')
!
!      OPEN(UNIT = 4997,File ='MainResults.Yr'//HistChar3//'Method1')
!      OPEN(UNIT = 4998,File ='MainResults.Yr'//HistChar3//'Method2')
!      OPEN(UNIT = 4999,File ='MainResults.Yr'//HistChar3//'Method3')
!      OPEN(UNIT = 5000,File ='MainResults.Yr'//HistChar3//'Method4')
!
!      OPEN(UNIT = 6992,File = 'SpatialResults.Yr'//HistChar1//'Method4')
!      OPEN(UNIT = 6996,File = 'SpatialResults.Yr'//HistChar2//'Method4')
!      OPEN(UNIT = 6000,File = 'SpatialResults.Yr'//HistChar3//'Method4')
!
!      OPEN(UNIT = 7992,File = 'SpatialResults.Yr'//HistChar1//'Method2')
!      OPEN(UNIT = 7996,File = 'SpatialResults.Yr'//HistChar2//'Method2')
!      OPEN(UNIT = 7000,File = 'SpatialResults.Yr'//HistChar3//'Method2')


	  OPEN(UNIT=51,FILE = "results/SUMM29.OUT")
	  OPEN(UNIT=52,FILE = "results/SUMM30.OUT")


C	Call a function here that within FUNK does the simulation NSim times, adds up the results from (1) a wide ranging profile and (2) the profile dialed in and then calls the FIT function to obtain the exact average best constant effort rate
C	The other option is to work within what's already built to do the same thing - use switches to turn off FindInitEffort and FindEffLev; set constant effort = 0 for these, make a new MinOpt in FUNK that is outside of isim and inside supersim = 1 and call the FIT function
C	FindMSYOnly = 1
	IF (DoMSYOnly.EQ.1) THEN
        !get Qval = from save.all
        READ(22,*) SEED1,SEED2,(QVAL(Ispec),ISpec=1,NSpec) !22 is save.all
	  CALL FindAvgMSY()
        STOP
	ENDIF

      !-----------------------------------------------------------------------------------------------------
C     Start simulations
      DO 10000 Isim = 1,Nsim 

      !Write out the simulation number and supersim each time
	WRITE(*,*)
	WRITE(*,*)
	WRITE(*,'(1x,A,I4,1x,A,I2)') 'Isim = ', Isim
	WRITE(*,*)
	WRITE(*,*)
C	Write(*,*) Isim



C      Reset the quotas and all areas to open
       Quota = 0
       Closed = 1


C      Read the seeds for this simulation 
       READ(2,*) SEED1,SEED2
       ISEEDP = SEED1 !ISEEDP = seed for Process error
       ISEEDO = SEED2 !ISEEDO = seed for Observation error
C
C      Generate the process errors for this simulation
!	 IF (SuperSim.EQ.1) THEN
         CALL GenProcErr() 
!	 ENDIF
C
!	QValOne = 0 !Added as an option to MPAMODEL.DAT
C	Carey used Qval = 1 for density ratio project
C==============================================
C      Find the catchability coefficients
       IF (USePrevVals.NE.1.AND.QValOne.NE.1) THEN
          CALL FindQvals()
          WRITE(22,'(2(I10,1x),F20.15,1x,F15.10)') SEED1,SEED2,
     +      (QVAL(Ispec),ISpec=1,NSpec),
     +      SpBioTot(1,BaseYr+1)/SpBioTot(1,MinPYr)*100
	    WRITE(39,'(1x,A)') 'Qval(Ispec),ISpec=1,NSpec)'
	    WRITE(39,*) (QVAL(Ispec),Ispec=1,NSpec)
	    WRITE(*,'(1x,A)') 'Qval(1) = '
	    WRITE(*,*) QVAL(1)
       ELSEIF (UsePrevVals.EQ.1.AND.QValOne.NE.1) THEN
          ! DO 44556 Ispec = 1,Nspec
          !   Qval(Ispec) = Qvalseed(Ispec)
!44556     ! CONTINUE
          READ(22,*) SEED1,SEED2,(QVAL(Ispec),ISpec=1,NSpec)
          DO 5151 Ispec = 1,Nspec
          !  X(Ispec) = Qval(Ispec)*NNSBlock*NEWBlock
             X(Ispec) = QVAL(Ispec)
5151      CONTINUE
          MinOpt = 1
          SS = FUNK(X)
       ELSEIF (QValOne.EQ.1) THEN
         !	 Carey is changing this one part so that Qval is 1
         !===================================
	   DO 66667 Ispec = 1,Nspec
            QVal(Ispec) = 1
C	      X(Ispec) = QVal(Ispec)
	
	      WRITE(*,*)
	      WRITE(*,*)
	      WRITE(*,*)
	      WRITE(*,'(1x,A)') 'Finding initial fishing effort 
     +	(Year FirstFishYr to BaseYr) to achieve InitDep'
	      CALL FindInitEffort()
66667	  CONTINUE
        !==================================
        !READ(22,*) SEED1,SEED2,(X(Ispec),ISpec=1,NSpec)
          ISEEDP = SEED1
          ISEEDO = SEED2
 
          MinOpt = 8
          SS = FUNKSCALAR(BestInitEffort)
      ENDIF
!=========================================================

C	Carey sets QVal to 1 for density ratio project

C
C      Find the effort levels
       IF (SolFinDepl.EQ.1) THEN
	WRITE(*,*)
	WRITE(*,*)
	WRITE(*,*)
	WRITE(*,'(1x,A)') 'Finding Effort level for BaseYr to 
     +	FinDepYr to achieve FinDep'
        !itarg = 1
!     Run the model and find effort levels such that the model is at particular depletion levels at particular times
!     Allows for the program to find multiple constant effort levels for different time periods, each ending in a particular depletion value
!     NTargetDepl is the number of target depletion levels that are read in from MPAMODEL.DAT
        DO 66695 itarg = 1,NTargetDepl
          IF (itarg.EQ.1) THEN
           CALL GetEffLev(BaseYr+1,DeplYr(itarg),TargetDeplVec(itarg,:))
          !Write out ConstantEff levels to ConstantEff.All
         WRITE(78,'(F26.15,1x,3(I3,1x),10(F10.6,1x))') ConstantEff,Isim,
     +     BaseYr+1,DeplYr(itarg),FinDepRecord,
     +       (TargetDeplVec(itarg,Ispec),Ispec=1,NSpec)
          ELSE
           mystyr = DeplYr(itarg-1)+1
           myendyr = DeplYr(itarg)
           mytarg = TargetDeplVec(itarg,:)
           CALL GetEffLev(mystyr,myendyr,mytarg)
           !Write out ConstantEff levels to ConstantEff.All
           WRITE(78,'(F26.15,1x,3(I3,1x),10(F10.6,1x))') ConstantEff,
     +    Isim,mystyr,myendyr,FinDepRecord,(mytarg(Ispec),Ispec=1,Nspec)
          ENDIF
!         Carey changed this from FUNKSCALAR to FUNK on 7/12/11
!	    SS = FUNKSCALAR(ConstantEff) !Here I think you need FUNKSCALAR: X doesn't mean anything here.
          SS = FUNK(ConstantEffVec)
          


66695  CONTINUE
       ENDIF 

       !Use previous values for the effort level that produces the depletion level that's desired.
       IF (SolFinDepl.NE.1) THEN     
       DO 66696 itarg = 1,NTargetDepl
           READ(78,*) ConstantEff,TheSim,StartYr,EndYr,DeplFound,
     +        (TheTarg(Ispec),Ispec=1,Nspec)
           ConstantEffVec(1) = ConstantEff
           MinOpt = 3
          DO 54991 Ispec = 1,Nspec
            FinalDep(Ispec) = TheTarg(Ispec) !change this for multiple species!
54991     CONTINUE 
           SS = FUNK(ConstantEffVec)
66696    CONTINUE
       ENDIF


  
C----------------------------------------------------------------------------------------
C----------------------------------------------------------------------------------------      
C      Generate the historical data: Loop from MinPYr to FinDepYr
C      ----------------------------
       Closed = 1  !Reset all of the areas to open
C       DO 11000 Yr = -50,FinDepYr
       DO 11000 Yr = MinPYr,FinDepYr
C
C       Impose MPA
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C	See if changing the following will work:
        IF (Yr.EQ.MPAYear) THEN
         DO 11100 Iarea = 1,NClosedAreas
          DO 11110 IblkNS = MinCNS(Iarea),MaxCNS(Iarea)
           DO 11110 IBlkEW = MinCEW(Iarea),MaxCEW(Iarea)
            Closed(IBlkNS,IBlkEW) = 0
11110     CONTINUE
11100    CONTINUE          
        ENDIF
C        
        CALL GenData(Yr)
        CALL Summary2(Yr) 
C
C       Check whether additional runs of the estimator
        DO 11200 Iass = 1,NhistAss
C
C        Conduct stock assessments
         IF (HistAssYr(Iass).EQ.Yr) THEN !6/4/11 start here on stepping through code

!---------------------------------------------------------------------
C          WRITE(*,*) Iass,Yr
          IF (OMonly.EQ.0) THEN
            CALL DoEstAll(Yr)
          ENDIF
         ENDIF !histassyr.EQ.yr
C        
11200   CONTINUE   !end loop over number of assessments     
C
11000  CONTINUE  !End loop over time up to FinDepYr.




C      Project forward
       DO 12000 TheProjYear = FinDepYr+1,FinDepYr+NprojYr
        WRITE(*,*) TheProjYear
        GlobalProjYr = TheProjYear !globally defines which year of the MSE loop is the current year; trickery for passing the ProjYr to FUNK and FIT and subroutines used within FUNK and FIT.
C       Impose MPA
        IF (Yr.EQ.MPAYear) THEN
         DO 12100 Iarea = 1,NClosedAreas
          DO 12110 IblkNS = MinCNS(Iarea),MaxCNS(Iarea)
           DO 12110 IBlkEW = MinCEW(Iarea),MaxCEW(Iarea)
            Closed(IBlkNS,IBlkEW) = 0
12110     CONTINUE
12100    CONTINUE          
        ENDIF
C
C       Conduct stock assessments
        !DoEstAll is effectively happening in the previous year (data aren't generated for TheProjYr yet)
        CALL DoEstAll(TheProjYear) !includes quota calculations. (change this later to make things more transparent?)
        CALL AnnualProj()
        CALL GenData(TheProjYear) !TheProjYear is now the current year
        CALL Summary2(TheProjYear)
        CALL SummaryPstar(TheProjYear)
        CALL SummaryPstar_Old(TheProjYear)

12000  CONTINUE 
   
10000 CONTINUE !Ends the Isim = 1,Nsim main loop

C==========================================================

C==================================================================

	CLOSE(14)
      CLOSE(15)
      CLOSE(16)
      CLOSE(17)
      CLOSE(18)
	CLOSE(19)
      CLOSE(22)
	CLOSE(26)
	CLOSE(27)
	CLOSE(28)
      CLOSE(29) 
	CLOSE(32)
	CLOSE(33)
	CLOSE(34)
	CLOSE(35)
	CLOSE(38)
	CLOSE(40)
	CLOSE(41)
      CLOSE(78)
      CLOSE(79)
      CLOSE(8000)
      CLOSE(8001)
      CLOSE(260)
C This CONTINUE refers to SuperSim = 1,2, which does not optimize, but performs Nsim simulations for an optimized control rule 
!10001 CONTINUE



	CLOSE(2)
	CLOSE(24)
	CLOSE(25)
	CLOSE(39)
	CLOSE(3)

      !Close Pstar files:
      CLOSE(400) !Pstar1.out
      CLOSE(401) !Pstar2.out
      CLOSE(402) !Pstar3.out
      CLOSE(4000)
      CLOSE(4001)
      CLOSE(4002)

      total = etime(elapsed)
	WRITE(*,*)
	WRITE(*,*)
      WRITE(*,'(1x,A)') 'The end total= user =  system ='
	WRITE(*,*) total, elapsed(1), elapsed(2)
	WRITE(*,*)
	WRITE(*,*)
C      


      CALL mydate('End time and date        ')
      CLOSE(2001)

      STOP
      END PROGRAM MPAMODEL
C


C      
      INCLUDE 'COMMON.FOR'
                              

C ===========================================================================
      SUBROUTINE DoEstAll(TheYr)
C
C     This subroutine handles calling the estimator 
C
	USE ShareDataRatioSubs
      IMPLICIT NONE
C
C     Global variables
      INTEGER TheYr 
C
C     Local variables 
      INTEGER Ispec,Iarea,II,Yr
C
C     Clear summary array
      PassEstVar = 0

C
C     Conduct stock assessments
       DO 10000 Ispec = 1,Nspec
        QuotaAct(Ispec) = 0

!Run assessments (either in cab or in ss3)
!---------------------------------------------------------------------------

          CALL DOCAB2(Ispec,0,0,MinPYr,TheYr-1,0,"A")

!--------------------------------------------------------------------------
         CLOSE(957)
10000   CONTINUE 
      

      
801   FORMAT(1x,I3,1x,2(F9.1,1x),2(F7.4,1x),2(F9.1,1x),2(F9.1,1x),
     +      2(F7.4,1x))
      RETURN
      END SUBROUTINE DoEstAll

C
C ===========================================================================
C    
	  
      SUBROUTINE ReadSpec()
C
C This subroutine reads in the specifications for the run
C      
	USE ShareDataRatioSubs
      IMPLICIT NONE
C
C     Local variables
      INTEGER Ispec,Isex,Itime,Iage,IblkNS,IblkEW,JblkNS,JblkEW,itarg
      INTEGER Yr,Iyr,Iarea,II,i
C took out Linf variables
      REAL*8 Vtime,Weight,Temp,TempTot,WeightBeg

      REAL(8),DIMENSION(MaxSpec,0:Maxage) :: Length  !Mid-year
      REAL(8),DIMENSION(MaxSpec,0:Maxage) :: LengthBeg  !Beginning of the year
      REAL*8 CentrNS,CentrEW
      REAL*8 TmpY1,TmpY2,TmpX11,TmpX12,TmpX21,TmpX22
      REAL*8 Apar1,Apar2,Bpar1,Bpar2,Cpar1,Cpar2
      REAL*8 Cent2NL,Cent2EL,ValNF,ValEF,ValNL,ValEL
C
C     Read in the specifications
      OPEN(UNIT=13,FILE='OM.DAT')
      READ(13,'(40X,I4)') Nsim 
      READ(13,'(40X,I4)') OMonly
      READ(13,'(40X,I4)') UsePrevVals
      READ(13,'(40X,I4)') NprojYr
      READ(13,'(40X,I4)') BaseYr
      READ(13,'(40X,I4)') FinDepYr
      READ(13,'(40X,I4)') FirstFishYr
      READ(13,'(40X,I4)') SolFinDepl
	READ(13,'(40X,F25.10)') FinDeplEffort
	READ(13,'(40X,I4)') MPAYear
C	Carey commented this out!!
C      IF (SolFinDepl.NE.1) FinDepYr = BaseYr
      READ(13,'(40X,I2)') Nspec
      READ(13,'(40X,I2)') Ntime
      READ(13,'(40X,I3)') NNSBlock,NEWBlock
      Rtime = 1.0/FLOAT(Ntime)

	TotYrs=NProjYr+FinDepYr + 1 + (-MinPYr)
      
      READ(13,'(40X,5(F10.2,1x))') (Kvec(Ispec),Ispec=1,Nspec)
      READ(13,'(40X,5(F4.2,1x))') (M(Ispec),Ispec=1,Nspec)
      READ(13,'(40X,5(I4,1x))') (Amax(Ispec),Ispec=1,Nspec)
      READ(13,'(40X,I4)') PlusGroup
      READ(13,'(40X,5(F6.3,1x))') (InitDep(Ispec),Ispec=1,Nspec)
      READ(13,'(40X,5(F6.3,1x))') (FinalDep(Ispec),Ispec=1,Nspec)
      READ(13,'(40X,I4)') QvalOne
      READ(13,'(40X,5(F5.2,1x))') (Qvalseed(Ispec),Ispec=1,Nspec)
      READ(13,*)
      READ(13,'(40X,I4)') DoProcessError
      READ(13,'(40X,5(F8.5,1x))') (SigmaR1(Ispec),Ispec=1,Nspec)
      READ(13,'(40X,5(F8.5,1x))') (SigmaR2(Ispec),Ispec=1,Nspec)
      READ(13,'(40X,5(F4.3,1x))') (Steep(ISpec),Ispec=1,Nspec)
      READ(13,'(40X,5(  I4,1x))') (SrType(ISpec),Ispec=1,Nspec)
	READ(13,'(40X,5(  I4,1x))') (RecruitOpt(ISpec),Ispec=1,Nspec)
      READ(13,'(40X,5(F7.3,1x))') (RecX(Ispec),Ispec=1,Nspec)
      READ(13,'(40X,5(F7.3,1x))') (RecY(Ispec),Ispec=1,Nspec)
      READ(13,'(40X,5(F7.3,1x))') (VRecXL(Ispec),Ispec=1,Nspec)
      READ(13,'(40X,5(F7.3,1x))') (VRecYL(Ispec),Ispec=1,Nspec)
      READ(13,'(40X,5(F7.3,1x))') (VRecXH(Ispec),Ispec=1,Nspec)
      READ(13,'(40X,5(F7.3,1x))') (VRecYH(Ispec),Ispec=1,Nspec)
      READ(13,*)
      READ(13,'(40X,5(F5.1,1x))') (Linf(Ispec,1),Ispec=1,Nspec)
      READ(13,'(40X,5(F5.1,1x))') (Linf(Ispec,2),Ispec=1,Nspec)
      READ(13,'(40X,5(F5.3,1x))') (kpar(Ispec,1),Ispec=1,Nspec)
      READ(13,'(40X,5(F5.3,1x))') (kpar(Ispec,2),Ispec=1,Nspec)
      READ(13,'(40X,5(F5.2,1x))') (t0(Ispec,1),Ispec=1,Nspec)
      READ(13,'(40X,5(F5.2,1x))') (t0(Ispec,2),Ispec=1,Nspec)
      READ(13,'(40X,5(F5.3,1x))') (Apar(Ispec,1),Ispec=1,Nspec)
      READ(13,'(40X,5(F5.3,1x))') (Apar(Ispec,2),Ispec=1,Nspec)
      READ(13,'(40X,5(F5.3,1x))') (Bpar(Ispec,1),Ispec=1,Nspec)
      READ(13,'(40X,5(F5.3,1x))') (Bpar(Ispec,2),Ispec=1,Nspec)
	READ(13,'(40X,5(  I4,1x))') (LengthBinIntervals(ISpec),Ispec=1,Nspec)
	READ(13,'(40X,5(  I4,1x))') (NLengthBins(ISpec),Ispec=1,Nspec)
      READ(13,*)
C     Fecundity 
      READ(13,'(40X,5(F5.2,1x))') (L50(Ispec),Ispec=1,Nspec)
      READ(13,'(40X,5(F5.2,1x))') (L95(Ispec),Ispec=1,Nspec)
      READ(13,*)
      READ(13,'(40X,5(F7.3,1x))') (MixN(Ispec),Ispec=1,Nspec)
      READ(13,'(40X,5(F7.3,1x))') (MixE(Ispec),Ispec=1,Nspec)
      READ(13,'(40X,15(F7.3,1x))') (X1N(Ispec),V1N(Ispec),Ispec=1,Nspec)
      READ(13,'(40X,15(F7.3,1x))') (X2N(Ispec),V2N(Ispec),Ispec=1,Nspec)
      READ(13,'(40X,15(F7.3,1x))') (X1E(Ispec),V1E(Ispec),Ispec=1,Nspec)
      READ(13,'(40X,15(F7.3,1x))') (X2E(Ispec),V2E(Ispec),Ispec=1,Nspec)
      READ(13,*)
C     Selectivity for the fishery 
      READ(13,'(40X,5(F7.3,1x))') (SCL50(ISpec),Ispec=1,Nspec)
      READ(13,'(40X,5(F7.3,1x))') (SCL95(ISpec),Ispec=1,Nspec)
      READ(13,*)
C     Selectivity for the survey 
      READ(13,'(40X,5(F7.3,1x))') (SSL50(ISpec),Ispec=1,Nspec)
      READ(13,'(40X,5(F7.3,1x))') (SSL95(ISpec),Ispec=1,Nspec)
      READ(13,*)
      READ(13,'(40X,5(F8.5,1x))') (SigmaSpat(ISpec),Ispec=1,Nspec)
      READ(13,'(40X,I1)') Hyperstability
!      READ(13,*)
!      READ(13,'(40X,F15.10)') TotEffortSeed
!      DO 11000 Yr= -18,BaseYr  !Why was this -18 and all other refs are for -19?
!      DO 11000 Yr= FirstFishYr,BaseYr
!       TotEffort(Yr) = TotEffortSeed
!11000 CONTINUE  
      READ(13,*)
      READ(13,*)     
      DO 11000 Yr= MinPyr,BaseYr
       READ(13,*) Iyr,TotEffort(Yr)
11000 CONTINUE 
      
      DO 11112 Yr = MinPyr+1,BaseYr
         IF (TotEffort(Yr-1).LE.0.0.AND.TotEffort(Yr).
     +   GT.0.0) THEN
           StartFishing = Yr
         ENDIF
11112 CONTINUE


C     Specifications for closed areas
      READ(13,*)
      READ(13,*)
      READ(13,'(40X,I2)') NClosedAreas
      READ(13,'(40X,100(I2,1x))') (MinCNS(Iarea),Iarea=1,NClosedAreas)
      READ(13,'(40X,100(I2,1x))') (MaxCNS(Iarea),Iarea=1,NClosedAreas)
      READ(13,'(40X,100(I2,1x))') (MinCEW(Iarea),Iarea=1,NClosedAreas)
      READ(13,'(40X,100(I2,1x))') (MaxCEW(Iarea),Iarea=1,NClosedAreas)

C	Specifications for open areas
      READ(13,*)
      READ(13,*)	
      READ(13,'(40X,100(I2,1x))') (MinOstrataNS(Iarea),
     +	Iarea=1,NClosedAreas+2)
      READ(13,'(40X,100(I2,1x))') (MaxOstrataNS(Iarea),
     +	Iarea=1,NClosedAreas+2)	
      READ(13,'(40X,100(I2,1x))') (MinOstrataEW(Iarea),
     +	Iarea=1,NClosedAreas+2)
      READ(13,'(40X,100(I2,1x))') (MaxOstrataEW(Iarea),
     +	Iarea=1,NClosedAreas+2)

C     Over-ride input (0th block is always the whole area) 
      MinENS(0) = 1
      MaxENS(0) = NNSBlock
      MinEEW(0) = 1
      MAXEEW(0) = NEWBlock
    
C     Data-related inputs
      READ(13,*)
      READ(13,*)
      READ(13,'(40X,I1)') IsDeterm
      READ(13,'(40X,F4.2)') CpueCV
      READ(13,'(40X,I5)') EffectSampCatch1
      READ(13,'(40X,I5)') EffectSampCatch2
      READ(13,'(40X,I4)') FirstSurvYr
      READ(13,'(40X,I3)') FreqSurv
      READ(13,'(40X,F4.2)') SurvCV
      READ(13,'(40X,I5)') EffectSampSurvey1
      READ(13,'(40X,I5)') EffectSampSurvey2
      READ(13,'(40X,I5)') NoSurveyClosed
      READ(13,*)
      READ(13,'(40X,I5)') NHistAss
      READ(13,'(40X,100(I4,1x))') (HistAssYr(II),II=1,NHistAss)
      READ(13,*)
      READ(13,'(40X,I5)') ImplOpt
	READ(13,*)
      READ(13,*)
      READ(13,'(40X,15(F7.3,1x))') (MixNLarval(Ispec),Ispec=1,Nspec)
      READ(13,'(40X,15(F7.3,1x))') (MixELarval(Ispec),Ispec=1,Nspec)
      READ(13,'(40X,15(F7.3,1x))') (X1NLarval(Ispec),V1NLarval(Ispec),
     +	Ispec=1,Nspec)
      READ(13,'(40X,15(F7.3,1x))') (X2NLarval(Ispec),V2NLarval(Ispec),
     +	Ispec=1,Nspec)
      READ(13,'(40X,15(F7.3,1x))') (X1ELarval(Ispec),V1ELarval(Ispec),
     +	Ispec=1,Nspec)
      READ(13,'(40X,15(F7.3,1x))') (X2ELarval(Ispec),V2ELarval(Ispec),
     +	Ispec=1,Nspec)
      READ(13,*)
      READ(13,'(40X,I5)') NTargetDepl
      READ(13,*)
      DO 3459 itarg = 1,NTargetDepl
        READ(13,'(I4,5(1x,F5.2))') DeplYr(itarg), 
     +    (TargetDeplVec(itarg,Ispec),Ispec = 1,Nspec)
3459  CONTINUE
      CLOSE(13)
C
!     Calculate length bins
      DO 101 Ispec = 1,Nspec
      LengthBins(Ispec,0) = 0
      DO 101 i=1,NLengthBins(Ispec)
        LengthBins(Ispec,i) = LengthBins(Ispec,i-1) + 
     +     LengthBinIntervals(Ispec)
101   CONTINUE    
  
        !Find params for SS von-bert parameterization: L(amin) instead of t0
      DO 4357 Ispec = 1,Nspec
       DO 4357 Isex = 1,Nsex
        Lamin(Ispec,Isex) = MAX(Linf(Ispec,Isex)*(1-exp(-kpar(Ispec,
     +   Isex)* + (0.0-t0(Ispec,Isex)))),0.001)
4357  CONTINUE 

C     Calculate weight-at-age
      DO 80000 Ispec = 1,Nspec
       DO 80000 Isex = 1,NSex
        DO 81000 Itime = 0,Ntime
C
C        Adjustment factor for the time - makes it mid-year calculations.
         IF (Itime.EQ.0) THEN 
          Vtime = 0
         ELSE
          Vtime = (Itime-0.5)*Rtime
         ENDIF  
C

C        Now calculate weight-at-age
         DO 81100 Iage = 0,Amax(Ispec)   
    !-------------------------------------------------------------------
    ! This was the regular von-bert way; changed from this to match SS3                    
    !      Length(Ispec,Iage) = Linf(Ispec,Isex)*
    ! +      (1-exp(-Kpar(Ispec,ISex)*(Iage+Vtime-T0(ISpec,Isex)))) !Mid-year
    !      LengthBeg(Ispec,Iage) = Max(Linf(Ispec,Isex)*
    ! +      (1-exp(-Kpar(Ispec,ISex)*(Iage-T0(ISpec,Isex)))),0.001) !Beginning of year
          !-------------------------------------------------------
    ! SS3 way of calculating length at age:
          Length(Ispec,Iage) = Linf(Ispec,Isex) + (Lamin(Ispec,Isex)-
     + Linf(Ispec,Isex))*exp(-Kpar(Ispec,ISex)*(Iage+Vtime-0.0)) !Mid-year
         LengthBeg(Ispec,Iage) = Linf(Ispec,Isex) + (Lamin(Ispec,Isex)-
     + Linf(Ispec,Isex))*exp(-Kpar(Ispec,ISex)*(Iage-0.0)) !Beginning of year
    !---------------------------------------------------------------------------
          !Bin lengths
           DO 87100 i = 1,NLengthBins(Ispec)
             IF (Length(Ispec,Iage).GT.LengthBins(Ispec,i-1).AND.
     +Length(Ispec,Iage).LE.LengthBins(Ispec,i)) THEN
               BinnedLength(Ispec,Iage) = (LengthBins(Ispec,i)+
     +          LengthBins(Ispec,i-1))/LengthBinIntervals(Ispec)
            ENDIF
             IF (LengthBeg(Ispec,Iage).GT.LengthBins(Ispec,i-1).AND.
     +LengthBeg(Ispec,Iage).LE.LengthBins(Ispec,i)) THEN
               BinnedLengthBeg(Ispec,Iage) = (LengthBins(Ispec,i)+
     +          LengthBins(Ispec,i-1))/LengthBinIntervals(Ispec)
            ENDIF
87100      CONTINUE
           !Done binning lengths
           !-----------------------------------------------------
           !Making an adjustment here to use BinnedLength instead of Length to calculate weight for SS3
!          Weight = Apar(Ispec,Isex)*(Length(Ispec,Iage)**Bpar
!     +      (ISpec,Isex))
          Weight = Apar(Ispec,Isex)*(BinnedLength(Ispec,Iage)**Bpar
     +      (ISpec,Isex))
          WeightBeg =Apar(Ispec,Isex)*(BinnedLengthBeg(Ispec,Iage)**Bpar
     +      (ISpec,Isex))
           !-------------------------------------------------------
          Wght(Ispec,Isex,iage,Itime) = Weight
          WghtBeg(Ispec,Isex,iage,Itime) = WeightBeg
81100    CONTINUE        
81000   CONTINUE        
80000 CONTINUE      
C
C     Calculate fecundity
      DO 90000 Ispec = 1,Nspec          
       DO 91000 Iage = 1,Amax(ISpec)
        Temp = (Iage-L50(Ispec))/(L95(Ispec)-L50(ISpec))
        Mature(Ispec,Iage) = 1.0/(1.0+EXP(-Log(19.0)*Temp))
        Fec(Ispec,Iage) = Mature(Ispec,Iage) * WghtBeg(Ispec,1,Iage,0) !This changed from mid-year weight (Wght) to beg. year weight (WghtBeg) for SS3
91000  CONTINUE
90000 CONTINUE       
C
C     Calculate the mixing matrices
      DO 92000 Ispec = 1,Nspec
       TmpY1 = V1N(Ispec) - 1
       TmpY2 = V2N(Ispec) - 1
       TmpX11 = X1N(Ispec)**2.0-0.25
       TmpX12 = X1N(Ispec)-0.5
       TmpX21 = X2N(Ispec)**2.0-0.25
       TmpX22 = X2N(Ispec)-0.5
       Bpar1=(TmpY2-TmpY1*TmpX21/TmpX11)/(TmpX22-TmpX12*TmpX21/TmpX11)
       Apar1 = (TmpY1-Bpar1*TmpX12)/TmpX11
       Cpar1=1-0.25*Apar1-0.5*Bpar1
       TmpY1 = V1E(Ispec) - 1
       TmpY2 = V2E(Ispec) - 1
       TmpX11 = X1E(Ispec)**2.0-0.25
       TmpX12 = X1E(Ispec)-0.5
       TmpX21 = X2E(Ispec)**2.0-0.25
       TmpX22 = X2E(Ispec)-0.5
       Bpar2=(TmpY2-TmpY1*TmpX21/TmpX11)/(TmpX22-TmpX12*TmpX21/TmpX11)
       Apar2 = (TmpY1-Bpar2*TmpX12)/TmpX11
       Cpar2=1-0.25*Apar2-0.5*Bpar2
C       WRITE(*,*) Apar1,Bpar1,Cpar1
C       WRITE(*,*) Apar2,Bpar2,Cpar2
       DO 92100 IBlkNS = 1,NNSBlock
        DO 92100 IBlkEW = 1,NEWBlock
         CentrNS = IblkNS - 0.5
         CentrEW = IblkEW - 0.5
         ValNF = MAX(0.01,Apar1*CentrNS**2.0+Bpar1*CentrNS+Cpar1)
         ValEF = MAX(0.01,Apar2*CentrEW**2.0+Bpar2*CentrEW+Cpar2)
C         WRITE(*,*) "A",ISpec,IBlkNS,IBlkEW,ValNF,ValEF
         TempTot = 0
         DO 92110 JBlkNS = 1,NNSBlock
          DO 92110 JBlkEW = 1,NEWBlock
           Cent2NL = JBlkNS-0.5
           Cent2EL = JBlkEW-0.5
           ValNL = MAX(0.01,Apar1*Cent2NL**2.0+Bpar1*Cent2NL+Cpar1)
           ValEL = MAX(0.01,Apar2*Cent2EL**2.0+Bpar2*Cent2EL+Cpar2)
           Temp = (CentrNS-Cent2NL)**2.0/(2.0*MixN(ISpec)**2.0) + 
     +            (CentrEW-Cent2EL)**2.0/(2.0*MixE(Ispec)**2.0)
           Temp = Exp(-Temp) * (ValNL/ValNF) * (ValEL/ValEF)
           MixMat(Ispec,IBlkNS,IBlkEW,JBlkNS,JBlkEW) = Temp
C           IF (ISpec.EQ.2) 
C     +       WRITE(*,*) JBlkNS,JBlkEW,ValNL/ValNF,ValEL/ValEF
           TempTot = TempTot + Temp
92110    CONTINUE          
         DO 92120 JBlkNS = 1,NNSBlock
          DO 92120 JBlkEW = 1,NEWBlock
           MixMat(Ispec,IBlkNS,IBlkEW,JBlkNS,JBlkEW) =
     +        MixMat(Ispec,IBlkNS,IBlkEW,JBlkNS,JBlkEW) / TempTot
92120    CONTINUE          
C
92100  CONTINUE        
92000 CONTINUE  

C     Calculate the larval mixing matrices
      DO 92002 Ispec = 1,Nspec
       TmpY1 = V1NLarval(Ispec) - 1
       TmpY2 = V2NLarval(Ispec) - 1
       TmpX11 = X1NLarval(Ispec)**2.0-0.25
       TmpX12 = X1NLarval(Ispec)-0.5
       TmpX21 = X2NLarval(Ispec)**2.0-0.25
       TmpX22 = X2NLarval(Ispec)-0.5
       Bpar1=(TmpY2-TmpY1*TmpX21/TmpX11)/(TmpX22-TmpX12*TmpX21/TmpX11)
       Apar1 = (TmpY1-Bpar1*TmpX12)/TmpX11
       Cpar1=1-0.25*Apar1-0.5*Bpar1
       TmpY1 = V1ELarval(Ispec) - 1
       TmpY2 = V2ELarval(Ispec) - 1
       TmpX11 = X1ELarval(Ispec)**2.0-0.25
       TmpX12 = X1ELarval(Ispec)-0.5
       TmpX21 = X2ELarval(Ispec)**2.0-0.25
       TmpX22 = X2ELarval(Ispec)-0.5
       Bpar2=(TmpY2-TmpY1*TmpX21/TmpX11)/(TmpX22-TmpX12*TmpX21/TmpX11)
       Apar2 = (TmpY1-Bpar2*TmpX12)/TmpX11
       Cpar2=1-0.25*Apar2-0.5*Bpar2
C       WRITE(*,*) Apar1,Bpar1,Cpar1
C       WRITE(*,*) Apar2,Bpar2,Cpar2
       DO 92102 IBlkNS = 1,NNSBlock
        DO 92102 IBlkEW = 1,NEWBlock
         CentrNS = IblkNS - 0.5
         CentrEW = IblkEW - 0.5
         ValNF = MAX(0.01,Apar1*CentrNS**2.0+Bpar1*CentrNS+Cpar1)
         ValEF = MAX(0.01,Apar2*CentrEW**2.0+Bpar2*CentrEW+Cpar2)
C         WRITE(*,*) "A",ISpec,IBlkNS,IBlkEW,ValNF,ValEF
         TempTot = 0
         DO 92112 JBlkNS = 1,NNSBlock
          DO 92112 JBlkEW = 1,NEWBlock
           Cent2NL = JBlkNS-0.5
           Cent2EL = JBlkEW-0.5
           ValNL = MAX(0.01,Apar1*Cent2NL**2.0+Bpar1*Cent2NL+Cpar1)
           ValEL = MAX(0.01,Apar2*Cent2EL**2.0+Bpar2*Cent2EL+Cpar2)
           Temp = (CentrNS-Cent2NL)**2.0/(2.0*MixNLarval(ISpec)**2.0) + 
     +            (CentrEW-Cent2EL)**2.0/(2.0*MixELarval(Ispec)**2.0)
           Temp = Exp(-Temp) * (ValNL/ValNF) * (ValEL/ValEF)
           MixMatLarval(Ispec,IBlkNS,IBlkEW,JBlkNS,JBlkEW) = Temp
C           IF (ISpec.EQ.2) 
C     +       WRITE(*,*) JBlkNS,JBlkEW,ValNL/ValNF,ValEL/ValEF
           TempTot = TempTot + Temp
92112    CONTINUE          
         DO 92122 JBlkNS = 1,NNSBlock
          DO 92122 JBlkEW = 1,NEWBlock
           MixMatLarval(Ispec,IBlkNS,IBlkEW,JBlkNS,JBlkEW) =
     +        MixMatLarval(Ispec,IBlkNS,IBlkEW,JBlkNS,JBlkEW) / TempTot
92122    CONTINUE          
C
92102  CONTINUE        
92002 CONTINUE       

     
C
C     Set up selectivity (commercial)
      DO 93000 Ispec = 1,Nspec
       DO 93100 Iage = 0,Amax(Ispec)
        Temp=-log(19.0)*(Iage-SCL50(Ispec))/(SCL95(Ispec)-SCL50(ISpec))
        Sel(Ispec,Iage) = 1.0/(1.0+exp(Temp))
93100  CONTINUE        
93000 CONTINUE      
C
C     Set up selectivity (survey)
      DO 94000 Ispec = 1,Nspec
       DO 94100 Iage = 0,Amax(Ispec)
        Temp=-log(19.0)*(Iage-SSL50(Ispec))/(SSL95(Ispec)-SSL50(ISpec))
        SSel(Ispec,Iage) = 1.0/(1.0+exp(Temp))
94100  CONTINUE        
94000 CONTINUE  

   

      RETURN
      END SUBROUTINE ReadSpec

C
C ===========================================================================
C    
	  
      SUBROUTINE ReadCntl()
C
C This subroutine reads in the specifications of the estimator
C      
	USE ShareDataRatioSubs
      IMPLICIT NONE
C
C     Local variables
      INTEGER Iarea,Ispec,imove,jmove
C
C     Read in the area specifications
      OPEN(UNIT=13,FILE='EM.DAT')
      READ(13,'(40X,I2)') AssOpt
      READ(13,'(40X,I2)') Narea
      READ(13,'(40X,100(I2,1x))') (MinNS(Iarea),Iarea=1,Narea)
      READ(13,'(40X,100(I2,1x))') (MaxNS(Iarea),Iarea=1,Narea)
      READ(13,'(40X,100(I2,1x))') (MinEW(Iarea),Iarea=1,Narea)
      READ(13,'(40X,100(I2,1x))') (MaxEW(Iarea),Iarea=1,Narea)
C
C     Adjust the effort areas to match the options
      IF (AssOpt.EQ.-1.OR.AssOpt.EQ.0.OR.AssOpt.EQ.1.OR.AssOpt.EQ.2.
     +	OR.AssOpt.EQ.3.OR.AssOpt.EQ.4) THEN
       NEffortAreas = 1
       DO 21000 Iarea = 1,NEffortAreas
        MinENS(Iarea) = 1
        MaxENS(Iarea) = NNSBlock
        MinEEW(Iarea) = 1
        MaxEEW(Iarea) = NEWBlock
21000  CONTINUE        
      ELSE
       NEffortAreas = Narea
       DO 22000 Iarea = 1,NEffortAreas
        MinENS(Iarea) = MinNS(Iarea)
        MaxENS(Iarea) = MaxNS(Iarea)
        MinEEW(Iarea) = MinEW(Iarea)
        MaxEEW(Iarea) = MAXEW(Iarea)
22000  CONTINUE        
      ENDIF
C
C     Read in the biological parameters
      DO 10000 Ispec = 1,Nspec
       IestIgR(Ispec) = 0
10000 CONTINUE       
      READ(13,*)
      READ(13,'(40X,10(I9,1x))') (AmaxEst(Ispec),Ispec=1,Nspec)
      READ(13,'(40X,10(F9.2,1x))') (EstM(Ispec,0),Ispec=1,Nspec)
      READ(13,'(40X,10(F9.2,1x))') (IestSigmaR(Ispec),Ispec=1,Nspec)
      READ(13,'(40X,10(F9.2,1x))') (IestSteep(Ispec,0),ISpec=1,Nspec)
      READ(13,'(40X,10(F9.6,1x))') (Iestln_S0(Ispec),Ispec=1,Nspec)
      READ(13,*)
      READ(13,'(40X,I2,1x,F9.4)') PhaseM,BoundM
      READ(13,'(40X,I2,1x,F9.4)') PhaseS0,BoundS0
      READ(13,'(40X,I2,1x,F9.4)') PhaseH,BoundH 
      READ(13,'(40X,I2)') PhaseRecDevs
      READ(13,'(40X,I2,2(1x,F5.1))') PhaseSigmaR,loBoundSigmaR,
     +hiBoundSigmaR
C------------------------------------------------------------------
C     Selectivity options and parameters
      READ(13,'(40X,I2,2(1x,F5.3))') PhaseSel,BoundSela50,BoundSelSlope
      READ(13,*)
      READ(13,*)
C------------------------------------------------------------
   
      READ(13,*)
      READ(13,'(40X,F8.4)') EstCpueCV
      READ(13,'(40X,F8.1)') EstEffCatchSamp
      READ(13,'(40X,F8.4)') EstSurveyCV
      READ(13,'(40X,F8.1)') EstEffSurveySamp
      READ(13,'(40X,I4)') StartYrRecDevs
      READ(13,'(40X,I4)') DebugAssessment
      READ(13,'(40X,I5)') NoCPUEMPA
      READ(13,'(40X,I5)') UseFisheryCPUE


      IF (EstCpueCV.LT.0) EstCpueCV = CpueCV
      IF (EstEffCatchSamp.LT.0) EstEffCatchSamp = EffectSampCatch1
      IF (EstSurveyCV.LT.0) EstSurveyCV = SurvCV
      IF (EstEffSurveySamp.LT.0) EstEffSurveySamp = EffectSampSurvey1
      IbiasType = 0

	READ(13,'(40X,I5)') DoMSYOnly	    

      CLOSE(13)
C
      RETURN
      END SUBROUTINE ReadCntl

C
C ===========================================================================
C    
	  
      SUBROUTINE Equn()
C
C This subroutine calculates the equiliibrium numbers-at-age spatially -
C essentially a spatial yield-per-recruit analysis.
C	This calculates the unfished equilibrium number at age per recruit (spatially if the spatial model is activated)
C	Steps:
C		1. the first step is to calculate the number at age per recruit
C		2. then to convert it to spawning number per recruit (multiply number per recruit with fecundity)
C			 The notation SPBR is used here but it is number per recruit
C		3. Determine the proportion of recruitment into areas. If =0, then it is a unique value. 
C			 If >0, then the recruitment is distributed to different area based on the specification in 
C			 OM.dat (which specifies the center of distribution AND the variance around it (Gaussian distribution)
     
C      
	USE ShareDataRatioSubs
      IMPLICIT NONE
C
C     Local variables
      INTEGER Ispec,Iage,isex,IBlkNS,IBLkEW
      REAL*8 Tmp1,Tmp2
C
C     Calculate the equilibrium numbers-at-age; start of plus group programming
      DO 10000 Ispec = 1,Nspec
       DO 11000 Isex = 1,Nsex
        Nequ(Ispec,Isex,0) = 1
        DO 11100 Iage = 1,Amax(Ispec)
         Nequ(Ispec,Isex,Iage) = Nequ(Ispec,isex,Iage-1)*EXP(-M(Ispec))
11100   CONTINUE

      !-----------------------------------------------------------------------
!      PlusGroup = 0 removes the plus group from the initial conditions. Do PlusGroup=0 if doing local larval dynamics.
        IF (PlusGroup.EQ.1) THEN
        Nequ(Ispec,ISex,Amax(Ispec)) = Nequ(Ispec,Isex,Amax(Ispec))
     +         / (1-EXP(-M(Ispec)))
        ENDIF
       !------------------------------------------------------------------

11000  CONTINUE
       SBPR(Ispec) = 0
       DO 12000 Iage = 1,Amax(Ispec)
        SBPR(ISpec) = SBPR(Ispec) + Nequ(Ispec,1,Iage) * Fec(Ispec,Iage)
12000  CONTINUE        
10000 CONTINUE         
C
C     Find the recruitment profile
C ***************Carey is eliminating the recruitment movement profile if recruits are sedentary!!***********
C Add a new RecVirg option where cells are randomly assigned a RecVirg - some are just lower in productivity than others by a rote vector

      DO 20001 ISpec = 1,NSpec
	 IF (RecruitOpt(Ispec).EQ.0) THEN
	   RecVirg(Ispec,:,:) = 1
	 ELSEIF (RecruitOpt(Ispec).EQ.2) THEN
C         CALL PatchyRecVirg(Ispec)
	 ELSE
        DO 20000 IBlkNS = 1,NNSBlock
         DO 20000 IblkEW = 1,NEWBlock

         IF (IblkNS.LT.RecX(Ispec)) THEN
          Tmp1 = EXP(-1*(IblkNS-RecX(Ispec))**2.0/VrecXL(ISpec))
         ELSE
          Tmp1 = EXP(-1*(IblkNS-RecX(Ispec))**2.0/VrecXH(ISpec))
         ENDIF
         IF (IblkEW.LT.RecY(Ispec)) THEN
          Tmp2 = EXP(-1*(IblkEW-RecY(Ispec))**2.0/VrecYL(ISpec))
         ELSE
          Tmp2 = EXP(-1*(IblkEW-RecY(Ispec))**2.0/VrecYH(ISpec))
         ENDIF
         RecVirg(Ispec,IBlkNS,IBlkEW) = Tmp1*Tmp2
20000   CONTINUE   
	 ENDIF
20001	CONTINUE
C**************END Carey's modification**************************    
C
      RETURN
      END SUBROUTINE Equn

C
C ===========================================================================
C    
	  
      SUBROUTINE InitProj()
C
C	This subroutine calculates the spawning stock number and the number at age per recruit
C	Steps:
C		1. calculates the number at age for 1 recruit because RecVirg = 1
C		2. using this, calculate the survival rate per animal at age
C		3. if the species is moving, calculate this survival rate with the movement matrix
C		4. there is a weird thing for dealing with the plus group but basically the code copies the plus group 
C		   value to both the plus group and the Amax age. Aplus = Amax-1
C		5. then calculate the Spawning stock number per 1 unit recruit (the above survival rate per animal at age multiplied by fecundity)
C      
	USE ShareDataRatioSubs
      IMPLICIT NONE
C
C     Local variables 
      INTEGER IBlkNS,IBlkEW,JBlkNS,JBlkEW,AmaxM
      INTEGER ISpec,Isex,Yr,Iage,Itime
      REAL*8 Temp1(MaxNSBlk,MaxEWBlk,0:MaxAge)
      REAL*8 Temp2(MaxNSBlk,MaxEWBlk,0:MaxAge)
      REAL*8 CareyTmp
      REAL*8 Surv,SSBT,SSBT2
      INTEGER :: Aplus
C
C     Project forward to find equilibrium conditions
      DO 10000 Ispec = 1,Nspec
       DO 10000 Isex = 1,Nsex
C
C       Set the maximum age
        AmaxM = Amax(ISpec)
C
C       Set up Temporary storage
        DO 11000 IBlkNS = 1,NNSBlock
         DO 11000 IBlkEW = 1,NEWBlock
          DO 11000 Iage = 0,MaxAge
           Temp1(IBlkNS,IBlkEW,Iage) = 
     +        Nequ(Ispec,Isex,Iage)*RecVirg(ISpec,IBlkNS,IBlKEW) !Why does this line make sense? Apparently doesn't matter what distribution of age 1 and up.
11000   CONTINUE
C
C       Project forward    
!        DO 12000 Yr = 1,(Amax(Ispec))
        DO 12000 Yr = 1,50
C	  Projects forward by year: burns in distribution of adults after movement (needs to run over years to do this; 100 years will burn in perfectly, but 50 is also quite precise)
	   DO 12000 Itime = 1,Ntime
C
C         Project forwards in time period within one year (if Ntime = 1, dynamics are just yearly)
          Temp2 = 0
          DO 12100 IBlkNS = 1,NNSBlock
           DO 12100 IBlkEW = 1,NEWBlock
            DO 12100 Iage = 0,AmaxM
             Surv = Temp1(IBlkNS,IBlkEW,Iage)*EXP(-M(Ispec)*Rtime)     ! survival per animal at age
             DO 12110 JBlkNS = 1,NNSBlock
              DO 12110 JBlkEW = 1,NEWBlock
               Temp2(JBlkNS,JBlkEW,IAge) = Temp2(JBlkNs,JBlkEW,IAge) +  
     +           Surv*MixMat(Ispec,IBlkNS,IBlkEW,JBlkNS,JBlkEW)
12110        CONTINUE              
12100     CONTINUE           
          Temp1 = Temp2                                                ! it is simply to save the temporary object
C
C         Update age-structure
          IF (Itime.EQ.NTime) THEN
           DO 12200 IBlkNS = 1,NNSBlock
            DO 12200 IBlkEW = 1,NEWBlock
C commenting out the next two lines should remove the existence of the plus group within this burn-in loop.
              IF (PlusGroup.EQ.1) THEN
                Aplus = AmaxM-1
              Temp1(IBlkNS,IBlkEW,AmaxM) = 
     +         Temp1(IBlkNS,IBlkEW,AmaxM) + Temp1(IBlkNS,IBlkEW,AmaxM-1)
               ELSE
                Aplus = AmaxM
               ENDIF
C             DO 12210 Iage = AmaxM-1,1,-1
			DO 12210 Iage = Aplus,1,-1                               !Equals AmaxM if no plus group and AmaxM-1 otherwise (because plus group is updated before this line)
              Temp1(IBlkNS,IBlkEW,Iage) = Temp1(IBlkNS,IBlkEW,Iage-1) 
12210        CONTINUE  !Small loop over ages          
             Temp1(IBlkNS,IBlkEW,0) = RecVirg(Ispec,IBlkNS,IBlkEW)
12200      CONTINUE !Small loop over NS and EW space
C

C          Calculate SSB (Note: Fec(ISpec,Iage) = Maturity(Ispec,Iage)*BegWtAtAge(Ispec,Iage)
           IF (ISex.EQ.1) THEN
            SSBT = 0
            DO 12220 IBlkNS = 1,NNSBlock
              DO 12220 IBlkEW = 1,NEWBlock
              SSBT2 = 0
              DO 12221 Iage = 1,AmaxM
               SSBT2 = SSBT2 + Temp1(IBlkNS,IBlkEW,Iage)*Fec(ISpec,Iage)
12221         CONTINUE
              SSBT = SSBT + SSBT2
12220       CONTINUE
            SSBEq(Ispec) = SSBT
           ENDIF 
C             
          ENDIF
C
12000   CONTINUE         
       TestNeq = 0
C       Copy the results into the N-matrix
        DO 13000 IBlkNS = 1,NNSBlock
         DO 13000 IBlkEW = 1,NEWBlock
          DO 13000 Iage = 0,AmaxM
           Neq(ISpec,Isex,IBlkNS,IBlkEW,Iage)=Temp1(IBlkNS,IBlkEW,Iage)
C           TestNeq(ISpec,Isex,Iage) = TestNeq(Ispec,Isex,Iage) +
C     +       Neq(Ispec,Isex,IBlkNS,IBlkEW,Iage)
            CareyTmp=TestNeq(Ispec,Isex,IAge)
            TestNeq(Ispec,Isex,IAge) = CareyTmp + 
     +            Neq(ISpec,Isex,IBlkNS,IBlkEW,Iage)


13000   CONTINUE      
10000 CONTINUE      
C       
      RETURN
600   FORMAT(1x,I2,1x,I2,1x,I2,1x,I2,1x,I2,1x,F9.2)
      END SUBROUTINE InitProj

C
C ===========================================================================
C    
	  
      SUBROUTINE ProjForward(Kvals)
C
C	This subroutine basically uses the results from InitProj and calculates the spawning stock number at age based on the number at age per recruit 
C	Steps:
C		1. first time it calculates the recruitment value			 
C		2. multiplies the above with the spawning stock number at age (calculated in InitPorj)
C		3. the stock recruit relationship
C		4. then project forward the SSN and recruitment using the ProjYr subroutine
C      
	USE ShareDataRatioSubs
      IMPLICIT NONE
C
C     Global variables
      REAL*8 KVals(MaxSpec) 
C
C     Local variables 
      INTEGER ISpec,IBlkNS,IBlkEW,Isex,Iage,Year,Iarea
      REAL*8 TheEffort(10)
      REAL*8 Recruits,SSBT,SSBT2,Rec2
      REAL(8),DIMENSION(Nspec,Nsex,0:Amax(1)) :: Nsums
C
C	Reset larval abundance
C	Larvae = 0
      Nsums = 0           ! sum of N over all ages
      NatAgeCompact(:,:,MinPYr,:,:)=0
      R0forSS3 = 0
C     Set up initial conditions given the KVals
      DO 10000 Ispec = 1,Nspec
C
C      Calculate recruitment (SSBEq is SSB per recruit, Kvals is SSB total 
C         => Recruits*SSBEq = Kvals => Recruits = Kvals/SSBEq)
       Recruits = Kvals(Ispec) / SSBEq(ISpec)

C      Copy from N to Temp1
       DO 11000 Isex = 1,Nsex
        DO 11000 IBlkNS = 1,NNSBlock
         DO 11000 IBlkEW = 1,NEWBlock
          DO 11100 Iage = 0,Amax(Ispec)
           N(Ispec,Isex,IBlkNS,IBlkEW,MinPyr,Iage) = Recruits *
     +        Neq(Ispec,Isex,IBlkNs,IBlkEW,Iage)
      !Carey adds in a non-spatial calculation of numbers at age at the beg. of year
           NatAgeCompact(Ispec,Isex,MinPYr,Iage,0) = NatAgeCompact
     +(Ispec,Isex,MinPYr,Iage,0)+N(Ispec,Isex,IBlkNS,IBlkEW,MinPyr,Iage)
           Nsums(Ispec,Isex,Iage) = Nsums(Ispec,Isex,Iage) + 
     +       N(Ispec,Isex,IBlkNS,IBlkEW,MinPYr,Iage)
11100     CONTINUE     
          RecSum(Ispec,MinPyr,IblkNS,IblkEW) =
     +      N(ISpec,1,IBlkNS,IBlkEW,MinPYr,0)  ! This just records the number of recruits when the model is in unfished equilibrium, but doesn't get used anywhere, except to write the info to a file.
          R0forSS3(Ispec) = R0forSS3(Ispec) + 
     +      N(Ispec,Isex,IBlkNS,IBlkEW,MinPyr,0)
11000  CONTINUE        
C
!Right here add a loop for figuring out R0forSS3 by assessment area. Use MinNS,MinEW,...
      R0forSS3B = 0
        DO 11144 Isex = 1,Nsex
         DO 11144 Iarea = 1,Narea
          DO 11133 IBlkNS = MinNS(Iarea),MaxNS(Iarea)
            DO 11133 IBlkEW = MinEW(Iarea),MaxEW(Iarea)
           R0forSS3B(Ispec,Iarea) = R0forSS3B(Ispec,Iarea) +
     +     N(Ispec,Isex,IBlkNS,IBlkEW,MinPYr,0)
11133   CONTINUE
11144   CONTINUE

C      Calculate initial mature
       SSBT = 0
       DO 12000 Iarea = 1,Narea
        SpBioArea(Ispec,Iarea,MinPYr) = 0
        DO 12100 IBlkNS = MinNS(Iarea),MaxNS(Iarea)
         DO 12100 IBlkEW = MinEW(Iarea),MaxEW(Iarea)
          SSBT2 = 0
          DO 12110 Iage = 1,Amax(Ispec)    
           SSBT2=SSBT2+N(Ispec,1,IBlkNS,IBlkEW,MinPYr,Iage) * 
     +       Fec(Ispec,Iage)
12110     CONTINUE 
          SSBT = SSBT + SSBT2
          SpBio(Ispec,IBlkNS,IBlkEW,MinPYr) = SSBT2
          SpBioArea(Ispec,Iarea,MinPyr) = SpBioArea(Ispec,Iarea,MinPyr) 
     +      + SSBT2
12100   CONTINUE
12000  CONTINUE           
       SPBioTot(Ispec,MinPYr) = SSBT
C
C      Stock-recruitment parameters - 1 GLOBAL density dependence
       IF (SrType(Ispec).EQ.1) THEN
        DO 13100 IBlkNS = 1,NNSBlock
         DO 13100 IBlkEW = 1,NEWBlock
          Rec2 = Recruits*Neq(Ispec,1,IBlkNs,IBlkEW,0)
          SRR0(Ispec,IBlkNS,IblkEW) = Rec2
          SRAlpha(Ispec,IBlkNS,IblkEW) = 
     +        (1-Steep(Ispec))/(4*Steep(Ispec)*Rec2)
          SRBeta(Ispec,IBlkNS,IblkEW) = 
     +        (5*Steep(Ispec)-1)/(4*Steep(ISpec)*Rec2)     
13100   CONTINUE         
       ENDIF
C
C      Stock-recruitment parameters - 2 LOCAL density dependence
       IF (SrType(Ispec).EQ.2) THEN
        DO 13200 IBlkNS = 1,NNSBlock
         DO 13200 IBlkEW = 1,NEWBlock
          Rec2 = Recruits*Neq(Ispec,1,IBlkNs,IBlkEW,0)
          SRR0(Ispec,IBlkNS,IblkEW) = Rec2
          SRAlpha(Ispec,IBlkNS,IblkEW) = 
     +        (1-Steep(Ispec))/(4*Steep(Ispec)*Rec2)
          SRBeta(Ispec,IBlkNS,IblkEW) = 
     +        (5*Steep(Ispec)-1)/(4*Steep(ISpec)*Rec2)     
13200   CONTINUE         
       ENDIF
C
C      Stock-recruitment parameters - 3 LOCAL density dependence
       IF (SrType(Ispec).EQ.3) THEN
        DO 13300 IBlkNS = 1,NNSBlock
         DO 13300 IBlkEW = 1,NEWBlock
          Rec2 = Recruits*Neq(Ispec,1,IBlkNs,IBlkEW,0)
          SRR0(Ispec,IBlkNS,IblkEW) = Rec2
          SRAlpha(Ispec,IBlkNS,IblkEW) = 
     +        (1-Steep(Ispec))/(4*Steep(Ispec)*Rec2)
          SRBeta(Ispec,IBlkNS,IblkEW) = 
     +        (5*Steep(Ispec)-1)/(4*Steep(ISpec)*Rec2)     
13300   CONTINUE         
       ENDIF
C
C      Stock-recruitment parameters - 4 LOCAL and GLOBAL density dependence
       IF (SrType(Ispec).EQ.4) THEN
        DO 13400 IBlkNS = 1,NNSBlock
         DO 13400 IBlkEW = 1,NEWBlock
          Rec2 = Recruits*Neq(Ispec,1,IBlkNs,IBlkEW,0)
          SRR0(Ispec,IBlkNS,IblkEW) = Rec2
          SRAlpha(Ispec,IBlkNS,IblkEW) = 
     +        (1-Steep(Ispec))/(4*Steep(Ispec)*Rec2)
          SRBeta(Ispec,IBlkNS,IblkEW) = 
     +        (5*Steep(Ispec)-1)/(4*Steep(ISpec)*Rec2)     
13400   CONTINUE         
       ENDIF
C
10000 CONTINUE
C
C     Project forward
      Closed = 1
      DO 20000 Year = MinPYr,BaseYr
C      
C      Impose MPA
       IF (Year.EQ.MPAYear) THEN
        DO 21100 Iarea = 1,NClosedAreas
         DO 21110 IblkNS = MinCNS(Iarea),MaxCNS(Iarea)
          DO 21110 IBlkEW = MinCEW(Iarea),MaxCEW(Iarea)
           Closed(IBlkNS,IBlkEW) = 0
21110    CONTINUE
21100   CONTINUE          
       ENDIF
C       
       TheEffort(1) = TotEffort(Year)
       CALL ProjYr(Year,TheEffort)

	IF (WriteOpt.EQ.1) THEN  !Write out information during an MSY run
        CALL Summary2(Year)
	ENDIF

20000 CONTINUE       

      RETURN
      END SUBROUTINE ProjForward

C
C ===========================================================================
C    
	  
      SUBROUTINE ProjYr(Year,TheEffort)
C
C	This subroutine projects the population forward (the SSN, recruits and also adds on the SSB, total biomass
C	Steps
C		1. for each year, it calculates the total biomass at age
C		2. it calculates the exploitable biomass ExpBio based on biomass*catchability and the biomass = Temp1 is based on mid year = 1/2M removed
C		3. Then effort (input) is distributed based on gravity rule (or uniformly depending on the option) if spatial model
C		4. F is then determined as F = Qval*Effort
C		5. catch is then N*exp(-M/2)*Sel*q*E (where N is the mid-year N)
C		6. survey catch at age is then based on after the catch
C		7. After removing catch and all M, then movement (if present)
C		8. calculate the biomass and the spawning biomass
C		9. then determine the Recruit for next year (including or not the process error)
C
C      
	USE ShareDataRatioSubs
      IMPLICIT NONE
C
C     Global variables
      INTEGER Year
      REAL*8 TheEffort(10),TheQ
C
C     Local variables
      INTEGER Ispec,Isex,IBlkNS,IBlkEW,Iage,Itime,AmaxM,Iarea,Jarea
      INTEGER JBlkNS,JBlkEW,MinEffA,MaxEffA,Ipnt,IeffA
      REAL*8 Temp1(MaxSpec,Nsex,MaxNsBlk,MaxEWBlk,0:MaxAge)
      REAL*8 Temp2(MaxNsBlk,MaxEWBlk,0:MaxAge)
      REAL*8 Eff(MaxNsBlk,MaxEWBlk)
      REAL*8 Surv,SSBT,SSBT2,Tmp
      REAL*8 ExpBioT(MaxNSBlk,MaxEWBlk),TotalExpBio,EffBlock
      REAL*8 CatN,Depl,Depl1,Alpha,Beta,TotRec

	REAL(8), DIMENSION(NSpec,MinPYr:MaxPYr) :: SumSPBio,SumLarvae
      INTEGER :: Aplus,TotOpen

C	Carey wants to record the total effort (NOT BY BLOCK) that is used in ProjYr!
C only if there is no spatial effort zoning!!!
      !This makes no sense if TheEffort is unrealistically high from finding optimal effort levels
! Moving this!	
!      GlobalEffort(Year) = MIN(MAX(TheEffort(1),0.0),9E+12)

      NatAgeCompact(:,:,Year,:,0) = 0
      TotalBiomass(:,Year) = 0
C     Copy from N to Temp1
      DO 10000 Ispec = 1,Nspec
       DO 10000 Isex = 1,Nsex
        DO 10000 IBlkNS = 1,NNSBlock
         DO 10000 IBlkEW = 1,NEWBlock
          DO 10000 Iage = 0,Amax(Ispec)
           Temp1(Ispec,Isex,IBlkNs,IBlkEW,Iage) = 
     +      N(Ispec,Isex,IBlkNS,IBlkEW,Year,Iage)      
       NatAgeCompact(Ispec,Isex,Year,Iage,0) = NatAgeCompact
     +(Ispec,Isex,Year,Iage,0)+N(Ispec,Isex,IBlkNS,IBlkEW,Year,Iage)
      TotalBiomass(Ispec,Year) = TotalBiomass(Ispec,Year) + 
     +  NatAgeCompact(Ispec,Isex,Year,Iage,0)*WghtBeg(Ispec,Isex,Iage,0)
10000 CONTINUE  
        
        TotalRecruits(:,Year) = 0
        DO 201 Ispec = 1,Nspec
         DO 201 Isex = 1,Nsex
          TotalRecruits(Ispec,Year) = TotalRecruits(Ispec,Year) +
     +     NatAgeCompact(Ispec,Isex,Year,0,0)
201    CONTINUE

C     Reset the total catch
      DO 11000 Ispec = 1,Nspec
       DO 11000 Iarea = 1,Narea
        CatBio(Ispec,Iarea,Year) = 0
        TruSuv(Ispec,Iarea,Year) = 0
        TruExp(Ispec,Iarea,Year) = 0
        Effort(Iarea,Year) = 0
        SpBioArea(Ispec,Iarea,Year+1) = 0
        DO 11100 Isex = 1,Nsex 
         DO 11100 Iage = 0,Amax(Ispec)
          CatNum(Ispec,Iarea,Year,Isex,Iage) = 0
          SuvNum(Ispec,Iarea,Year,Isex,Iage) = 0
11100    CONTINUE           
11000 CONTINUE         
C
C     Move ahead by months
      DO 20000 Itime = 1,Ntime
  
C      Remove first-half of the natural mortality
       DO 21000 Ispec = 1,Nspec
        Surv = EXP(-M(Ispec)*Rtime/2.0) !Survival halfway through the year if Rtime = 1 (because we assume that fishing happened at the middle of the year)
        DO 21100 ISex = 1,Nsex
         DO 21100 IBlkNS = 1,NNSBlock
          DO 21100 IBlkEW = 1,NEWBlock
           DO 21100 Iage = 0,Amax(Ispec)
C	Temp1 is N; below it's N but with the first half of the nat. mort. removed
            Temp1(ISpec,ISex,IBlkNs,IBlkEW,Iage) = Surv *
     +        Temp1(ISpec,Isex,IBlkNS,IBlkEW,Iage)       
21100   CONTINUE     
21000  CONTINUE       
C
C      Find the effort by spatial cell
       TotalExpBio = 0
       DO 22100 IBlkNS = 1,NNSBlock
        Do 22100 IBlkEW = 1,NEWBlock
         ExpBioT(IBlkNS,IBlkEW) = 0
         DO 22110 ISpec = 1,Nspec
          Tmp= 0
          DO 22111 Isex = 1,Nsex
           DO 22111 Iage = 0,Amax(Ispec)
            Tmp = Tmp + Temp1(ISpec,Isex,IBlkNs,IBlkEW,Iage)*
     +           Closed(IbLKNS,IBlkEW)*Sel(Ispec,Iage)*
     +           Wght(Ispec,Isex,Iage,Itime)
          IF (Year.GT.0) THEN
C           WRITE(*,*) Isex,Iage,Temp1(Ispec,Isex,IBlkNs,IblkEW,Iage)
          ENDIF
22111     CONTINUE       
!          IF (DoProcessError.EQ.1) THEN   ! Qval is the globla catchability and theQ is the catchability by space (defined by RecQ)
            TheQ = Qval(Ispec)*RecQ(Isim,ISpec,IBlkNS,IblkEW)
!          ELSE
!            TheQ = Qval(Ispec)
!          ENDIF
          ExpBio(Ispec,IBlkNS,IBlkEW,Year,Itime) = Tmp*TheQ       
          ExpBioT(IBlkNS,IBlkEW) = ExpBioT(IBlkNS,IBlkEW) + Tmp*TheQ
22110    CONTINUE
         TotalExpBio = TotalExpBio + ExpBioT(IBlkNS,IBlkEW)
22100  CONTINUE          
C

C	Write out TotalExpBio, ExpBioT(IBlkNS,IBlkEW)
C	IF (Year.EQ.-18) THEN
C	  WRITE(*,'(1x,A)') 'line 1891 TotalExpBio Temp1age0 Nage0'
C	  WRITE(*,*) TotalExpBio, Temp1(1,1,1,1,0), N(1,1,1,1,-18,0) 
C	ENDIF

C      Effort allocation module
       IF (Year.LE.FinDepYr) THEN
        MinEffA = 0
        MaxEffA = 0
       ELSE
        MinEffA = 1
        MaxEffA = NeffortAreas 
       ENDIF 
       Ipnt = 0
       DO 45000 IeffA = MinEffA,MaxEffA
C
C       Update Effort Pointer
        Ipnt = Ipnt + 1
C
C       Find the total exploitable biomass in this area
        TotalExpBio = 0
        DO 45100 IBlkNS = MinENS(IeffA),MaxENS(IeffA)
         Do 45100 IBlkEW = MinEEW(IeffA),MaxEEW(IeffA)
          TotalExpBio = TotalExpBio + ExpBioT(IBlkNS,IBlkEW)
45100   CONTINUE          
C
C       If no hyperstability, need this to allocate effort to blocks:
        TotOpen = 0
        DO 45201 IBlkNS = 1,NNSBlock
         DO 45201 IBlkEW = 1,NEWBlock
           TotOpen = TotOpen + Closed(IBlkNS,IBlkEW)
45201  CONTINUE
C       Allocate effort to blocks
        DO 45200 IBlkNS = MinENS(IeffA),MaxENS(IeffA)
         Do 45200 IBlkEW = MinEEW(IeffA),MaxEEW(IeffA)
          IF (TheEffort(Ipnt).GT.0) THEN
            IF (Hyperstability.EQ.0) THEN
            EffBlock = 1.0/REAL(TotOpen)*TheEffort(Ipnt)
            ELSE
           EffBlock = ExpBioT(IBlkNS,IBlkEW)/TotalExpBio*TheEffort(Ipnt)
            ENDIF
          ELSE
           EffBlock = 0
          ENDIF 
          Eff(IblkNS,IBlkEW) = EffBlock
          NetEffort(IblkNS,IBlkEW,Year) = EffBlock
          DO 45200 Ispec = 1,NSpec
!          IF (DoProcessError.EQ.1) THEN
            FullFish(Ispec,IBlkNS,IBlkEW,Year,Itime) = 
     +      MIN(0.99d0,Qval(Ispec)*RecQ(Isim,Ispec,IBlkNs,IBlkEW)*
     +     EffBlock)
 !         ELSE
 !         FullFish(Ispec,IBlkNS,IBlkEW,Year,Itime) = 
 !    +      MIN(0.99d0,Qval(Ispec)*EffBlock)
 !         ENDIF
           IF (FullFish(Ispec,IBlkNS,IBlkEW,Year,Itime).EQ.0.99d0.AND.
     +       ASSOPT.NE.-1) THEN
             Eff(IblkNS,IBlkEW) = 0.99d0/Qval(Ispec)
             NetEffort(IblkNS,IBlkEW,Year) = Eff(IblkNS,IBlkEW)
           ENDIF
45200   CONTINUE     
C
45000  CONTINUE       
C
!     Moved from above: record the global effort
      GlobalEffort(Year) = 0.0
      DO 65110 IBlkNS = 1,NNSBlock
        DO 65110 IBlkEW = 1,NEWBlock     
      GlobalEffort(Year)=MIN(MAX(GlobalEffort(Year)+NetEffort(IBlkNS,
     + IBlkEW,Year),0.0),9E+12)
65110 CONTINUE
C      Record the effort
       DO 45300 Iarea = 1,Narea
        DO 45300 IBlkNS = MinNS(Iarea),MaxNS(Iarea)
         DO 45300 IBlkEW = MinEW(Iarea),MaxEW(Iarea)
          Effort(Iarea,Year) = Effort(Iarea,Year) + Eff(IBlkNS,IBlkEW)
45300  CONTINUE
C
C	Write out TotalExpBio, ExpBioT(IBlkNS,IBlkEW)
C	IF (Year.EQ.-18) THEN
C	  WRITE(*,'(1x,A)') 'line 1943 TotalExpBio Nage0 Effort FullFish'
C	  WRITE(*,*) TotalExpBio, N(1,1,1,1,-18,0), Effort(1,-18), 
C     +	  FullFish(1,1,1,-18,1) 
C	ENDIF

	SpatialSuvNumSum(:,Year,:,:) = 0
	SpatialSuvNum(:,:,:,:,Year,:) = 0
C      Store the expected value of the survey and CPUE index as well
C      as the survey numbers-at-age
       IF (Itime.EQ.(Ntime-1)/2+1) THEN
        DO 25100 Ispec = 1,nspec
         DO 25110 Iarea = 1,Narea
          DO 25110 IBlkNS = MinNS(Iarea),MaxNS(Iarea)
           DO 25110 IBlkEW = MinEW(Iarea),MaxEW(Iarea)
            DO 25110 Isex = 1,Nsex
             DO 25110 Iage = 0,Amax(Ispec)
              CatN = Temp1(Ispec,Isex,IBlkNS,IBlkEw,Iage)*
     +          Closed(IBlkNS,IBlkEW)*Sel(Ispec,Iage)*
     +          FullFish(Ispec,IBlkNS,IBlkEW,Year,Itime)
C   
C	Tmp is N with 1/2 nat. mort. removed but don't save to Temp1...  
             Tmp = Temp1(Ispec,Isex,IBlkNS,IBlkEW,Iage) - CatN/2.0
C	Spatially explicit survey numbers: SpatialSuvNum (added by Carey)
	       SpatialSuvNum(Ispec,Isex,IBlkNS,IBlkEW,Year,Iage) = 
     +		   (Temp1(Ispec,Isex,IBlkNS,IBlkEW,Iage)-CatN/2.0)*
     +            SSel(Ispec,Iage)
		   
	       SpatialSuvNumSum(Ispec,Year,IBlkNS,IBlkEW) = 
     +		   SpatialSuvNumSum(Ispec,Year,IBlkNS,IBlkEW)
     +           +SpatialSuvNum(Ispec,Isex,IBlkNS,IBlkEW,Year,Iage)

             IF (NoSurveyClosed.EQ.1) Tmp = Tmp * CLosed(IblkNS,IBlkEW)
             TruSuv(Ispec,Iarea,Year) = TruSuv(Ispec,Iarea,Year) + 
     +         SSel(Ispec,Iage)*Tmp*Wght(Ispec,Isex,Iage,Itime)        
             TruExp(Ispec,Iarea,Year) = TruExp(Ispec,Iarea,Year) +
     +         CLosed(IBlkNS,IBlkEW)*Sel(Ispec,Iage) 
     +                         *Tmp*Wght(Ispec,Isex,Iage,Itime)        
             SuvNum(Ispec,Iarea,Year,Isex,Iage) =
     +         SuvNum(Ispec,Iarea,Year,Isex,Iage) + SSel(Ispec,Iage)*Tmp
25110    CONTINUE         
25100   CONTINUE            
       ENDIF
	 IF (Year.EQ.3) THEN
C	  WRITE(*,'(1x,A)') 'Year,Isex,IBlkNS,IBlkEW,IAGE,SpatialSuvNum'
C	  WRITE(*,*) Year,Isex,IBlkNS,IBlkEW,IAGE,
C     +	  SpatialSuvNum(1,1,1:NNSBlock,1:NEWBlock,Year,0:Amax(1))
	 END IF
C

C	Write out TotalExpBio, ExpBioT(IBlkNS,IBlkEW)
C	IF (Year.EQ.-18) THEN
C	  WRITE(*,'(1x,A)') 'line 1994 Temp1age0 Tmp CatN Nage0'
C	  WRITE(*,*) Temp1(1,1,1,1,0),Tmp,CatN, N(1,1,1,1,-18,0)
C	ENDIF


C      Store the catch data
       DO 25200 Ispec = 1,Nspec
        DO 25210 Iarea = 1,Narea
         DO 25210 IBlkNS = MinNS(Iarea),MaxNS(Iarea)
          DO 25210 IBlkEW = MinEW(Iarea),MaxEW(Iarea)
           CatArea(ISpec,IBlkNS,IBlkEW,Year) = 0
           DO 25210 Isex = 1,Nsex
            DO 25210 Iage = 0,Amax(Ispec)
             CatN = Temp1(Ispec,Isex,IBlkNS,IBlkEw,Iage)*
     +          Closed(IBlkNS,IBLkEW)*Sel(Ispec,Iage)*
     +          FullFish(Ispec,IBlkNS,IBlkEW,Year,Itime)
C            
             CatNum(Ispec,Iarea,Year,Isex,Iage) =
     +         CatNum(Ispec,Iarea,Year,Isex,Iage) + CatN
             CatBio(Ispec,Iarea,Year) = CatBio(Ispec,Iarea,Year) +
     +         CatN*Wght(Ispec,ISex,Iage,Itime)       
              CatArea(Ispec,IBlkNS,IBlkEW,Year) = 
     +         CatArea(Ispec,IBlkNS,IBlkEW,Year) +
     +         CatN*Wght(Ispec,ISex,Iage,Itime)
25210   CONTINUE        
25200  CONTINUE       
C
C      Remove catch and the 2nd half of natural mortality
       DO 22300 Ispec = 1,Nspec
        Surv = EXP(-M(Ispec)*Rtime/2.0)
        DO 22300 IBlkNS = 1,NNSBlock
         DO 22300 IBlkEW = 1,NEWBLock
          DO 22310 Isex = 1,Nsex
           DO 22310 Iage = 0,Amax(Ispec)
            CatN = Temp1(Ispec,Isex,IBlkNS,IBlkEw,Iage)*
     +          Closed(IBlkNS,IBlkEW)*Sel(Ispec,Iage)*
     +          FullFish(Ispec,IBlkNS,IBlkEW,Year,Itime)
C          Temp1 is N removing all of the catch and nat mort, but before movement  
            Temp1(Ispec,ISex,IBlkNS,IBLkEW,Iage) = 
     +          (Temp1(Ispec,Isex,IBLkNS,IBlkEW,Iage) - CatN)*Surv
22310     CONTINUE             
22300  CONTINUE         
C	Write out TotalExpBio, ExpBioT(IBlkNS,IBlkEW)
C	IF (Year.EQ.-18) THEN
C	  WRITE(*,'(1x,A)') 'line 2038 Temp1age0 CatN Nage0 FullFish'
C	  WRITE(*,*) Temp1(1,1,1,1,0),CatN, N(1,1,1,1,-18,0), 
C     +	  FullFish(1,1,1,-18,1)
C	ENDIF

C      Move animals around spatially
      CALL Movement2ForSS3(Temp1)
      IF (Isim.eQ.1.AND.Year.EQ.MinPYr) THEN
        WRITE(79,*) 'Nsim,Year, FromArea, ToArea, Age, ProportionMoved'
        WRITE(79,'(4(I4,1x))') Nsim,FinDepYr,Narea,Amax(1)
      ENDIF

      IF (Isim.EQ.1) THEN
      DO 56 Iarea = 1,NArea
        DO 56 Jarea = 1,NArea
          DO 56 Iage = 0,Amax(1)
           !Write out SS-style movement
           WRITE(79,'(5(I4,1x),F20.15)') Isim,Year,Iarea,Jarea,Iage,
     +       PropN(Iarea,Jarea,Iage)
56      CONTINUE
      ENDIF

       DO 24000 Ispec = 1,Nspec
        DO 24000 Isex = 1,Nsex
         Temp2 = 0
         DO 24100 IBlkNS = 1,NNSBlock
          DO 24100 IBlkEW = 1,NEWBlock
           DO 24100 Iage = 0,Amax(Ispec)
            Tmp = Temp1(Ispec,Isex,IBlkNS,IBlkEW,IAge)
            DO 24110 JBlkNS = 1,NNSBlock
             DO 24110 JBlkEW = 1,NEWBlock
              Temp2(JBlkNS,JBlkEW,IAge) = Temp2(JBlkNs,JBlkEW,IAge) +
     +          MixMat(Ispec,IBlkNS,IBlkEW,JBlkNS,JBlkEW) * Tmp
24110       CONTINUE              
24100    CONTINUE 
         DO 24200 IBlkNS = 1,NNSBlock
          DO 24200 IBlkEW = 1,NEWBlock
           DO 24200 Iage = 0,Amax(Ispec)
C	Temp1 is now N with all nat mort and catch gone executing movement
            Temp1(ISpec,Isex,IBlkNS,IBlkEW,Iage) = 
     +           Temp2(IBlkNS,IBlkEW,Iage)
24200     CONTINUE            
24000  CONTINUE    
C
C	IF (Year.EQ.-18) THEN
C	  WRITE(*,'(1x,A)') 'line 2067 Temp1age0 Tmp Temp2(1,1,0) CatN'
C	  WRITE(*,*) Temp1(1,1,1,1,0),Tmp, Temp2(1,1,0), CatN
C	ENDIF


20000 CONTINUE      
C
C     Copy from Temp1 to N
      DO 31000 Ispec = 1,Nspec
       DO 31000 Isex = 1,Nsex
        DO 31000 IBlkNS = 1,NNSBlock
         DO 31000 IBlkEW = 1,NEWBlock
C	N as of the beginning of next year, updating age structure (this was already a comment before removing the plus group)
        ! Adding a switch to do the plus group (Read in during READSPEC)
          AmaxM = Amax(Ispec)
          IF (PlusGroup.EQ.1) THEN
            Aplus = AmaxM-1
C           Calculate the plus group if PlugGroup = 1
            N(Ispec,Isex,IBlkNS,IBlkEW,Year+1,AmaxM) =
     +       Temp1(Ispec,Isex,IBlkNs,IBlkEW,AmaxM) +
     +       Temp1(Ispec,Isex,IBlkNs,IBlkEW,AmaxM-1)
          ELSE
            Aplus = AmaxM
          ENDIF
C          DO 31100 Iage = 1,AmaxM-1 if there is a plus group and Iage=1,AmaxM with no plus group
          DO 31100 Iage = 1,Aplus !Aplus accounts for whether calculating the plus group or not
           N(Ispec,Isex,IBlkNS,IBlkEW,Year+1,Iage) =
     +      Temp1(Ispec,Isex,IBlkNs,IBlkEW,Iage-1)
31100     CONTINUE     
31000 CONTINUE        
C
C	IF (Year.EQ.-18) THEN
C	  WRITE(*,'(1x,A)') 'line 2093 Temp1age0 Nage0 CatN'
C	  WRITE(*,*) Temp1(1,1,1,1,0),N(1,1,1,1,Year,0),CatN
C	ENDIF

C     Calculate the spawning biomass
      DO 32000 Ispec = 1,Nspec
       SSBT = 0
       DO 32100 Iarea = 1,Narea
        DO 32100 IBlkNS = MinNS(Iarea),MaxNS(Iarea)
         DO 32100 IBlkEW = MinEW(Iarea),MaxEW(Iarea)
          SSBT2 = 0
          DO 32110 Iage = 1,Amax(Ispec)    
           SSBT2=SSBT2+N(Ispec,1,IBlkNS,IBlkEW,Year+1,Iage) * 
     +        Fec(Ispec,Iage)
32110     CONTINUE 
          SSBT = SSBT + SSBT2
          SpBio(Ispec,IBlkNS,IBlkEW,Year+1) = SSBT2
          SpBioArea(Ispec,Iarea,Year+1) = SpBioArea(Ispec,Iarea,Year+1) 
     +      + SSBT2
32100  CONTINUE           
       SPBioTot(Ispec,Year+1) = SSBT


32000 CONTINUE 
C	IF (Year.EQ.-18) THEN !2nd WRITE
C	  WRITE(*,'(1x,A)') 'line 2116 SSBT N(1,1,1,1,-18,0) N(1,1,1,1,-17,1)'
C	  WRITE(*,*) SSBT, N(1,1,1,1,-18,0), N(1,1,1,1,-17,1)      
C	ENDIF
C     Add in some recruitment
      DO 33000 Ispec = 1,Nspec
C
C      How to proceed depends on the form of density dependence
C      ========================================================
C
C      Type = 1; Global density-dependence
       IF (SrType(Ispec).EQ.1) THEN
        Depl =  SPBioTot(Ispec,Year+1)/SPBioTot(Ispec,MinPYr)
        DO 33100 IBlkNS = 1,NNSBlock
         DO 33100 IBlkEW = 1,NEWBlock
          Alpha = SRalpha(Ispec,IBlkNS,IBlkEW)
          Beta = SRbeta(Ispec,IBlkNS,IBlkEW)
          DO 33110 Isex = 1,Nsex
           N(Ispec,Isex,IBlkNs,IBlkEW,Year+1,0) = 
     +        Depl / (Alpha + Beta * Depl)
33110     CONTINUE
33100   CONTINUE         
       ENDIF
C       
C      Type = 2; Local density-dependence (at settlement)
C       IF (SrType(Ispec).EQ.2) THEN
C        DO 33200 IBlkNS = 1,NNSBlock
C         DO 33200 IBlkEW = 1,NEWBlock
C          Depl =  MAX(0.001,SPBio(Ispec,IBLkNS,IBlkEW,Year+1)/
C     +                             SPBio(Ispec,IBLkNS,IBlkEW,MinPYr))
C          Alpha = SRalpha(Ispec,IBlkNS,IBlkEW)
C          Beta = SRbeta(Ispec,IBlkNS,IBlkEW)
C          DO 33200 Isex = 1,Nsex
C           N(Ispec,Isex,IBlkNs,IBlkEW,Year+1,0) = 
C     +        Depl / (Alpha + Beta * Depl)
C33200   CONTINUE
C       ENDIF
C    
C      Type = 2; Local density-dependence (at settlement) with Carey's larval movement
	

       IF (SrType(Ispec).EQ.2) THEN
C	right here apply MixMat to Depl (make depl in a separate IBlk and JBlk loop)  Depl is basically serving as # of eggs.
	Larvae(:,:,:,Year+1) = 0
	IF (Year.EQ.MinPYr) THEN
        Larvae(:,:,:,MinPYr) = 0
	ENDIF
        DO 33202 IBlkNS = 1,NNSBlock
         DO 33202 IBlkEW = 1,NEWBlock

	     DO 33201 JBlkNS = 1,NNSBlock
	      DO 33201 JBlkEW = 1,NEWBlock
             Larvae(Ispec,JBlkNS,JBlkEW,Year+1)=
     +	   Larvae(Ispec,JBlkNS,JBlkEW,Year+1) +
     +       MixMatLarval(Ispec,IBlkNS,IBlkEW,JBlkNS,JBlkEW) * 
     +       SPBio(Ispec,IBlkNS,IBlkEW,Year)    
	       IF (Year.EQ.MinPYr) THEN
             Larvae(Ispec,JBlkNS,JBlkEW,MinPYr)=
     +	   Larvae(Ispec,JBlkNS,JBlkEW,MinPYr) +
     +       MixMatLarval(Ispec,IBlkNS,IBlkEW,JBlkNS,JBlkEW) * 
     +       SPBio(Ispec,IBlkNS,IBlkEW,MinPYr)
	       ENDIF
33201      CONTINUE
33202   CONTINUE
        DO 33200 IBlkNS = 1,NNSBlock
         DO 33200 IBlkEW = 1,NEWBlock
          Depl =  MAX(0.001,Larvae(Ispec,IBLkNS,IBlkEW,Year+1)/
     +                             Larvae(Ispec,IBLkNS,IBlkEW,MinPYr))
          Alpha = SRalpha(Ispec,IBlkNS,IBlkEW)
          Beta = SRbeta(Ispec,IBlkNS,IBlkEW)
          DO 33200 Isex = 1,Nsex
           N(Ispec,Isex,IBlkNs,IBlkEW,Year+1,0) = 
     +        Depl / (Alpha + Beta * Depl)
33200   CONTINUE
       ENDIF

C	IF (Year.EQ.-18) THEN !2nd WRITE
C	  WRITE(*,'(1x,A)') 'line 2190 Depl,Alpha,Beta,Larvae(1,1,1,-18),
C     +	  Larvae(1,1,1,-17),
C     +	  Larvae(1,1,1,MinPYr),SPBio(1,1,1,-18),SPBio(1,1,1,-50)'
C	  WRITE(*,*) Depl,Alpha,Beta,Larvae(1,1,1,-18),Larvae(1,1,1,-17),
C     +	  Larvae(1,1,1,MinPYr),SPBio(1,1,1,-18),SPBio(1,1,1,-50) 
C	ENDIF

C Print some output information to see if there are differences between total "SPBio" and total "Larvae"
	DO 44445 IBlkNS = 1,NNSBlock
	  DO 44445 IBlkEW = 1,NEWBlock
         SumSPBio(Ispec,Year+1) = SumSPBio(Ispec,Year+1) + 
     +	   SPBio(Ispec,IBlkNS,IBlkEW,Year+1)
         SumLarvae(Ispec,Year+1) = SumLarvae(Ispec,Year+1) + 
     +	   Larvae(Ispec,IBlkNS,IBlkEW,Year+1)	  
44445	CONTINUE 

C      Type = 3; Local density-dependence (at dispersal)
       IF (SrType(Ispec).EQ.3) THEN
        TotRec = 0
        DO 33310 IBlkNS = 1,NNSBlock
         DO 33310 IBlkEW = 1,NEWBlock
          Depl =  SPBio(Ispec,IBLkNS,IBlkEW,Year+1)/
     +                             SPBio(Ispec,IBLkNS,IBlkEW,MinPYr)
          Alpha = SRalpha(Ispec,IBlkNS,IBlkEW)
          Beta = SRbeta(Ispec,IBlkNS,IBlkEW)
          TotRec = TotRec + Depl / (Alpha + Beta * Depl)
33310   CONTINUE
        DO 33320 IBlkNS = 1,NNSBlock
         DO 33320 IBlkEW = 1,NEWBlock
          DO 33320 Isex = 1,Nsex
           N(Ispec,Isex,IBlkNs,IBlkEW,Year+1,0) =        !This line needs work.
     +      TotRec*SPBio(ISpec,IBlkNS,IBlkEW,Year+1)/
     +      SPBioTot(Ispec,Year+1)
33320   CONTINUE
       ENDIF
C       
C      Type = 4; Local and global density-dependence (at settlement)
       IF (SrType(Ispec).EQ.4) THEN
        Depl =  SPBioTot(Ispec,Year+1)/SPBioTot(Ispec,MinPYr)
        DO 33400 IBlkNS = 1,NNSBlock
         DO 33400 IBlkEW = 1,NEWBlock
          Depl1 =  SPBio(Ispec,IBLkNS,IBlkEW,Year+1)/
     +                             SPBio(Ispec,IBLkNS,IBlkEW,MinPYr)
          Alpha = SRalpha(Ispec,IBlkNS,IBlkEW)
          Beta = SRbeta(Ispec,IBlkNS,IBlkEW)
          DO 33400 Isex = 1,Nsex
           N(Ispec,Isex,IBlkNs,IBlkEW,Year+1,0) = 
     +        Depl / (Alpha + Beta * Depl1)
33400   CONTINUE
       ENDIF
C
C      Add noise
       DO 33500 IBlkNS = 1,NNSBlock
        DO 33500 IBlkEW = 1,NEWBlock
!-------------------------------------------------------------------------------
!     Adding a switch to do or exclude process error for the purpose of debugging. 7/19/11
         IF (DoProcessError.EQ.1) THEN
           DO 33510 Isex = 1,Nsex
          N(Ispec,Isex,IBlkNs,IBlkEW,Year+1,0) = 
     +      N(Ispec,Isex,IBlkNs,IBlkEW,Year+1,0) *
     +       RecRes(Isim,Ispec,Year,IBlkNS,IBlkEW)
33510      CONTINUE     
         ENDIF
!-------------------------------------------------------------------------------
         RecSum(Ispec,Year+1,IblkNS,IblkEW) = 
     +      N(ISpec,1,IBlkNS,IBlkEW,Year+1,0)
33500  CONTINUE
C 
33000 CONTINUE         
C

      RETURN
      END SUBROUTINE ProjYr

C
C ===========================================================================
C    
	  
      SUBROUTINE FindQvals()
C
C	This subroutine finds the value of Q that produces the specified value for depletion
C	this calls a function called "FIT" which uses as input (X=Qval,SS=output,Nspec,CareyGRD)
C
	USE ShareDataRatioSubs
      IMPLICIT NONE
C
C     Local variables
      REAL*8 X(100),SS
      INTEGER Ispec
      REAL*8 FUNK
      EXTERNAL FUNK
C
      MinOpt = 1
	CareyGRD = 1.2
      DO 10000 Ispec = 1,Nspec
       X(Ispec) = Qvalseed(Ispec)*NNSBlock*NEWBlock
10000 CONTINUE
      CALL FIT(X,SS,Nspec,CareyGRD)  
C     Carey: FIT finds the Qval that produces the specified value for depletion; 
C     FIT calls FUNK(X) to do this; X = QVal, so in the end SS = FUNK(X) just finds the value of the objective function (SS) that produced the
C     correct depletion.
      SS = FUNK(X)

      RETURN
      END SUBROUTINE FindQvals

C
C ===========================================================================
C    
	  
      SUBROUTINE AnnualProj
C
C	This is an important subroutine because it projects the population forward one year. 
C	This function is used within FUNK and it is the main function that advances the pop by one
C	year in the MSE
C      
	USE ShareDataRatioSubs
      IMPLICIT NONE
C
C     Local variables
      REAL*8 X(100),SS
	INTEGER :: II,III
      REAL*8 FUNK
      EXTERNAL FUNK
C
C	MinOpt = 2 to find the effort level that produces the quota.
      MinOpt = 2


C	ELSE
      DO 1000 II = 1,NeffortAreas
       X(II) = -10
1000  CONTINUE    
	CareyGRD = 1.2   

C     Problem with this is that Fit is not getting at all close to the quota sometimes. The seed has to be close to start with. CRM.
      CALL Fit(X,SS,NeffortAreas,CareyGRD)
C	ENDIF
      SS = FUNK(X)
C      WRITE(*,*) 'Minimum',SS,EXP(X(1))*qval(1)
C
      RETURN
      END SUBROUTINE AnnualProj

C
C ===========================================================================
C    
	  
      REAL*8 FUNCTION FUNK(X)
C
C 	This function is the objective function that projects forward to see how close depletion is
C 	to what is specified as initial depletion or final depletion
C 	This subroutine has different behavior depending on the value of MinOpt
C	MinOpt=1 --> find q for matching depletion value OUTPUTS depletion, FUNK
C	MinOpt=2 --> minimize catch to quota by changing effort value OUTPUTS FUNK
C	MinOpt=3 --> calculate the depletion value for a specific effort OUTPUTS FinDepRecord, FUNK
C	MinOpt=6 --> calls findMSY and calculate the corresponding OUTPUTS TotSpatialSuvNum, FUNK
C	MinOpt=8 --> weird. again find effort value that minimizes difference in initial or final depletion value compared to the input OUTPUTS depletion, FUNK
C     MinOpt=9 --> finds the equilibrium condition for each Isim and OUTPUTS AvgDepletion, AvgSSB, AvgCatch accross Nsim, FUNK
C     MinOpt=10 --> to calculate an estimated F(t) for the P* control rule from a given Nvector,estimated Selex, and specified M to yield a particular Catch biomass (a trial ABC or the OFL for instance) (I think deprecated as FUNK2 is used)
C     MinOpt=12 --> finds the F required for estimated SPR = 0.35     	   !!!!!!!!!!!!!!!! That is the one I need to use for my work
C     MinOpt=13 --> finds the F required for estimated SPR = 0.40 	       !!!!!!!!!!!!!!!! That is the one I need to use for my work
C     MinOpt=14 --> exactly the same as option 12					           !!!!!!!!!!!!!!!! That is the one I need to use for my work
C     MinOpt=15 --> exactly the same as option 13					           !!!!!!!!!!!!!!!! That is the one I need to use for my work
 


	USE ShareDataRatioSubs
      IMPLICIT NONE
      
C     Global variables
      REAL*8 X(100)
      REAL(8) :: Depl
      REAL(8), DIMENSION(Nspec) :: Depl2
C
C     Local variables
      REAL*8 Kvals(MaxSpec),CatPred(MaxSpec),Tmp,TheEffort(10),Penalty1
	REAL(8)::temp,temp1,temp2
      INTEGER :: Ispec,Iarea,II,Yr,IblkNS,IblkEW,IYr,Spec,Year,Ivec
      INTEGER :: Isex,Iage,Iconfig
	REAL(8) :: SumDepl,SumSSB,SumCatch,ISEED1,ISEED2,SEED1,SEED2

C      FOR MinOpt.EQ.10: (most of these defined globally and called outside of funk now because Ispec,Iconfig,Ivec loops are outside of funk)
!      REAL :: TheM,TheX
!      REAL, DIMENSION(1:Nsex,0:Maxage) :: TheNs
!      REAL, DIMENSION(1:Nsex,0:Maxage) :: TheSelex
      INTEGER :: TheNrands
      
C     Reset FUNK
      FUNK = 1.0E20 
C
C     How to proceed depends on MinOpt

C=======================================================================
C     MinOpt=1 means to find q value that leads to specific depletion level
      IF (MinOpt.EQ.1) THEN 
C
C      Copy values across
       DO 10000 Ispec = 1,Nspec
        IF (X(Ispec).LT.0) RETURN
        Qval(Ispec) = X(Ispec)
10000  CONTINUE       
C
C      Project forward 
C      Kvals = 10000.0
       Kvals = Kvec
       CALL ProjForward(Kvals)
C
C      Calculate the objective function
       FUNK = 0
       DO 11000 ISpec = 1,Nspec
!       Carey is changing this line:
!        Depl = SPBIoTot(Ispec,1)/SpBioTot(Ispec,MinPYr)*100
        Depl = SPBIoTot(Ispec,BaseYr+1)/SpBioTot(Ispec,MinPYr)*100
C        WRITE(*,*) Depl,InitDep(Ispec)
        IF (ABS(Depl-InitDep(Ispec)).GT.0.01) Funk = 
     +		Funk+(Depl-InitDep(Ispec))**2.0
11000  CONTINUE
C       WRITE(*,*) FUNK,Qval
C
      ENDIF

C=======================================================================
C     MinOpt=2 means to find E value that minimizes catch to quota
      IF (MinOpt.EQ.2) THEN 
      
       DO 22000 II = 1,NeffortAreas
        TheEffort(II) = EXP(X(II))
22000  CONTINUE         
       CALL projYr(TheProjYear,TheEffort)
C
C      Calculate the predicted catch - for MinOpt = 2, FUNK minimizes diffs between Catch and quota by changing total effort
       Funk = 0
       DO 20000 Ispec = 1,Nspec
        CatPred(Ispec) = 0
        DO 21000 Iarea = 1,Narea
         CatPred(Ispec) = CatPred(Ispec)+CatBio(Ispec,Iarea,TheProjYear)
21000   CONTINUE
        Tmp = CatPred(Ispec) - QuotaAct(Ispec)
        IF (ImplOpt.EQ.0) THEN
         Funk = Funk + Tmp**2.0
        ELSEIF (ImplOpt.EQ.1) THEN 
         IF (Tmp.GT.0) THEN
          Funk = Funk + 100*Tmp**4.0
         ELSE
          Funk = Funk + Tmp**2.0
         ENDIF
        ELSE
         IF (Tmp.LT.0) THEN
          Funk = Funk + 100*Tmp**4.0
         ELSE
          Funk = Funk + Tmp**2.0
         ENDIF
        ENDIF 
        
20000  CONTINUE
C
      ENDIF

C=======================================================================   
C    MinOpt=3 is just to calculate the depletion value for a specific effort level (TheEffort(1) = EXP(X(1) before now, TheEffort(1) = X(1))
      IF (MinOpt.EQ.3) THEN 
C      
       DO 31000 Yr = StartYr,EndYr
C      
C       Impose MPA
        IF (Yr.EQ.MPAYear) THEN
         DO 21100 Iarea = 1,NClosedAreas
          DO 21110 IblkNS = MinCNS(Iarea),MaxCNS(Iarea)
           DO 21110 IBlkEW = MinCEW(Iarea),MaxCEW(Iarea)
            Closed(IBlkNS,IBlkEW) = 0
21110     CONTINUE
21100    CONTINUE          
        ENDIF
C        TheEffort(1) = EXP(X(1))  !CHANGETHIS
        DO 31024 Ispec = 1,Nspec
          TheEffort(Ispec) = X(Ispec)
31024   CONTINUE
        CALL projYr(Yr,TheEffort)
31000  CONTINUE
       FUNK = 0
       DO 32024 Ispec = 1,Nspec       
       Depl2(Ispec) = SPBIoTot(Ispec,EndYr+1)/SpBioTot(Ispec,MinPYr)*100

       IF (ABS(Depl2(Ispec)-FinalDep(Ispec)).GT.0.01) Funk = 
     +      (Depl2(Ispec)-FinalDep(Ispec))**2.0
C       WRITE(*,*) FUNK,TheEffort(1),Depl
	FinDepRecord = Depl2(1)
32024  CONTINUE
      ENDIF

	
C Carey removed MinOpt.EQ.4 as it was related to the density ratio project
C Carey removed MinOpt.EQ.5 as it was to find optimal effort for DR project


C=======================================================================
C     MinOpt=6 is to call FindMSY
C     for funk assuming optimum already achieved:
	IF (MinOpt.EQ.6) THEN

	Steep(1) = X(1)
	CALL FindMSY()

	DO 37000 IYr = MinPYr,(MinPYr+200)
	DO 36000 Ispec = 1,NSpec
	  DO 36000 IBlkNS =1,NNSBlock
	    DO 36000 IBlkEW=1,NEWBlock
	      TotSpatialSuvNum(IYr) = TotSpatialSuvNum(IYr) + 
     +	    SpatialSuvNumSum(Ispec,IYr,IBlkNS,IBlkEW)
36000	CONTINUE
37000 CONTINUE

C	FUNK = (TotSpatialSuvNum(MinPYr+199)-(TotSpatialSuvNum(-50))/2)**2
	FUNK = (TotSpatialSuvNum(MinPYr+199)-(TotSpatialSuvNum(MinPYr))/2)**2
	ENDIF
C have to add up the exploitable numbers in the regular program (Prior to fishing and in most current years).
C	FUNK = (exploitable numbers - K(exploitable numbers)/2)*(exploitable numbers - K(exploitable numbers)/2)

C=======================================================================
C     MinOpt=8 finds effort value that minimizes difference in initial or final depletion value compared to the input
      IF (MinOpt.EQ.8) THEN 
C
C      Copy values across
        IF (X(1).LT.0) RETURN
C	  DO 33333 Yr = -19,BaseYr
    	  DO 33333 Yr = FirstFishYr,BaseYr
         TotEffort(Yr) = X(1)
33333  CONTINUE       
C
C      Project forward 
C       Kvals = 10000.0
      Kvals = Kvec
C     Project forward from MinPYr (-50 in original projects) to BaseYr with TotEffort(-19:BaseYr) = X(1)
       CALL ProjForward(Kvals)
C
C      Calculate the objective function
       FUNK = 0
       DO 11001 ISpec = 1,Nspec
        Depl = SPBioTot(Ispec,BaseYr)/SpBioTot(Ispec,MinPYr)*100
C        WRITE(*,*) Depl,InitDep(Ispec)
        IF (ABS(Depl-InitDep(Ispec)).GT.0.01) Funk = 
     +		Funk+(Depl-InitDep(Ispec))**2.0
11001  CONTINUE
C       WRITE(*,*) FUNK,Qval
	InitDeplRecord = Depl
C
      ENDIF

C=======================================================================
C     MinOpt=9 finds the equilibrium condition for each Isim and calculates AvgDepletion, AvgSSB, AvgCatch accross Nsim

	IF (MinOpt.EQ.9) THEN
	
	SumDepl = 0.0
	SumSSB = 0.0
	SumCatch = 0.0
	FUNK = 0.0

	AvgDepletion = 0.0
	AvgSSB = 0.0
	AvgCatch = 0.0

	DO 1000 Isim = 1,Nsim
         WRITE(*,*) "Isim = "
         WRITE(*,*) Isim


	  DO 1001 Yr = MinPYr,(MinPYr+150)
          TotEffort(Yr) = 0.00000000001 
1001	  CONTINUE !End Yr loop for initializing the model with no fishing effort (wipe out OM.dat entries)

C      Project forward to get the population dynamics started and moving
C          Kvals = 10000.0
           Kvals = Kvec
          CALL ProjForward(Kvals) !MinPYr to BaseYr (make BaseYr = -48)

!	  DO 1002 Year = (BaseYr + 1),(BaseYr + NProjYr)
        DO 1002 Year = MinPYr,(MinPYr+150)
          TheEffort(1) = EXP(X(1))
          CALL projYr(Year,TheEffort) 
		IF (WriteOpt.EQ.1) THEN
		  CALL Summary2(Year)
		ENDIF       
1002	  CONTINUE
C========================================================
	  FUNK = FUNK - CatBio(1,1,(MinPYr+149))
C=======================================================

C	Some extras:
       Depl = SPBIoTot(1,(MinPYr+149))/SpBioTot(1,MinPYr)*100  
	 SumDepl = SumDepl + Depl
       SumSSB = SumSSB + SPBioTot(1,(MinPYr+149))
	 SumCatch  = SumCatch + CatBio(1,1,(MinPYr+149))
1000	CONTINUE !end Isim loop

	AvgDepletion = SumDepl/Nsim
	AvgSSB = SumSSB/Nsim
	AvgCatch = SumCatch/Nsim

	ENDIF !MinOpt.EQ.9


C=======================================================
C     MinOpt=10 is to calculate an estimated F(t) for the P* control rule (I think deprecated as FUNK2 is used)
      IF (MinOpt.EQ.10) THEN
          !This MinOpt is to calculate an estimated F(t) for the P* control rule from
          !a given Nvector,estimated Selex, and specified M
          !to yield a particular Catch biomass (a trial ABC or the OFL for instance)

          FUNK = 0.0
 !         X(1) = 0.1 !This will be specified elsewhere later: Trial fishing mortality to produce the trial ABC value.
 !         Ivec = 1 !This will need to be changed later to a loop of some sort
         ! Ispec = 1 !This will need to be changed into a loop later as well: Ivec nested in Ispec loop.
 !         Iconfig = 0

 !         DO 4343 Ispec = 1,Nspec
 !          DO 4343 Iconfig = 0,Nconfigs
 !            TheM = EstM(Ispec,Iconfig) !Not vector specific now or with MCMC.
!             DO 4343 Ivec = 1,Nmvn
              !Let's do this in a loop because it's not going quite right here:
 !             DO 2117 Isex = 1,Nsex
 !!               DO 2117 Ivec = 1,Nmvn
 !                DO 2117 Iage = 0,Amax(Ispec)
 !                   TheNs(Isex,Iage) = exp(Nvec(Ispec,Isex,Ivec,
 !    +                                     Iage,Iconfig)) !Vector specific
!2117          CONTINUE !End Isex and Iage loop
!          
!              TheSelex = EstS(Ispec,1,:,:,Iconfig) !not vector specific at the moment, but with MCMC will be.
              TheX = X(1)
                !CalcCatchFromF(Ispec,Fval,Mval,TheNs,SelexVec,Ivec)
              CALL CalcCatchFromF(TheX,TheM,
     +        ThoseNs(:,:)
     +        ,TheSelex(:,:),MyWeight(:,:)) !The 1 will eventually be the ivec in a loop.

              FUNK = FUNK + (pCatBio-TrialQuota)**2.0 !Trial ABC = SeedQuota(Ispec) = CatBio from the most recent year of catch and assumes no implementation error, so is catch limit AND realized catch.
!4343      CONTINUE !End Ispec,Iconfig,Ivec loop
          
          
          ENDIF !MinOpt.EQ.10

          IF (MinOpt.EQ.11) THEN
           !This is just like MinOpt.EQ.10 for figuring out an F for a particular trial catch
              !but applied to the OM, not applied to a bunch of random, simulated vectors.
              FUNK = 0.0
!              Ispec = 1
!              TheM = M(Ispec)
!              TheNs = MyNs(Ispec,1:Nsex,0:Amax(Ispec))
!              TheSelex(1,:) =  Sel(Ispec,:)
!              TheSelex(2,:) = Sel(Ispec,:)
              TheX = X(1)
    
          CALL CalcCatchFromF(TheX,TheM,ThoseNs,
     +TheSelex(:,:),MyWeight(:,:)) !The 1 will eventually be the ispec in a loop.
         
!          IF (((pCatBio(Ispec)-500.0)**2).GT.0.01) THEN
             FUNK = FUNK + (pCatBio-TrialQuota)**2.0 !TrialABC = SeedQuota(Ispec) = CatBio from the most recent year of catch and assumes no implementation error, so is catch limit AND realized catch.
!          ENDIF
              
              ENDIF

C==========================================================
C     MinOpt=12 finds the F required for estimated SPR = 0.35
       IF (MinOpt.EQ.12) THEN
         !Find the estimated F required for estimated SPR = 0.35
         !MinOpt = 12 and MinOpt = 13 are sisters; 13 the same, but for SPR = 0.40.
         !Likewise, MinOpt 14 and 15 are for finding true F35, B35, F40, and B40
           FUNK = 0.0
         
         !TheM = M(Ispec)
         !TheX = X
         !RecVec = EstRec(1,1,1:81) !Define RecVec in FUNK
           
           !would have to make a change to Fec if Fecundity were estimated.
!         Ispec = 1 !This changed to be defined at a much higher level
  !     CALL AFSC_SPR(X,EstM(Ispec,Iconfig),Fec(1,0:Amax(Ispec)),
  !   +EstRec(1,1,1:GlobalProjYr),1,1) !do you have Iconfig right here in EstM and EstRec? what is 1,1 at the end?
        CALL AFSC_SPR(X,TheM,TheFecVec,TheRecVec,TheSelex(1,:),
     +ThePlusAge)
   
           
         FUNK = FUNK + ((EstSPR- 0.35)**2)*10000  !Check into what the dimensions of EstSPR are (??)  Maybe need to get rid of the dimensions.
       
         ENDIF

C==========================================================
C     MinOpt=13 finds the F required for estimated SPR = 0.40
      IF (MinOpt.EQ.13) THEN
        !Find the estimated F required for estimated SPR = 0.40
        FUNK = 0.0
        !Ispec = 1
        CALL AFSC_SPR(X,TheM,TheFecVec,
     +  TheRecVec,TheSelex(1,:),ThePlusAge)   !This is halfway converted! Haven't converted the outer shell.
        
        FUNK = FUNK + ((EstSPR- 0.40)**2)*10000       
      ENDIF

C==========================================================
C     MinOpt=14 is the same as 12
      IF (MinOpt.EQ.14) THEN
        FUNK = 0.0
 
        CALL AFSC_SPR(X,TheM,TheFecVec,
     +  TheRecVec,TheSelex(1,:),ThePlusAge)  
        
        FUNK = FUNK + ((EstSPR- 0.35)**2)*10000       
      ENDIF            

C==========================================================
C     MinOpt=15 is the same as 13 
      IF (MinOpt.EQ.15) THEN
        FUNK = 0.0

        CALL AFSC_SPR(X,TheM,TheFecVec,
     +  TheRecVec,TheSelex(1,:),ThePlusAge)   
        
        FUNK = FUNK + ((EstSPR- 0.40)**2)*10000       
      ENDIF              
       
      RETURN
      END FUNCTION FUNK

C
C ===========================================================================
C    
	  
      SUBROUTINE GenProcErr()
C
C This subroutine generates the process errors that will be needed for the
C current simulation
C
	USE ShareDataRatioSubs
      IMPLICIT NONE
C
C     Local variables
      INTEGER Year,Ispec,IBlkNS,IBlkEW
      REAL*8 Error1,Error2
      REAL*8 XNORM
      EXTERNAL XNORM
C
C     Generate the recruitment error terms
      IF (DoProcessError.EQ.1) THEN
      DO 10000 Ispec = 1,NSpec
       DO 10000 Year = MinPYr,MaxPYr
        !global process error on recruitment (Error1)
        Error1 = XNORM(1,0.0d0,SigmaR1(Ispec),ISEEDP)
     +                             -SigmaR1(ISpec)**2.0/2.0
        DO 11000 IBlKNS = 1,NNSBlock
         DO 11000 IBlkEW = 1,NEWBlock
         !Local process error on recruitment (Error2)
          Error2 = XNORM(1,0.0d0,SigmaR2(Ispec),ISEEDP)
     +                             -SigmaR2(ISpec)**2.0/2.0

           RecRes(ISim,ISpec,Year,IBlkNS,IBlkEW) = EXP(Error1+Error2)

11000   CONTINUE         
10000 CONTINUE       
         ELSE
           RecRes = 1
         ENDIF
C
C     Generate the catchability error terms
      !Carey did this 7/19/11
      !Carey 1/14/14: does this need to be fixed such that there is an ISEED for spatial catchability process error
      ! and this is Error3, not Error2 (as Error2 already refers to local process error on recruitment and ISEEDP is being used multiple times)
      IF (DoProcessError.EQ.1) THEN
       DO 20000 Ispec = 1,NSpec
        DO 21000 IBlKNS = 1,NNSBlock
         DO 21000 IBlkEW = 1,NEWBlock
          Error2 = XNORM(1,0.0d0,SigmaSpat(Ispec),ISEEDP)
     +                             -SigmaSpat(ISpec)**2.0/2.0
         RecQ(ISim,ISpec,IBlkNS,IBlkEW) = EXP(Error2)
21000   CONTINUE         
20000  CONTINUE       
      ELSE
        RecQ = 1
      ENDIF
      RETURN
      END SUBROUTINE GenProcErr

C
C ===========================================================================
C    
	  
      SUBROUTINE Summary()
C
C This subroutine dumps various outputs
C	1. Age, Fecundity, and weight_mid_year
C
      USE ShareDataRatioSubs
      IMPLICIT NONE
C
C     Local variables
      INTEGER Ispec,Isex,Itime,Iage,IBlkNS,IBlkEW
C
C     Open file
      OPEN(UNIT=14,FILE="results/SUMM.OUT")
      DO 10000 Ispec = 1,Nspec
       DO 10000 Isex = 1,Nsex
        WRITE(14,601) Ispec,Isex
        WRITE(14,602)
        DO 11000 iage = 0,Amax(Ispec)
         WRITE(14,603) Iage,Fec(Ispec,Iage),
     +         (Wght(Ispec,Isex,Iage,Itime),Itime=0,Ntime)
11000   CONTINUE        
10000 CONTINUE        
      CLOSE(14)   
      RETURN
601   FORMAT(1x\,1x,I2,1x,I2)
602   FORMAT(1x,"Age Fecundity   Weight")
603   FORMAT(1x,I3,1x,F9.4,1x,100(F8.3,1x))

      END SUBROUTINE Summary

C
C ===========================================================================
C    
	  
      SUBROUTINE Summary2(Year)
C
C This subroutine dumps results to a file
C	1. Year, Effort
C 	2. Species, Year, Relative Spbio (to the first year)
C 	3. Species, Year, Relative SpBio (to the max value)
C 	4. Species, Year, SpBio
C 	5. Species, Year, Catch
C 	6. Species, Year, Recruitment 
C 	7. Species, Year, SpBio in MPA 
C 	8. Species, Year, Age, N at age (for both sex)
C 	9. Species, Year, N at age 0 (for both sex)
C 	10. Species, Year, SpBio in MPA 
C 	11. Species, Year, SpBio in MPA 
C 	12. Year, Total catch, maximum effort, quota per species
C 	13. age, selectivity at age by species
C 	14. some spatial model input variables
C 	15. species, block, N at equilibrium (not useful for non spatial model)
C 	16. species, block, N at equilibrium (not useful for non spatial model)
C 	17. Year, SpBio, Relative SpBio, Total catch, Quota, Total effort
C 	18. Year, net effort / maximum effort
C 	19. Year, Nsum, Catch
C 	20. Year, species, block, larvae, recruitments
C 	21. Isim, year, species, CPUE
C 	22. Isim, year, species, area, Survey estimates
C
      USE ShareDataRatioSubs
      IMPLICIT NONE
C
C     Local variables 
C
C     Local variables
      INTEGER Ispec,Year,IBlkNS,IBlkEW,JBlkNS,JBlkEW,Iarea,II,Isex,Iage
      REAL*8 TotCat,Val(100),MaxRec,MaxSSB,MaxEff,TotEff,IMPA

C	Write out the total effort not divided by blocks (according to Carey)
C GlobalEffort def'n has to be changed if spatial effort blocks become activated.
	WRITE(38,'(I5,1x,F30.7)') Year, GlobalEffort(Year)
C
C     Total effort
      TotEff = 0
      DO 9000 IBLkEW = 1,NEWBlock
       DO 9000 IBlkNS = 1,NNSBlock
        TotEff = TotEff + NetEffort(IBlkNS,IBlkEW,Year)
9000  CONTINUE        
      IF (TotEff.Gt.90000) TotEff = 90000
C
C     Open file
      DO 10000 Ispec = 1,Nspec
       TotCat = 0
       DO 11000 Iarea = 1,Narea
         Totcat = Totcat + CatBio(Ispec,Iarea,Year)
11000  CONTINUE          
       Val(Ispec*4-3) = SpBioTot(ISpec,Year)
       Val(Ispec*4-2) = SPBioTot(Ispec,Year)/SPBioTot(Ispec,MinPYr)*100
       IF (totCat.GT.90000) Totcat = 90000
       Val(ispec*4-1) = TotCat
       Val(ispec*4  ) = QuotaAct(ISpec)
       IF (Val(ispec*4).GT.90000) val(ispec*4) = 90000
       
       MaxRec = 0
       MaxSSB = 0
       DO 12000 IBLkEW = 1,NEWBlock
        DO 12000 IBlkNS = 1,NNSBlock
         IF (SRR0(Ispec,IBlkNS,IBlkEW).GT.MaxRec) 
     +    MaxRec = SRR0(Ispec,IBlkNS,IBlkEW)
         IF (SPBio(Ispec,IBlkNs,IBlkEW,MinPyr).GT.MaxSSB)
     +    MaxSSB = SPBio(Ispec,IBlkNs,IBlkEW,MinPyr) 
12000  CONTINUE     
       WRITE(15,601) ISpec,Year,
     +      ((SpBio(Ispec,IBlkNS,IBlkEW,Year)
     +       /SPBio(Ispec,IBlkNs,IBlkEW,MinPyr)*100,IBlkEW=1,NEWBlock),
     +        IBlkNS=1,NNSBlock)
       WRITE(17,601) ISpec,Year,
     +      ((SpBio(Ispec,IBlkNS,IBlkEW,Year)/MaxSSB*100,
     +        IBlkEW=1,NEWBlock),IBlkNS=1,NNSBlock)
       WRITE(34,601) ISpec,Year,
     +      ((SpBio(Ispec,IBlkNS,IBlkEW,Year),
     +        IBlkEW=1,NEWBlock),IBlkNS=1,NNSBlock)
       WRITE(35,601) Ispec,Year,
     +      ((CatArea(Ispec,IBlkNS,IBlkEW,Year),
     +        IBlkEW=1,NEWBlock),IBlkNS=1,NNSBlock)
       WRITE(16,601) Ispec,Year,
     +      ((RecSum(Ispec,Year,IBlkNS,IBlkEW)/MaxRec*100,
     +        IBlkEW=1,NEWBlock),IBlkNS=1,NNSBlock)
C Carey added this one
	 SpBioMPA(Ispec,:,Year) = 0
	DO 95000 IBlkNS = 1,NNSBlock
	  DO 95000 IBlkEW = 1,NEWBlock
	IF (Closed(IBlkNS,IBlkEW).EQ.0) THEN
	SpBioMPA(Ispec,2,Year) = SpBioMPA(Ispec,2,Year) + 
     +	SpBio(Ispec,IBlkNS,IBlkEW,Year)
	ENDIF
	IF (Closed(IBlkNS,IBlkEW).NE.0) THEN
	SpBioMPA(Ispec,1,Year) = SpBioMPA(Ispec,1,Year) + 
     +	SpBio(Ispec,IBlkNS,IBlkEW,Year)
	ENDIF
95000	CONTINUE
      WRITE(41,618) ISpec,Year,
     +      (SpBioMPA(Ispec,IMPA,Year),IMPA=1,2)


C	 IF (Year.EQ.MinPYr) THEN  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C       IF (Year.EQ.0.OR.Year.EQ.15.OR.Year.EQ.30.OR.Year.EQ.45.
C     +	 OR.Year.EQ.60.OR.Year.EQ.-50) THEN !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        DO 12001 Iage = 0,Amax(ISpec)
         WRITE(26,613) ISpec,Year,Iage,
     +		(((N(Ispec,Isex,IBlkNS,IBlkEW,Year,Iage),
     +		 IBlkEW=1,NEWBlock),IBlkNS=1,NNSBlock),Isex = 1,Nsex)
12001   CONTINUE
         WRITE(260,624) ISpec,Year,
     +		(((N(Ispec,Isex,IBlkNS,IBlkEW,Year,0),
     +		 IBlkEW=1,NEWBlock),IBlkNS=1,NNSBlock),Isex = 1,Nsex)
C       ENDIF  

       WRITE(27,614) Year,TotCat,MaxEff,QuotaAct(Ispec)

C      Carey's additional output
C      Write out selectivity at age to a file to compare to output from cab (= selA.dat).
       IF (Year.EQ.MinPYr) THEN
        DO 94012 Iage = 0,Amax(Ispec)
         WRITE(28,699) 1,Iage, Sel(Ispec,Iage)
94012   CONTINUE
        DO 94013 Iage = 0,Amax(Ispec)
         WRITE(28,699) 2,Iage, Sel(Ispec,Iage)
94013   CONTINUE
       ENDIF
Cfs
       IF (Year.EQ.MinPYr) THEN
        DO 16000 JBlkNS = 1,NNSBlock
         DO 16100 JBlkEW = 1,NEWBlock 
           WRITE(19,612) ISpec,JBlkNS,JBlkEW,
     +				((MixMat(Ispec,IBlkNS,IBlkEW,JBlkNS,JBlkEW)
     +              ,IBlkEW=1,NEWBlock),IBlkNS=1,NNSBlock)
16100    CONTINUE 
16000   CONTINUE
        DO 17000 IBlkNS = 1,NNSBlock
           WRITE(33,615) ISpec,1,IBlkNS,
     +          (Neq(Ispec,1,IBlkNS,IBlkEW,0),IBlkEW=1,NEWBlock)                   
           WRITE(33,615) ISpec,2,IBlkNS,
     +          (Neq(Ispec,2,IBlkNS,IBlkEW,0),IBlkEW=1,NEWBlock)                   
17000    CONTINUE 
       ENDIF

10000 CONTINUE
      WRITE(14,600) Year,(Val(II),II=1,Nspec*4),TotEff
      MaxEff = -1
      DO 21000 IBLkEW = 1,NEWBlock
       DO 21000 IBlkNS = 1,NNSBlock
        IF (NetEffort(IBlkNS,IBlkEW,Year).GT.MaxEff) 
     +   MaxEff = NetEffort(IBlkNS,IBlkEW,Year)
21000 CONTINUE     
      IF (MaxEff.LE.0) MaxEff = 0.001
      WRITE(18,601) -1,Year,
     +      ((NetEffort(IBlkNS,IBlkEW,Year)/MaxEff*100,
     +        IBlkEW=1,NEWBlock),IBlkNS=1,NNSBlock)
      WRITE(32,602) -1,Year,
     +      ((NetEffort(IBlkNS,IBlkEW,Year),
     +        IBlkEW=1,NEWBlock),IBlkNS=1,NNSBlock)
C
C	Carey is adding a single N and Catch in numbers file
	DO 22000 Ispec = 1,NSpec
	  DO 22000 Iage = 0,MaxAge
	    DO 22000 Isex = 1,2
	      CatNumSum(Year) = CatNumSum(Year) + 
     +		  CatNum(Ispec,1,Year,Isex,Iage)
	      DO 22100 IBlkNS = 1,NNSBlock
	        DO 22100 IBlkEW = 1,NEWBlock
                NSum(Year) = NSum(Year)+
     +			  N(Ispec,Isex,IBlkNS,IBlkEW,Year,Iage)
22100 CONTINUE
22000	CONTINUE



	WRITE(40,620) Year,NSum(Year),CatNumSum(Year)

C	Carey is adding a file for larvae and recruits by cell
	DO 45678 ISpec = 1,NSpec
	  DO 45678 IBlkNS = 1,NNSBlock
	    DO 45678 IBlkEW = 1,NEWBlock
           TotRecruits(ISpec,IBlkNS,IBlkEW,Year) = 
     +		 N(ISpec,1,IBlkNS,IBlkEW,Year,0) +
     +	     N(ISpec,2,IBlkNS,IBlkEW,Year,0)

	     WRITE(49,621) Year,ISpec,IBlkNS,IBlkEW,
     +	     Larvae(Ispec,IBlkNS,IBlkEW,Year), 
     +         TotRecruits(ISpec,IBlkNS,IBlkEW,Year)
45678	CONTINUE

C     Carey is adding 2 files, one for catch cpue and one for survey cpue
      DO 9878 Ispec = 1,Nspec
   !     DO 9878 Iarea = 0,Narea
      WRITE(1500,623) Isim,Year,Ispec,(CPUE(Ispec,Iarea,Year),Iarea=
     +   0,Narea)
      WRITE(1600,623) Isim,Year,Ispec,Iarea,(SurveyEst(Ispec,Iarea,
     +  Year),Iarea=0,Narea)
9878  CONTINUE

      RETURN
600   FORMAT(1x,I4,10(3(1x,F8.2)))
601   FORMAT(1x,I2,1x,I4,1x,200(F8.2,1x))
602   FORMAT(1x,I2,1x,I4,1x,200(F8.4,1x))
612   FORMAT(1x,I2,1x,I2,1x,I2,1x,200(F8.6,1x))
613   FORMAT(1x,I2,1x,I4,1x,I4,1x,200(F20.12,1x))
614   FORMAT(1x,I4,1x,F8.2,1x,F6.3,1x,F8.2,1x)
615   FORMAT(1x,I2,1x,I2,1x,I2,1x,200(F8.6,1x))
616   FORMAT(1x,I2,1x,I4,1x,3(F10.6),1x,2(F10.3),1x,2(I4))
699   FORMAT(I2,1x,I2,1x,F8.4)
617   FORMAT(1x,I4,1x,2(F10.4))
618   FORMAT(1x,I2,1x,I4,1x,2(F12.2,1x))
619	FORMAT(2(1x,I4),1x,F8.6,1x,F10.3,1x,F8.6,1x,2(F8.6,1x),F10.3)
620   FORMAT(1x,I4,1x,2(F15.4))
621	FORMAT(1x,I4,1x,3(I2,1x),2(F25.7,1x))
622   FORMAT(1x,I2,1x,I4,1x,F10.3,1x,I4)
!Carey added to output CPUE data (Catch and Survey CPUE)
623   FORMAT(2(I4,1x),2(I2,1x),F18.8)
624   FORMAT(1x,I2,1x,I4,1x,200(F20.12,1x))
!Sim,Yr,Ispec,Iarea,CPUE
      END SUBROUTINE Summary2

C
C ===========================================================================
C    
	  
      SUBROUTINE Summary3(Year)
C
C This subroutine dumps results to a file
C 	1. Species, Year, Relative Spbio (to the first year)
C 	2. Species, Year, Relative SpBio (to the max value)
C 	3. Species, Year, SpBio
C 	4. Species, Year, relative Recruitment to the max recruits
C
	USE ShareDataRatioSubs
      IMPLICIT NONE
C
C     Local variables 
C
C     Local variables
      INTEGER Ispec,Year,IBlkNS,IBlkEW,Iarea,II
      REAL*8 TotCat,Val(0:100),MaxRec,MaxSSB
C
C     Open file
      DO 10000 Ispec = 1,Nspec
       TotCat = 0
       DO 11000 Iarea = 1,Narea
         Totcat = Totcat + CatBio(Ispec,Iarea,Year)
11000  CONTINUE          
       Val(Ispec*4-3) = SpBioTot(ISpec,Year)
       Val(Ispec*4-2) = SPBioTot(Ispec,Year)/SPBioTot(Ispec,MinPYr)*100
       Val(ispec*4-1) = -1
       Val(ispec*4  ) = -1
       
       MaxRec = 0
       MaxSSB = 0
       DO 12000 IBLkEW = 1,NEWBlock
        DO 12000 IBlkNS = 1,NNSBlock
         IF (SRR0(Ispec,IBlkNS,IBlkEW).GT.MaxRec) 
     +    MaxRec = SRR0(Ispec,IBlkNS,IBlkEW)
         IF (SPBio(Ispec,IBlkNs,IBlkEW,MinPyr).GT.MaxSSB)
     +    MaxSSB = SPBio(Ispec,IBlkNs,IBlkEW,MinPyr) 
12000  CONTINUE     
       WRITE(15,601) ISpec,Year,
     +      ((SpBio(Ispec,IBlkNS,IBlkEW,Year)
     +       /SPBio(Ispec,IBlkNs,IBlkEW,MinPyr)*100,IBlkEW=1,NEWBlock),
     +        IBlkNS=1,NNSBlock)
       WRITE(17,601) ISpec,Year,
     +      ((SpBio(Ispec,IBlkNS,IBlkEW,Year)/MaxSSB*100,
     +        IBlkEW=1,NEWBlock),IBlkNS=1,NNSBlock)
       WRITE(34,601) ISpec,Year,
     +      ((SpBio(Ispec,IBlkNS,IBlkEW,Year),
     +        IBlkEW=1,NEWBlock),IBlkNS=1,NNSBlock)
       WRITE(16,601) ISpec,Year,
     +      ((RecSum(Ispec,Year,IBlkNS,IBlkEW)/MaxRec*100,
     +        IBlkEW=1,NEWBlock),IBlkNS=1,NNSBlock)

10000 CONTINUE
      WRITE(14,600) Year,(Val(II),II=1,Nspec*4),0.0d0
      RETURN
600   FORMAT(1x,I4,10(3(1x,F8.2)))
601   FORMAT(1x,I2,1x,I4,1x,200(F8.2,1x))
602   FORMAT(1x,I2,1x,I4,1x,200(F8.4,1x))
616   FORMAT(1x,I2,1x,I4,1x,3(F10.6),1x,2(F10.3),1x,2(I4))

      END SUBROUTINE Summary3

C
C ===========================================================================
C    
	  
      SUBROUTINE GenData(Yr)
C
C  This subroutine creates the data used in the assessment model = CAB.exe
C  Steps:
C 	1. Catch is observed without errors
C 	2. Fishery catch rate index with lognormal errors (Kotaro modified it so that is is lognormally distributed)
C 	3. Survey index of abundance with lognormal errors 
C 	4. Extract the commercial catch at age (fishery age composition data)
C 	5. Extract the survey catch at age (survey age composition data)
C      
	USE ShareDataRatioSubs
      IMPLICIT NONE
C
C     Global variables
      INTEGER Yr
C
C     Local variables
      INTEGER Ispec,Iarea,ISex,Iage,EffSS
      REAL*8 Total,Error,C(2000),Prop(2000),Top,Bot
      REAL*8 XNORM,RAN2
      EXTERNAL XNORM,RAN2
C
C     Extract the catch by assessment area
      DO 10000 Ispec = 1,Nspec
       Obscatch(Ispec,0,Yr) = 0
       DO 10000 Iarea = 1,Narea
        ObsCatch(Ispec,Iarea,Yr) = CatBio(Ispec,Iarea,Yr)
        ObsCatch(Ispec,0,Yr) = ObsCatch(Ispec,0,Yr) + 
     +    CatBio(Ispec,Iarea,Yr)
10000 CONTINUE       
C
      !------------------------------------------------------------
C      
C     Extract the catch-rate index    
C     Kotaro modified a few things so that the CPUE (catch) follows a lognormal distribution in EACH area
C     There was a controversy on whether to add the error on the catch (numerator) or the effort (denomiator) (discussion between CRM, AEP, RM)
C     However, we decided to add error on catch because the technical interaction model will be essentially used with ONE strata (area)
C 	  and the EM assumes that fishery CPUE is lognormally distributed
C     Need to modify this code if one want to use multiple area with different assumption (right now lognormality is assumed for each strata)
C
      DO 20000 Ispec = 1,Nspec
       Top = 0
       Bot = 0
       DO 21000 Iarea = 1,Narea
        Error = XNORM(2,0.0d0,CPUECV,ISEEDO)-CPUECV**2.0/2.0
        IF (CatBIo(Ispec,Iarea,Yr).GT.0.AND.
     +     NoCPUEMPA.EQ.0.AND.Effort(Iarea,Yr).GT.0) THEN
         IF (IsDeterm.EQ.1) THEN 
	     Error = 0
	   ENDIF
    !     CPUE(Ispec,Iarea,Yr) = TruExp(Ispec,Iarea,Yr)
         CPUE(Ispec,Iarea,Yr) = CatBio(Ispec,Iarea,Yr)/Effort(Iarea,Yr)
     +    * EXP(Error)
    !    CPUE(Ispec,Iarea,Yr) = CPUE(ispec,Iarea,Yr) * EXP(Error)
         CVCPUE(Ispec,Iarea,Yr) = CPUECV

       !  Top = Top + CatBio(Ispec,Iarea,Yr)
       !  Bot = Bot + EXP(Error)/Effort(Iarea,Yr)
       !  Top = Top + CatBio(Ispec,Iarea,Yr)
       !  Bot = Bot + EXP(Error)*Effort(Iarea,Yr)
         Top = Top + CPUE(Ispec,Iarea,Yr)
        ELSE
         CPUE(ISpec,Iarea,Yr) = -1
         CVCPUE(Ispec,Iarea,YR) = -1
        ENDIF  
21000  CONTINUE        
       IF (Top.Gt.0) THEN
  !      CPUE(ISpec,0,Yr) = Top*Bot
  !      CPUE(ISpec,0,Yr) = Top/Bot
         CPUE(ISpec,0,Yr) = Top/Narea       ! this calculates the average CPUE accross areas
         CVCPUE(Ispec,0,YR) = CPUECV
       ELSE
         CPUE(ISpec,0,Yr) = -1
         CVCPUE(Ispec,0,YR) = -1
       ENDIF
20000    CONTINUE        
C
      !------------------------------------------------------------
C        
C     Extract the survey index
      DO 30000 Ispec = 1,Nspec
       Top = 0
       DO 31000 Iarea = 1,Narea
        Error = XNORM(2,0.0d0,SurvCV,ISEEDO)-SurvCV**2.0/2.0
        IF (Yr.GE.FirstSurvYr.AND.MOD(Yr-FirstSurvYr,FreqSurv).EQ.0) 
     +   THEN
         IF (IsDeterm.EQ.1) Error = 0
         SurveyEst(Ispec,Iarea,Yr) = TruSuv(Ispec,Iarea,Yr)*EXP(Error)
         SurveyCV(Ispec,Iarea,Yr) = SurvCV
         Top = Top + SurveyEst(Ispec,Iarea,Yr)
        ELSE
         SurveyEst(Ispec,Iarea,Yr) = -1
         SurveyCV(ISpec,Iarea,Yr) = -1
        ENDIF 
31000  CONTINUE        
       IF (Top.GT.0) THEN
    !    SurveyEst(Ispec,0,Yr) = Top !Here Top is the sum of SurveyEst over areas
         SurveyEst(Ispec,0,Yr) = Top/Narea  ! this calculates the average CPUE accross areas
         SurveyCV(ISpec,0,Yr) = SurvCV
       ELSE
        SurveyEst(Ispec,0,Yr) = -1
        SurveyCV(ISpec,0,Yr) = -1
       ENDIF
30000 CONTINUE        
C
      !------------------------------------------------------------
C           
C     Extract the commercial catch-at-age data
      DO 40000 Ispec = 1,Nspec
       DO 40000 Isex = 1,Nsex
        ComSS(ispec,0,Yr) = 0
        DO 40100 Iage = 0,AmaxEst(Ispec)
         CNObs(Ispec,0,Yr,Isex,Iage) = 0
40100   CONTINUE          
        DO 41000 Iarea = 1,Narea
C   
C        Renormalize 
         Total = 0
         DO 41100 Iage = 0,AmaxEst(Ispec)
          Total = Total + CatNum(Ispec,Iarea,Yr,Isex,Iage)
41100    CONTINUE          
         IF (Total.EQ.0) THEN
          ComSS(Ispec,Iarea,Yr) = -1
          Error = RAN2(ISEEDO)
          DO 41200 Iage = 0,AmaxEst(Ispec)
           CNObs(Ispec,Iarea,Yr,Isex,Iage) = 0
41200     CONTINUE          
         ELSE
C          IF (Yr.LE.0) THEN !Carey changed this to MPAYear.
          IF (Yr.LE.MPAYear) THEN 
           EffSS = EffectSampCatch1
          ELSE
           EffSS = EffectSampCatch2
          ENDIF 
          ComSS(Ispec,Iarea,Yr) = EffSS
          ComSS(Ispec,0,Yr) = ComSS(Ispec,0,Yr) + EffSS
          Prop = 0
          C = 0
          DO 41310 Iage = 0,AmaxEst(Ispec)
           C(Iage+1) = CatNum(Ispec,Iarea,Yr,Isex,Iage) / Total !Total is the catch in numbers summed over ages (specific to species, area, yr, sex)
           Prop(Iage+1) = C(Iage+1) !Proportion at age
41310     CONTINUE          
          IF (ISDeterm.EQ.1) THEN
           Error = Ran2(ISEEDO)
          ELSE
           C = 0
           Call GenMul(Prop,ComSS(Ispec,Iarea,Yr),C,ISEEDO,.FALSE.) !This takes the proportions and the effective sample size and makes C contain effective catch at age (sampled from the multinomial) that adds to the effective sample size (if adding over ages)
           DO 41320 Iage = 0,AmaxEst(Ispec)
            C(Iage+1) = C(Iage+1)/FLOAT(ComSS(Ispec,Iarea,Yr)) !Now C contains the multinomially sampled proportions at age (adding to 1)
41320      CONTINUE
          ENDIF  
          DO 41330 Iage = 0,AmaxEst(Ispec) !This loop makes the aggregated CNObs from the CNObs made above for each assessment area (Iarea = 0)
           CNObs(Ispec,Iarea,Yr,Isex,Iage) = C(Iage+1)
           CNObs(Ispec,0,Yr,Isex,Iage) = CNObs(Ispec,0,Yr,Isex,Iage) + 
     +      C(Iage+1)*ComSS(Ispec,Iarea,Yr)
41330     CONTINUE          
         ENDIF
C
41000   CONTINUE
        IF (ComSS(Ispec,0,Yr).EQ.0) THEN
         ComSS(Ispec,0,Yr) = -1
        ELSE
         DO 42000 Iage = 0,AmaxEst(Ispec) !And this switches everything to proportions for the aggregated CNObs (Iarea = 0)
          CNObs(Ispec,0,Yr,Isex,Iage) = CNObs(ispec,0,Yr,Isex,Iage) / 
     +     ComSS(Ispec,0,Yr)
42000    CONTINUE        
        ENDIF
40000 CONTINUE         
C
      !------------------------------------------------------------
C
C     Extract the survey catch-at-age data (Carey's addition to this comment: this looks like a repeat of the process for effective catch at age.)
      DO 50000 Ispec = 1,Nspec
       DO 50000 Isex = 1,Nsex
        SuvSS(ispec,0,Yr) = 0
        DO 50100 Iage = 0,AmaxEst(Ispec)
         SNObs(Ispec,0,Yr,Isex,Iage) = 0
50100   CONTINUE          
        DO 51000 Iarea = 1,Narea
C
C        Renormalize 
         Total = 0
         DO 51100 Iage = 0,AmaxEst(Ispec)
          Total = Total + SuvNum(Ispec,Iarea,Yr,Isex,Iage)
51100    CONTINUE          
         IF (Yr.LT.FirstSurvYr.OR.MOD(Yr-FirstSurvYr,FreqSurv).NE.0.OR.
     +     SurveyEst(Ispec,Iarea,Yr).LE.0) THEN
          SuvSS(Ispec,Iarea,Yr) = -1
          Error = RAN2(ISEEDO)
          DO 51200 Iage = 0,AmaxEst(Ispec)
           SNObs(Ispec,Iarea,Yr,Ispec,Iage) = 0
51200     CONTINUE
         ELSE  
C          IF (Yr.LE.0) THEN
          IF (Yr.LE.MPAYear) THEN
           EffSS = EffectSampSurvey1
          ELSE
           EffSS = EffectSampSurvey2
          ENDIF 
          SuvSS(Ispec,Iarea,Yr) = EffSS
          SuvSS(Ispec,0,Yr) = SuvSS(Ispec,0,Yr) + EffSS
          Prop = 0
          C = 0
          DO 51310 Iage = 0,AmaxEst(ispec)
           C(Iage+1) = SuvNum(Ispec,Iarea,Yr,Isex,Iage) / Total
           Prop(Iage+1) = C(Iage+1)
51310     CONTINUE
          IF (IsDeterm.EQ.1) THEN
           Error = Ran2(ISEEDO)
          ELSE
           C = 0
           Call GenMul(Prop,SuvSS(Ispec,Iarea,Yr),C,ISEEDO,.FALSE.)
           DO 51320 Iage = 0,AmaxEst(Ispec)
            C(Iage+1) = C(Iage+1)/FLOAT(SuvSS(Ispec,Iarea,Yr)) 
51320      CONTINUE
          ENDIF  
          DO 51330 Iage = 0,AmaxEst(Ispec)
           SNObs(Ispec,Iarea,Yr,Isex,Iage) = C(Iage+1)
           SNObs(Ispec,0,Yr,Isex,Iage) = SNObs(Ispec,0,Yr,Isex,Iage) + 
     +      C(Iage+1)*SuvSS(Ispec,Iarea,Yr)
51330     CONTINUE          
         ENDIF
C
51000   CONTINUE
        IF (SuvSS(Ispec,0,Yr).EQ.0) THEN
         SuvSS(Ispec,0,Yr) = -1
        ELSE
         DO 52000 Iage = 0,AmaxEst(Ispec)
          SNObs(Ispec,0,Yr,Isex,Iage) = SNObs(ispec,0,Yr,Isex,Iage) / 
     +     SuvSS(Ispec,0,Yr)
52000    CONTINUE        

        ENDIF
50000 CONTINUE         
C
      RETURN
      END SUBROUTINE GenData

C
C ===========================================================================
C    
	  
      SUBROUTINE DoCab2(Ispec,Iarea1,Iarea2,Year1,Year2,MinA,INDC)
C
C This subroutine sets up a run using CAB.TPL and conducts the assessment for a single configuration run
C Steps:
C 	1. Set the CAB.dat file
C 	2. Set the CAB.ctl file --> this is influenced by the model specification
C 	3. Call CAB.exe
C 	4. Read the cab report file 
C 	5. Read the cab par file
C 	6. then if the cab.cor file exists, then run extracts the number at age based on MVN distribution for the projection
C 	7. Use of "system" command to run other program. Here ".bat" file
C      
	USE ShareDataRatioSubs
      IMPLICIT NONE
C
C     Global variables
      INTEGER ISpec,Iarea1,Iarea2,Year1,Year2,MinA,Istuff,Iconfig
      CHARACTER*1 INDC
C
C     Local variables
      REAL*8 Pred(6),L50a,LDiff,SpTru,B0Tru,ExpTru
      CHARACTER*8 TsCode1,TsCode2
      CHARACTER*21 TstOut
      INTEGER Nrel1(10),Nrel2(10),Nage1(10),Nage2(10),EstFleet,Iarea
      INTEGER Iage,Isex,Itime,Yr,Ilen,II,Ioff,Iarea1a,Iarea2a,OutUnit
	INTEGER Addone,CArea1
      CHARACTER*2 cIsim,cIYr,cIconfig
      CHARACTER*21 cstring

      
C      WRITE(*,*) INDC
C	Carey added the following: this adds 1 to Narea1 in array subscripts to stay in bounds if Narea1 = 0 (this happens when the single assessment method is being performed)
C	Under the old program, not adding the 1 stops the program.
 	Addone = 0
	IF (Iarea1.LT.1.) Addone = 1
	Carea1 = Iarea1 + Addone
C	Carey is done adding new code, Andre's code continues from here

C     Set up weight and fecundity
C     Carey: Itime puts you in the middle of the year
      Itime = (Ntime-1)/2 + 1
 
C-----------------------------------------------------------
      !CRM commented out this section in 2014 - outdated, I think
      !DO 91000 Iage = 0,Amax(Ispec)
      ! FecEst(Iage) = Fec(Ispec,Iage)/1000.0
      ! DO 91000 Isex = 1,Nsex
      !  WghtEst(Isex,Iage) = Wght(Ispec,Isex,Iage,Itime)/1000.0
!91000 CONTINUE         
C---------------------------------------------------------
C     How many fleets
      Estfleet = Iarea2-Iarea1+1
C      
      OPEN(UNIT=99,FILE='CAB.DAT')
C     Pre-set inputs (N fleets, 2 sexes, dummy lengths (44,2), age-specific
C     selectivity (1), simple weight-at-age calculate (0) 
      WRITE(99,601) 1,Year2-Year1+1,1,EstFleet,EstFleet,1,1,
     + 2,AmaxEst(Ispec),
     + 44,6,2,1,0 !return to this line and fill in new entries that make logical sense
C     Dummy growth parameters - not used 
      WRITE(99,602) 13.2d0,20.665d0,48.2d0,52.2d0,
     +              kpar(ISpec,1),kpar(Ispec,2),
     +              0.02d0,0.02d0,0.09d0,0.09d0,
     +              Apar(Ispec,1),Bpar(Ispec,1)
      WRITE(99,603)
      WRITE(99,6031) (WghtBeg(Ispec,1,Iage,0),Iage=0,AmaxEst(Ispec))
      WRITE(99,6031) (WghtBeg(Ispec,2,Iage,0),Iage=0,AmaxEst(Ispec))
      WRITE(99,604)
      WRITE(99,6041) (Wght(Ispec,1,Iage,1),Iage=0,AmaxEst(Ispec))
      WRITE(99,6041) (Wght(Ispec,2,Iage,1),Iage=0,AmaxEst(Ispec))
      WRITE(99,605)
      WRITE(99,6051) (Mature(Ispec,Iage),Iage=0,AmaxEst(Ispec))
C It's no longer fecundity that is written to cab, it's maturity.
C     WRITE(99,6051) (FecEst(Iage),Iage=0,AmaxEst(Ispec))
C
C     Catch data 
      WRITE(99,606)
      DO 10100 Yr = Year1,Year2
       WRITE(99,6061) Yr-Year1+1,(ObsCatch(Ispec,Iarea,Yr),Iarea=Iarea1,
     +    Iarea2)
10100 CONTINUE       
C
C     Survey and CPUE data
      NRel1 = 0
      NRel2 = 0
      DO 10210 Yr = Year1,Year2
       DO 10210 Iarea = Iarea1,Iarea2
        Ioff = Iarea-Iarea1+1
        IF (Cpue(Ispec,Iarea,Yr).GT.0) Nrel1(Ioff) = Nrel1(Ioff) + 1
        IF (SurveyEst(Ispec,Iarea,Yr).GT.0) Nrel2(Ioff) = Nrel2(Ioff)+1       
10210 CONTINUE
C     2,1,-3 (two series, one CPUE, one survey
      IF (UseFisheryCPUE.EQ.1) THEN
        WRITE(99,6071) 2*EstFleet !EstFleet refers to whether there are fleets by area or not, but is multiplied by 2 here to account for survey and fishery. This is something that needs to be updated to reflect current options.
        WRITE(99,7071) 1 2 !!1 is fishery, 2 is survey (links to fleetpops) Carey is adding this line to match new cab - not compatible with "areas as fleets" option.
        WRITE(99,6072) 1 1 !1+ = CPUE index, -1 = SSB, -2 = Recruitment, -3- means skip, do for each index (fishery and survey)
        WRITE(99,6073) (Nrel1(Iarea),Iarea=1,EstFleet),
     +  (Nrel2(Iarea),Iarea=1,EstFleet)
      ELSE
        WRITE(99,6071) EstFleet !CRM: Here there is 1 survey fleet and no abundance index for the fishery.
        WRITE(99,7071) 2 !2 is the index for the survey; Carey is adding this line to match new cab
        WRITE(99,6072) 1 !!1+ = CPUE index, -1 = SSB, -2 = Recruitment, -3- means skip (for survey index only)
        WRITE(99,6073) (Nrel2(Iarea),Iarea=1,EstFleet)
      ENDIF

C Carey is re-writing this line in the "IF" statement above and eliminating the "fleets" option.
C      WRITE(99,6072) (Iarea,Iarea=1,EstFleet),
C     +  (-2-Iarea,Iarea=1,EstFleet)
C      WRITE(99,6073) (Nrel1(Iarea),Iarea=1,EstFleet),
C     +  (Nrel2(Iarea),Iarea=1,EstFleet)
      WRITE(99,6074)
      IF (UseFisheryCPUE.EQ.1) THEN
      DO 10230 Iarea = Iarea1,Iarea2
       DO 10230 Yr = Year1,Year2
        IF (Cpue(Ispec,Iarea,Yr).GT.0) WRITE(99,6081) Iarea-Iarea1+1,
     +    Yr-Year1+1,Cpue(Ispec,Iarea,Yr),EstCpueCV
10230 CONTINUE   
      ENDIF
      
      WRITE(99,*)
      DO 10240 Iarea = Iarea1,Iarea2
       DO 10240 Yr = Year1,Year2
        IF (SurveyEst(Ispec,Iarea,Yr).GT.0) WRITE(99,6081) 
     +   Iarea-Iarea1+1,Yr-Year1+1,
     +   SurveyEst(Ispec,Iarea,Yr),EstSurveyCV
10240 CONTINUE       
C
C     Dummy length data (not used): may want to change this for tech interact project at some point.
      WRITE(99,609) 1,1,1,1
      WRITE(99,6091) 1,1,0,0.0d0,0.0d0,(1.0d0,ILen=1,44)
C
      Nage1 = 0
      Nage2 = 0
      DO 10310 ISex=1,2
       DO 10310 Iarea = Iarea1,Iarea2
        DO 10310 Yr = Year1,Year2
         Ioff = Iarea-Iarea1+1
         IF (ComSS(Ispec,Iarea,Yr).GT.0) Nage1(Ioff) = Nage1(Ioff) + 1
         IF (SuvSS(Ispec,Iarea,Yr).GT.0) Nage2(Ioff) = Nage2(Ioff) + 1
10310 CONTINUE         
      WRITE(99,6101) 2*EstFleet
      WRITE(99,7101) 1,2
      WRITE(99,6102) (Iarea,Iarea=1,EstFleet),
     +  (-2-Iarea,Iarea=1,EstFleet)
      WRITE(99,6103) (Nage1(Iarea),Iarea=1,EstFleet),
     +  (Nage2(Iarea),Iarea=1,EstFleet)
      WRITE(99,6104)
      DO 10320 Isex = 1,2
       DO 10320 Iarea = Iarea1,Iarea2
        DO 10320 Yr= Year1,Year2
         IF (ComSS(Ispec,Iarea,Yr).GT.0) THEN
          WRITE(99,6105) 1,Yr-Year1+1,Isex,
     +        ComSS(Ispec,Iarea,Yr)*1.0d0,ComSS(Ispec,Iarea,Yr)*1.0d0,
     +       (CNObs(Ispec,Iarea,Yr,ISex,IAge),IAge=0,AmaxEst(Ispec))
         ENDIF !currently there is no option to leave out the fishery age data
10320 CONTINUE         
      DO 10330 Isex = 1,2
       DO 10330 Iarea = Iarea1,Iarea2
        DO 10330 Yr= Year1,Year2
         IF (SuvSS(Ispec,Iarea,Yr).GT.0) THEN
          WRITE(99,6105) 2,Yr-Year1+1,Isex,
     +        SuvSS(Ispec,Iarea,Yr)*1.0d0,SuvSS(Ispec,Iarea,Yr)*1.0d0,
     +       (SNObs(Ispec,Iarea,Yr,ISex,IAge),IAge=0,AmaxEst(Ispec))
        ENDIF !currently there is no option to leave out the survey age data
10330 CONTINUE
C
      WRITE(99,611)
      CLOSE(99)
C******************************************************************************************
C     Iterate over configurations (configs only influence ctl file): loop should run through cab.exe, then call a batch file to rename output within the loop.
      DO 1298 Iconfig = 0,Nconfigs
C******************************************************************************************
C       The control file
        OPEN(UNIT=99,FILE="CAB.CTL")
        WRITE(99,701)  DebugAssessment
        WRITE(99,901)  1.0 !Not doing an an area-specific or population-specific assessment.
        WRITE(99,7021) EstM(Ispec,Iconfig)-BoundM,
     +EstM(Ispec,Iconfig)+BoundM,PhaseM !Careful not to use Bounds that are too big when using P* configurations and not actually estimating M or H
        WRITE(99,7022) Iestln_S0(Ispec)-BoundS0,Iestln_S0(Ispec)+
     +               BoundS0,PhaseS0
        WRITE(99,7023) IestSteep(Ispec,Iconfig)-BoundH,
     +IestSteep(Ispec,Iconfig)+BoundH,
     +               PhaseH
      !------------------------------------------------------------------------
      ! Carey does not want to estimate sigmaR so is changing the following lines:
      !      WRITE(99,7024) IestSigmaR(Ispec)**2.0-0.05,
      !     +      IestSigmaR(Ispec)**2.0+0.05
        WRITE(99,7024) IestSigmaR(Ispec)**2.0-0.000001,
     +      IestSigmaR(Ispec)**2.0+0.000001
! End of Carey's change.
!------------------------------------------------------------------------
      WRITE(99,703) StartYrRecDevs,(Year2-Year1+1)-IestIgR(Ispec),
     +  -10.0d0,10.0d0,2
C     The line below makes the program run for just one year to id problems (note by Carey)
C     WRITE(99,703) 2,2,-10.0d0,10.0d0,2
      
        WRITE(99,704) (1,Iarea=1,EstFleet),(1,Iarea=1,EstFleet)
 
      WRITE(99,904) 4 !Number of selex parameters (sum number of parameters separately for each fleet - this currently has to be counted manually.
      L50a = SCL50(Ispec)
      LDiff = SCL95(Ispec) - SCL50(Ispec)

C CRM commented out the following write statement to update to new cab:
C      WRITE(99,705) 1,L50a-BoundSel,L50a+BoundSel,PhaseSel,
C     +                LDiff-BoundSel,LDiff+BoundSel,PhaseSel

C CRM: the following write line is for new cab for fleet 1 = fishery data
      WRITE(99,905) SCL50(Ispec)-BoundSela50,SCL50(Ispec)+BoundSela50,
     +              SCL50(Ispec),PhaseSel,
     +              (SCL95(Ispec) - SCL50(Ispec))-BoundSelSlope,
     +              (SCL95(Ispec) - SCL50(Ispec))+BoundSelSlope,
     +              (SCL95(Ispec) - SCL50(Ispec)),PhaseSel 
 
      WRITE(99,9051) SSL50(Ispec)-BoundSela50,SSL50(Ispec)+BoundSela50,
     +              SSL50(Ispec),PhaseSel,
     +              (SSL95(Ispec) - SSL50(Ispec))-BoundSelSlope,
     +              (SSL95(Ispec) - SSL50(Ispec))+BoundSelSlope,
     +              (SSL95(Ispec) - SSL50(Ispec)),PhaseSel       

C This line is no longer in CAB (CRM):
C      WRITE(99,707) -1
C
C     The following line needs to be fixed - why double the number of fleets twice???
C      WRITE(99,708) (1,Iarea=1,2*EstFleet),1,(1,Iarea=1,2*EstFleet)
      IF (UseFisheryCPUE.EQ.1) THEN
       WRITE(99,708) (1,Iarea=1,2*EstFleet),1,(1,Iarea=1,2*EstFleet)
      ELSE
        WRITE(99,708) (1,Iarea=1,EstFleet),1,(1,Iarea=1,2*EstFleet)         
      ENDIF
      
      WRITE(99,7091) 
      WRITE(99,7092) 0
      WRITE(99,7093) (NINT(EstEffCatchSamp),Iarea=1,EstFleet),
     +      (NINT(EstEffSurveySamp),Iarea=1,EstFleet)
      WRITE(99,710) IBiasTYpe
      WRITE(99,711)
      CLOSE(99)
C
C     Copy files (for testing) 
      WRITE(TsCode1,8051) INDC,Iarea1
      WRITE(TsCode2,8052) INDC,Iarea1
      WRITE(tstOut,8053) TsCode1
      CALL SYSTEM(TstOut)
      WRITE(tstOut,8054) TsCode2
      CALL SYSTEM(TstOut)
C
C     Summation limits for aggregating results
      IF (Iarea1.EQ.0) THEN
       Iarea1a = 1
       Iarea2a = Narea
      ELSE
       Iarea1a = Iarea1
       Iarea2a = Iarea2
      ENDIF
     	Write(*,'(1x,A)') 'About to call CAB.EXE'
      Write(*,*) 'Sim: ', Isim, 'Year: ', Year2, 'Config: ', Iconfig
      !CALL SYSTEM("CAB.EXE -est -nox >output.log")
      CALL SYSTEM("CAB.EXE -nox >output.log")
      Write(*,*) 'Finished cab assessment'
      !Option 2: Potentially rename important output files to be config-specific here.
      !Option 2: Potentially end Iconfig loop here.


      CALL ReadCabRep(Year1,Year2,INDC,Iarea1,Iarea2,Carea1,Ispec,
     +Iarea1a,Iarea2a,Iconfig)  ! The commented out section below was moved to ReadCabRep
     
      CALL ReadPar(Ispec,Iarea1+Addone,Iconfig)
      
      !Prepare to copy cab.* to sim, year, and config specific files
      Write(cIsim,'(I2)') Isim
      Write(cIYr,'(I2)') GlobalProjYr
      Write(cIconfig,'(I2)') Iconfig
      cstring = 'SaveCabStuff '//cIYr//' '//cIsim//' '//
     +cIconfig
      
      !put a stop to further processing if the .cor file doesn't actually exist.
      INQUIRE(FILE = "cab.cor",EXIST=cor_exists)
      IF (cor_exists) THEN
          Write(*,*) 'cab.cor exists for Sim ',Isim,', Year ', Year2, 
     +' Config ', Iconfig
          CALL IterateNvecsMVN(Ispec,Iarea1+Addone,Iconfig,Year2+1) !Might need to move this to avoid iterating for iconfig = 0 (which will be the true value and not the assessment author's specified trial values)
      ENDIF
      !Make a flag
      IF (cor_exists.EQ..FALSE..OR.TestMVN.EQ.'Correlation') THEN
        WRITE(402,*) Isim,GlobalProjYr,Ispec,Iconfig,TestMVN
        !return without re-calculating the quota (auto-sets quota to last year's)
        QuotaAct(Ispec) = CatBio(Ispec,1,GlobalProjYr-1)
        Write(*,*) 'Using last years catch as the quota:'
        Write(*,*) '.cor file DNE or prob with MVN process for Sim ',
     +Isim,', Year ',Year2,', Config ', Iconfig
        !cleanup from this config 
        !Call a batch file to copy cab.rep to cab.rep_YrXX_configYY and delete extra cab files
        CALL SYSTEM(cstring)
        CALL SYSTEM("DeleteCabStuff")
        
        !1/6/15: need to delete mvninput files and EstNatAge,EstStdNatAge,EstNatCor here??
        
        !Return and run through annual proj, etc with last year's quota. Don't bother with control rule calculations
        RETURN !no use in continuing with the making of Nvecs as plan is to set quota to last year's, make a flag, and continue on.
      ENDIF
      
        !cleanup from this config 
        !Call a batch file to copy cab.rep to cab.rep_YrXX_configYY and delete extra cab files
        CALL SYSTEM(cstring)
        CALL SYSTEM("DeleteCabStuff")
      
 !       CLOSE(957)  !not used outside of ReadCabRep,so closed at the end of that function
 !       CLOSE(958)  !not used outside of ReadCabRep,so closed at the end of that function
 !       CLOSE(959)  !needed to be closed within ReadCabRep or else weird mistakes when something doesn't converge.
        
1298  CONTINUE !End Iconfig loop.


C
C----------------------------------------------------
!From the original project: Not sure these actually work anymore - there's a mysterious input file required. Better to re-program if needed (these are west coast rules)
C     Try a "complicated" control rule 
 !     CALL CalcQuota(Ispec,Iarea1+Addone,Year2)
 !     CALL TruQuota(Ispec,Iarea1,Year2,Year1)
C---------------------------------------------------

C---------------------------------------------------
!This is the p* control rule (Created by CRM)
      !This function sets the quota:
        !Will only get to Pstar calcs if all configs converged:
      IF (UsePstar.EQ.1) THEN 
         CALL FindPstarQuota(Year2,Ispec) !The above two subroutines are called within FindPstarQuota.
 !      FindFforCatchReg(Year2), FindFofl() both done within FindPstarQuota
      ENDIF
     
C-----------------------------------------------------

C---------------------------------------------------
!This is the AFSC tier system quota
      !This function sets the quota:
        !Will only get to Pstar calcs if all configs converged:
    !  IF (UsePstar.EQ.0) THEN 
    ! +     CALL FindPstarQuota(Year2,Ispec) !The above two subroutines are called within FindPstarQuota.
 !      FindFforCatchReg(Year2), FindFofl() both done within FindPstarQuota
C-----------------------------------------------------
        RETURN

      
8051  FORMAT("CABD.",A1,"0",I1)
8052  FORMAT("CABC.",A1,"0",I1)
8053  FORMAT("COPY CAB.DAT ",A8)
8054  FORMAT("COPY CAB.CTL ",A8)
701   FORMAT("# Debugging: turn off parameters, 1- Yes; 0 - No"/I4/
     +       "# Controls: Estimate the additional variation"/" 0")
901   FORMAT("# Proportion of recruitment among populations"/F7.4)
7021  FORMAT("# Prior for M"/
     +        F5.3,1x,F6.3,1x,I3)
7022  FORMAT("# Prior for ln(s0)"/
     +        F9.6,1x,F9.6,1x,I3)
7023  FORMAT("# Steepness (h)"/
     +        F5.3,1x,F6.3,1x,I3)
7024  FORMAT("# SigmaR"/
     +        F10.8,1x,F10.6,1x,"-2")             
703   FORMAT("# Start_rec_yr, End_rec_yr, Low_lmt, Upp_lmt, Phase"/
     +       I4,1x,I4,1x,F5.1,1x,F5.1,1x,I2)
704   FORMAT("# Sel patrns (3=Dome;2=DblLog;1=Log"/
     +       100(I2,1x))
904   FORMAT("#NselexCheck: Number of selex params"/I4)
705   FORMAT("# Logistic Selectivity pattern fleet ",I1/
     +       2(F8.4,1x,F8.4,1x,I2/))
905   FORMAT("# Logistic Selectivity pattern fishery (fleet 1, 
     + a50 first line, a95-a50 second line) "/,
     +      2(3(F8.4,1x),I2,/))
9051  FORMAT("# Logistic Selectivity pattern survey (fleet 2,
     + a50 first line, a95-a50 second line) "/,
     +      2(3(F8.4,1x),I2/))
706   FORMAT("# Dbl Logistic Selectivity"/
     +       7("-1 -1 -1"/))
707   FORMAT("# Dbl normal Selectivity"/
     +       (I2,1x,I2,1x,I2))
708   FORMAT("# Switches: use all data"/
     +       100(I2,1x))
7091  FORMAT("# Lambda for recruitment"/"1.0")
7092  FORMAT("# Effective sample size (length data)"/I5)
7093  FORMAT("# Effective sample size (age data)"/100(I5,1x))
710   FORMAT("# Bias correction"/I3)
711   FORMAT("# End-of-File"/" 999")
601   FORMAT('# Start Year'/I4/
     +       '# End Year'/I4/
     +       '# Number of pops'/I2/
     +       '# Number of fleets'/I2/
     +       '# Number of surveys'/I2/
     +       '# Which pop caught by each fleet (fisheries then surveys, 
     + each fleet catches only one pop)'/I2,1x,I2/
     +       '# Number of sexes'/I2/
     +       '# Number of age-classes'/I2/
     +       '# Number of length-classes'/I2/
     +       '# Start bin for length classes'/I2/
     +       '# Size of length bins'/I2/
     +       '# Age-specific selectivity (1=Yes;0=No)?'/I2/
     +       '# Weight-at-age account for selectivity?'/I2)
602   FORMAT('# Biological information (males then females)'/
     +        '#  Lmin, Lmax, K CVLmin, CVLMax: females then males'/
     +        5(1x,F8.3,1x,F8.3/)
     +        '# Weight-length relationship: intercept then slope'/
     +        F10.8,1x,F10.8)
603   FORMAT('# Weight-at-age (Start-year) females then males')
6031  FORMAT(100(F20.9,1x))      
604   FORMAT('# Weight-at-age (Mid-year) females then males')
6041  FORMAT(100(F20.9,1x))      
!605   FORMAT('# Fecundity-at-age (Maturity-at-age x Wt-at-age-beg-yr')
605   FORMAT('# Maturity-at-age (only maturity)')
6051  FORMAT(100(F12.9,1x))      
606   FORMAT('# Catch data')
6061  FORMAT(I4,5(1x,F15.7))
6071  FORMAT(/"# Abundance index data"/
     +       "# Number of indices"/I2)
7071  FORMAT(/"# Indx_fleet: which pop matches with which fleet??"/
     +       5(I2,1x))
6072  FORMAT("# Index type determines the nature of the index
     +: 1+ CPUE; -1 SSB; -2 Recruitment"/100(I2,1x))
6073  FORMAT("# Number of years for each index"/100(I3,1x))
6074  FORMAT("# The abundance data (series, year)")
6081  FORMAT(I2,1x,I4,1x,F12.4,1x,F5.3)     
609   FORMAT(/"# Length-frequency data"/
     +       "# Number of fleets"/I2/
     +       "# lenf_fleets: 1=fishery, 2 = survey "/I2/
     +       "# Length-types: outdated - just fill the space"/I2/
     +       "# Years of data"/I3/
     +       "# The actual data")
6091  FORMAT(I2,I4,1x,I2,1x,2(F8.4,1x),100(F7.5,1x))
6101  FORMAT(/"# Age-frequency data"/
     +       "# Number of fleets"/I2)
7101  FORMAT(/"# agef_fleet: 1= fishery, 2 = survey "/100(I2,1x))
6102  FORMAT("# Age-types: outdated, just fill the space"/100(I2,1x))
6103  FORMAT("# Years of data"/100(I3,1x))
6104  FORMAT("# The actual data")
6105  FORMAT(I2,I4,1x,I2,1x,2(F12.4,1x),100(F12.8,1x))
611   FORMAT(//"# End of file"/"999")

      END SUBROUTINE DoCab2

C
C ===========================================================================
C    
	  
      SUBROUTINE TruQuota(Ispec,Iarea,Year2,Year1)
C
C This subroutine computes the quota using the 40:10 rule
C      
	USE ShareDataRatioSubs
      IMPLICIT NONE
C
C     Global variables
      INTEGER Ispec,Iarea,Year2,Year1
C
C     Local variables
      REAL*8 NN(0:MaxAge),VSPR,Fmin,Fmax,FF,FSPR,Surv,Tmp,Depl
      REAL*8 OrigQuota,ExpBio2,Quota2(MaxSpec)
      REAL*8 Fpass,Finc,Clast,Cpass,Cat1,Cat2,Deriv,Quota3
      INTEGER Iage,MaxA,II,Isex,Itime,IBlKEW,IBlkNS
      REAL*8 GetYield
      EXTERNAL GetYield
C
C     Extract maximum age
      MaxA = Amax(Ispec) 
C
C     Find Virgin SPR
      NN(0) = 1
      DO 10100 Iage = 1,MaxA
       NN(Iage) = NN(Iage-1)*EXP(-M(Ispec))
10100 CONTINUE
      NN(MaxA) = NN(MaxA)/(1-EXP(-M(Ispec)))
      VSPR = 0
      DO 10200 Iage = 1,MaxA
       VSPR = VSPR + Fec(Ispec,Iage)*NN(Iage)
10200 CONTINUE       
C
C     Find the 40:10 exploitation rate
C	Carey: First find FSPR = Fished SPR and VSPR = Virgin SPR; then find the fishing mortality (FF) that satisfied FSPR = 0.4*VSPR using the bisection method
      Fmin = 0
      Fmax = 1 
      DO 20000 II = 1,50
C
C      Bisect
       FF = (Fmin+Fmax)/2.0

        NN(0) = 1
        DO 20100 Iage = 1,MaxA
         Surv =EXP(-M(Ispec))*(1-FF*Sel(Ispec,Iage-1))
         NN(Iage) = NN(Iage-1)*Surv
20100  CONTINUE
       Surv =EXP(-M(Ispec))*(1-FF*Sel(Ispec,MaxA))
       NN(MaxA) = NN(MaxA)/(1-Surv)
       FSPR = 0
       DO 20200 Iage = 1,MaxA
        FSPR = FSPR + Fec(Ispec,Iage)*NN(Iage)
20200  CONTINUE       
       IF (ABS(FSPR/VSPR-0.45).LT.0.0001) GOTO 20999
C       
       IF (FSPR.GT.0.45*VSPR) THEN
        Fmin = FF
       ELSE
        Fmax = FF 
       ENDIF 
C
20000 CONTINUE      
20999 CONTINUE
C
C     Calculate the ABC
      Quota2(Ispec) = 0
      Itime = (Ntime-1)/2+1
      DO 30000 Isex = 1,Nsex
       DO 30000 Iage = 0,MaxA
        DO 30000 IBlkNS = 1,NNSBlock
         DO 30000 IBlkEW = 1,NEWBlock
          Tmp = N(Ispec,Isex,IBlkNS,IBlkEW,Year2+1,Iage)*
     +        EXP(-M(Ispec)/2.0)
          Tmp = Tmp * Sel(Ispec,Iage)*FF
          Tmp = Tmp * Wght(Ispec,Isex,Iage,Itime)
          Quota2(Ispec) = Quota2(Ispec) + Tmp
30000 CONTINUE       
      ExpBio2 = Quota2(Ispec)/FF
C
C     Do the 40-10 adjustment
      Depl =SpBioTot(Ispec,Year2+1)/SpBioTot(Ispec,Year1)
      OrigQuota = Quota2(Ispec)
      IF (Depl.LT.0.4) THEN
       Quota2(Ispec) = (Depl-0.1)/0.3*Quota2(Ispec)
       IF (Quota2(Ispec).LT.0) Quota2(ISpec) = 0
      ENDIF
      WRITE(29,601) Ispec,Iarea,Year2+1,FSPR/VSPR,FF,
     +     OrigQuota,Depl,Quota2(Ispec),ExpBio2,SpBioTot(Ispec,Year2+1)
C
      Fpass = 0
      Finc = 0.001
      Clast = -1
      DO 40000 II = 1,1000
       Fpass = Fpass + Finc
       Cpass = GetYield(Fpass,Ispec)
       IF (Cpass.LT.Clast) THEN
        GOTO 40001
       ELSE
        Clast = Cpass
       ENDIF
40000 CONTINUE
40001 CONTINUE
C
C     Bisection
      Fmin = Fpass - Finc
      Fmax = Fpass
      DO 50000 II = 1,50
       FF = (Fmin+Fmax)/2.0
       Cat1 = GetYield(FF-0.0005,Ispec)
       Cat2 = GetYield(FF+0.0005,Ispec)
       Deriv = (Cat1+Cat2)/2
       IF (Deriv.LT.0) THEN
        FF = Fmax
       ELSE
        FF = Fmin
       ENDIF
50000 CONTINUE
      Cpass = GetYield(FF,Ispec)   
C
C     Calculate the Quota
      Quota3 = 0
      DO 60000 Isex = 1,Nsex
       DO 60000 Iage = 0,MaxA
        DO 60000 IBlkNS = 1,NNSBlock
         DO 60000 IBlkEW = 1,NEWBlock
          Tmp = N(Ispec,Isex,IBlkNS,IBlkEW,Year2+1,Iage)*
     +        EXP(-M(Ispec)/2.0)
          Tmp = Tmp * Sel(Ispec,Iage)*FF
          Tmp = Tmp * Wght(Ispec,Isex,Iage,Itime)
          Quota3 = Quota3 + Tmp
60000 CONTINUE       
      WRITE(29,602) Ispec,Iarea,Year2+1,-1.0d0,FF,
     +     Quota3,Depl,Quota3,ExpBio2,SpBioTot(Ispec,Year2+1)
C
      RETURN
601   FORMAT("T45",1x,I2,1x,I2,1x,I4,1x,F6.3,1x,F7.5,1x,F8.2,1x,F5.3,
     +       1x,F8.2,2(1x,F10.1))
602   FORMAT("TDM",1x,I2,1x,I2,1x,I4,1x,F6.3,1x,F7.5,1x,F8.2,1x,F5.3,
     +       1x,F8.2,2(1x,F10.1))
      END SUBROUTINE TruQuota

C
C ===========================================================================
C    
	  
      REAL*8 FUNCTION getYield(Fpass,Ispec)
C
C This function calculates Yields given F. ONLY used with the TruQuota (40:10 rule) (But there seems to be an error)
C Steps:
C 	1. Calculates the unfished SPR
C 	2. Calculates the fished SPR for a specific F level	
C 	3. Calculates the yield per recruit
C 	4. The calculation of recruitment if erroneous. It is not based on Spawning biomass
C      

	USE ShareDataRatioSubs
      IMPLICIT NONE
C
C     Global variables
      REAL*8 Fpass
      INTEGER Ispec
C
C     Local variables  
      REAL*8 NN(0:MaxAge),VSPR,FSPR,Surv,YPR,Top,Bot,Recr
      INTEGER Iage,MaxA,Itime
C
C     Extract maximum age
      MaxA = Amax(Ispec) 
      Itime = (Ntime-1)/2+1
C
C     Find Virgin SPR
      NN(0) = 1
      DO 10100 Iage = 1,MaxA
       NN(Iage) = NN(Iage-1)*EXP(-M(Ispec))
10100 CONTINUE
      NN(MaxA) = NN(MaxA)/(1-EXP(-M(Ispec)))
      VSPR = 0
      DO 10200 Iage = 1,MaxA
       VSPR = VSPR + Fec(Ispec,Iage)*NN(Iage)
10200 CONTINUE       
C
C     Age-stucture
      NN(0) = 1
      DO 20100 Iage = 1,MaxA
       Surv =EXP(-M(Ispec))*(1-Fpass*Sel(Ispec,Iage-1))
       NN(Iage) = NN(Iage-1)*Surv
20100 CONTINUE
      Surv = EXP(-M(Ispec))*(1-Fpass*Sel(Ispec,MaxA))
      NN(MaxA) = NN(MaxA)/(1-Surv)
      FSPR = 0
      DO 20200 Iage = 1,MaxA
       FSPR = FSPR + Fec(Ispec,Iage)*NN(Iage)
20200 CONTINUE       
      YPR = 0
      DO 20300 Iage = 1,MaxA
       YPR = YPR + NN(Iage)*Fpass*EXP(-M(Ispec)/2)*Sel(Ispec,Iage)*
     +     Wght(Ispec,1,Iage,Itime)
20300 CONTINUE
C
C     Recruitment 
      Top = 4*Steep(Ispec)*FSPR - (1-Steep(ISpec))*VSPR
      Bot = 5*Steep(ISpec)-1
      Recr = Top/Bot/FSPR
      GetYield = YPR *Recr
C
      RETURN
      END FUNCTION getYield
      
C
C===============================================================================
C
      
	SUBROUTINE FindMSY()
C
C This subroutine is supposed to find the MSY value
C Steps:
C 	1. It calls FIT2 
C 	2. then FUNK2
C	
        USE ShareDataRatioSubs
        IMPLICIT NONE
        
	  REAL(8), DIMENSION(100):: Y
	  REAL(8) :: FUNK2,Z
        EXTERNAL FUNK2
        Funk2_MinOpt = 1 !CRM 7/30/14 I don't know if this optimization still works, but now FUNK2 is used for the P* approach, so different minimization options are necessary.
	  CareyGRD = 0.2
        Y(1) = LOG(1.5)
        CALL FIT2(Y,Z,1,CareyGRD)
        
	  Z = FUNK2(Y)
	  MSYEffort1 = Y(1)

      END SUBROUTINE FindMSY

C
C===============================================================================
C

	REAL*8 FUNCTION FUNK2(Y)
C This function is a totally different version that the FUNK function above but produces the TrialP value
C	Steps:
C 	1. For Funk2_minOpt=1, the function calculates the catch per projected year
C 		1.1 It projects the population forward
C 		1.2 It generates the data that goes along and the catch per year
C 	2. For Funk2_minOpt=2, the function optimizes the F level that leads to a desired P* value
C 		2.1 It finds F for a proposed ABC = proposed catch levels
C 		2.2 It then find the FOFL for each of the randomly drawn population number at age
C 		2.3 It then calculates the number of times the proposed ABC lead to F>FOFL = P*
C
      
      USE ShareDataRatioSubs
      IMPLICIT NONE

	  REAL(8), DIMENSION(100) :: Y
	  
        !LocalVariables
        INTEGER :: Yr,Ispec,Isex,Iage,Iconfig,Ivec
        REAL(8) :: WeightedOverfishNum,WeightedAllNum
        INTEGER,DIMENSION(1:Nconfigs) :: OverfishNum
        REAL(8) :: TrialP
        
    ! Case when Funk2_MinOpt = 1 --> this calculates the catch per projection year   
        IF (Funk2_MinOpt.EQ.1) THEN
C  	    DO 35000 Yr = -50,250
	    DO 35000 Yr = MinPYr,250
            CALL ProjYr(Yr,Y(1))
	      CALL GenData(Yr)

	      CatNumSum(Yr) = 0
C           Carey is calculating Catch in numbers
	      DO 22000 Ispec = 1,NSpec
	        DO 22000 Iage = 0,MaxAge
	          DO 22000 Isex = 1,2
	            CatNumSum(Yr) = CatNumSum(Yr) + 
     +	    	CatNum(Ispec,1,Yr,Isex,Iage)
22000  	  CONTINUE
35000	     CONTINUE

          !FUNK2 = -CatNumSum(250)
	    !FUNK2 = -CatNumSum(MinPYr+200)
          FUNK2 = 0
      ENDIF !end Funk2_MinOpt

    ! Case when Funk2_MinOpt = 2 --> this     
       IF (Funk2_MinOpt.EQ.2) THEN
        !This is for the outer layer of the nested optimization developed for the p* approach    
           !Below is cut and pasted into FUNK2 with global year and global ispec.
          !Find F(t) that would produce SeedQuota
          TrialQuota = EXP(Y(1))
          CALL FindFforCatchReg(GlobalProjYr-1,GlobalIspec)
   
          !Find FOFL for all of the random vectors
          CALL FindFofl(GlobalIspec,GlobalProjYr-1) !Enter Ispec as an input! CRM 7/28/14
       
          !Integrate over (F(t),FOFL) distribution
          !Want the proportion of time F(t) > FOFL among the pairs (F(t),FOFL) sampled by way of the N-at-age vectors.
          OverfishNum = 0
          WeightedOverfishNum = 0.0
          WeightedAllNum = 0.0
          DO 303 Iconfig = 1,Nconfigs
              !ADD a line: skip this loop if cor_exists is false
              !set Fdistn to -1.0 and reconfigure Nmvn to take into account the missing configuration.
              DO 304 Ivec = 1,Nmvn
                  IF (FDistn(GlobalIspec,GlobalProjYr-1,Iconfig,Ivec).
     +GT.EstFOFL(GlobalISpec,Iconfig,Ivec)) THEN
                    OverfishNum(Iconfig) = OverFishNum(Iconfig)+1
                    ENDIF
304           CONTINUE !end Ivec loop
              WeightedOverfishNum = WeightedOverfishNum + 
     +OverfishNum(Iconfig)*ConfigWts(GlobalIspec,Iconfig) 
              WeightedAllNum = WeightedAllNum + 
     +Nmvn*ConfigWts(GlobalIspec,Iconfig)
303       CONTINUE !end Iconfig loop
      
          TrialP = WeightedOverfishNum/WeightedAllNum 

 !         FUNK2 = ((TrialP-Pstar)**2)*10000
           FUNK2 = TrialP
          
      ENDIF
      
           END FUNCTION FUNK2
C
C==============================================================================
C
           
      SUBROUTINE FIT2(X,SSS,NDIM,GRD)
C This subroutine calls amoeba2 to find the proposed ABC level that leads to a proposed P* value
C Set up the parameters for a fit

C     IMPLICIT NONE
      REAL*8 P(101,100),Y(101),X(100),SSS,FUNK,TOL
	REAL(8), INTENT(IN) :: GRD
      INTEGER*4 NDIM,I,J,ITER
      EXTERNAL FUNK2

C     SET UP TOLERANCES AND GRIDDING
      TOL = 0.00000001
C      GRD = 1.2

      DO 1 I=1,NDIM+1
         DO 2 J=1,NDIM
            P(I,J)=X(J)
            IF((I-1).EQ.J) P(I,J) = GRD*P(I,J)
2        CONTINUE
1     CONTINUE

      DO 3 I=1,NDIM+1
         DO 4 J=1,NDIM
            X(J)=P(I,J)
4        CONTINUE
         Y(I)=FUNK2(X)
3     CONTINUE

      CALL AMOEBA2(P,Y,101,100,NDIM,TOL,ITER)

      DO 5 J = 1,NDIM
         X(J) = P(1,J)
5     CONTINUE
      SSS = Y(1)

      RETURN
      END

C
C=============================================================
C
      
      SUBROUTINE AMOEBA2(P,Y,MP,NP,NDIM,FTOL,ITER)
C
C     MULTIDIMENSIONAL MINIMISATION OF THE FUNCTION FUNK(X) WHERE X IS
C     AN NDIM-DIMENSIONAL VECTOR, BY THE DOWNHILL SIMPLEX METHOD OF
C     NELDER AND MEAD. INPUT IS A MATRIX P WHOSE NDIM+1 ROWS ARE THE
C     NDIM-DIMENSIONAL VECTORS WHICH ARE THE VERTICES OF THE STARTING
C     SIMPLEX. [LOGICAL DIMENSIONS OF P ARE P(NDIM+1,NDIM); PHYSICAL
C     DIMENSIONS ARE INPUT AS P(MP,NP).] ALSO INPUT IS THE VECTOR Y
C     OF LENGTH NDIM+1, WHOSE COMPONENTS MUST BE PRE-INITIALISED TO
C     THE VALUES OF FUNK EVALUATED AT THE NDIM+1 VERTICES (ROWS) OF P;
C     AND FTOL IS THE FRACTIONAL CONVERGENCE TOLERANCE TO BE ACHIEVED
C     IN THE FUNCTION VALUE (N.B.!). ON OUTPUT, P AND Y WILL HAVE BEEN
C     RESET TO NDIM+1 NEW POINTS ALL WITHIN FTOL OF A MINIMUM FUNCTION
C     VALUE, AND ITER GIVES THE NUMBER OF ITERATIONS TAKEN.
C
C     FROM: NUMERICAL RECIPES - THE ART OF SCIENTIFIC COMPUTING
C           BY W. H. PRESS ET AL, CAMBRIDGE UNIVERSITY PRESS
C           ISBN 0-251-30811-9
C
C     ********************************************************************
C

C     SPECIFY THE MAXIMUM NUMBER OF DIMENSIONS, THREE PARAMETERS WHICH
C     DEFINE THE EXPANSIONS AND CONTRACTIONS, AND THE MAXIMUM NUMBER OF
C     ITERATIONS ALLOWED

C     IMPLICIT NONE
      REAL*8 ALPHA,BETA,GAMMA,EPS
      INTEGER*4 ITMAX
      PARAMETER (ALPHA=1.0,BETA=0.5,GAMMA=2.0,ITMAX=5000,EPS=1.0E-10)

C     Global Data
      INTEGER*4 MP,NP,NDIM,ITER
      REAL*8 P(MP,NP),Y(MP),PR(100),PRR(100),PBAR(100),FTOL,FUNK2
      EXTERNAL FUNK2

C     Local Data
      REAL*8 YPR,YPRR
      INTEGER*4 I,ILO,IHI,INHI,J,MPTS

C     NOTE THAT MP IS THE PHYSICAL DIMENSION CORRESPONDING TO THE LOGICAL
C     DIMENSION MPTS, NP TO NDIM.

      MPTS = NDIM+1
      ITER=0

C     FIRST WE MUST DETERMINE WHICH POINT IS THE HIGHEST (WORST), NEXT
C     HIGHEST, AND LOWEST (BEST).

1     ILO=1
      IF(Y(1).GT.Y(2)) THEN
         IHI=1
         INHI=2
      ELSE
         IHI=2
         INHI=1
      ENDIF
      DO 11 I=1,MPTS
         IF(Y(I).LT.Y(ILO)) ILO=I
         IF(Y(I).GT.Y(IHI)) THEN
            INHI=IHI
            IHI=I
         ELSE IF(Y(I).GT.Y(INHI)) THEN
            IF(I.NE.IHI) INHI=I
         ENDIF
11    CONTINUE

C     COMPUTE THE FRACTIONAL RANGE FROM THE HIGHEST TO THE LOWEST AND
C     RETURN IF SATISFACTORY

      IF(2.0*ABS(Y(IHI)-Y(ILO)).LT.FTOL*(ABS(Y(IHI))+ABS(Y(ILO))+EPS)) 
     +   RETURN
      IF(ITER.EQ.ITMAX) THEN
c        WRITE(6,200)
200      FORMAT(1H ,'AMOEBA EXCEEDING MAXIMUM ITERATIONS')
         RETURN
      ENDIF
      ITER=ITER+1

      DO 12 J=1,NDIM
         PBAR(J)=0.
12    CONTINUE

C     BEGIN A NEW ITERATION. COMPUTE THE VECTOR AVERAGE OF ALL POINTS
C     EXCEPT THE HIGHEST, I.E. THE CENTRE OF THE "FACE" OF THE SIMPLEX
C     ACROSS FROM THE HIGH POINT. WE WILL SUBSEQUENTLY EXPLORE ALONG
C     THE RAY FROM THE HIGH POINT THROUGH THE CENTRE.

      DO 14 I=1,MPTS
         IF(I.NE.IHI) THEN
            DO 13 J=1,NDIM
               PBAR(J)=PBAR(J)+P(I,J)
13          CONTINUE
         ENDIF
14    CONTINUE

C     EXTRAPOLATE BY A FACTOR ALPHA THROUGH THE FACE, I.E. REFLECT THE
C     SIMPLEX FROM THE HIGH POINT

      DO 15 J=1,NDIM
         PBAR(J)=PBAR(J)/NDIM
         PR(J)=(1.+ALPHA)*PBAR(J)-ALPHA*P(IHI,J)
15    CONTINUE

C     EVALUATE THE FUNCTION AT THE REFLECTED POINT

      YPR=FUNK2(PR)

C     GIVES A BETTER RESULT THAN THE BEST POINT, SO TRY AN ADDITIONAL
C     EXTRAPOLATION BY A FACTOR GAMMA

      IF(YPR.LE.Y(ILO)) THEN
         DO 16 J=1,NDIM
            PRR(J)=GAMMA*PR(J)+(1.-GAMMA)*PBAR(J)
16       CONTINUE

C        CHECK THE FUNCTION THERE

         YPRR=FUNK2(PRR)

C        THE ADDITIONAL EXTRAPOLATION SUCCEEDED, AND REPLACES THE
C        HIGHEST POINT

         IF(YPRR.LT.Y(ILO)) THEN
            DO 17 J=1,NDIM
               P(IHI,J)=PRR(J)
17          CONTINUE
            Y(IHI)=YPRR
         ELSE

C        THE ADDITIONAL EXTRAPOLATION FAILED, BUT WE CAN STILL USE THE
C        REFLECTED POINT

         DO 18 J=1,NDIM
            P(IHI,J)=PR(J)
18       CONTINUE
         Y(IHI)=YPR
      ENDIF

C     THE REFLECTED POINT IS WORSE THAN THE SECOND HIGHEST

      ELSE IF(YPR.GE.Y(INHI)) THEN

C        IF IT'S BETTER THAN THE HIGHEST, THEN REPLACE THE HIGHEST

         IF(YPR.LT.Y(IHI)) THEN
            DO 19 J=1,NDIM
               P(IHI,J)=PR(J)
19          CONTINUE
            Y(IHI)=YPR
         ENDIF

C        BUT LOOK FOR AN INTERMEDIATE LOWER POINT; IN OTHER WORDS
C        PERFORM A CONTRACTION OF THE SIMPLEX ALONG ONE DIMENSION
C        AND THEN EVALUATE THE FUNCTION

         DO 21 J=1,NDIM
            PRR(J)=BETA*P(IHI,J)+(1.-BETA)*PBAR(J)
21       CONTINUE
         YPRR=FUNK2(PRR)

C        CONTRACTION GIVES AN IMPROVEMENT, SO ACCEPT IT

         IF(YPRR.LT.Y(IHI)) THEN
            DO 22 J=1,NDIM
               P(IHI,J)=PRR(J)
22          CONTINUE
            Y(IHI)=YPRR
         ELSE

C           CAN'T SEEM TO GET RID OF THAT HIGH POINT. BETTER CONTRACT
C           AROUND THE LOWEST (BEST) POINT

            DO 24 I=1,MPTS
               IF(I.NE.ILO) THEN
                  DO 23 J=1,NDIM
                     PR(J)=0.5*(P(I,J)+P(ILO,J))
                     P(I,J)=PR(J)
23                CONTINUE
                  Y(I)=FUNK2(PR)
               ENDIF
24          CONTINUE
         ENDIF

      ELSE

C        WE ARRIVE HERE IF THE ORIGINAL REFLECTION GIVES A MIDDLING
C        POINT. REPLACE THE OLD HIGH POINT AND CONTINUE

         DO 25 J=1,NDIM
            P(IHI,J)=PR(J)
25       CONTINUE
         Y(IHI)=YPR

      ENDIF

      GO TO 1

      END

C      
C=========================================================================
C

      SUBROUTINE RPERM(P,N,RANDOMSEED)
!===Generate a random permutation, P, of the first N integers.
!   (Shuffling, equivalent to sampling WITHOUT REPLACEMENT).
!   Adaptation of Knuth Volume 2, Algorithm 3.4.2P.
      INTEGER N,P(N), K,J,I,IPJ,ITEMP,M
	INTEGER, INTENT(IN), DIMENSION(3) :: RANDOMSEED
      REAL U(100)
      DO I=1,N
        P(I)=I
      END DO
!---Generate up to 100 U(0,1) numbers at a time.
      DO I=1,N,100
        M=MIN(N-I+1,100)
        CALL RAND(U,M,RANDOMSEED)
        DO J=1,M
          IPJ=I+J-1
          K=INT(U(J)*(N-IPJ+1))+IPJ
          ITEMP=P(IPJ)
          P(IPJ)=P(K)
          P(K)=ITEMP
        END DO
      END DO
      RETURN
      END

C
C=======================================================================
C

      SUBROUTINE RAND(U,N,RANDOMSEED)
      INTEGER N,X,Y,Z
      REAL U(N),V
	INTEGER, INTENT(IN), DIMENSION(3) :: RANDOMSEED

	X = RANDOMSEED(1)
	Y = RANDOMSEED(2)
	Z = RANDOMSEED(3)
C      COMMON/RANDOM/X,Y,Z
      IF(N.LE.0) RETURN
      DO  I=1,N
        X=171*MOD(X,177)-2*(X/177)
        IF(X.LT.0) X=X+30269
        Y=172*MOD(Y,176)-35*(Y/176)
        IF(Y.LT.0) Y=Y+30307
        Z=170*MOD(Z,178)-63*(Z/178)
        IF(Z.LT.0) Z=Z+30323
        V=X/30269.0 + Y/30307.0 + Z/30323.0
        U(I)=V-INT(V)
      END DO
      RETURN
      END

C      
C====================================================================
C

      SUBROUTINE WriteSpecs(UNUM)
C
C This subroutine writes the specs for the run to a file for tracking what the simulation did
C The file (SUMM1.OUT) can be copied and pasted into an excel spreadsheet which tracks all of the simulations that have been completed.     
C This should be called just after READSPEC() and READCTL(), Outside of Supersim and NSim loops
	USE ShareDataRatioSubs
      IMPLICIT NONE

C	Local Variables
	INTEGER Iarea,Ispec
	INTEGER, INTENT(IN) :: UNUM

	
	WRITE(UNUM,'(I8)') Nsim			!Number of Simulations
	WRITE(UNUM,'(I8)') NNSBlock		!Number of North-South Cells
	WRITE(UNUM,'(I8)') NEWBlock		!Number of East-West Cells
      WRITE(UNUM,'(5(F6.3,1x))') (InitDep(Ispec),Ispec=1,Nspec)	!Initial depletion level (year 0)
      WRITE(UNUM,'(5(F6.3,1x))') (FinalDep(Ispec),Ispec=1,Nspec)	!Final depletion level (prior to implementing the control rule)
	WRITE(UNUM,'(I4)') BaseYr		!Base Year for depletion
	WRITE(UNUM,'(I4)') FinDepYr		!Final Year for depletion
      WRITE(UNUM,'(100(I8,1x))') (MinCNS(Iarea),Iarea=1,NClosedAreas)  !Limits of MPAs
      WRITE(UNUM,'(100(I8,1x))') (MaxCNS(Iarea),Iarea=1,NClosedAreas)
      WRITE(UNUM,'(100(I8,1x))') (MinCEW(Iarea),Iarea=1,NClosedAreas)
      WRITE(UNUM,'(100(I8,1x))') (MaxCEW(Iarea),Iarea=1,NClosedAreas)
	WRITE(UNUM,'(I8)') OptCRrecord			!Optimize the control rule?
	WRITE(UNUM,'(F8.4)') TotEffort(BaseYr) !Year 0 Effort (used as a seed)
	WRITE(UNUM,*)
	WRITE(UNUM,'(I8)') IsDeterm		!Is the sampling deterministic?
	WRITE(UNUM,'(5(F8.5,1x))') (SigmaR1(Ispec),Ispec=1,Nspec)	!Global recruitment variation
	WRITE(UNUM,'(5(F8.5,1x))') (SigmaR2(Ispec),Ispec=1,Nspec)	!Local recruitment variation
	WRITE(UNUM,'(F8.2)') CpueCV		!CV of CPUE data	
	WRITE(UNUM,'(I8)') FreqSurv		!Frequency of surveys 1 = every year; 3 = every third year
	WRITE(UNUM,'(F8.2)') SurvCV		!CV of survey data
	WRITE(UNUM,'(I8)') EffectSampSurvey2 !Effective survey sample size POST MPA
	WRITE(UNUM,'(I8)') ImplOpt		!Is there implementation error (1 or 0)	
	WRITE(UNUM,'(5(F8.5,1x))') (SigmaSpat(ISpec),Ispec=1,Nspec)  !SigmaQ = CV on catchability.
	WRITE(UNUM,*)
	WRITE(UNUM,'(5(F8.2,1x))') (M(Ispec),Ispec=1,Nspec)
	WRITE(UNUM,'(5(F8.3,1x))') (Steep(ISpec),Ispec=1,Nspec)	!Steepness
	WRITE(UNUM,'(5(I8,1x))') (SrType(ISpec),Ispec=1,Nspec)   !Type of Density 
	WRITE(UNUM,*)
      WRITE(UNUM,'(5(F8.1,1x))') (Linf(Ispec,1),Ispec=1,Nspec)
      WRITE(UNUM,'(5(F8.1,1x))') (Linf(Ispec,2),Ispec=1,Nspec)
      WRITE(UNUM,'(5(F8.3,1x))') (kpar(Ispec,1),Ispec=1,Nspec)
      WRITE(UNUM,'(5(F8.3,1x))') (kpar(Ispec,2),Ispec=1,Nspec)
      WRITE(UNUM,'(5(F8.2,1x))') (t0(Ispec,1),Ispec=1,Nspec)
      WRITE(UNUM,'(5(F8.2,1x))') (t0(Ispec,2),Ispec=1,Nspec)
      WRITE(UNUM,'(5(F8.3,1x))') (Apar(Ispec,1),Ispec=1,Nspec)
      WRITE(UNUM,'(5(F8.3,1x))') (Apar(Ispec,2),Ispec=1,Nspec)
      WRITE(UNUM,'(5(F8.3,1x))') (Bpar(Ispec,1),Ispec=1,Nspec)
      WRITE(UNUM,'(5(F8.3,1x))') (Bpar(Ispec,2),Ispec=1,Nspec)	
	WRITE(UNUM,*)
      WRITE(UNUM,'(5(F8.2,1x))') (L50(Ispec),Ispec=1,Nspec)
      WRITE(UNUM,'(5(F8.2,1x))') (L95(Ispec),Ispec=1,Nspec)
	WRITE(UNUM,*)
	WRITE(UNUM,'(5(I8,1x))') (RecruitOpt(ISpec),Ispec=1,Nspec)
      WRITE(UNUM,'(5(F8.3,1x))') (RecX(Ispec),Ispec=1,Nspec)
      WRITE(UNUM,'(5(F8.3,1x))') (RecY(Ispec),Ispec=1,Nspec)
      WRITE(UNUM,'(5(F8.3,1x))') (VRecXL(Ispec),Ispec=1,Nspec)
      WRITE(UNUM,'(5(F8.3,1x))') (VRecYL(Ispec),Ispec=1,Nspec)
      WRITE(UNUM,'(5(F8.3,1x))') (VRecXH(Ispec),Ispec=1,Nspec)
      WRITE(UNUM,'(5(F8.3,1x))') (VRecYH(Ispec),Ispec=1,Nspec)
      WRITE(UNUM,*)
      WRITE(UNUM,'(5(F8.3,1x))') (MixN(Ispec),Ispec=1,Nspec)
      WRITE(UNUM,'(5(F8.3,1x))') (MixE(Ispec),Ispec=1,Nspec)
      WRITE(UNUM,'(5(F8.3,1x))') (X1N(Ispec),Ispec=1,Nspec)
      WRITE(UNUM,'(5(F8.3,1x))') (X2N(Ispec),Ispec=1,Nspec)
      WRITE(UNUM,'(5(F8.3,1x))') (X1E(Ispec),Ispec=1,Nspec)
      WRITE(UNUM,'(5(F8.3,1x))') (X2E(Ispec),Ispec=1,Nspec)
	WRITE(UNUM,'(5(F8.3,1x))') (V1N(Ispec),Ispec=1,Nspec)
	WRITE(UNUM,'(5(F8.3,1x))') (V2N(Ispec),Ispec=1,Nspec)
	WRITE(UNUM,'(5(F8.3,1x))') (V1E(Ispec),Ispec=1,Nspec)
	WRITE(UNUM,'(5(F8.3,1x))') (V2E(Ispec),Ispec=1,Nspec)
C     Selectivity for the fishery 
	WRITE(UNUM,*)
      WRITE(UNUM,'(5(F8.3,1x))') (SCL50(ISpec),Ispec=1,Nspec)
      WRITE(UNUM,'(5(F8.3,1x))') (SCL95(ISpec),Ispec=1,Nspec)
      WRITE(UNUM,*)
C     Selectivity for the survey 
      WRITE(UNUM,'(5(F8.3,1x))') (SSL50(ISpec),Ispec=1,Nspec)
      WRITE(UNUM,'(5(F8.3,1x))') (SSL95(ISpec),Ispec=1,Nspec)
C	Optimal policy information
C	Extras
	WRITE(UNUM,'(I4)') MPAYear  !Year that MPA is implemented
	WRITE(UNUM,'(5(I4,1x))') (Amax(Ispec),Ispec=1,Nspec)  !Oldest age in the model
            WRITE(UNUM,*)
      WRITE(UNUM,'(5(F8.3,1x))') (MixNLarval(Ispec),Ispec=1,Nspec)
      WRITE(UNUM,'(5(F8.3,1x))') (MixELarval(Ispec),Ispec=1,Nspec)
      WRITE(UNUM,'(5(F8.3,1x))') (X1NLarval(Ispec),Ispec=1,Nspec)
      WRITE(UNUM,'(5(F8.3,1x))') (X2NLarval(Ispec),Ispec=1,Nspec)
      WRITE(UNUM,'(5(F8.3,1x))') (X1ELarval(Ispec),Ispec=1,Nspec)
      WRITE(UNUM,'(5(F8.3,1x))') (X2ELarval(Ispec),Ispec=1,Nspec)
	WRITE(UNUM,'(5(F8.3,1x))') (V1NLarval(Ispec),Ispec=1,Nspec)
	WRITE(UNUM,'(5(F8.3,1x))') (V2NLarval(Ispec),Ispec=1,Nspec)
	WRITE(UNUM,'(5(F8.3,1x))') (V1ELarval(Ispec),Ispec=1,Nspec)
	WRITE(UNUM,'(5(F8.3,1x))') (V2ELarval(Ispec),Ispec=1,Nspec)

	RETURN
	END
C===================================================================

C
C ===========================================================================
C
	
      SUBROUTINE FindInitEffort()
C
C This subroutine finds the q-values that leads to the specified initial depletion level
C
	USE ShareDataRatioSubs
      IMPLICIT NONE
C
C     Local variables
      REAL*8 X(100),SS,ii,jj,IEffort,BestRoughEffort
	REAL*8 i,ObjFun,ObjFunRough
      INTEGER Ispec,kk,ll
      REAL*8 FUNK
      REAL*8 FUNKSCALAR
      EXTERNAL FUNK
      EXTERNAL FUNKSCALAR
C
      MinOpt = 8
C	don't need CareyGRD if not using FIT.
C	CareyGRD = 1.2
	ObjFun = 1000000
	DO 10001 ii = 0,30
C	  DO 10000 Ispec = 1,Nspec
	   i = ii*3 + 0.00001 !to avoid divide by zeros
         X(1) = i    
C        X(Ispec) = TotEffortSeed
C10000   CONTINUE
        SS = FUNK(X)	
	    IF (SS.LT.ObjFun) THEN
	      ObjFun = SS
	      BestRoughEffort = X(1)
          ENDIF

10001 CONTINUE

Ctest
	ObjFunRough = ObjFun

C	WRITE(*,'(1x,A)') 'Multiple trials with BestRoughEffort'
C	SS = FUNK(BestRoughEffort)
C	WRITE(*,*)
C	WRITE(*,'(1x,A)') 'SpBioTot(MinPYr),SpBioTot(BaseYr),
C     +	N(1,1,1,1,-16,1),SSBEq(1),Neq(1,1,1,1,0)'
C	WRITE(*,*) SPBioTot(1,BaseYr),SPBioTot(1,MinPYr),
C     +	N(1,1,1,1,-16,1),SSBEq(1),Neq(1,2,3,2,5)
C	X(1) = BestRoughEffort
C	SS = FUNK(X)
C	WRITE(*,*) SPBioTot(1,BaseYr),SPBioTot(1,MinPYr),
C     +	N(1,1,1,1,-16,1),SSBEq(1),Neq(1,2,3,2,5)
C	X(1) = 1.01

C=====================================================
C	Take the rough estimate and do a finer estimate to get the exact effort level that gets you to InitDep.
	IEffort = MAX(0.0001,(BestRoughEffort-3))
	DO 10002 jj = 1,60
	  IEffort = IEffort + 0.1 !figure out the increment.         
        X(1) = IEffort
	  SS = FUNK(X)
	  IF (SS.LT.ObjFun) THEN
	    ObjFun = SS
	    BestInitEffort = X(1)
        ENDIF

10002 CONTINUE
C=======================================================

C	Profile, don't fit.  Stupid FORTRAN.
C      CALL FIT(X,SS,Nspec,CareyGRD)  
C     Carey: FIT finds the Qval that produces the specified value for depletion; 
C     FIT calls FUNK(X) to do this; X = QVal, so in the end SS = FUNK(X) just finds the value of the objective function (SS) that produced the
C     correct depletion.
	IF (ObjFunRough.LE.ObjFun) THEN
        BestInitEffort = BestRoughEffort
	ENDIF
	
      SS = FUNKSCALAR(BestInitEffort)
	WRITE(*,'(1x,A)') 'Actual BaseYr InitDep (InitDeplRecord), ObjFun'
	WRITE(*,*) InitDeplRecord, SS

	X(1) = BestInitEffort
C	After profiling you have an excellent seed:
      CALL FIT(X,SS,Nspec,CareyGRD)  
C     Carey: FIT finds the Qval that produces the specified value for depletion; 
C     FIT calls FUNK(X) to do this; X = QVal, so in the end SS = FUNK(X) just finds the value of the objective function (SS) that produced the
C     correct depletion.

	SS = FUNK(X)
	BestInitEffort = X(1)
      RETURN
      END SUBROUTINE FindInitEffort
      
C
C ===========================================================================
C
	
      SUBROUTINE GetEffLev(LastDepYear,DepYear,DepTarg)
C
C This subroutine finds the constant effort level 
C Carey is changing this for TechInteractions 8/6/14
	USE ShareDataRatioSubs
      IMPLICIT NONE
C
C     Local variables
      REAL*8 X(100),SS
      INTEGER Ispec
      INTEGER,INTENT(IN) :: LastDepYear, DepYear
      REAL(8),DIMENSION(Nspec),INTENT(IN) :: DepTarg
      REAL*8 FUNK,ObjFun,i,ii,jj,IEffort
      REAL(8) :: BestRoughEffort
	REAL(8) :: ObjFunRough, BestRoughDepl
      REAL*8 FUNKSCALAR
      EXTERNAL FUNK
      EXTERNAL FUNKSCALAR
!      StartYr = BaseYr+1
!      EndYr = FinDepYr
       DO 54990 Ispec = 1,Nspec
         FinalDep(Ispec) = DepTarg(Ispec) !change this for multiple species!
54990  CONTINUE       
       StartYr = LastDepYear
       EndYr = DepYear
       MinOpt = 3
	 CareyGRD = 1.2
C      X(1) = LOG(TotEffort(0)) !before SS3MPA stuff
	X(1) = TotEffort(LastDepYear-1) !Carey changed this from 1
	ObjFun = 1000000
      CALL FIT(X,SS,1,CareyGRD)
C     FIT finsds the constant effort from year 1 to year FinDepYr+1 which
C     minics the (pre-specified) final depletion. Note that this only works for ONE
C     species
 !     SS = FUNKSCALAR(X(1))
      SS = FUNK(X)

	ConstantEff = X(1)
      ConstantEffVec = X
	Write(*,'(1x,A)') 'Final Depletion, Constant effort for yrs , 
     + 1 to FinDepYr'
	Write(*,*) FinDepRecord, X(1)
	Write(51,'(1x,A)') 'Constant effort for yrs 1 to FinDepYr, 
     +	FinDepRecord'
	Write(51,'(2(F25.10))') X(1), FinDepRecord
C
      RETURN
      END SUBROUTINE GetEffLev

C
C ===========================================================================
C

	SUBROUTINE FindAvgMSY()
C This subroutine finds the effort that leads to MSY, MSY  
C
      USE ShareDataRatioSubs
      IMPLICIT NONE

	REAL(8), DIMENSION(100) :: X
	REAL(8) :: i,ObjFun
	REAL(8) :: IEffort,BMSY,SEED1,SEED2
      REAL(8) :: FUNK,SS
	EXTERNAL FUNK

	MinOpt = 9
	ObjFun = 1000000.0
	WriteOpt = 0

!	QVal = 1.0 !delete this and do the qval measurement. Qval is set from save.all to match what's being used elsewhere.
	Quota = 0.0
	Closed = 1

C      Read the seeds for this simulation 
	DO 22222 Isim = 1,Nsim
        READ(2,*) SEED1,SEED2
        ISEEDP = SEED1
        ISEEDO = SEED2
	  CALL GenProcErr()
22222 CONTINUE

      
	X(1) = LOG(TotEffort(StartFishing)) !change to log in a minute.
	CareyGRD = 1.2
      CALL FIT(X,SS,1,CareyGRD)
C     FIT finsds the constant effort from year 1 to year FinDepYr+1 which
C     minics the (pre-specified) final depletion. Note that this only works for ONE
C     species

!	WriteOpt = 1  !Write out the results of doing simulations with the average optimal effort


      SS = FUNK(X)
	EffortMSY = X(1)
      
      !Can't remember what this is about anymore (leftover from DR project)
	CALL WriteSpecs(14)
	CALL WriteSpecs(26)
	CALL WriteSpecs(38)
	CALL WriteSpecs(41)
	!CALL WriteSpecs(52)
      
      !Find the F that produces the MSY using new function:
      !AvgCatch = trial quota.

	Write(*,'(1x,A)') 'M,MaxAge,EffortMSY,AvgSSB,AvgDepletion,
     +	AvgCatch'
	Write(*,*) M(1),Amax(1),EffortMSY,AvgSSB,AvgDepletion,AvgCatch

	WRITE(52,'(1x,A)') 'M,MaxAge,EffortMSY,AvgSSB,AvgDepletion,
     +	AvgCatch'
	WRITE(52,700) M(1),Amax(1),EffortMSY,AvgSSB,AvgDepletion,AvgCatch

700	FORMAT(1x,F5.3,1x,I4,1x,F15.13,1x,F15.5,1x,F15.10,1x,F15.5)

	RETURN
      END SUBROUTINE FindAvgMSY
      
C      
C================================================================
C

      SUBROUTINE GenObsErr()
C
C This subroutine generates the density ratio survey errors that will be needed for the
C  simulation. NOT NEEDED for the techinteraction model

C	Call this for the years for which density ratio sampling happens and only
C	when SuperSim = 1
C
	USE ShareDataRatioSubs
      IMPLICIT NONE
C
C     Local variables
      INTEGER Year,Ispec
      REAL*8 OpenError,ClosedError
      REAL*8 XNORM
      EXTERNAL XNORM
C
C     Generate the recruitment error terms
      DO 10000 Ispec = 1,NSpec
       DO 10000 Year = MPAYear,(FinDepYr+NProjYr)
        OpenError = XNORM(2,0.0d0,SurvCV,ISEEDO)-SurvCV**2.0/2.0
	  ClosedError = XNORM(2,0.0d0,SurvCV,ISEEDO)-SurvCV**2.0/2.0 

         ObsResOpen(ISim,ISpec,Year) = EXP(OpenError)
	   ObsResClosed(ISim,ISpec,Year) = EXP(ClosedError)
       
10000 CONTINUE       
C
      RETURN
      END SUBROUTINE GenObsErr

C
C ===========================================================================
C
	
      REAL*8 FUNCTION FUNKSCALAR(X)
C
C This function is the objective function associated with MinOpt = 8. It was made so to deal with actual arguments that were scalars and not vectors.
C Regular FUNK deals with X as a vector.  You can't use actual (read-in arguments) that are scalars when the dummy argument (X) is not a scalar.
C http://software.intel.com/en-us/blogs/2009/03/31/doctor-fortran-in-ive-come-here-for-an-argument/
C to what is specified as initial depletion or final depletion
C Options:
C 	1. If MinOpt = 3, it calls ProjYr and loops over start to end year to calculate the final depletion level for an input E level for ONE species. Then it minimizes the SS difference between the observed depletion level and the specified depletion level. 
C 	2. If MinOpt = 8, it calls ProjForward and calculates the final depletion level for an input E level for ALL species. Then it minimizes the SS difference between the observed depletion level and the specified depletion level SUMMED across ALL species. 
C
	USE ShareDataRatioSubs
      IMPLICIT NONE
      
C     Global variables
      REAL*8 X,Depl !The difference between FUNK and FUNKSCALAR is here: 
C      X is a scalar and not a vector.  Thus, the actual argument matches the dimensions of the dummy argument (X).
C
C     Local variables
      REAL*8 Kvals(MaxSpec),CatPred(MaxSpec),Tmp,TheEffort(10),Penalty1
	REAL(8)::temp,temp1,temp2
      INTEGER Ispec,Iarea,II,Yr,IblkNS,IblkEW,IYr,Spec,Year
	REAL(8) :: SumDepl,SumSSB,SumCatch,ISEED1,ISEED2,SEED1,SEED2

C
C     Reset FUNKSCALAR
      FUNKSCALAR = 1.0E20 
C
C     How to proceed depends on MinOpt
      
      !--------------------------------------------------------      
C     If MinOpt.EQ.3 --> this objective function is the SS difference in final depletion for a specific E level for 1 species
      IF (MinOpt.EQ.3) THEN 
C      
       DO 31000 Yr = StartYr,EndYr  !Define in GetEffLev
C      
C       Impose MPA
        IF (Yr.GE.MPAYear) THEN
         DO 21100 Iarea = 1,NClosedAreas
          DO 21110 IblkNS = MinCNS(Iarea),MaxCNS(Iarea)
           DO 21110 IBlkEW = MinCEW(Iarea),MaxCEW(Iarea)
            Closed(IBlkNS,IBlkEW) = 0
21110     CONTINUE
 
21100    CONTINUE          
        ENDIF
C        TheEffort(1) = EXP(X(1))  !CHANGETHIS
        TheEffort(1) = X
        CALL projYr(Yr,TheEffort)
31000  CONTINUE       
       Depl = SPBIoTot(1,EndYr+1)/SpBioTot(1,MinPYr)*100
       FUNKSCALAR = 0
       IF (ABS(Depl-FinalDep(1)).GT.0.01) FUNKSCALAR = 
     +      (Depl-FInalDep(1))**2.0
C       WRITE(*,*) FUNKSCALAR,TheEffort(1),Depl
      FinDepRecord = Depl 
      ENDIF

      !--------------------------------------------------------      
C     If MinOpt.EQ.8 --> this objective function is the SS difference in final depletion for a specific E level SUMMED accross species
      IF (MinOpt.EQ.8) THEN 
C
C      Copy values across
        IF (X.LT.0) RETURN
C	  DO 33333 Yr = -19,BaseYr
      	  DO 33333 Yr = FirstFishYr,BaseYr
         TotEffort(Yr) = X
33333  CONTINUE       
C
C      Project forward 
C       Kvals = 10000.0
        Kvals = Kvec
       CALL ProjForward(Kvals)
C
C      Calculate the objective function
       FUNKSCALAR = 0
       DO 11001 ISpec = 1,Nspec
        Depl = SPBioTot(Ispec,BaseYr)/SpBioTot(Ispec,MinPYr)*100
C        WRITE(*,*) Depl,InitDep(Ispec)
        IF (ABS(Depl-InitDep(Ispec)).GT.0.01) FUNKSCALAR = 
     +		FUNKSCALAR+(Depl-InitDep(Ispec))**2.0
11001  CONTINUE
C       WRITE(*,*) FUNKSCALAR,Qval
	InitDeplRecord = Depl
C
      ENDIF

C========================================================================

      RETURN
      END FUNCTION FUNKSCALAR

C
C ===========================================================================
C

      SUBROUTINE MovementForSS3
      USE ShareDataRatioSubs
      IMPLICIT NONE

      !Local Variables
      INTEGER :: Ispec,Iarea,Jarea,IBlkNS,IBlkEW,JBlkNS,JBlkEW
      REAL(8),DIMENSION(Nspec,Narea,Narea) :: AreaMove1
      REAL(8),DIMENSION(Nspec,Narea) :: SumAreaMove

      AreaMove1 = 0
      
      DO 99 Ispec = 1,Nspec
      
      DO 100 Iarea = 1,Narea !These will be loops
      DO 100 Jarea = 1,Narea  !These will be loops
      DO 100 IBlkNS = MinNS(Iarea),MaxNS(Iarea)
       DO 100 IBlkEW = MinEW(Iarea),MaxEW(Iarea)
        DO 100 JBlkNS = MinNS(Jarea),MaxNS(Jarea)
         DO 100 JBlkEW = MinEW(Jarea),MaxEW(Jarea)
          AreaMove1(Ispec,Iarea,Jarea) = AreaMove1(Ispec,Iarea,Jarea) + 
     +      MixMat(Ispec,IBlkNS,IBlkEW,JBlkNS,JBlkEW)
100   CONTINUE

      SumAreaMove = 0
      DO 200 Iarea = 1,Narea
       DO 200 Jarea = 1,Narea
        SumAreaMove(Ispec,Iarea) = SumAreaMove(Ispec,Iarea) + 
     +    AreaMove1(Ispec,Iarea,Jarea)
200   CONTINUE

      !Normalize Area movement so that prob(moving from Iarea) sums to 1 (and includes pr(Iarea to itself)
      DO 300 Iarea = 1,Narea
       DO 300 Jarea = 1,Narea
        AreaMove(Ispec,Iarea,Jarea) = AreaMove1(Ispec,Iarea,Jarea)/
     +   SumAreaMove(Ispec,Iarea)
300   CONTINUE

99    CONTINUE


      END SUBROUTINE MovementForSS3

C
C ===========================================================================
C

       SUBROUTINE mydate(mystring)
C    This subroutine outputs the date and time of the simulation
C
       IMPLICIT NONE

      ! this is based on a website example: http://download.oracle.com/docs/cd/E19205-01/819-5259/aetcf/index.html
       integer date_time(8)
       character*10 b(3)
       character*25 mystring

       call date_and_time(b(1), b(2), b(3), date_time)

      ! WRITE(2001,*) 'End time and date'
        WRITE(2001,*) mystring
    !   WRITE(2001,300) date_time(5),date_time(6),date_time(7)
       WRITE(2001,301) date_time(2),date_time(3),date_time(1),
     +  date_time(5),date_time(6),date_time(7)
       WRITE(*,301) date_time(2),date_time(3),date_time(1),
     +  date_time(5),date_time(6),date_time(7)

     
!       print *,'date_time    array values:'
!       print *,'year=',date_time(1)
!       print *,'month_of_year=',date_time(2)
!       print *,'day_of_month=',date_time(3)
!       print *,'time difference in minutes=',date_time(4)
!       print *,'hour of day=',date_time(5)
!       print *,'minutes of hour=',date_time(6)
!       print *,'seconds of minute=',date_time(7)
!       print *,'milliseconds of second=',date_time(8)
!       print *, 'DATE=',b(1)
!       print *, 'TIME=',b(2)
!       print *, 'ZONE=',b(3)

 300   FORMAT('Time: ',I2,':',I2,':',I2)
 301   FORMAT('Date:',I2,'/',I2,'/',I4,', Time: ',I2,':',I2,':',I2)

      END SUBROUTINE mydate

C
C ===========================================================================
C

      SUBROUTINE Movement2ForSS3(TempN)
      USE ShareDataRatioSubs
      IMPLICIT NONE

      !This is put in the model at the same time as movement to measure the proportions of animals
      !of age a moving from one assessment area to another; hence it's an idea of what the movement
      !function looks like for the operating model.
      !Local Variables

      REAL(8),DIMENSION(Nspec,Nsex,NNSBlock,NEWBlock,0:Maxage) :: TempN
      INTEGER :: Ispec,Isex,Iage,Iarea,Jarea,IBlkNS,IBlkEW,JBlkNS,JBlkEW
      !REAL(8),DIMENSION(Nspec,Narea,Narea) :: AreaMove1
      !REAL(8),DIMENSION(Nspec,Narea) :: SumAreaMove
      REAL(8),DIMENSION(Narea,Narea,0:Maxage)::MyN
      Real(8),DIMENSION(0:Maxage,Narea)::SumN

      Ispec = 1
      Isex = 1
      MyN = 0
      
  !    DO 99 Ispec = 1,Nspec
  !    DO 99 Isex = 1,Nsex
      DO 99 Iage = 0,Amax(Ispec)

      DO 100 Iarea = 1,Narea !These will be loops
      DO 100 Jarea = 1,Narea  !These will be loops
      DO 100 IBlkNS = MinNS(Iarea),MaxNS(Iarea)
       DO 100 IBlkEW = MinEW(Iarea),MaxEW(Iarea)
        DO 100 JBlkNS = MinNS(Jarea),MaxNS(Jarea)
         DO 100 JBlkEW = MinEW(Jarea),MaxEW(Jarea)
    !      AreaMove1(Ispec,Iarea,Jarea) = AreaMove1(Ispec,Iarea,Jarea) + 
    ! +      MixMat(Ispec,IBlkNS,IBlkEW,JBlkNS,JBlkEW)
           MyN(Iarea,Jarea,Iage) = MyN(
     +      Iarea,Jarea,Iage) + 
     +      MixMat(Ispec,IBlkNS,IBlkEW,JBlkNS,JBlkEW)*
     +      TempN(Ispec,Isex,IBlkNS,IBlkEW,Iage)

100   CONTINUE

      SumN(Iage,:) = 0
      DO 200 Iarea = 1,Narea
       DO 200 Jarea = 1,Narea
        SumN(Iage,Iarea) = SumN(Iage,Iarea) + MyN(Iarea,Jarea,Iage)
200   CONTINUE

      !Normalize Area movement so that prob(moving from Iarea) sums to 1 (and includes pr(Iarea to itself)
      DO 300 Iarea = 1,Narea
       DO 300 Jarea = 1,Narea
 !       AreaMove(Ispec,Iarea,Jarea) = AreaMove1(Ispec,Iarea,Jarea)/
 !    +   SumAreaMove(Ispec,Iarea)
      
        PropN(Iarea,Jarea,Iage) = 
     +   MyN(Iarea,Jarea,Iage)/SumN(Iage,
     +   Iarea)

300   CONTINUE
99    CONTINUE


      END SUBROUTINE Movement2ForSS3

C
C ===========================================================================
C

      SUBROUTINE ReadCabRep(Year1,Year2,INDC,Iarea1,Iarea2,Carea1,Ispec,
     +Iarea1a,Iarea2a,Icon)
      USE ShareDataRatioSubs
      IMPLICIT NONE

      !Attn: this would have to be updated to accomodate multiple species
C      Read the output from each .rep file from one of each of multiple configs (but right now there's just one):
C      Create a loop over running cab and reading in rep files later on.
C      Right now just read in a cab file, as before.      

      !Global Variables
      INTEGER :: Icon
      !Local Variables
      INTEGER::Istuff,Yr,Year1,Year2,Isex,II,Ispec,OutUnit,Iarea1,Iarea2
      INTEGER :: Iarea1a,Iarea2a,Iarea,Iage,Carea1
      CHARACTER :: INDC
      CHARACTER(4) :: FirstStuff
      CHARACTER(18) :: SecondStuff
      CHARACTER(6) :: w1
      CHARACTER(4) :: w2,w3,w4
      REAL(8),DIMENSION(6) :: Pred
      REAL(8) :: SpTru,ExpTru,B0tru
      INTEGER::TheYr,t1,t2,t3
      
      OPEN(UNIT=957,FILE='results/SUMM1.EST',ACCESS='APPEND')
      OPEN(UNIT=958,FILE='results/SUMM2.EST',ACCESS='APPEND')
      OPEN(UNIT=959,FILE='results/REC.EST',ACCESS='APPEND')
      IF (INDC.EQ."B") THEN !Carey changed from "C" to "B" to make consistent with WSGReport and numbers
       OutUnit = 958
      ELSE 
       OutUnit = 957
      ENDIF  
      
      IF (Year2.EQ.FinDepYr.AND.Icon.EQ.0) THEN
        !WRITE(OutUnit,*) Year2-Year1+1,Ispec,Iarea1,Iarea2
         WRITE(OutUnit,*) 'Yr,Ispec,Icon,TruSSB,EstSSB,TruDepl,EstDepl'
         WRITE(OutUnit,*) Year2,Nspec,Nconfigs
      ENDIF
      B0tru = 0
      
      DO 50300 Iarea = Iarea1a,Iarea2a
       B0Tru = B0Tru + SpBioArea(Ispec,Iarea,Year1)
50300 CONTINUE    

 !     DO 60300 Yr = Year1,Year2+1
 !       Sptru = 0
 !       ExpTru = 0
 !       DO 50110 Iarea = Iarea1a,Iarea2a
 !         SpTru = SpTru + SpBioArea(Ispec,Iarea,Yr)
 !         ExpTru = ExpTru + TruExp(ISpec,Iarea,Yr)
!50110   CONTINUE  
!60300 CONTINUE
      
      
      OPEN(UNIT=99,FILE='CAB.REP')
!      DO 50110 Yr=Year1,Year2+1
!        Sptru = 0
!        ExpTru = 0
!        DO 50110 Iarea = Iarea1a,Iarea2a
!          SpTru = SpTru + SpBioArea(Ispec,Iarea,Yr)
!          ExpTru = ExpTru + TruExp(ISpec,Iarea,Yr)
!50110   CONTINUE  
          
        READ(99,*) 
        READ(99,*) cab_obj_fun(1)
       
       
         Do 440 Istuff = 1,10000
           READ(99,*) FirstStuff
C           IF (FirstStuff.EQ.'Year, SpBio(Year), SpBio(Year)/SpBio(styr)
C     +,BioOut, Recruit, Hrate ') THEN
           IF (FirstStuff.EQ.'Year') THEN
!            PRINT 'Exiting loop to find the first line of the rep file'
             EXIT
          ENDIF
440      CONTINUE
         
C      New cab meaning of Pred is as follows (June 10, 2014)
C     Pred(1) = Year, Pred(2) =estimated SpBio(Year), Pred(3) = Est Depletion(Year)
C     Pred(4) = Est TotBio(Year), Pred(5) = EstRecruitment(Year), Pred(6) = Harvest rate(Year)
        DO 50100 Yr=Year1,Year2+1
         READ(99,*) (Pred(II),II=1,6) 
         IF (Pred(1).GT.9000000) Pred(1) = 9000000
         IF (Pred(2).GT.9000000) Pred(2) = 9000000
         IF (Pred(4).GT.9000000) Pred(4) = 9000000
         IF (Pred(6).GT.9000000) Pred(6) = 9000000
C       
C     This next "if" statement will need to be modified when testing "B" to match teh current meaning of Pred
         IF (INDC.EQ."B") THEN !Carey changed "C" to "B" to be consistent with WSG report and numbers
           PassEstVar(Ispec,Yr,1) = PassEstvar(Ispec,Yr,1) + SpTru
           PassEstVar(Ispec,Yr,2) = PassEstvar(Ispec,Yr,2) + Pred(1)
           PassEstVar(Ispec,Yr,5) = PassEstvar(Ispec,Yr,5) + ExpTru
           PassEstVar(Ispec,Yr,6) = PassEstvar(Ispec,Yr,6) + Pred(6)
         ENDIF 
         !Set up for old CAB
C         WRITE(OutUnit,801) Yr,SpTru,Pred(2),
C     +    SpTru/B0Tru,Pred(3),ExpTru,Pred(6) !what is ExpTru?
C         EstSB(Ispec,Carea1,Yr) = Pred(2)
C         EstRec(Ispec,Carea1,yr) = Pred(5)
         
        Sptru = 0
        ExpTru = 0
        DO 50110 Iarea = Iarea1a,Iarea2a
          SpTru = SpTru + SpBioArea(Ispec,Iarea,Yr)
          ExpTru = ExpTru + TruExp(ISpec,Iarea,Yr) !ExpTru not used anymore
50110   CONTINUE

         WRITE(OutUnit,801) Yr,Ispec,Icon,SpTru,Pred(2),
     +    SpTru/B0Tru,Pred(3) !removed ExpTru (no longer a pred column to match); add back in later if needed.
         EstSB(Ispec,Carea1,Yr) = Pred(2)
         EstRec(Ispec,Carea1,Icon,yr) = Pred(5)

         WRITE(959,802) Yr,Ispec,Icon,Pred(5) !Write out to Rec.est: recruitment
50100   CONTINUE   !End loop over years
      
      !Need some lines to skip down to selectivity-at-age reporting in cab.rep
         Do 541 Istuff = 1,10000
           READ(99,*) SecondStuff

C           IF (FirstStuff.EQ.'Year, SpBio(Year), SpBio(Year)/SpBio(styr)
C     +,BioOut, Recruit, Hrate ') THEN
           
           
           !This needs to be split up into reading in fleet sex age selA as separate dummy variables'
           !just need SecondStuff to equal selA!
           
           IF (SecondStuff.EQ.'Selectivity-at-age') THEN
!           PRINT 'Exiting loop to find the first line of the rep file' 
           READ(99,*) w1,w2,w3,w4
             IF (w4.EQ.'selA') THEN
               EXIT
             ENDIF
          ENDIF
541      CONTINUE   
           
       !read in selectivity-at-age (selA)
      !nope, not gonna work like this.
      DO 542 IAge=0,AmaxEst(Ispec)
        READ(99,*) t1,t2,t3,EstS(Ispec,Carea1,1,IAge,Icon)
542   CONTINUE
      DO 543 IAge=0,AmaxEst(Ispec)     
        READ(99,*) t1,t2,t3,EstS(ISpec,Carea1,2,IAge,Icon)
543   CONTINUE              

      !Need some lines to skip down to numbers at age section
         Do 441 Istuff = 1,10000
           READ(99,*) FirstStuff
C           IF (FirstStuff.EQ.'Year, SpBio(Year), SpBio(Year)/SpBio(styr)
C     +,BioOut, Recruit, Hrate ') THEN
           IF (FirstStuff.EQ.'Age ') THEN
!            PRINT 'Exiting loop to find the first line of the rep file'
             EXIT
          ENDIF
441      CONTINUE        
        
      DO 50200 Yr=Year1,Year2
       READ(99,*) TheYr,(EstN(ISpec,Carea1,1,Yr,IAge),
     +IAge=0,AmaxEst(Ispec)) !EstN is not currently used. EstNatAge is read in in the MVN routine from cab.std (same info) and that is used in the MVN routine.  (CRM 7/28/14)
50200 CONTINUE
       READ(99,*)
       READ(99,*)
      DO 50201 Yr=Year1,Year2
       READ(99,*) TheYr,(EstN(Ispec,Carea1,2,Yr,IAge),IAge=0, 
     +AmaxEst(Ispec)) !EstN is not currently used. EstNatAge is read in in the MVN routine from cab.std (same info) and that is used in the MVN routine (CRM 7/28/14)
50201 CONTINUE
      
! This info is not part of the rep file anymore!
! SteepCAEst, SigmaRCAest, R0CAEst,S0CAEst would need to be read in from elsewhere.
!Only S0CAEst is used.      
 !     READ(99,*) SteepCAEst(ISpec,Carea1),SigmaRCAest(ISpec,Carea1),
 !    +           R0CAEst(ISpec,Carea1),S0CAEst(ISpec,CArea1)
 !     SigmaRCAest(Ispec,Carea1) = SQRT(SigmaRCAEst(ISpec,Carea1))
      
      CLOSE(99)
      CLOSE(957)
      CLOSE(958)
      CLOSE(959)
 
801   FORMAT(3(1x,I5),1x,2(F15.5,1x),2(F15.6,1x))
802   FORMAT(3(1x,I5),1x,F20.9,1x)
      END SUBROUTINE ReadCabRep

C
C ===========================================================================
C

      SUBROUTINE ReadPar(Ispec,Iarea,Icon)
      USE ShareDataRatioSubs
      IMPLICIT NONE
      
      !This subroutine reads in a single par file
      
      !Global variables
      INTEGER :: Ispec,Iarea,Icon

      !Local Variables
      CHARACTER(6) :: parstuff
      CHARACTER :: hash
      REAL(8), DIMENSION(Nspec) :: ln_S0
      INTEGER :: Ipar
      
      
      !Open the par file
      OPEN(Unit=50,FILE = 'cab.par')
      
      !read in the parameters of interest
      DO 100 Ipar = 1,200
          READ(50,1000) hash,parstuff
          
          IF (parstuff.EQ.'ln_S0:') THEN
            READ(50,*) ln_S0(Ispec)
            EXIT
          ENDIF
100   CONTINUE
          !formats
      
      S0CAEst(ISpec,Iarea,Icon) = EXP(ln_S0(Ispec))/Nsex
      
      CLOSE(50) !Close the par file.
1000  FORMAT(A1,1x,A6)      
      END SUBROUTINE ReadPar

C
C ===========================================================================
C

      SUBROUTINE IterateNvecsMVN(Ispec,Iarea,Iconfig,Iyear)
C     This function was written to sample vectors of log(N-at-age) from a multivariate normal (MVN) distribution
C     As a quick way to generate a distribution for numerically calculating a P* control rule
C     An alternative to this method would be to do MCMC to generate parameter vectors (including N-at-age)
C     MVN is a program written by John Uebersax to generate multivariate normal random numbers, given means, stds, and a correlation matrix.
      
      USE ShareDataRatioSubs
      IMPLICIT NONE
      !Global variables
      INTEGER :: Ispec,Iarea,Iconfig,Iyear
      
      !Local variables
      INTEGER :: varnum,FirstVarNum
      REAL(8) :: value,std
      CHARACTER(12) :: varname
      REAL(8),DIMENSION(1:Nsex,0:Amax(Ispec)) :: EstNatAge,EstStdNatage
      REAL(8),DIMENSION(1:Nsex,0:Amax(Ispec),0:200) :: EstNatAgeCor
      CHARACTER(12) :: stuff1
      CHARACTER(11) :: stuff2
      
      INTEGER :: Isims,Isex
      INTEGER :: Iage,Ipar,Iagevec,Icolvec
      !These to be added to the rest of the code later (assessment config and assessment year+1):
      INTEGER :: NumAges
      
      !This probably needs to be global?? Yep
 !     REAL(8),DIMENSION(1:Nspec,1:Nmvn,0:Amax(Ispec)) :: Nvec
      EstNatAge = 0.0
      EstStdNatage = 0.0
      EstNatAgeCor= 0.0
      
      NumAges = AINT((Amax(Ispec) + 1)/12.0)
      IF (MOD((Amax(Ispec) + 1),12).GT.0) THEN
        NumAges = NumAges + 1
      ENDIF
      !------------------------------------------------------
      !Get mean and std information from cab.std
      OPEN(Unit = 551,FILE = 'cab.std')
      
      READ(551,*) varname
      DO 100 Ipar = 1,200
          READ(551,1001) varnum,varname,EstNatAge(1,0),
     +EstStdNatage(1,0)
          IF (varname.EQ.'LogTheNatAge') THEN
              BACKSPACE(551)
              FirstVarNum = varnum
              DO 101 Isex = 1,Nsex
              DO 101 Iage = 0,Amax(Ispec)
                  READ(551,1001) varnum,varname,EstNatAge(Isex,Iage),
     +EstStdNatage(Isex,Iage)
101          CONTINUE !end Iage loop
              EXIT
          ENDIF !varname.EQ.LogT
100   CONTINUE !end loop reading through std file
      
      CLOSE(551) !Close cab.std
      
      !-------------------------------------------------------
      !Get correlation information
      
      EstNatAgeCor = 0.0
      OPEN(Unit = 552,FILE = 'cab.cor')
      READ(552,*) varname
      READ(552,*) varname
      DO 102 Ipar = 1,200
          READ(552,1003) varnum,varname !,value,std,
!     +(EstNatAgeCor(0,Iagevec),Iagevec=0,varnum-1) !this should be 200, not 2, but 200 doesn't work.
 
          IF (varname.EQ.'LogTheNatAge') THEN
              BACKSPACE(UNIT = 552)
              DO 103 Isex = 1,Nsex
              DO 103 Iage = 0,Amax(Ispec)
!                  READ(552,1002) varnum,varname,value,std,
!     +(EstNatAgeCor(Isex,Iage,Iagevec),Iagevec=0,varnum-1)
                  READ(552,*) varnum,varname,stuff1,stuff2,
     +(EstNatAgeCor(Isex,Iage,Iagevec),Iagevec=0,varnum-1)
103          CONTINUE !end Iage loop
              EXIT
          ENDIF !varname.EQ.LogT
102   CONTINUE !end loop reading through std file
      
      CLOSE(552)
      
      !------------------------------------------------------
      !Generate vectors of N-at-age for this configuration
      !Write the text file for MVN.exe and then call MVN.exe:
      OPEN(Unit = 53,File = 'mvninputFem.txt')
        WRITE(53,1004) Amax(Ispec)+1 !Number of variables
        WRITE(53,1005) Nmvn !Number of random vectors (need to add to input file)
        WRITE(53,1009) 500 !Random number seed - we want this to be the same for each scenario (the means, stds, cor will change)
        WRITE(53,1006) 1    !means are supplied
        WRITE(53,1006) 1    !stds are supplied
        WRITE(53,1006) 1    !correlations are supplied
        
        !write out mean std of log-N-at-age
        WRITE(53,1007) (EstNatAge(1,Iage),Iage=0,Amax(Ispec))
        WRITE(53,1007) (EstStdNatAge(1,Iage),Iage=0,Amax(Ispec))
      
        !DO 104 Isex = 1,Nsex
        DO 104 Iage = 1,Amax(Ispec)
            WRITE(53,1007) (EstNatAgeCor(1,Iage,Icolvec),Icolvec = 
     +(FirstVarNum-1),((FirstVarNum-1)+(Iage-1)))
104     CONTINUE
        
!        !DO 1041 Isex = 1,Nsex
!            DO 1041 Iage = 0,Amax(Ispec)
!              WRITE(53,1007) (EstNatAgeCor(2,Iage,Icolvec),Icolvec = 
!     +(FirstVarNum-1),((FirstVarNum-1)+Amax(Ispec)+1 +(Iage-1)))
!1041    CONTINUE
      CLOSE(53) !Close input file cabmvn.txt (input for MVN program)
      CALL SYSTEM("(echo mvninputFem.txt& echo mvnoutputFem.txt&
     + echo .& echo .) | mvn") !Draw random log(N-at-age) vectors from the multivariate normal distribution.
      !------------------------------------------------------
      
      !Generate male vector of random numbers:
      !------------------------------------------------------
      !Generate vectors of N-at-age for MALES for this configuration
      !Write the text file for MVN.exe and then call MVN.exe:
      OPEN(Unit = 53,File = 'mvninputMale.txt')
        WRITE(53,1004) Amax(Ispec)+1 !Number of variables
        WRITE(53,1005) Nmvn !Number of random vectors (need to add to input file)
        WRITE(53,1009) 500 !Random number seed - we want this to be the same for each scenario (the means, stds, cor will change)
        WRITE(53,1006) 1    !means are supplied
        WRITE(53,1006) 1    !stds are supplied
        WRITE(53,1006) 1    !correlations are supplied
        
        !write out mean std of log-N-at-age
        WRITE(53,1007) (EstNatAge(2,Iage),Iage=0,Amax(Ispec))
        WRITE(53,1007) (EstStdNatAge(2,Iage),Iage=0,Amax(Ispec))
      
!        !DO 104 Isex = 1,Nsex
!        DO 104 Iage = 1,Amax(Ispec)
!            WRITE(53,1007) (EstNatAgeCor(1,Iage,Icolvec),Icolvec = 
!     +(FirstVarNum-1),((FirstVarNum-1)+(Iage-1)))
!104     CONTINUE
        
        !Test this line carefully
            DO 1041 Iage = 1,Amax(Ispec)
              WRITE(53,1007) (EstNatAgeCor(2,Iage,Icolvec),Icolvec = 
     +(FirstVarNum+Amax(Ispec)),((FirstVarNum+Amax(Ispec))+(Iage-1)))
1041    CONTINUE
      CLOSE(53) !Close input file cabmvn.txt (input for MVN program)
      CALL SYSTEM("(echo mvninputMale.txt& echo mvnoutputMale.txt&
     + echo .& echo .) | mvn") !Draw random log(N-at-age) vectors from the multivariate normal distribution.
      !------------------------------------------------------      

      !Test that MVN created the needed files ok
      OPEN(Unit = 54,File = "mvnoutputFem.txt")
      READ(54,*)
      READ(54,*) TestMVN
      CLOSE(54)
      IF (TestMVN.EQ.'Correlation') THEN
        RETURN    
      ENDIF
  
      !Test that MVN created the needed files ok
      OPEN(Unit = 54,File = "mvnoutputMale.txt")
      READ(54,*)
      READ(54,*) TestMVN
      CLOSE(54)
      IF (TestMVN.EQ.'Correlation') THEN
        RETURN    
      ENDIF
      
      !Needs to be updated to read in female then male information into Nvec that accepts females and males separately.
      !Read in the vectors of N values
      OPEN(Unit = 54,File = 'mvnoutputFem.txt')
      DO 105 Isims = 1,Nmvn !program in this number - this has to stop now!
        READ(54,1008) (Nvec(Ispec,1,Isims,Iage,Iconfig),
     +Iage=0,Amax(Ispec)) !this is a prob bc only reads in to the end of a line, not all 31 obs for each Nvec.
105   CONTINUE
      CLOSE(54)

      OPEN(Unit = 54,File = 'mvnoutputMale.txt')
      DO 1051 Isims = 1,Nmvn 
        READ(54,1008) (Nvec(Ispec,2,Isims,Iage,Iconfig),
     +Iage=0,Amax(Ispec)) !this is a prob bc only reads in to the end of a line, not all 31 obs for each Nvec.
1051  CONTINUE
      
      CLOSE(54) !close mvnoutput.txt file
      
      !Calculate current SSB for each vector
      !Iconfig = 1 is the last dimension of SSB_MVNvec
      DO 1061 Isims = 1,Nmvn
          DO 1061 Iage = 0,Amax(Ispec)
             SSB_MVNvec(Ispec,Isims,Iconfig) = SSB_MVNvec(Ispec,Isims,
     +Iconfig) +
     +exp(Nvec(Ispec,1,Isims,Iage,Iconfig))*Fec(Ispec,Iage)
1061  CONTINUE
      !Run minopt for finding the effort required to meet a trial quota, but run it on the dynamics with each Nvec with the output being estF(Ivec)
      !make a subroutine that does this, taking in Est Selex, Nvec, M, and a trial ABC
      
      !Likewise run minopt for finding the FOFL associated with each vector (do you mean FOFL or FMSY?) FOFL for Tier 3; ow could use a pdf of FOFL for this?
      
      
      
      !Shertzer 2008 numerical integration to find Pr(Ft > FOFL) = p* for estimated Ft and FOFL
      

      
      !------------------------------------------------------
1001  FORMAT(4x,I2,3x,A12,1x,2(ES11.4E4,1x))
1002  FORMAT(I6,3x,A12,1x,ES12.4E4,ES12.4E4,2x,200(F7.4,1x))
1003  FORMAT(I6,3x,A12)
      
1004  FORMAT(I2,15x,"#Number of variables")
1005  FORMAT(I5,12x,"#Number of random vectors")
1006  FORMAT(I2,15x, "#means, stds, and correlations are supplied")
1007  FORMAT(100(F30.27,1x)) !might need to back down on decimals; we'll see.
1008  FORMAT(<NumAges>(1x,12(F8.4,1x),/))
1009  FORMAT(I5,12x,"#Random Number Seed")
        END SUBROUTINE IterateNvecsMVN 

C
C ===========================================================================
C

      SUBROUTINE ReadCR()
      USE ShareDataRatioSubs
      IMPLICIT NONE
      
      !Local variables:
      INTEGER :: Ispec,Iconfig
      !Call this up near ReadSpec and ReadCntl
      !Read in the harvest control rule specs
      OPEN(Unit = 540,File = "CR.dat")
      READ(540,'(40X,I5)') UsePstar
      
      READ(540,*)
      READ(540,*)
      READ(540,('(40X,F4.2)')) Pstar
      READ(540,('(40X,I5)')) Nmvn
      READ(540,('(40X,I5)')) Nconfigs

      READ(540,*)
      DO 1278 Ispec = 1,Nspec
        READ(540,*)
        READ(540,('(40X,15(F5.3,1X))')) 
     +   (EstM(Ispec,Iconfig),Iconfig = 1,Nconfigs)
        READ(540,('(40X,15(F5.3,1X))')) 
     +   (IEstSteep(Ispec,Iconfig),Iconfig = 1,Nconfigs)
        READ(540,('(40X,15(F5.3,1X))')) 
     +   (ConfigWts(Ispec,Iconfig),
     +Iconfig = 1,Nconfigs)
1278  CONTINUE
      
      READ(540,*)
      READ(540,*)
      READ(540,('(40X,F4.2)')) AFSC_alpha
      READ(540,('(40X,I5)')) AFSC_MinusYrs !Not being used yet. minus this many rs for avgrec calc in FOFL control rule to find B40.
      CLOSE(540)
      RETURN
      END SUBROUTINE ReadCR

C
C ===========================================================================
C

      SUBROUTINE DummyCR(Ispec)
C This subroutine is to assign a dummy catch that is fed in to the OM for the following year
C It's a dummy control rule for testing purposes.      
      USE ShareDataRatioSubs
      IMPLICIT NONE
      
      INTEGER :: Ispec

      Quota(Ispec) = 500
      
      END SUBROUTINE DummyCR

C
C ===========================================================================
C

      SUBROUTINE Tier3CR()
C AFSC Tier 3 control rule:
      USE ShareDataRatioSubs
      IMPLICIT NONE

      !Write out the input files for Jim I's projection model
      !Run the projection model
      !read in and grab the ABC, OFL, B40, F40, whatever else
      !continue with MSE.
      
      
      
      
      END SUBROUTINE Tier3CR

C
C ===========================================================================
C

      SUBROUTINE CalcCatchFromF(Fval,Mval,TheNs,SelexVec,WeightA)
C Subroutine to calculate F given a trial ABC to be used for the P* approach HCR
C Attention that CRM is only making this work with a catch limit applied to one aggregated area.
C This function is called in the FUNK2 function (for MinOpt = 2)
C
      USE ShareDataRatioSubs
      IMPLICIT NONE

      !CHANGE NEEDS TO BE MADE SO THIS IS JUST FOR 1 RANDOM VECTOR - separate optimization for each vector.
      !GlobalVariables
      REAL(8), DIMENSION(1:Nsex,0:Maxage) :: TheNs
      REAL(8), DIMENSION(1:Nsex,0:Maxage) :: Selexvec
      REAL(8), DIMENSION(1:Nsex,0:Maxage) :: WeightA
      REAL(8) :: Mval,Fval
      !Local variables - dimensions need to be expanded to include configurations (M's h's etc).
!      INTEGER :: Nrands
      INTEGER :: Iage,Itime,Isex
      REAL(8), DIMENSION(1:Nsex,0:Maxage) :: pCatAge
!      REAL, DIMENSION(1:Nspec,1:Nmvn) :: pCatBio  !This needs to be global to access in FUNK.
      
! Might abandon doing this within the same function.
!      !----------------------------------------
!      Base calculations on true (OM) or estimated (EM) values?
!       IF (TrueOrEst.EQ.0) THEN
!         NAge =       
!       ENDIF
!      !----------------------------------------
      
      !get rid of Ivec dimension in catch statistics?
      pCatAge(:,:) = 0.0
      pCatBio = 0.0
 !     DO 100 Ivec = 1,Nmvn
          DO 100 Isex = 1,Nsex
             DO 100 Iage = 0,Maxage
          pCatAge(Isex,Iage) = TheNs(Isex,Iage)*
     +    *exp(-Mval/2)
     +    *(1-exp(-SelexVec(Isex,Iage)*Fval))
          
          pCatBio = pCatBio + pCatAge(Isex,Iage)*
     +              WeightA(Isex,Iage)!Assumes Itime = 1 for Wght (I've never run the model under other settings (CRM)) 
          
100   CONTINUE

      END SUBROUTINE CalcCatchFromF

C
C ===========================================================================
C

      SUBROUTINE FindFforCatchReg(IYear,Ispec)
C Uses CalcCatchFromF (called from MinOpt.EQ.2 in FUNK2) to find the F that results in a particular catch regulation
C given a vector of Numbers-at-age and values for M, selex, and weight-at-age.
C Attention: this was called previously in FUNK but CRM modified it to use FUNK2
C      
      USE ShareDataRatioSubs
      IMPLICIT NONE
      
       !Global variables
       INTEGER :: IYear,Ispec
       
       
       !Local variables
       INTEGER :: Iconfig,Ivec,Isex,Iage
       REAL(8) :: SS
       REAL(8),DIMENSION(100) :: X
       REAL(8) :: FUNK
       EXTERNAL FUNK
       
       MinOpt = 10
       CareyGRD = 1.2
       X(1) = 0.1
       
       !Loop over Ispec,Iconfig,Ivec and save the output from each individual optimization
!       DO 100 Ispec = 1,Nspec !Ispec loop is now at the P* subroutine level and is entered into this subroutine 7/28/14
           MyWeight(:,:) = Wght(Ispec,:,:,1) 
           DO 100 Iconfig = 0,Nconfigs
               !ADD A LINE HERE TO SKIP iCONFIG IF NO COR FILE; SET fIDSTN = -1.0, ssfdISTN = -1.0
               TheM = EstM(Ispec,Iconfig)
               DO 101 Ivec = 1,Nmvn
       
       !Define the inputs that will be used within FUNK (for CalcCatchFromF) as global variables that are reset for each iteration of this big loop
                 DO 2117 Isex = 1,Nsex
                   DO 2117 Iage = 0,Amax(Ispec)
                    ThoseNs(Isex,Iage) = exp(Nvec(Ispec,Isex,Ivec,
     +                                     Iage,Iconfig)) !Vector specific
2117             CONTINUE !End Isex and Iage loop
                 TheSelex = EstS(Ispec,1,:,:,Iconfig) !Not Ivec specific right now, but would be under MCMC
                   
       !Here you should cal FIT.
       CALL Fit(X,SS,1,CareyGRD)
       SS = FUNK(X) !X is the initial guess for F that produces a pre-specified amount of catch for the p* approach; 
      
       !make a variable to save F(t)'s and SS's for each part of the loop (so that you can check that they all worked)
       FDistn(Ispec,Iyear,Iconfig,Ivec) = X(1)
       SSFDistn(Ispec,IYear,Iconfig,Ivec) = SS
      
101            CONTINUE !End Ivec loop
100     CONTINUE !end loop over Ispec,Iconfig
       !Get the true F(t) that matches the trial ABC.
       MinOpt = 11

!       DO 201 Ispec = 1,Nspec !CRM 7/28/14 taken out
         X(1) = 0.1 
         MyWeight(:,:) = Wght(Ispec,:,:,1) 
         TheM = M(Ispec)
         MyNs = NatAgeCompact(:,:,IYear,:,0)
         ThoseNs = MyNs(Ispec,1:Nsex,0:Maxage)
         TheSelex(1,:) = Sel(Ispec,:)
         TheSelex(2,:) = Sel(Ispec,:)
              
         CALL Fit(X,SS,1,CareyGRD)
         SS = FUNK(X)
         
         !Save true F(t) information (pCatBio and F(t) = X(1)
         TrueF(Ispec,IYear) = X(1)
         SSTrueF(Ispec,IYear) = SS !As a check.
!201   CONTINUE !End Ispec loop        

      RETURN       
      END SUBROUTINE FindFforCatchReg

C
C ===========================================================================
C

!      SUBROUTINE AFSC_SPR(Fval,Mval,Fecvec,Recvec,Ispec,Iconfig)
      SUBROUTINE AFSC_SPR(Fval,Mval,Fecvec,Recvec,SelexVec,PlusAge)
C Find the SPR rate and SSB_SPR given an F (Fval) and other population numbers
C Using the average recruitment from 1978 to IYear - 2 to calc B_SPR is what makes this function specific to AFSC
C This is used within a FIT and FUNK function to calculate F_35% and B_40% - 2 diff. calls to this.
C Attention: Need to be careful because CRM uses FUNK2 for some of the operations instead of FUNK so need to be sure that we change this for the tech int model
C   
      USE ShareDataRatioSubs
      IMPLICIT NONE
      
      !Global variables
      REAL(8),DIMENSION(100) :: Fval  !Because this is X, but only 1 dim of X is used X(1)
      REAL(8) :: Mval
      REAL(8), DIMENSION(0:Maxage) :: Fecvec
      REAL(8), DIMENSION(1:GlobalProjYr-1) :: Recvec !Things starting in year 1 needs to change - things need to start in 1978.
!      INTEGER :: Ispec,Iconfig !Can't enter these directly when AFSC_SPR is called from FUNK, but it will be both Ispec, Iconfig-specific and year specific and maybe Ivec specific.
      REAL(8), DIMENSION(0:Maxage) :: SelexVec
      INTEGER :: PlusAge
      
      !Local variables
      INTEGER :: Iage,ii,Iyr
      REAL(8), DIMENSION(0:Maxage) :: NperRec0,NperRecF
      REAL(8) :: EstAvgRec,SBPR0,SBPR_F
      
      SBPR0 = 0.0
      SBPR_F = 0.0
      EstSPR = 0.0
      EstB_SPR = 0.0
      
!Calculate SBPR0
      NperRec0(0) = 1
      DO 1000 Iage = 1,PlusAge
        NperRec0(Iage) = NperRec0(Iage-1)*exp(-Mval)
1000  CONTINUE
      
      NperRec0(PlusAge) = NperRec0(PlusAge)/(1-exp(-Mval))
      
      DO 1001 Iage = 0,PlusAge
        SBPR0 = SBPR0 + NperRec0(Iage)*Fecvec(Iage)
1001    CONTINUE

!Calculate SBPR(F)
      NperRecF(0) = 1
      DO 1002 Iage = 1,PlusAge
          NperRecF(Iage) = NperRecF(Iage-1)*exp(-Mval-SelexVec(Iage-1)*
     +Fval(1))
1002  CONTINUE
          NperRecF(PlusAge) = NperRecF(PlusAge)/
     +(1-exp(-Mval-SelexVec(Iage)*Fval(1)))
         
      DO 1003 Iage = 0,PlusAge
          SBPR_F = SBPR_F + NperRecF(Iage)*Fecvec(Iage)
1003  CONTINUE

!Calculate SPR(Fval):
      EstSPR = SBPR_F/SBPR0 !define globally
      
      !Calculate B
      ii = 0
      EstAvgRec = 0.0
      DO 1004 Iyr = 1,(GlobalProjYr-1) !This will change when using actual years.
        EstAvgRec = EstAvgRec + Recvec(Iyr)
        ii = ii+1
1004  CONTINUE 
      EstAvgRec = EstAvgRec/ii
      
      EstB_SPR = EstAvgRec*SBPR_F !Define globally
      
      RETURN
      END SUBROUTINE AFSC_SPR

C
C ===========================================================================
C

      SUBROUTINE FindFofl(Ispec,Iyear)
C This is the subroutine that set up the MinOpt rule for the FUNK function in order to determine the F35%, F40%, B35%, and B40% = one of them is the Fofl depending in the tier system (not the true Fofl)
      USE ShareDataRatioSubs
      IMPLICIT NONE
      
      !Global variables
      INTEGER :: Ispec,Iyear
      
      !Local variables
      INTEGER :: Iconfig,Ivec
      REAL(8) :: Depl
      REAL(8),DIMENSION(100) :: X    
      REAL(8) :: SS
      REAL(8) :: FUNK
      EXTERNAL FUNK
      
!      Ispec = 1  !eventually a loop
        MinOpt = 12 !Find F35% and B35%
      DO 1033 Iconfig = 1,Nconfigs !loop over Iconfig; no Ivec here because EstRec is only config-specific, but would be iVec loop if switch to MCMC; don't include true (iconfig = 0)
        !ADD A LINE HERE TO SKIP TO NEXT ICONFIG IF COR_EXISTS(ICONFIG) = .FALSE.
        X = 0.0
        X(1) = 0.05
        TheM = EstM(Ispec,Iconfig)
        TheFecVec = Fec(1,0:Maxage)
        TheRecVec(1:IYear) = EstRec(Ispec,1,Iconfig,
     +1:IYear) !1 is area 1 (maybe get rid of this dimension later to clean things up)
        TheSelex = EstS(Ispec,1,:,:,Iconfig) 
        ThePlusAge = Amax(Ispec)
        EstB_SPR = 0.0
        EstSPR = 0.0
!Info on what goes into AFSC_SPR:      
!      AFSC_SPR(X,EstM(Ispec,Iconfig),Fec(1,0:Amax(Ispec)),
!     +EstRec(1,1,1:GlobalProjYr),1,1) 
      
        CALL Fit(X,SS,1,CareyGRD)
        SS = FUNK(X)
      
        EstF35(Ispec,Iconfig) = X(1)
        EstB35(Ispec,Iconfig) = EstB_SPR
        
        !Save info to check out whether optimization worked (written to file in PstarSummary) 
        EstSPR35(Ispec,Iconfig) = EstSPR
        SSEstSPR35(Ispec,Iconfig) = SS
        
 
        MinOpt = 13 !Find F40% and B40%
        X = 0.0
        X(1) = EstF35(Ispec,Iconfig) !Might as well use a good seed to speed things up
        !Continue to use same defn's of AFSC_SPR() arguments.
        CALL Fit(X,SS,1,CareyGRD)
        SS = FUNK(X)
      
        EstF40(Ispec,Iconfig) = X(1)
        EstB40(Ispec,Iconfig) = EstB_SPR
 
        !Save info to check out whether optimization worked (written to file in PstarSummary) 
        EstSPR40(Ispec,Iconfig) = EstSPR
        SSEstSPR40(Ispec,Iconfig) = SS
     
        !OFL control rule equation, based on above somersaults:
        !Need an estimate of current SSB for each of NMVN numbers-at-age vectors generated
        !Need an estimate SSB0 that is configuration-specific (with MCMC it would be vector-specific)
      
        !THE AFSC TIER 3 CONTROL RULE:
        Depl = 0.0
      
        DO 200 Ivec = 1,Nmvn
          Depl = SSB_MVNvec(Ispec,Ivec,Iconfig)/S0CAEst(Ispec,1,Iconfig)

          IF ((Depl.GT.0.40)) THEN
            EstFOFL(Ispec,Iconfig,Ivec) = EstF35(Ispec,Iconfig)

          ELSEIF (Depl.LE.0.40.AND.Depl.GT.AFSC_alpha) THEN
            EstFOFL(Ispec,Iconfig,Ivec) = EstF35(Ispec,Iconfig)*
     +      ((SSB_MVNvec(Ispec,Ivec,Iconfig)/EstB40(Ispec,Iconfig))
     +      -AFSC_alpha)*(1-AFSC_alpha)

          ELSE
            EstFOFL(Ispec,Iconfig,Ivec) = 0.0
          ENDIF
        
        !Calculate estimated OFL for each random vector based on FOFL
          !This info is not really needed, but keep code just in case??
          !Need to calculate FOFL and OFL based on the single base assessment for comparison to AFSC control rule.
          !This here is calculated because need a distribution of FOFL for the p* control rule.
          !call subroutine
          !CALL CalcCatchFromF(Ispec,EstFOFL(Ispec,Iconfig,Ivec),M(Ispec),
!     +    Nvec(Ispec,:,Ivec,0:Amax(Ispec)),
!     +    EstS(Ispec,1,:,:))!subroutine doesn't work because EstFOFL is a different type of argument than needed (??)
          !Est_OFL(Ispec,Ivec,Iconfig) = pCatBio(Ispec) !Still need to define Est_OFL (maybe use a different name if this one gets used for something else).
200     CONTINUE !End loop over Ivec (Ivec loop will need to be expanded to include calls to FIT and FUNK above if switch to MCMC method.
1033  CONTINUE !End loop over Iconfig  
       ! TheNs(Isex,Iage) = exp(Nvec(Ispec,Isex,Ivec,Iage))
      

      
      
      !----------------------------------------------------------------------------------------------      
      !Find the true FOFL for each species and year - can always move this to another subroutine if it's too much here.
      !No Iconfigs and no vecs for the true FOFL.
      MinOpt = 14 !Find F35% and B35%
      X = 0.0
      X(1) = EstF35(Ispec,1) !Might as well use a good seed to speed things up
      TheRecVec = 0.0 !Re-initialize after the steps above.
      TheM = M(Ispec)
      TheFecVec = Fec(1,0:Maxage)
      TheRecVec(1:IYear) = TotalRecruits(Ispec,1:IYear)
      TheSelex(1,:) = Sel(Ispec,:)
      TheSelex(2,:) = Sel(Ispec,:)
      
      EstB_SPR = 0.0
      EstSPR = 0.0
      
        CALL Fit(X,SS,1,CareyGRD)
        SS = FUNK(X)
      
      TrueF35(Ispec) = X(1)
      TrueB35(Ispec,IYear) = EstB_SPR  
      
      TrueSPR35(Ispec) = EstSPR
      SSTrueSPR35(Ispec) = SS
      
       !----------------------------------------------------------------------------------------------      
      MinOpt = 15 !Find F40% and B40%
      X = 0.0
      X(1) = TrueF35(Ispec) !Might as well use a good seed to speed things up
      !Continue to use other args to AFSC_SPR() from MinOpt = 15.
      CALL Fit(X,SS,1,CareyGRD)
      SS = FUNK(X)
      
      TrueF40(Ispec) = X(1)
      TrueB40(Ispec,IYear) = EstB_SPR   
      
      TrueSPR40(Ispec) = EstSPR
      SSTrueSPR40(Ispec) = SS
      
      !TrueFOFL calculations:
        Depl = SpBioTot(Ispec,IYear)/Kvec(Ispec)

        IF ((Depl.GT.0.40)) THEN
          TrueFOFL(Ispec,IYear) = TrueF35(Ispec)

        ELSEIF (Depl.LE.0.40.AND.Depl.GT.AFSC_alpha) THEN
          TrueFOFL(Ispec,IYear) = TrueF35(Ispec)*
     +    ((SpBioTot(Ispec,IYear)/TrueB40(Ispec,IYear))
     +    -AFSC_alpha)*(1-AFSC_alpha)

        ELSE
          TrueFOFL(Ispec,IYear) = 0.0
        ENDIF      
      
      !Calculate true OFL
      !Don't need this right now, but might need it in the future
        !use CalcCatchFromF() for this calculation - straightforward (except for inability to enter in certain arguments for some reason).
      
      MinOpt = 16 !Calculate the estimated FOFL and OFL based on the base assessment only using the correct M and h (?)
      !Fill this in later - not sure what exactly makes sense for comparison, as this rule would have more info on M and h inherently than the p* rule (which accounts for uncertainty in M and h).
      !Could be that estimating M and h and doing the AFSC OFL control rule would be the comparison?

      RETURN
      END SUBROUTINE FindFofl

C
C ===========================================================================
C

      SUBROUTINE FindPstarQuota(IYear,Ispec) !IYear = Year2 within DoCab2 (can change this to something more clear)    
C  This subroutine finds the ABC that results in a probability of p* of overfishing
C  The method is to use numerical integration to calculate Pr(estF(t)>estOFL) = p
C  and uses the FIT and FUNK functions to find the ABC that results in p = p* (p* is specified in the CR.dat input file)
C  There is a nested optimization using FIT and FUNK to find the F(t) resulting in a particular trial ABC value
C  for each of Nmvn randomly drawn vectors of Numbers-at-age (drawn in the IterateNvecsMVN subroutine; F(t) found in FindFforCatchReg)
C  
      USE ShareDataRatioSubs
      IMPLICIT NONE
       !Global variables
       INTEGER :: IYear,Ispec

       !Local variables
       REAL(8) :: SSS,CareyGRD2
       REAL(8),DIMENSION(100) :: Y
       REAL(8) :: FUNK2
       EXTERNAL :: FUNK2
       
       !configuration calculations will happen within subroutines for finding F(t) and FOFL. 
       !Used in FUNK and FUNK does not know species (Ispec), so need to change where this is defined and make not dependent on Ispec
!       DO 114 Ispec = 1,Nspec
          SSS = 1.0E+20
          SeedQuota = CatBio(Ispec,1,IYear) !not config specific, of course.
          GlobalIspec = Ispec
          Funk2_MinOpt = 2
          CareyGRD2 = 1.2
          
          Y(1) = LOG(SeedQuota)
           CALL Bisect(Y,SSS,0.01,3.0,1000)
!          CALL FIT2(Y,SSS,1,CareyGRD2) !This used to be used (AMOEBA optimization)
          SSS = FUNK2(Y) !Y might not be the best variable either, as it appears within other functions, but probably ok in that it's a local variable.
          
!          !Below is cut and pasted into FUNK2 with global year and global ispec.
!          !Find F(t) that would produce SeedQuota
!          CALL FindFforCatchReg(IYear,Ispec)
!   
!          !Find FOFL for all of the random vectors
!          CALL FindFofl(Ispec,IYear) !Enter Ispec as an input! CRM 7/28/14
!       
!          !Integrate over (F(t),FOFL) distribution
!          !Want the proportion of time F(t) > FOFL among the pairs (F(t),FOFL) sampled by way of the N-at-age vectors.
          
          
          !FIT function (or new fit function) for achieving p(overfishing) = p*
          !by way of changing SeedQuota
       
      
          !QuotaAct(Ispec) = X-type results from new FIT function
          QuotaAct(Ispec) = EXP(Y(1))  
          EstimatedP = SSS
!114    CONTINUE
       
       RETURN !What does RETURN get you anyway?
      END SUBROUTINE FindPstarQuota
      
C
C ===========================================================================
C
    
!!===================================================================
!      SUBROUTINE FindMSY(IYear,Ispec) !IYear = Year2 within DoCab2 (can change this to something more clear)
!      USE ShareDataRatioSubs
!      IMPLICIT NONE
!      
!      !Just run this at teh end or after find Qvals as long as that's still there.
!      !Global variables
!      INTEGER :: IYear,Ispec
!      
!      !Local variables
!      
!      !Find the F that, when applied, results in the greatest catch over 100 years
!      
!      !X(1) = F
!      !FUNK = 1E+20
!      
!      !Within funk, loop over years and call ProjYr and some combo thing. Double check how it worked to FUNK for a particular depletion level
!      
!      !Within funk: FUNK = FUNK + -(CatBio(Iyear,Ispec,...)**2) looped over IYear    
!     END SUBROUTINE FindMSY
      
C
C ===========================================================================
C

      SUBROUTINE SummaryPstar(IYear) !IYear = Year2 within DoCab2 (can change this to something more clear)
C This produces the summary output from the Pstar control rule along with SPR35%, SPR40% for example. 	
C
      USE ShareDataRatioSubs
      IMPLICIT NONE
  
      !Global variables
      INTEGER :: IYear,Ispec,Iconfig,Ivec
      
      !----------------------------------------------------------------
      !Pstar1.OUT (UNIT = 400):
      IF (IYear.EQ.(FinDepYr+1).AND.Isim.EQ.1) THEN
      WRITE(400,*) 'Isim IYearMinus1appliedinIYear Ispec Iconfig Ivec 
     + FDistn SSFDistn
     + EstFOFL TrueF SSTrueF TrueFOFL EstimatedP Pstar QuotaAct'
      WRITE(400,*) Nsim,NProjYr,Nspec,Nconfigs,Nmvn
      ENDIF
      
      DO 100 Ispec = 1,Nspec
          DO 100 Iconfig = 1,Nconfigs
              DO 100 Ivec = 1,Nmvn
      WRITE(400,1) Isim,IYear-1,Ispec,Iconfig,Ivec,
     +             FDistn(Ispec,IYear-1,Iconfig,Ivec),
     +             SSFDistn(Ispec,IYear-1,Iconfig,Ivec),
     +             EstFOFL(ISpec,Iconfig,Ivec),
     +             TrueF(Ispec,IYear-1),SSTrueF(Ispec,IYear-1),
     +             TrueFOFL(Ispec,IYear-1),
     +             EstimatedP,Pstar,QuotaAct(Ispec)
100   CONTINUE
      !----------------------------------------------------------------
!         TrueF(Ispec,IYear) = X(1)
!         SSTrueF(Ispec,IYear) = SS !As a check.      

      !----------------------------------------------------------------
      !Pstar2.OUT (UNIT 401): Write out optimization info for estimated F35% and F40% (all configuration-specific)
      IF (IYear.EQ.(FinDepYr+1).AND.Isim.EQ.1) THEN
      WRITE(401,*) 'Isim,IYr,Ispec,Icon,EstF35,EstB35,EstSPR35,
     +SSEstSPR35,EstF40,EstB40,EstSPR40,SSEstSPR40,TruSPR35,
     +SSTruSPR35,TruSPR40,SSTruSPR40'
      WRITE(401,*) Nsim,NProjYr,Nspec,Nconfigs
      ENDIF
      
      DO 101 Ispec = 1,Nspec
          DO 101 Iconfig = 1,Nconfigs
      WRITE(401,2) Isim,IYear-1,Ispec,Iconfig,EstF35(Ispec,Iconfig),
     +EstB35(Ispec,Iconfig),EstSPR35(Ispec,Iconfig),
     +SSEstSPR35(Ispec,Iconfig),EstF40(Ispec,Iconfig),
     +EstB40(Ispec,Iconfig),EstSPR40(Ispec,Iconfig),
     +SSEstSPR40(Ispec,Iconfig),TrueSPR35(Ispec),SSTrueSPR35(Ispec),
     +TrueSPR40(Ispec),SSTrueSPR40(Ispec)

101   CONTINUE

      
      
      !----------------------------------------------------------------
      
1     FORMAT(5(1x,I6),1x,F17.10,1x,F20.5,1x,F17.10,1x,F17.10,1x,F20.5,
     +1x,F17.10,2(1x,F15.12),1x,F15.5)
2     FORMAT(4(1x,I4),8(1x,F10.8,1x,F10.4))
!3     FORMAT(20(25x,1x))

      END SUBROUTINE SummaryPstar

C
C ===========================================================================
C

      SUBROUTINE SummaryPstar_Old(IYear) !IYear = Year2 within DoCab2 (can change this to something more clear)
C This produces another type of summary output from the Pstar control rule in additio to the SummaryPstar	
C
      USE ShareDataRatioSubs
      IMPLICIT NONE
  
      !Global variables
      INTEGER :: IYear,Ispec,Iconfig,Ivec
      
      !----------------------------------------------------------------
      !Pstar1.OUT (UNIT = 400):
      IF (IYear.EQ.(FinDepYr+1)) THEN
      WRITE(4000,*) 'Isim IYearMinus1appliedinIYear Ispec Iconfig Ivec 
     + FDistn SSFDistn
     + EstFOFL TrueF SSTrueF TrueFOFL EstimatedP Pstar QuotaAct'
      WRITE(4000,*) Nsim,NProjYr,Nspec,Nconfigs,Nmvn
      ENDIF
      
      DO 100 Ispec = 1,Nspec
          DO 100 Iconfig = 1,Nconfigs
              DO 100 Ivec = 1,Nmvn
      WRITE(4000,1) Isim,IYear-1,Ispec,Iconfig,Ivec,
     +             FDistn(Ispec,IYear-1,Iconfig,Ivec),
     +             SSFDistn(Ispec,IYear-1,Iconfig,Ivec),
     +             EstFOFL(ISpec,Iconfig,Ivec),
     +             TrueF(Ispec,IYear-1),SSTrueF(Ispec,IYear-1),
     +             TrueFOFL(Ispec,IYear-1)
100   CONTINUE
      !----------------------------------------------------------------
!         TrueF(Ispec,IYear) = X(1)
!         SSTrueF(Ispec,IYear) = SS !As a check.      

      !----------------------------------------------------------------
      !Pstar2.OUT (UNIT 401): Write out optimization info for estimated F35% and F40% (all configuration-specific)
      IF (IYear.EQ.(FinDepYr+1)) THEN
      WRITE(4001,*) 'Isim,IYr,Ispec,Icon,EstF35,EstB35,EstSPR35,
     +SSEstSPR35,EstF40,EstB40,EstSPR40,SSEstSPR40,TruSPR35,
     +SSTruSPR35,TruSPR40,SSTruSPR40'
      WRITE(4001,*) Nsim,NProjYr,Nspec,Nconfigs
      ENDIF
      
      DO 101 Ispec = 1,Nspec
          DO 101 Iconfig = 1,Nconfigs
      WRITE(4001,2) Isim,IYear-1,Ispec,Iconfig,EstF35(Ispec,Iconfig),
     +EstB35(Ispec,Iconfig),EstSPR35(Ispec,Iconfig),
     +SSEstSPR35(Ispec,Iconfig),EstF40(Ispec,Iconfig),
     +EstB40(Ispec,Iconfig),EstSPR40(Ispec,Iconfig),
     +SSEstSPR40(Ispec,Iconfig),TrueSPR35(Ispec),SSTrueSPR35(Ispec),
     +TrueSPR40(Ispec),SSTrueSPR40(Ispec)

101   CONTINUE

      
      
      !----------------------------------------------------------------
      
1     FORMAT(5(1x,I6),1x,F17.10,1x,F20.5,1x,F17.10,1x,F17.10,1x,F20.5,
     +1x,F17.10)
2     FORMAT(4(1x,I4),8(1x,F10.8,1x,F10.4))
!3     FORMAT(20(25x,1x))

      END SUBROUTINE SummaryPstar_Old
!======================================================================

!!======================================================================
 !     SUBROUTINE FindFMSY(Ispec) !IYear = Year2 within DoCab2 (can change this to something more clear)
!      USE ShareDataRatioSubs
!      IMPLICIT NONE
!
!      !Global variables
!      INTEGER :: Ispec
!      REAL(8),DIMENSION(100) :: X
!      
!      !First find the effort that produces highest long-term catch, then find the F that produces that catch using above functions
!      MinOpt = 16 !is that right?
!      Qval(Ispec) = 1.0
!      X(1) = LOG(TotEffort(StartFishing))
!
!      
!      !CALL FIT(X,... 
!      END SUBROUTINE FindFMSY
!======================================================================

C
C ===========================================================================
C

      SUBROUTINE Bisect(Yinit,SSS,Tol,GridSize,NIter)
C Uses the bisection method to find the P* value (take two points, take the mid point as the next trial, and see if the objective function is improved and keep taking the midpoint until the convergence criteria is achieved	
C
      USE ShareDataRatioSubs
      IMPLICIT NONE      
      
      !works for monotonic relatioships.
      !Global variables
      REAL(8),DIMENSION(100) :: Yinit
      REAL(8) :: FUNK2
      REAL(8) :: SSS
      REAL :: Tol,GridSize
      INTEGER :: NIter
      
      EXTERNAL FUNK2
      
      !Local variables
      REAL(8),DIMENSION(100):: Ysmall,Ybig
      REAL(8) :: SSSsmall, SSSbig
      INTEGER :: i
      REAL(8) :: SSSfinal
 

      !Calculate smaller and larger catch guesses based on GridSize
 !     Ysmall(1) = (1/GridSize)*Yinit(1)
      Ysmall(1) = 0.0
      Ybig(1) = GridSize*Yinit(1)
 
      !Loop over Niter times or until tolerance is reached:
      i = 0
      DO WHILE ((i.LE.Niter).OR.(SSSfinal.LE.Tol))
         SSSsmall = 1.0E+20
         SSSbig = 1.0E+20
         SSS = 1.0E+20
         i = i+1
        !Calculate objective function values for difference between p and p* under each catch guess (small, init, large)
        SSS = FUNK2(Yinit)
 !       SSSsmall = FUNK2(Ysmall)
 !       SSSbig = FUNK2(Ybig)   

        IF (SSS.LT.Pstar) THEN
            !search the upper half of the range: find the mid-point of upper half of range
            Ysmall(1) = Yinit(1)
            Yinit(1) = Ysmall(1) + (Ybig(1)-Ysmall(1))/2   
        ELSEIF (SSS.GT.Pstar) THEN
            !search the lower half of the range: find the mid-point of the lower half of range
            Ybig(1) = Yinit(1)
            Yinit(1) = Ybig(1) - (Ybig(1)-Ysmall(1))/2   
        ELSE
            !you already found p = p*; return SSS and Yinit
                
           RETURN  
       ENDIF
        
      END DO
      
      
      
      END SUBROUTINE Bisect

C
C ===========================================================================
C

C NOTES on the tech int model 
C
C Modifications:
C 1. line 1745: changed Aplus to Amax because it doesn't make sense if not
C 2. line 2492: included the statement about process error = recruit deviation
C 3. subroutine getYield is WRONG (the calculation of recruitment)
C 4. line 4775: refined the increment value for E
C 5. CRM: line 1745: changed modification 1 back to be Aplus (which has a different meaning depending on whether there's a plus group or not) - Kotaro and Carey discussed.
C 6. CRM: Uncommented lines 3884 and 3885 to correct inconsistency between Carey and Kotaro's code and make sure old simulation files get cleaned up properly
C 7. CRM: added back in lines 679-683 to delete old files (again to avoid mix-ups between old and new simulation results)
C 8. CRM: 
C
C TO DO:
C 0. Be sure to use AFSC_SPR in the MPA_model (the main program) instead of the P* controle rule for the multispecies MSE
C 1. change the fishery age composition data that goes into CAB to be sampled from a overdispersed distribution (right now it is multinomial)  
C 2. Change the RECDEVS simulation. We need the possibility to include user-specified rec-devs or a multivariate time series model
C 3. change the MPAmodel program to choose the effort level for each species that leads to the initial or final depletion level      
C 4. check that the model is adapted for N species
C 5. check the run and see if the model makes sense (the three species should have the same dynamics)
C 6. make the model more realistic
C 7. include the tech interaction into the model (CALL Jim Ianelli's ADMB .exe. Need to modify this and compile it first)
C 	7.1. RUN the fishing strategy to the observer data
C 	7.2. use the above information and build the capacity to have a time varying fishing strategy 
C 	7.3. include the harvest cap into the model (that is about 1.7 millions based on 2 million - the three species TAC)
C 8. include the control rule (how to calculate TAC) into the model
C 9. include an option to generate recdevs based on multivariate time series model

C MERGE to cbisection  