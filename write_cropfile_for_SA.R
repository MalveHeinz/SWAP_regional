 ### write crop input file for the combined SA


  #-------------------------  05 : write crop text file --------------------------#
# SWAP expects floating numbers (=1.0, 0.0) not integers (1, 0), which is why I use the function format() here

  
  tx=c("***********************************************************************************************",
      "* Filename: PotatoD.CRP",
      "* Contents: SWAP 4 - input Data for detailed crop model for LHS",
      "***********************************************************************************************",
      "** Input data based on:",
      "* Potato (Solanum tuberosum L.)",
      "* Input data set were based on :",
      "** $Id: pot701.cab 1.3 1997/09/25 14:07:00 LEM release $",
      "** File POT701.CAB",
      "** CROP DATA FILE for use with WOFOST Version 5.4, June 1992",
      "**",
      "** default potato file",
      "**********************************************************************************",
      "",
      "***********************************************************************************************",
      "",
      "*** PLANT GROWTH SECTION ***",
      "",
      "***********************************************************************************************",
      "* Part 1: Crop factor or crop height                             ",
      "",
      "* Choose between crop factor and crop height",
      "* Choose crop factor if ETref is used, either from meteo input file (SWETR = 1) or with Penman-Monteithv",
      "* Choose crop height if Penman-Monteith should be used with actual crop height, albedo and canopy resistance",
      "SWCF = 1 ! 1 = crop factor ",
      "! 2 = crop height",
      "",
      "* If SWCF = 1, list crop factor CF [0..2 -, R],   as function of dev. stage DVS [0..2 -, R]:",
      "* If SWCF = 2, list crop height CH [0..1.d4 cm, R], as function of dev. stage DVS [0..2 -, R]:",
      "",
      "DVS       CH     CF   ! ( maximum MAGRS records) ",
      paste0("0.0      1.0    ", format(lhs_paras$CF[i], nsmall=1, digits=1),""),
      paste0("1.0      40.0    ",format(lhs_paras$CF[i]*1.1, nsmall=1, digits=1),""),
      paste0("2.0      50.0    ",format(lhs_paras$CF[i]*1.1, nsmall=1, digits=1),""),
      "* End of table",
      "",
      "* If SWCF = 2, in addition to crop height list crop specific values for:",
      "ALBEDO =   0.23 ! crop reflection coefficient [0..1.0 -, R]  ",                  
      paste0("RSC    =  ",format(lhs_paras$RSC[i], nsmall=1, digits=1)," ! Minimum canopy resistance [0..1d6 s/m, R]  "),                  
      "RSW    =    0.0 ! Canopy resistance of intercepted water [0..1d6 s/m, R]  ",     
      "***********************************************************************************************",
      "",
      "***********************************************************************************************",
      "* Part 2 : Crop development",
      "",
      "* Switch for crop development:",
      "IDSL   = 0 ! 0 = Crop development before anthesis depends on temperature only",
      "! 1 = Crop development before anthesis depends on daylength only",
      "! 2 = Crop development before anthesis depends on temperature and daylength",
      "",
      "* If IDSL = 1 or 2, specify:",
      "DLO    = 14.0     ! Optimum day length for crop development [0..24 h, R]",
      "DLC    =  8.0     ! Minimum day length [0..24 h, R]",
      "",
      "* If IDSL = 0 or 2 specify:",
      "TSUMEA =   150.0 ! Temperature sum from emergence to anthesis [0..1d4 ºC, R]",
      "TSUMAM =  1550.0 ! Temperature sum from anthesis to maturity  [0..1d4 ºC, R]",
      "",
      "* List increase in temperature sum [0..60 ºC, R] as function of daily average temperature [0..100 ºC, R]",
      "",
      "*         TAV  DTSM    (maximum 15 records)",
      "DTSMTB =",
      "0.00   0.00",
      "2.00   0.00",
      "13.00  11.00",
      "30.00  28.00",
      "* End of Table",
      "",
      "DVSEND =      2.00 ! Development stage at harvest [0..3 -, R]",
      "*",
      "* Germination defined with INITCRP in .swp-file : ",
      "*  INITCRP=2: CROPSTART defines emergence (default), INITCRP=2: CROPSTART defines sowing",
      "* IF INITCRP = 2 specify",
      "TSUMEMEOPT  =  170.   ! temperature sum needed for crop emergence     [0..1000 C d, R]",
      "TBASEM      =  3.0    ! minimum temperature, used for germination trajectory  [0..40 C, R]  ",
      "TEFFMX      =  18.0   ! maximum temperature, used for germination trajectory  [0..40 C, R]  ",
      "HDRYGERM    =  -500.0 ! pressure head rootzone for dry germination trajectory [-1000..-0.01 cm, R]",
      "HWETGERM    =  -100.0 ! pressure head rootzone for wet germination trajectory [-100..-0.01 cm, R]",
      "AGERM       =  203.   ! a-coefficient Eq. 24/25 Feddes & Van Wijk     [1..1000, R]",
      "CGERM       = -432.   ! c-coefficient Eq. 24    Feddes & Van Wijk     [1..1000, R]",
      "BGERM       =  522.   ! b-coefficient Eq. 25    Feddes & Van Wijk     [1..1000, R]  ",
      "***********************************************************************************************",
      "",
      "***********************************************************************************************",
      "* Part 3: Initial values",
      "",
      paste0("TDWI   =    ",format(lhs_paras$TDWI[i], nsmall=1)," ! Initial total crop dry weight [0..10000 kg/ha, R]"),
      "LAIEM  =  0.0589 ! Leaf area index at emergence [0..10 m2/m2, R]",
      paste0("RGRLAI = ",format(lhs_paras$RGRLAI[i], nsmall=4)," ! Maximum relative increase in LAI [0..1 m2/m2/d, R]"),
      "***********************************************************************************************",
      "",
      "***********************************************************************************************",
      "* Part 4: Green surface area",
      "",
      "SPA    =  0.0000 ! Specific pod area  [0..1 ha/kg, R]",
      "SSA    =  0.0000 ! Specific stem area [0..1 ha/kg, R]",
      paste0("SPAN   =   ",format(lhs_paras$SPAN[i], nsmall=1)," ! Life span under leaves under optimum conditions [0..366 d, R]"),
      "TBASE  =    2.00 ! Lower threshold temperature for ageing of leaves [-10..30 ºC, R]",
      "",
      "* List specific leaf area [0..1 ha/kg, R] as function of crop development stage [0..2 -, R]",
      "",
      "*           DVS    SLA    (maximum 15 records)",
      "SLATB =",
      paste0("0.00 ",format(lhs_paras$SLATB[i], nsmall=4),""),
      paste0("1.10 ",format(lhs_paras$SLATB[i], nsmall=4),""),
      paste0("2.00 ",format(lhs_paras$SLATB[i]*0.5, nsmall=4),""),
      "* End of Table ",
      "***********************************************************************************************",
      "",
      "***********************************************************************************************",
      "* Part 5: Assimilation",
      "",
      "KDIF   =    1.00 ! Extinction coefficient for diffuse visible light [0..2 -, R]",
      "KDIR   =    0.75 ! Extinction coefficient for direct visible light [0..2 -, R]",
      "EFF    =    0.45 ! Light use efficiency [0..10 kg/ha/hr/(Jm2s), R]",
      "",
      "* List max CO2 assimilation rate [0..100 kg/ha/hr, R] as function of development stage [0..2 -, R]",
      "",
      "*          DVS    AMAX   (maximum 15 records)",
      "AMAXTB =",
      paste0("0.00 ",format(lhs_paras$AMAXTB[i],nsmall=2),""),
      paste0("1.57 ",format(lhs_paras$AMAXTB[i],nsmall=2),""),
      "2.00   0.00 ",
      "* End of table ",
      "",
      "* List reduction factor of AMAX [-, R] as function of average day temperature [-10..50 ºC, R]",
      "",
      "*          TAVD   TMPF  (maximum 15 records)",
      "TMPFTB =",
      "0.00  0.010",
      "3.00  0.010",
      paste0("10.00  ",format(lhs_paras$TMPFTB[i], nsmall=3),""),
      "15.00  1.000",
      "24.00  1.000",
      paste0("29.00  ",format(lhs_paras$TMPFTB[i], nsmall=3),""),
      "36.00  0.010",
      "* End of table ",
      "",
      "* List reduction factor of AMAX [-, R] as function of minimum day temp. [-10..50 ºC, R]",
      "",
      "*          TMNR    TMNF  (maximum 15 records)",
      "TMNFTB = ",
      "0.00  0.000",
      "3.00  1.000",
      "* End of table ",
      "***********************************************************************************************",
      "",
      "***********************************************************************************************",
      "* Part 6: Conversion of assimilates into biomass",
      "",
      "CVL    =  0.7200 ! Efficiency of conversion into leaves         [0..1 kg/kg, R]",
      "CVO    =  0.8500 ! Efficiency of conversion into storage organs [0..1 kg/kg, R]",
      "CVR    =  0.7200 ! Efficiency of conversion into roots          [0..1 kg/kg, R]",
      "CVS    =  0.6900 ! Efficiency of conversion into stems          [0..1 kg/kg, R]",
      "***********************************************************************************************",
      "",
      "***********************************************************************************************",
      "* Part 7: Maintenance respiration",
      "",
      paste0("Q10    =  ",format(lhs_paras$Q10[i], nsmall=1, digits=2)," ! Increase in respiration rate with temperature  [0..5 -/10 C, R]"),
      "RML    =  0.0300 ! Maintenance respiration rate of leaves         [0..1 kgCH2O/kg/d, R]",
      "RMO    =  0.0045 ! Maintenance respiration rate of storage organs [0..1 kgCH2O/kg/d, R]",
      "RMR    =  0.0100 ! Maintenance respiration rate of roots          [0..1 kgCH2O/kg/d, R]",
      "RMS    =  0.0150 ! Maintenance respiration rate of stems          [0..1 kgCH2O/kg/d, R]",
      "",
      "* List reduction factor of senescence [-, R] as function of development stage [0..3 -, R]",
      "",
      "*          DVS    RFSE  (maximum 15 records)",
      "RFSETB = ",
      "0.00   1.00",
      "2.00   1.00",
      "* End of table ",
      "***********************************************************************************************",
      "",
      "***********************************************************************************************",
      "* Part 8: Partitioning",
      "",
      "* List fraction of total dry matter increase partitioned to the roots [kg/kg, R]",
      "* as function of development stage [0..3 -, R]",
      "",
      "*          DVS     FR    (maximum 15 records)",
      "FRTB = ",
      "0.00   0.2",
      "1.00   0.2",
      "1.36   0.00",
      "2.00   0.00",
      "* End of table ",
      "",
      "* List fraction of total above ground dry matter increase partitioned to the leaves [kg/kg, R]",
      "* as function of development stage [0..3 -, R]",
      "",
      "*          DVS     FL   (maximum 15 records)",
      "FLTB = ",
      " 0.00  0.83", 
      " 1.00  0.586617",
      " 1.30  0.0100666",
      " 1.40  0.00",
      " 2.00  0.00",
      "* End of table ",
      "",
      "* List fraction of total above ground dry matter increase partitioned to the stems [kg/kg, R]",
      "* as function of development stage [0..3 -, R]",
      "",
      "*          DVS    FS   (maximum 15 records)",
      "FSTB = ",
      " 0.00  0.17", 
      " 1.00  0.00",
      " 1.30  0.1664687",
      " 1.40  0.102477",
      " 2.00  0.00",
      "* End of table ",
      "",
      "* List fraction of total above ground dry matter increase partitioned to the storage organs [kg/kg, R]",
      "* as function of development stage [0..3 -, R]",
      "",
      "*          DVS    FO    (maximum 15 records)",
      "FOTB = ",
      " 0.00  0.00", 
      " 1.00  0.4133842",
      " 1.30  0.8234647",
      " 1.40  0.897523",
      " 2.00  1.00",
      "* End of table",
      "***********************************************************************************************",
      "",
      "***********************************************************************************************",
      "* Part 9: Death rates",
      "" ,
      "PERDL =   0.030 ! Maximum relative death rate of leaves due to water stress [0..3 /d, R]",
      "",
      "* List relative death rates of roots [kg/kg/d] as function of development stage [0..2 -, R]",
      "",
      "*          DVS    RDRR    (maximum 15 records)",
      "RDRRTB = ",
      "0.0000 0.0000",
      "1.5000 0.0000",
      "1.5001 0.0200",
      "2.0000 0.0200",
      "* End of table",
      "",
      "* List relative death rates of stems [kg/kg/d] as function of development stage [0..2 -, R]",
      "",
      "*          DVS     RDRS    (maximum 15 records)",
      "RDRSTB = ",
      "0.0000 0.0000",
      "1.5000 0.0000",
      "1.5001 0.0200",
      "2.0000 0.0200",
      "* End of table",
      "***********************************************************************************************",
      "",
      "***********************************************************************************************",
      "* Part 10: Soil water extraction by plant roots  ",                                          
      "",
      "* -- Part 10a: Oxygen stress -----------------------",
      "",
      "* Switch for oxygen stress:",
      "SwOxygen = 2      ! 0 = No oxygen stress",
      "! 1 = Oxygen stress according to Feddes et al. (1978)",
      "! 2 = Oxygen stress according to Bartholomeus et al. (2008)",
      "",
      "* If SwOxygen = 1, specify:",
      "HLIM1  =    -10.0    ! No water extraction at higher pressure heads [-100..100 cm, R]",
      "HLIM2U =    -25   ! h below which optimum water extr. starts for top layer [-1000..100 cm, R]",
      "HLIM2L =    -25    ! h below which optimum water extr. starts for sub layer [-1000..100 cm, R]",
      "",
      "* If SwOxygen = 2, specify:",
      "SwOxygenType = 1      ! Switch for physical processes or reproduction functions to calculate oxygen stress:",
      " ! 1 = Use physical processes",
      "! 2 = Use reproduction functions",
      "",
      "* In case of physical processes (SwOxygenType = 2), specify:",
      paste0("Q10_microbial       = ",format(lhs_paras$Q10_microbial[i],nsmall=2),"          ! Relative increase in microbial respiration at temperature rise of 10 ºC [1.0..4.0 -, R]"),
      paste0("Specific_resp_humus = ",format(lhs_paras$resp_humus[i],nsmall=2),"         ! Respiration rate of humus at 25 ºC [0.0..1.0 kg O2/kg ºC/d, R] "),
      paste0("SRL                 = ",format(lhs_paras$SRL[i],nsmall=1),"      ! Specific root length [0.d0..1d10 (m root)/(kg root), R]    "),  
      "SwRootRadius        = 2              ! Switch for calculation of root radius:",
      "! 1 = Calculate root radius",
      "! 2 = Root radius is given in input file",
      "* If SwRootRadius = 1, specify:",
      " Dry_mat_cont_roots      = 0.075d0    ! Dry matter content of roots [0..1 -, R]",
      "Air_filled_root_por     = 0.05d0     ! Air filled root porosity [0..1 -, R]",
      "Spec_weight_root_tissue = 1.0d3      ! Specific weight of non-airfilled root tissue [0..1d5 (kg root)/(m3 root), R]",
      "Var_a                   = 4.175d-10  ! Variance of root radius [0..1 -, R]",
      "",
      "*If SwRootRadius = 2, specify:",
      paste0("Root_radiusO2 = ",format(lhs_paras$root_radius[i],nsmall=5),"            ! Root radius (mind: in meter!) for oxygen stress module [1d-6..0.1 m, R]"),
      "",
      "* In case of reproduction functions (SwOxygenType = 2), specify:",
      "SwTopSub     = 2      ! Switch for topsoil or subsoil: 1 = topsoil, 2 = subsoil",
      "NrStaring    = 3      ! Number of soil type according to Staring series (Wosten et al., 2001), [1..18, I]",
      "",
      "* -- Part 10b: Drought stress -----------------------",
      "",
      "*Switch for drought stress:",
      "SwDrought = 1      ! 1 = Drought stress according to Feddes et al. (1978)",
      "! 2 = Drought stress according to De Jong van Lier et al. (2008)",
      "",
      "* If SwDrought = 1, or in case of irrigation scheduling (SCHEDULE = 1), specify:",
      paste0("HLIM3H =    ",format(lhs_paras$HLIM3H[i],nsmall=1),"     ! Pressure head below which water uptake reduction starts at high Tpot [-1d4..100 cm, R]"),
      paste0("HLIM3L =    ",format(lhs_paras$HLIM3L[i],nsmall=1),"    ! Pressure head below which water uptake reduction starts at low Tpot [-1d4..100 cm, R]"),
      paste0("HLIM4  =  ",format(lhs_paras$HLIM4[i],nsmall=1),"     ! No water extraction at lower soil water pressure heads [-2d4..100 cm, R]"),
      paste0("ADCRH  =       ",format(lhs_paras$ADCRH[i],nsmall=2),"     ! Level of high atmospheric demand, corresponding to HLIM3H [0..5 cm/d, R] "),    
      paste0("ADCRL  =       ",format(lhs_paras$ADCRL[i],nsmall=2),"     ! Level of low atmospheric demand, corresponding to HLIM3L [0..5 cm/d, R] "),    
      paste0("ALPHACRIT =    ",format(lhs_paras$ALPHACRIT[i],nsmall=2),"     ! Critical stress index (Jarvis, 1989) for compensation of root water uptake [0.2..1 -, R]  "),
      "",
      "* If SwDrought = 2, specify:",
      "WILTPOINT  = -20000.0 ! Minimum pressure head in leaves [-1d8..-1d2 cm, R]",
      "KSTEM =       1.03d-4 ! Hydraulic conductance between leaf and root xylem [1d-10..10 /d, R]",
      "RXYLEM =         0.02 ! Xylem radius [1d-4..1 cm, R]",
      "ROOTRADIUS =     0.05 ! Root radius [1d-4..1 cm, R]",
      "KROOT =        3.5d-5 ! Radial hydraulic conductivity of root tissue [1d-10..1d10 cm/d, R] ",
      "ROOTCOEFA  =     0.53 ! Defines relative distance between roots at which mean soil water content occurs [0..1 -, R]",
      "SWHYDRLIFT =        0 ! Switch for possibility hydraulic lift in root system [N=0, Y=1]",
      "ROOTEFF    =      1.0 ! Root system efficiency factor [0..1 -, R]",
      "STEPHR   =        1.0 ! Step between values of hroot and hxylem in iteration cycle [0..10 cm, R]",
      "CRITERHR =      0.001 ! Maximum difference of Hroot between iterations; convergence criterium [0...10 cm, R]",
      "TACCUR =        0.001 ! Maximum absolute difference between simulated and calculated potential transpiration rate (1d-5..1d-2 cm/d, R)",
      "***********************************************************************************************",
      "",
      "***********************************************************************************************",
      "* Part 11: salt stress ",                                           
      "",
      "* Switch salinity stress ",
      "SWSALINITY = 0  ! 0 = No salinity stress",
      " ! 1 = Maas and Hoffman reduction function",
      " ! 2 = Use osmotic head",
      "",
      "* If SWSALINITY = 1, specify threshold and slope of Maas and Hoffman",
      "SALTMAX   =  3.0 ! Threshold salt concentration in soil water  [0..100 mg/cm3, R] ",
      "SALTSLOPE =  0.1 ! Decline of root water uptake above threshold [0..1.0 cm3/mg, R] ",
      "",
      "* If SWSALINITY = 2, specify:",
      "  SALTHEAD  =  624.0 ! Conversion salt concentration (mg/cm3) into osmotic head (cm) [0..1000.0 cm/(mg/cm3), R]  ",
      "***********************************************************************************************",
      "",
      "***********************************************************************************************",
      "* Part 12: interception                                            ",
      "",
      "* For agricultural crops apply interception concept of Von Hoyningen-Hune and Braden",
      "SWINTER =  1  ! Switch for rainfall interception method:",
      "! 0 = No interception calculated",
      "! 1 = Agricultural crops (Von Hoyningen-Hune and Braden)",
      "! 2 = Trees and forests (Gash)",
      "COFAB =  0.25 ! Interception coefficient, corresponding to maximum interception amount [0..1 cm, R]",
      "***********************************************************************************************",
      "",
      "***********************************************************************************************",
      "* Part 13: Root growth and root density profile",
      "",
      "RDI    =   10.00 ! Initial rooting depth [0..1000 cm, R]",
      paste0("RRI    =   ",format(lhs_paras$RRI[i], nsmall=2)," ! Maximum daily increase in rooting depth [0..100 cm/d, R]"),
      paste0("RDC    =   ",format(lhs_paras$RDC[i], nsmall=1)," ! Maximum rooting depth of particular cultivar [0..1000 cm, R]"),
      "",
      "* List root density [0..100 cm/cm3, R] as function of relative rooting depth [0..1 -, R]:",
      "* In case of drought stress according to Feddes et al. (1978) (SWDROUGHT = 1), relative root density (-) is sufficient",
      "",
      "*    Rdepth Rdensity          ! (maximum 11 records)",
      "RDCTB =             ",  
      paste0("0.0,",(lhs_paras$RDCTB[i]*log(1) + 0.9966),""),
      paste0("0.1,",(lhs_paras$RDCTB[i]*log(2) + 0.9966),""),
      paste0("0.2,",(lhs_paras$RDCTB[i]*log(3) + 0.9966),""),
      paste0("0.3,",(lhs_paras$RDCTB[i]*log(4) + 0.9966),""),
      paste0("0.4,",(lhs_paras$RDCTB[i]*log(5) + 0.9966),""),
      paste0("0.5,",(lhs_paras$RDCTB[i]*log(6) + 0.9966),""),
      paste0("0.6,",(lhs_paras$RDCTB[i]*log(7) + 0.9966),""),
      paste0("0.7,",(lhs_paras$RDCTB[i]*log(8) + 0.9966),""),
      paste0("0.8,",(lhs_paras$RDCTB[i]*log(9) + 0.9966),""),
      paste0("0.9,",(lhs_paras$RDCTB[i]*log(10) + 0.9966),""),
      paste0("1.0,",(lhs_paras$RDCTB[i]*log(11) + 0.9966),""),
      "* End of table",
      "",
      "* Swith for calculation rooting depth:",
      " SWDMI2RD = 1  ! 0 = Rooting depth increase is related to availability assimilates for roots",
      "! 1 = Rooting depth increase is related to relative dry matter increase",
      "**********************************************************************************",
      "",
      "**********************************************************************************",
      "* Part 14: Management, other than irrigation, e.g. pests,diseases, or nutrients",
      "",
      "* Management factor",
      "flpotrelmf = .false. ! Flag indicating calculation of attainable yield instead of theoretical potential yield",
      "relmf = 0.90         ! Management factor to reduce theoretical potential yield to attainable yield [0..1.0 -, R]",
      "",
      "* Losses of organic matter ",
      "FraDeceasedLvToSoil = 0.3   ! Fraction of deceased leaves incorporated in soil  [0..1.0 kg/kg DM, R]",
      "************************************************************************************",
      "",
      "",
      "*** IRRIGATION SCHEDULING SECTION ***",
      "",
      "**********************************************************************************",
      "* Part 15: General",
      "",
      paste0("SCHEDULE = ",irr,"  ! Switch for application irrigation scheduling [Y=1, N=0] "),
      "",
      "* If SCHEDULE = 0, no more information is required in this input file!", 
      "* If SCHEDULE = 1, continue ....",
      "",
      "STARTIRR = 8 6 ! Specify day and month at which irrigation scheduling starts [dd mm]",
      "ENDIRR  = 1 8 ! Specify day and month at which irrigation scheduling stops [dd mm]",
      " CIRRS = 100.0     ! Solute concentration of irrigation water [0..100 mg/cm3, R]",
      "ISUAS = 0       ! Switch for type of irrigation method: ",
      "! 0 = sprinkling irrigation",
      "! 1 = surface irrigation",
      "",
      "* Specify pressure head at field capacity which will be used for irrigation timing options",
      "phFieldCapacity = -100.0   ! Soil water pressure head at field capacity [-1000..0 cm, R] ",
      "**********************************************************************************",
      "",
      "**********************************************************************************",
      "* Part 16: Irrigation time criteria",
      "",
      "*** Choose one of the following 5 timing options:",
      "TCS = 3  ! Switch for timing criterion [1..6 -, I]",
      "! 1 = Ratio actual/potential transpiration",
      "! 2 = Depletion of Readily Available Water",
      "! 3 = Depletion of Totally Available Water",
      "! 4 = Depletion of absolute Water Amount",
      "! 5 = Pressure head or moisture content",
      "! 6 = Fixed weekly irrigation, bring root zone back to field capacity",
      "",
      "",
      "* Ratio actual/potential transpiration (TCS = 1)",
      "* If TCS = 1, specify mimimum of ratio actual/potential transpiration Trel [0..1 -, R] as function of crop development stage",
      "DVS_tc1  Trel   ! (maximum 7 records) ",
      "0.0  0.95",
      "2.0  0.95",
      "* End of table",
      "",
      "",
      "* Depletion of Readily Available Water (TCS = 2) ",
      "* If TCS = 2, specify minimum fraction of readily available water RAW [0..1 -, R] as function of crop development stage",
      "DVS_tc2   RAW   ! (maximum 7 records)",
      "0.0  0.95",
      "2.0  0.95",
      "* End of table",
      "",
      "",
      "* Depletion of Totally Available Water (TCS = 3)",
      "* If TCS = 3, specify minimal fraction of totally available water TAW [0..1 -, R] as function of crop development stage",
      "DVS_tc3   TAW   ! (maximum 7 records)",
      "0.0  0.45",
      "2.0  0.45",
      "* End of table",
      "",
      "",
      "* Depletion of absolute Water Amount (TCS = 4)",
      "* If TCS = 4, specify maximum amount of water depleted below field capacity DWA [0..500 mm, R] as function of crop development stage",
      "DVS_tc4   DWA   ! (maximum 7 records)",
      "0.0  40.0",
      "2.0  40.0",
      "* End of table",
      "",
      "",
      "* Pressure head or Moisture content (TCS = 5), specify",
      "PHORMC = 0    ! Switch, use either pressure head (PHORMC = 0) or water content (PHORMC = 1)",
      "DCRIT = -30.0 ! Depth of the sensor [-100..0 cm, R]",
      "* Also specify critical pressure head [-1d6..-100 cm, R] or moisture content [0..1 cm3/cm3, R] as function of crop development stage",
      "DVS_tc5  Value_tc5",
      "0.0    -100.0",
      "2.0    -100.0",
      "* End of table",
      "",
      "* In case TCS = 5, over-irrigation can be applied if the salinity concentration exceeds a threshold salinity",
      "* Switch for over-irrigation:",
      "SWCIRRTHRES = 0    ! 0 = No over-irrigation",
      "! 1 = Apply over-irrigation",
      "* If SWCIRRTHRES = 1, specify:",
      "CIRRTHRES = 8.0    ! Threshold salinity concentration above which over-irrigation occurs [0..100 mg/cm3, R]",
      "PERIRRSURP = 10.0  ! Over-irrigation as percentage of the usually scheduled irrigation depth [0..100 %, R]",
      "",
      "* In case TCS = 6, specify: ",
      "* Fixed weekly irrigation, root zone back to field capacity (TCS = 6), specify",
      "* Threshold value for weekly irrigation; only irrigate when soil water deficit in root zone is larger than threshold",
      "IRGTHRESHOLD = 1.0       ! threshold value  [0..20 mm, R]",
      "",
      "",
      "* Switch for minimum time interval between irrigation applications",
      "TCSFIX = 0       ! 0 = no minimum time interval",
      "! 1 = define minimum time interval",
      "* If TCSFIX = 1, specify:",
      "IRGDAYFIX = 7    ! Minimum number of days between irrigation applications [1..366 d, I]",
      "",
      "**********************************************************************************",
      "",
      "",
      "**********************************************************************************",
      "* Part 17: Irrigation depth criteria",
      "",
      "* Choose one of the following two options for irrigation depth:",
      "DCS = 1      ! 1 = Back to Field Capacity",
      "! 2 = Fixed Irrigation Depth",
      "",
      "",
      "* Back to Field Capacity (DCS = 1)  ",
      "* If DCS = 1, specify amount of under (-) or over (+) irrigation dI [-100..100 mm, R],",
      "* as function of crop development stage [0..2, R]",
      "DVS_dc1    dI        ! (maximum 7 records)",
      "0.0  10.0",
      "2.0  10.0",
      "* End of table",
      "RAITHRESHOLD = 10.0 ! When rainfall exceeds RAITHRESHOLD, irrigation is reduced with rainfall [0..1000 cm, R]",
      "",
      "* Fixed Irrigation Depth (DCS = 2)",
      "* If DCS = 2, specify fixed irrigation depth FID [0..400 mm, R],",
      "* as function of crop development stage [0..2, R]",
      "DVS_dc2   FID         ! (maximum 7 records)",
      "0.0   0.0",
      "0.5  20.0",
      "1.0  25.0",
      "1.9  25.0",
      "2.0   0.0",
      "* End of table",
      "",
      "",
      "* Select minimum and maximum of irrigation depths:",
      "dcslim = 1         ! Switch, limit range irrigation depth  [Y=1, N=0]",
      "",
      "* If dcslim = 1, specify:",
      "irgdepmin = 15.0    !  Minimum irrigation depth [0..100 mm, I]",
      "irgdepmax = 37.0    !  Maximum irrigation depth [irgdepmin..1d7 mm, I]",
      "",
      "",
      "* End of .crp file ",
      "")

  writeLines(tx, con=pathy)
  
  