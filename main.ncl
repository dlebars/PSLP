;##############################################################################
; main.ncl: Compute future total sea level distribution.
;           Calls external functions for each contribution.
;            List of processes considered:
;            - Thermal expansion and SSH changes
;            - Glaciers and Ice Caps
;            - SMB Greenland and Antarctica
;            - Landwater
;            - Dynamics Antarctic Ice Sheet
;            - Dynamics Greenland Ice Sheet
;            - Iverse barometer effect
;##############################################################################
load "./func_ant.ncl"
load "./func_odyn.ncl"
load "./func_misc.ncl"
;##############################################################################
begin

;#### User defined parameters #################################################
VER    = 1               ; File version to write in outputs
N      = tointeger(1e4)  ; Number of sample in all distributions, normally: 1e4 (7*1e5)
MIN_IT = 50              ; Minimum of iterations (bipasses the convergence parameter)
er     = 0.1             ; Convergence parameter, cm difference for the 99th percentile 
                         ;of total ditribution
;####### End of user defined parameters ########################################
if (ANT_DYN.eq."DC16").or.(ANT_DYN.eq."DC16T") then
  ANT_SCE     = 1        ; Four scenarios are possible for DC16 and DC16T:
                         ; (see extended data table 1 from DeConto 2016)
                         ; 0: LIG=3.6-7.4m, Plio. 10-20m, no temp. adj.
                         ; 1: LIG=3.6-7.4m, Plio. 10-20m, with temp. adj., WORST CASE
                         ; 2: LIG=3.6-7.4m, Plio. 5-20m, no temp. adj.
                         ; 3: LIG=3.6-7.4m, Plio. 5-20m, with temp. adj.
end if

start_periods = ispan(2006,2100,1) ; start dates of periods
end_periods   = ispan(2006,2100,1) ; end dates of periods
;start_periods = (/2021,2041,2061,2081/)
;end_periods   = (/2040,2060,2080,2100/)
;start_periods = (/2021,2041,2061,2081,2101,2121/)
;end_periods   = (/2040,2060,2080,2100,2120,2140/)
period        = (/start_periods, \
                  end_periods /)

ROOT        = "/nobackup/users/bars/SeaLevelFromHylke/CMIP5_OCEAN/"
ROOTL       = "./InputData/"
DIR_T       = ROOTL+"Tglobal/"
DIR_IPCC    = ROOTL+"AR5calculations_Jonathan/project/"

if ODYN.eq."CMIP5"
  DIR_OCMIP5  = "/usr/people/bars/Project_ProbSLR/CMIP5_ThermalExp/"
end if

if LOC then
  DIR_F       = ROOT+"Fingerprints/"
  DIR_IBE     = ROOT+"InvertedBarometer/1x1_reg/"
  DIR_IBEmean = ROOT+"InvertedBarometer/globalmeans_from_1x1_glob/"

  ;Box of interest, used for local computations:
  lat_N = 60    ; Original box from Hylke: 51-60N, -3.5,7.5E
  lat_S = 51
  lon_W = -3.5
  lon_E = 7.5

  ;Point of interest (Netherlands), used for local computations:
  lat_Neth = 53
  lon_Neth = 5
  if ODYN.eq."KNMI" then
    DIR_O       = ROOT+"OCEANDATA_MARK/1x1_reg/" 
  else if ODYN.eq."IPCC"
    print("ERROR: This option of ocean dynamics can only be used for global computations")
    exit()
  end if
  end if

  else
    if ODYN.eq."KNMI" then
      DIR_OG      = ROOT+"OCEANDATA_MARK/globalmeans_from_1x1_glob/"
    end if
end if

ProcessNames = (/"Global steric", "Local ocean", "Inverse barometer", "Glaciers",    \
                 "Greenland SMB", "Antarctic SMB", "Landwater", "Antarctic dynamics",\ 
                 "Greenland dynamics", "sum anta.", "Total"/)

if Corr then
  if LOC then  ; Number of components, used to compute correlations efficently.
    nb_comp = 9+1          ;?  To adjust if needed
    NameComponents = (/"Glob. temp.", "Thermal exp.", "Barometer effect", "Glaciers", "Green. SMB", \
                       "Ant. SMB", "Land water", "Ant. dyn.", "Green dyn."/) 
    NameComponentsAll = ""
    do i=0,nb_comp-1
      NameComponentsAll = NameComponentsAll+", "+NameComponents(i)
    end do
    ; Used for correlation ouptput
    print("!!! Warning: This combination hasn't been tested")
    else
    nb_comp = 7+1 
    NameComponents = (/"Glob. temp.", "Thermal exp.", "Glaciers", "Green. SMB", "Ant. SMB", "Land water", \
                       "Ant. dyn.", "Green dyn."/)
    NameComponentsAll = ""
    do i=0,nb_comp-1
      NameComponentsAll = NameComponentsAll+", "+NameComponents(i)
    end do
  end if
end if

;List of model names
MOD = (/"ACCESS1-0","BCC-CSM1-1","CanESM2","CNRM-CM5","CSIRO-Mk3-6-0","GFDL-ESM2G", \
        "GFDL-ESM2M","GISS-E2-R","HadGEM2-CC","HadGEM2-ES","inmcm4","IPSL-CM5A-LR", \
        "IPSL-CM5A-MR","MIROC5","MIROC-ESM-CHEM","MIROC-ESM","MPI-ESM-LR","MPI-ESM-MR", \
        "MRI-CGCM3","NorESM1-ME","NorESM1-M"/)
nb_MOD_AR5 = dimsizes(MOD)

; For KNMI files the SSH has a different name for each scenario
if SCE.eq."rcp45" then
  SSH_VAR  = "ZOSH45"
  else if SCE.eq."rcp85" then
  SSH_VAR  = "ZOSH85"
  end if
end if

if dimsizes(dimsizes(period)).eq.2 then
  nb_PER = dimsizes(period(0,:))
  else if dimsizes(dimsizes(period)).eq.1 then
  nb_PER = 1    ; This option doesn't work beause of the "new" statements
  else
  print("ERROR: The period variable is wrongly defined")
  end if
end if

PER    = new(nb_PER,string) ; Dates of periods considered to plot
do i=0,nb_PER-1
  PER(i) = tostring(period(0,i))+"_"+tostring(period(1,i))
end do

;### General parameters
start_date  = 1986    ; 1862    ; Start reading data
ys   = 2006   ; Starting point for the integration, if this is changed problems in functions
ye   = 2100   ; End year for computation

nb_y  = ye-start_date+1 ; Number of years, 1860 to 2100, most models stop at 2100
                        ; Most models start in 1850 but the Hadley models start in
                        ; 1860 so we just read data from there.
nb_y2 = ye - ys +1

Aoc   = 3.6704e14 ; Ocean Area (m2)
rho_w = 1e3       ; Water density (kg.m-3)
fac   = -1e12/(Aoc*rho_w) ; Convert Giga tones to m sea level
alpha_98 = cdfnor_x(0.98,0,1)
alpha_95 = cdfnor_x(0.95,0,1)
alpha_05 = cdfnor_x(0.05,0,1)
PI       = 3.141592653589793

;### Specific parameters
;# GIC: Values of f and p, fitting parameter from formula B.1 (de Vries et al. 2014)
f    = (/3.02,4.96,5.45,3.44/)
p    = (/0.733,0.685,0.676,0.742/)
;# Greenland SMB:
;Values of the 3 parameters from formula B.2 (de Vries et al. 2014)
a = -71.5
b = -20.4
c = -2.8
;# Antarctic SMB
Cref     = 1923         ; Reference accumulation (Gt/year)
Delta_C  = 5.1/100      ; Increase in accumulation with local temperature (%/degC)
Delta_Ce = 1.5/100      ; Error bar
AmpA     = 1.1          ; Antarctic amplification of global temperature
AmpAe    = 0.2          ; Amplification uncertainty
;# Antarctic dynamics
a1_up_a                = 0.061    ; First order term (cm/y), equal to observations in 2006
a1_lo_a                = 0.021
Delta_ant_up_2100      = 18.5
Delta_ant_up_KNMI_2100 = 40.
Delta_ant_lo_2100      = 0.      ; IPCC takes -2
if ANT_DYN.eq."KNMI16" then
  if SCE.eq."rcp45" then
    Delta_ant_up_new_2100 = 49.
    else if SCE.eq."rcp85" then
      Delta_ant_up_new_2100 = 105.
    else
      print("Scenario not compatible with ANT_DYN is KNMI16")
    end if
  end if
end if

;### Refernce time period to take the reference global temperature
; Reference period for Glaciers and Ice Caps
ysr_gic  = 1986          ; Start of the reference period
yer_gic  = 2005          ; End of the reference period
; Reference period for Greenland SMB
ysr_g  = 1980          ; Start of the reference period
yer_g  = 1999          ; End of the reference period
; Reference period for Antarctic SMB
ysr_a  = 1985          ; Start of the reference period
yer_a  = 2005          ; End of the reference period
; Reference period for Antarctic Dynamics, DC16T option
ysr_ad = 2000          ; Start of the reference period
yer_ad = 2000          ; End of the reference period

;### Parameters to produce PDF
;Range to produce PDF of the distributions
opt                    = True
opt@bin_min            = -20.5
opt@bin_max            = 500.5     ;260.5
nbin                   = toint(opt@bin_max-opt@bin_min)

;###
nb_par     = dimsizes(f)   ; Number of parameter sets

TIME       = ispan(start_date,ye,1)
TIME2      = ispan(ys,ye,1)
ind_d      = ind(TIME2%10.eq.0) ;Select the indices of 2010, 2020... 2100
nb_yd      = dimsizes(ind_d)

i_ys       = closest_val(ys,TIME)
i_ye       = closest_val(ye,TIME)
i_ysr_gic  = closest_val(ysr_gic,TIME)
i_yer_gic  = closest_val(yer_gic,TIME)
i_ysr_g    = closest_val(ysr_g,TIME)
i_yer_g    = closest_val(yer_g,TIME)
i_ysr_a    = closest_val(ysr_a,TIME)
i_yer_a    = closest_val(yer_a,TIME)
i_ysr_ad   = closest_val(ysr_ad,TIME)
i_yer_ad   = closest_val(yer_ad,TIME)

i_ysep     = new((/2,nb_PER/),integer)      ; Start and end periods index
do i=0,nb_PER-1                             ; Loop periods
  do j=0,1                                  ; Start and end years (included)
    i_ysep(j,i)    = closest_val(period(j,i),TIME2)
  end do
end do

;#### Read finger prints, some are time dependent so make all of them  time 
; dependent for easy mutliplication at the end.
F_gic2  = new(nb_y2,float)
F_gsmb2 = new(nb_y2,float)
F_asmb2 = new(nb_y2,float)
F_gdyn2 = new(nb_y2,float)
F_adyn2 = new(nb_y2,float)
F_gw2   = new(nb_y2,float)

if LOC then
  FPs    = 1986 ; The fingerprints start in 1986 and end in 2100,
  FPe    = 2100 ; but do not have a useful time vector so we create one here
  TIME_F = ispan(FPs,FPe,1)
  ifs    = closest_val(ys,TIME_F) 
  ife    = closest_val(ye,TIME_F)
  jfs    = closest_val(FPs,TIME2)
  jfe    = closest_val(FPe,TIME2)

  ;# Read fingerprint for Glaciers and Ice Caps
  f_gic      = addfile(DIR_F+"Relative_GLACIERS_reg.nc","r")
  lat_gic    = tofloat(f_gic->latitude)
  lon_gic    = tofloat(f_gic->longitude)
  finger_gic = f_gic->RSL
  finger_gic@_FillValue = 0
  F_gic      = finger1D(lat_Neth, lon_Neth, lat_gic, lon_gic, finger_gic) 
  F_gic2(jfs:jfe) =  F_gic(ifs:ife)/100 ; Convert from % to fraction

  delete(f_gic)
  delete(lat_gic)
  delete(lon_gic)
  delete(finger_gic)
  delete(TIME_F)

  f_ic        = addfile(DIR_F+"Relative_icesheets_reg.nc","r")
  lat_ic      = tofloat(f_ic->latitude)
  lon_ic      = tofloat(f_ic->longitude)
  finger_gsmb = f_ic->SMB_GRE
  finger_gsmb@_FillValue = 0
  finger_asmb = f_ic->SMB_ANT
  finger_asmb@_FillValue = 0
  finger_gdyn = f_ic->DYN_GRE
  finger_gdyn@_FillValue = 0
  finger_adyn = f_ic->DYN_ANT
  finger_adyn@_FillValue = 0
  F_gsmb      = finger1D(lat_Neth, lon_Neth, lat_ic, lon_ic, finger_gsmb)
  F_gsmb2(jfs:jfe)  =  F_gsmb/100
  F_asmb      = finger1D(lat_Neth, lon_Neth, lat_ic, lon_ic, finger_asmb)
  F_asmb2(jfs:jfe)  =  F_asmb/100
  F_gdyn      = finger1D(lat_Neth, lon_Neth, lat_ic, lon_ic, finger_gdyn)
  F_gdyn2(jfs:jfe)  =  F_gdyn/100
  F_adyn      = finger1D(lat_Neth, lon_Neth, lat_ic, lon_ic, finger_adyn)
  F_adyn2(jfs:jfe)  =  F_adyn/100

  delete(f_ic)
  delete(lat_ic)
  delete(lon_ic)
  delete(finger_gsmb)
  delete(finger_asmb)
  delete(finger_gdyn)
  delete(finger_adyn)

  f_gw       = addfile(DIR_F+"Relative_GROUNDWATER_reg.nc","r")
  lat_gw     = tofloat(f_gw->latitude)
  lon_gw     = tofloat(f_gw->longitude)
  finger_gw  = f_gw->GROUND
  finger_gw@_FillValue = 0
  F_gw       = finger1D(lat_Neth, lon_Neth, lat_gw, lon_gw, finger_gw)
  F_gw2(jfs:jfe)  =  F_gw(ifs:ife)/100

  delete(f_gw)
  delete(lat_gw)
  delete(lon_gw)
  delete(finger_gw)

  ;Extend the values of the fingerprints after 2100:
  if jfe.lt.(nb_y2-1) then
    print("Extending fingerprints using last available value")
    F_gic2(jfe+1:) = F_gic2(jfe)
    F_gsmb2(jfe+1:) = F_gsmb2(jfe)
    F_asmb2(jfe+1:) = F_asmb2(jfe)
    F_gdyn2(jfe+1:) = F_gdyn2(jfe)
    F_adyn2(jfe+1:) = F_adyn2(jfe)
    F_gw2(jfe+1:) = F_gw2(jfe)
  end if
  delete(ifs)
  delete(ife)
  delete(jfs)
  delete(jfe)

  else  ;### For global average put all the fingerprints to 1
  F_gic2  = 1
  F_gsmb2 = 1
  F_asmb2 = 1
  F_gdyn2 = 1
  F_adyn2 = 1
  F_gw2   = 1

end if

;##############################################################################
if INFO then 
  print("### Read Tglob             #################")
end if

if TEMPf.eq."all" then
  files     = systemfunc("ls "+DIR_T+"global_tas_Amon_*_"+SCE+"_r1i1p1.dat")
  else if TEMPf.eq."AR5" then
    files_i     = new(nb_MOD_AR5,string)
    do m=0,nb_MOD_AR5-1
      if MOD(m).eq."BCC-CSM1-1" then
        loc_mod = "bcc-csm1-1"
        else
          loc_mod = MOD(m)
      end if
      files_i(m) = systemfunc("ls "+DIR_T+"global_tas_Amon_"+loc_mod+"_"+SCE+"_r1i1p1.dat")
    end do
    files = files_i(ind(.not.ismissing(files_i)))
    delete(files_i)    
  else
    print("Option TEMPf: "+TEMPf+" is not supported")
  end if
end if
nb_MOD    = dimsizes(files)

;Read Tglob
TGLOB    = new((/nb_MOD,nb_y/),float)
Tref_gic = new(nb_MOD,float)   ; Each process needs a different time anomaly
Tref_g   = new(nb_MOD,float)
Tref_a   = new(nb_MOD,float)
Tref_ad  = new(nb_MOD,float)   ; Antarctic dynamics for DC16T option

if INFO then
  print("Models used for "+SCE)
  print(files)
end if

do m=0,nb_MOD-1
  TEMP     = readAsciiTable(files(m),13,"float",3)
  time     = TEMP(:,0)
  dim_t    = dimsizes(time)
  i_start  = closest_val(start_date,time)
  i_end    = closest_val(ye,time)
  TGLOB(m,:i_end-i_start) = dim_avg_n(TEMP(i_start:i_end,1:12),1)   ;Data in degree Kelvin
  Tref_gic(m) = avg(TGLOB(m,i_ysr_gic:i_yer_gic))
  Tref_g(m)   = avg(TGLOB(m,i_ysr_g:i_yer_g))
  Tref_a(m)   = avg(TGLOB(m,i_ysr_a:i_yer_a))
  Tref_ad(m)  = avg(TGLOB(m,i_ysr_ad:i_yer_ad))
  ;### Issue of missing temperature value for rcp26 after 2100 for this scenario
  ; it is ok to assume it is constant after 2100
  if (SCE.eq."rcp26").and.(ye.gt.2100) then
    i2100 = closest_val(2099,time)
    print(i2100)
    TGLOB(m,i2100-i_start:) = TGLOB(m,i2100-i_start)
  end if
  delete(TEMP)
  delete(time)
end do
delete(files)

TGLOBs = TGLOB(:,i_ys:)
delete(TGLOB)

;#############################################################################
; Start loop 
CONV  = NewList("lifo")
nb_it = 0

X_O_G_pdf     = new((/nb_PER,nbin/),double)
X_O_A_pdf     = new((/nb_PER,nbin/),double)
X_B_pdf       = new((/nb_PER,nbin/),double)
X_gic_pdf     = new((/nb_PER,nbin/),double)
X_gsmb_pdf    = new((/nb_PER,nbin/),double)
X_asmb_pdf    = new((/nb_PER,nbin/),double)
X_landw_pdf   = new((/nb_PER,nbin/),double)
X_ant_pdf     = new((/nb_PER,nbin/),double)
X_gre_pdf     = new((/nb_PER,nbin/),double)
X_ant_tot_pdf = new((/nb_PER,nbin/),double)
X_tot_pdf     = new((/nb_PER,nbin/),double)

if Corr then
  X_all       = new((/nb_comp,N,nb_yd/),float)
  nb_el       = nb_comp*(nb_comp-1)/2
  M_Corr_P    = new((/nb_yd,nb_el/),float)    ; Pearson correlation matrice 
  M_Corr_S    = new((/nb_yd,nb_el/),float)    ; Rank correlation matrice
  M_Corr_P    = 0
  M_Corr_S    = 0
end if

X_O_G_pdf     = 0
X_O_A_pdf     = 0
X_B_pdf       = 0
X_gic_pdf     = 0
X_gsmb_pdf    = 0
X_asmb_pdf    = 0
X_landw_pdf   = 0
X_ant_pdf     = 0
X_gre_pdf     = 0
X_ant_tot_pdf = 0
X_tot_pdf     = 0

END           = False

do while .not.(END)
  nb_it = nb_it + 1
  comp  = 1    ; Re-initialize the component count for the Corr case

  ;### Set seeds for the random number generators
  rseed1 = toint(systemfunc(" date +%s"))
  rseed2 = toint((12345l*rseed1)%2147483398l)
  random_setallseed(rseed1, rseed2)

  ; Sample a normal distribution to use for temperature
  NormD  = random_normal(0,1,N)

  if COMB.eq."DEP" then
    ; Reorder contribution in ascending order
    ip   = dim_pqsort_n(NormD, 2, 0)
    delete(ip)
  end if

  if Corr then
    ; NormD is used in the correlation with temperature because all of the temperature
    ; distributions (different for each process) are based on a linear combination of it. 
    NormDc       = conform_dims((/N,nb_yd/),NormD,0)
    X_all(0,:,:) = NormDc
  end if

  ;##############################################################################
  if INFO then
    print("### Ocean dynamics ###################################################")
  end if
  if COMB.eq."IND" then
    CorrGT   = 0 ; Force a 0 correlation, even if the CorrGT coefficient has another value
    else if COMB.eq."DEP" then
      CorrGT = 1
    end if
  end if
  NormDT1 = random_normal(0,1,N)
  ; Build NormDT as a combination of NormD (the distribution of GMST) and an independent
  ; normal distribution.
  if CorrM.eq."Pearson" then
    rhoP  = CorrGT
    else if CorrM.eq."Spearman" then
      ; Convert correlation coefficient from Spearman to Pearson
      rhoP  = 2*sin(PI/6*CorrGT)
    end if
  end if
  NormDT          = NormD*rhoP + NormDT1*sqrt(1 - rhoP^2)

  if LOC then
    if ODYN.eq."KNMI" then
      X_Of = odyn_loc(SCE, MOD, nb_y, nb_y2, DIR_O, lat_N, lat_S, lon_W, \
                      lon_E, start_date, ye, SSH_VAR, N, i_ys, GAM, NormDT)
      else if ODYN.eq."CMIP5" then
        X_Of = odyn_cmip5(SCE, LOC, DIR_OCMIP5, N, ys, ye, GAM, NormDT)
    end if
    end if

    else
    if ODYN.eq."KNMI" then
      X_Of = odyn_glob_knmi(SCE, MOD, nb_y, nb_y2, DIR_O, DIR_OG, \
             start_date, ye, SSH_VAR, N, i_ys, GAM, NormDT)
    else if ODYN.eq."IPCC" then
      X_Of = odyn_glob_ipcc(SCE, DIR_IPCC, N, nb_y2, GAM, NormDT)
    else if ODYN.eq."CMIP5" then
      X_Of = odyn_cmip5(SCE, LOC, DIR_OCMIP5, N, ys, ye, GAM, NormDT)
    end if
    end if
    end if
  end if

  ;Compute the pdfs based on the chosen periods
  do per=0,nb_PER-1        ; Loop on the period
    X_O_G_pdf(per,:) = X_O_G_pdf(per,:) + \
                         pdfx(dim_avg_n(X_Of(1,:,i_ysep(0,per):i_ysep(1,per)),1),nbin,opt) 
    X_O_A_pdf(per,:) = X_O_A_pdf(per,:) + \ 
                         pdfx(dim_avg_n(X_Of(2,:,i_ysep(0,per):i_ysep(1,per)),1),nbin,opt)
  end do

  ;Update X_tot, the sum of all contributions
  if COMB.eq."DEP" then
    ; Reorder contribution in ascending order
    ip   = dim_pqsort_n(X_Of, 2, 1)
    delete(ip)
  end if

  X_tot = new((/N,nb_y2/),float)

  if NoU_O then
    do y=0,nb_y2-1        ; Loop on the period
      X_tot(:,y) = avg(X_Of(0,:,y))
    end do
  else
    X_tot = X_Of(0,:,:)
  end if

  if Corr then
    X_all(comp,:,:) = X_Of(0,:,ind_d)
    comp            = comp + 1
  end if

  delete(X_Of)
  delete(NormDT1)
  delete(NormDT)

  ;##############################################################################
  if INFO then
    print("### Inverse Barometer effect #########################################")
  end if

  ;Initialize the SSH matrix: (model, time (years))
  ;Use the AR5 model list and not the complete one
  ; because the fields should be downloaded for other models
  MAT = new((/nb_MOD_AR5,nb_y/),float)
  if LOC.and.IBarE then
    do m=0,nb_MOD_AR5-1
      if MOD(m).eq."BCC-CSM1-1" then
        loc_mod = "bcc-csm1-1"
        else
          loc_mod = MOD(m)
      end if
      fb    = addfile(DIR_IBE+"ssh_ibe_"+loc_mod+"_"+SCE+".nc","r")
      fmean = addfile(DIR_IBEmean+"ssh_ibe_"+loc_mod+"_"+SCE+".nc","r")
      lon = fb->longitude  ; Do not assume they all have the same grid, even though they DO
      lat = fb->latitude
      lat_Ni  = closest_val(lat_N,lat)
      lat_Si  = closest_val(lat_S,lat)
      lon_Wi  = closest_val(lon_W,lon)
      lon_Ei  = closest_val(lon_E,lon)
      TIMEt   = fb->TIME
      TIMEt2  = cd_calendar(TIMEt,0)
      i_start = closest_val(start_date,TIMEt2(:,0))
      i_end   = closest_val(ye,TIMEt2(:,0))
      IBE     = fb->ssh_ibe(i_start:i_end,lat_Si:lat_Ni,lon_Wi:lon_Ei)
      dimIBE  = dimsizes(IBE)
      IBEmean = fmean->ssh_ibe(i_start:i_end,0,0)
      nb_y_loop = i_end - i_start +1
      if nb_y_loop.eq.nb_y then
        MAT(m,:) = dim_avg_n(dim_avg_n(IBE,2),1) - IBEmean
		     ;RQ: No scaling depending on the area, gives more weight to the southern cell
      else
        MAT(m,:nb_y-2) = dim_avg_n(dim_avg_n(IBE,2),1) - IBEmean
        MAT(m,nb_y-1)  = MAT(m,nb_y-2)
      end if
      ; Remove the average SSH of the initial period from all models
      MAT(m,:) =  MAT(m,:) - avg(MAT(m,:20))
      delete(fb)
      delete(fmean)
      delete(lon)
      delete(lat)
      delete(IBE)
      delete(TIMEt)
      delete(TIMEt2)
    end do
  else ; Global case, no IBE
    MAT = 0
  end if

  MATs = MAT(:,i_ys:)*100 ; Convert from m to cm

  ;Build the distribution
  NormDi = random_normal(0,1,N) ; Independent Normal distribution
  X_B_m  = dim_avg_n(MATs,0)    ; Compute the inter-model mean for each time
  X_B_sd = dim_stddev_n(MATs,0) ; Compute the inter-model standard deviation

  X_B    = new((/N,nb_y2/),float)
  do t=0,nb_y2-1
    X_B(:,t)    = X_B_m(t) + GAM*NormDi(:)*X_B_sd(t)
  end do

  ;Compute the pdfs based on the chosen periods
  do per=0,nb_PER-1        ; Loop on the period
    X_B_pdf(per,:)  = X_B_pdf(per,:) + \
                        pdfx(dim_avg_n(X_B(:,i_ysep(0,per):i_ysep(1,per)),1),nbin,opt)
  end do

  ;Update X_tot, the sum of all contributions
  if COMB.eq."DEP" then
    ; Reorder contribution in ascending order
    ip   = dim_pqsort_n(X_B, 2, 0)
    delete(ip)
  end if

  X_tot = X_tot + X_B

  delete(MAT)
  delete(MATs)
  delete(NormDi)

  if Corr.and.LOC then
    X_all(comp,:,:) = X_B(:,ind_d)
    comp            = comp + 1
  end if
  delete(X_B)

  ;##############################################################################
  if INFO then
    print("### Glaciers and Ice Caps (Not including Antarctica) #################")
  end if

  if COMB.eq."IND" then
    ; Redefine NormD to loose correlation
    NormD  = random_normal(0,1,N)
  end if

  ;Build the distribution of global temperature for this process
  Td_gic = TempDist(TGLOBs, Tref_gic, GAM, NormD)

  ;Convert the global temperature to glaciers
  Icum = new((/N,nb_y2/),float)
  do r=0,N-1
    Icum(r,:) = cumsum(Td_gic(r,:),0)
  end do
  ; Following formula is not defined for negative values of Icum, which happens for
  ; a few extreme values of the distribution.
  Icum = where(Icum.lt.0,0,Icum)

  DEL_gic        = new((/N,nb_y2/),float)
  DEL_gic(:,:)   = (f(0)*Icum^p(0)+f(1)*Icum^p(1)+f(2)*Icum^p(2)+f(3)*Icum^p(3))/4
  DEL_gic        = DEL_gic*.1   ; Convert from mm to cm

  X_gic = new((/N,nb_y2/),float)
  NormDs  = random_normal(0,1,N)   ; This distribution is then kept for correlation
  NormDl  = 1 + NormDs*0.2
  do t=0,nb_y2-1
    X_gic(:,t) = DEL_gic(:,t)*NormDl
  end do

  X_gic = X_gic +0.95 ; Add 0.95 cm for the changes between 1996 to 2006

  do t=0,nb_y2-1
    X_gic(:,t) = X_gic(:,t)*F_gic2(t)
  end do

  ;Compute the pdfs based on the chosen periods
  do per=0,nb_PER-1        ; Loop on the period
    X_gic_pdf(per,:)  = X_gic_pdf(per,:) + \
                          pdfx(dim_avg_n(X_gic(:,i_ysep(0,per):i_ysep(1,per)),1),nbin,opt)
  end do

  ;Update X_tot, the sum of all contributions
  if COMB.eq."DEP" then
    ; Reorder contribution in ascending order
    ip   = dim_pqsort_n(X_gic, 2, 0)
    delete(ip)
  end if

  if NoU_Gl then
    do y=0,nb_y2-1        ; Loop on the period
      X_tot(:,y) = X_tot(:,y) + avg(X_gic(:,y))
    end do
  else
    X_tot = X_tot + X_gic
  end if

  if Corr then
    X_all(comp,:,:) = X_gic(:,ind_d)
    comp        = comp + 1
  end if
  delete(X_gic)
  delete(Icum)
  delete(DEL_gic)
  delete(NormDl)
  delete(Td_gic)

  ;##############################################################################
  if INFO then
    print("### SMB Greenland ####################################################")
  end if

  if COMB.eq."IND" then
    ; Redefine NormD to loose correlation
    NormD  = random_normal(0,1,N)
  end if

  ;Build the distribution of global temperature for this contributor
  Td_g  = TempDist(TGLOBs, Tref_g, GAM, NormD)

  ;Convert the global temperature to Greenland SMB changes
  G_SMB_cum = new((/N,nb_y2/),float)
  G_SMB     = fac*(a*Td_g + b*Td_g^2 +c*Td_g^3)

  do r=0,N-1
    G_SMB_cum(r,:) = cumsum(G_SMB(r,:),0)
  end do

  ;Build the distibution including other uncertainties
  Unif   = random_uniform(1,1.15,N)
  if CorrSMB then
    NormDl = NormDs
    else
      NormDl = random_normal(0,1,N)
  end if
  LogNor = exp(NormDl*0.4)
  X_gsmb = new((/N,nb_y2/),float)
  do t=0,nb_y2-1
    X_gsmb(:,t) = G_SMB_cum(:,t)*Unif*LogNor
  end do
  X_gsmb = X_gsmb*100 ; Convert from m to cm

  if GRE.eq."KNMI14" then
    X_gsmb = X_gsmb +0.291 ; 0.291 is used for KNMI14 because nothing is added to 
			   ; the dynamics while IPCC uses 0.15 cm as well there
    else if GRE.eq."IPCC" then
    X_gsmb = X_gsmb +0.15 ; Add 0.15 cm for contribution between 1995-2005, 
    else
    print("ERROR: Wrong GRE option")
    end if
  end if

  do t=0,nb_y2-1
    X_gsmb(:,t) = X_gsmb(:,t)*F_gsmb2(t)
  end do

  ;Compute the pdfs based on the chosen periods
  do per=0,nb_PER-1        ; Loop on the period
    X_gsmb_pdf(per,:)  = X_gsmb_pdf(per,:) + \
                         pdfx(dim_avg_n(X_gsmb(:,i_ysep(0,per):i_ysep(1,per)),1),nbin,opt)
  end do

  ;Update X_tot, the sum of all contributions
  if COMB.eq."DEP" then
    ; Reorder contribution in ascending order
    ip   = dim_pqsort_n(X_gsmb, 2, 0)
    delete(ip)
  end if

  if NoU_G then
    do y=0,nb_y2-1        ; Loop on the period
      X_tot(:,y) = X_tot(:,y) + avg(X_gsmb(:,y))
    end do
  else
    X_tot = X_tot + X_gsmb
  end if

  delete(Unif)
  delete(LogNor)
  delete(NormDl)
  delete(Td_g)
  delete(G_SMB_cum)
  delete(G_SMB)
  if Corr then
    X_all(comp,:,:) = X_gsmb(:,ind_d)
    comp            = comp + 1
  end if
  delete(X_gsmb)

  ;##############################################################################
  if INFO then
    print("### SMB Antarctica ###################################################")
  end if

  if COMB.eq."IND" then
    ; Redefine NormD to loose correlation
    NormD  = random_normal(0,1,N)
  end if

  if (ANT_DYN.eq."IPCC").or.(ANT_DYN.eq."KNMI16").or.(ANT_DYN.eq."LEV14") then
    ;Build the distribution of global temperature for this contributor
    Td_a = TempDist(TGLOBs, Tref_a, GAM, NormD)

    NormD2 = random_normal(AmpA,AmpAe,N)
    if CorrSMB then
      NormDl = NormDs
      else
        NormDl = random_normal(0,1,N)  ;
    end if
    NormDl = Delta_C + NormDl*Delta_Ce

    X_asmb = new((/N,nb_y2/),float)
    do t=0,nb_y2-1
      X_asmb(:,t) = fac*Cref*Td_a(:,t)*NormD2*NormDl
    end do
    
    ; Cumulate over time
    do r=0,N-1
      X_asmb(r,:) = cumsum(X_asmb(r,:),0)
    end do
    X_asmb = X_asmb*100 ; Convert from m to cm

    delete(Td_a)
    delete(NormD2)
    delete(NormDl)

  else
    ;In other cases DC16 and DC16T the SMB is included in the dynamics
    X_asmb = new((/N,nb_y2/),float)
    X_asmb = 0
  end if

  do t=0,nb_y2-1
    X_asmb(:,t) = X_asmb(:,t)*F_asmb2(t)
  end do

  ;Compute the pdfs for the the chosen periods
  do per=0,nb_PER-1        ; Loop on the period
    X_asmb_pdf(per,:)  = X_asmb_pdf(per,:) + \
                           pdfx(dim_avg_n(X_asmb(:,i_ysep(0,per):i_ysep(1,per)),1),nbin,opt)
  end do

  ;Update X_tot, the sum of all contributions
  if COMB.eq."DEP" then
    ; Reorder contribution in ascending order
    ip   = dim_pqsort_n(X_asmb, 2, 0)
    delete(ip)
  end if

  X_ant_tot  = new((/N,nb_y2/),float)
  if NoU_A then
    do y=0,nb_y2-1           ; Loop on the period
      X_tot(:,y)     = X_tot(:,y) + avg(X_asmb(:,y))
      X_ant_tot(:,y) = avg(X_asmb(:,y))
    end do
  else
    X_tot     = X_tot + X_asmb
    X_ant_tot = X_asmb
  end if

  if Corr then
    X_all(comp,:,:) = X_asmb(:,ind_d)
    comp            = comp + 1
  end if
  ;delete(X_asmb)

  ;##############################################################################
  if INFO then
    print("### Landwater changes ################################################")
  end if

  ;# Land water
  a1_up             = 0.049 ; First order term (cm/y), equal to observations in 2006
  a1_lo             = 0.026
  Delta_grw_up_2100 = 9
  Delta_grw_lo_2100 = -1

  ;Compute the second order coefficient of the equations:
  a2_up  = (Delta_grw_up_2100 - a1_up*(2100-ys))/(2100 - ys)^2
  a2_lo  = (Delta_grw_lo_2100 - a1_lo*(2100-ys))/(2100 - ys)^2

  Delta_grw_up = a1_up*(TIME2-ys) + a2_up*(TIME2-ys)^2
  Delta_grw_lo = a1_lo*(TIME2-ys) + a2_lo*(TIME2-ys)^2

  ;Build the distibution including other uncertainties
  sig2   = (Delta_grw_up - Delta_grw_lo)/(alpha_95 - alpha_05)
  NormD2 = random_normal(0,1,N)       ; The standard deviations sould be scaled.

  X_landw = new((/N,nb_y2/),float)  ; Independent of the scenario, size is to add up easily.
  do t=0,nb_y2-1
      X_landw(:,t) = (Delta_grw_up(t) - Delta_grw_lo(t))/2 + sig2(t)*NormD2
  end do

  do t=0,nb_y2-1
    X_landw(:,t) = X_landw(:,t)*F_gw2(t)
  end do

  ;Compute the pdfs based on the chosen periods
  do per=0,nb_PER-1        ; Loop on the period
    X_landw_pdf(per,:)  = X_landw_pdf(per,:) + \
                            pdfx(dim_avg_n(X_landw(:,i_ysep(0,per):i_ysep(1,per)),1),nbin,opt)
  end do

  ;Update X_tot, the sum of all contributions
  if COMB.eq."DEP" then
    ; Reorder contribution in ascending order
    ip   = dim_pqsort_n(X_landw, 2, 0)
    delete(ip)
  end if

  X_tot = X_tot + X_landw

  delete(Delta_grw_up)
  delete(Delta_grw_lo)
  delete(NormD2)
  if Corr then
    X_all(comp,:,:) = X_landw(:,ind_d)
    comp        = comp + 1
  end if
  delete(X_landw)

  ;##############################################################################
  if INFO then
    print("### Antarctic dynamics ###############################################")
  end if

  if COMB.eq."IND" then
    ; Redefine NormD to loose correlation
    NormD  = random_normal(0,1,N)
  end if

  if ANT_DYN.eq."KNMI16"
    X_ant = ant_dyn_knmi16(Delta_ant_up_2100, Delta_ant_up_KNMI_2100, Delta_ant_up_new_2100, \
			   Delta_ant_lo_2100, a1_up_a, a1_lo_a, ys, ye,      \
			   TIME2, N)
    else if ANT_DYN.eq."IPCC" then
      Unif_AA   = random_uniform(0,1,N)
      ;### The Delta_ant_lo_2100 is hardcoded to -2 because we changed it to 0 for KNMI2016
      X_ant = proj2order(TIME2, a1_up_a, a1_lo_a, Delta_ant_up_2100, -2, Unif_AA)

      ;### This is the influence of increased SMB on an increase in Antactic dynamics
      if COMB.eq."IND" then
        irg = generate_unique_indices(N)
        do t=0,nb_y2-1
          X_ant(:,t) = X_ant(:,t) - Unif_AA*0.35*X_asmb(irg,t)
        end do
      else
        do t=0,nb_y2-1
          X_ant(:,t) = X_ant(:,t) - Unif_AA*0.35*X_asmb(:,t)
        end do
      end if

    else if ANT_DYN.eq."DC16" then
      TempDep = new(N, float)
      TempDep = 0 ; No temperature dependence
      X_ant   = ant_dyn_dc16(SCE, N, TempDep, ANT_SCE, a1_up_a, a1_lo_a, TIME2) 

    else if ANT_DYN.eq."DC16T" then
      ;Build the distribution of global temperature for Antarctic dynamics, 
      ;focus on the last date (2100).

      ;Need to read temperature from another climate scenario
      if SCE.eq."rcp45" then
        SCE_loc = "rcp85"
        else if SCE.eq."rcp85" then
          SCE_loc = "rcp45"
        else
          print("ERROR: Scenario not supported by DC16T option")
        end if
      end if
      files     = systemfunc("ls "+DIR_T+"global_tas_Amon_*_"+SCE_loc+"_r1i1p1.dat")
      nb_MOD2   = dimsizes(files)

      ;Read Tglob
      TGLOB    = new((/nb_MOD2,nb_y/),float)
      Tref_ad2 = new(nb_MOD2,float)   ; Antarctic dynamics for DC16T option

      do m=0,nb_MOD2-1
        TEMP     = readAsciiTable(files(m),13,"float",3)
        time     = TEMP(:,0)
        dim_t    = dimsizes(time)
        i_start  = closest_val(start_date,time)
        i_end    = closest_val(ye,time)
        TGLOB(m,:i_end-i_start) = dim_avg_n(TEMP(i_start:i_end,1:12),1)   ;Data in degree Kelvin
        Tref_ad2(m)  = avg(TGLOB(m,i_ysr_ad:i_yer_ad))
        delete(TEMP)
        delete(time)
      end do
      delete(files)
      TGLOBs2 = TGLOB(:,i_ys:)
      delete(TGLOB)

      TGLOB_ad = new((/nb_MOD,nb_y2/),float)
      do m=0,nb_MOD-1
        TGLOB_ad(m,:)    = TGLOBs(m,:) - Tref_ad(m)
      end do
      TGLOB_ad_m    = dim_avg_n(TGLOB_ad,0)    ; Compute the inter-model mean for each time
      TGLOB_ad_sd   = dim_stddev_n(TGLOB_ad,0) ; Compute the inter-model standrad deviation

      TGLOB_ad2 = new((/nb_MOD2,nb_y2/),float)
      do m=0,nb_MOD2-1
        TGLOB_ad2(m,:)    = TGLOBs2(m,:) - Tref_ad2(m)
      end do
      delete(TGLOBs2)
      TGLOB_ad2_m    = dim_avg_n(TGLOB_ad2,0)    ; Compute the inter-model mean for each time

      ;Compute the temperature dependence
      TempDep    = GAM*NormD(:)*TGLOB_ad_sd(nb_y2-1)/abs(TGLOB_ad2_m(nb_y2-1)-TGLOB_ad_m(nb_y2-1))

      print("Difference of mean temperature")
      print(TGLOB_ad_m(1)-TGLOB_ad_m(0))
      print("Average TempDep:")
      print(dim_avg_n(TempDep,0))
      X_ant = ant_dyn_dc16(SCE, N, TempDep, ANT_SCE, a1_up_a, a1_lo_a, TIME2)

    else if ANT_DYN.eq."LEV14" then
      UnifDd = random_uniform(0,1,N)
      X_ant  = ant_dyn_lev14(ROOTL, SCE, MOD, ys, GAM, NormD, UnifDd, DIR_T)
    else
      print("ERROR: This parameterisation of Antarctic dynamical processes does not exist")
    end if
    end if
    end if  
    end if
  end if

  X_ant = X_ant + 0.25 ; Add 0.25cm, the conribution from 1995 to 2005

  do t=0,nb_y2-1
    X_ant(:,t) = X_ant(:,t)*F_adyn2(t)
  end do

  ;Compute the pdfs based on user defined periods of time
  do per=0,nb_PER-1        ; Loop on the period
    X_ant_pdf(per,:)  = X_ant_pdf(per,:) + \
                        pdfx(dim_avg_n(X_ant(:,i_ysep(0,per):i_ysep(1,per)),1),nbin,opt)
  end do

  ;Update X_tot, the sum of all contributions
  if COMB.eq."DEP" then
    ; Reorder contribution in ascending order
    ip   = dim_pqsort_n(X_ant, 2, 0)
    delete(ip)
  end if

  if NoU_A then
    do y=0,nb_y2-1           ; Loop on the period
      X_tot(:,y)     = X_tot(:,y) + avg(X_ant(:,y))
      X_ant_tot(:,y) = X_ant_tot(:,y) + avg(X_ant(:,y))
    end do
  else
    X_tot     = X_tot + X_ant
    X_ant_tot = X_ant_tot + X_ant
  end if

  ;Compute the pdfs based on user defined periods of time
  do per=0,nb_PER-1        ; Loop on the period
    X_ant_tot_pdf(per,:)  = X_ant_tot_pdf(per,:) + \
                            pdfx(dim_avg_n(X_ant_tot(:,i_ysep(0,per):i_ysep(1,per)),1),nbin,opt)
  end do

  if Corr then
    X_all(comp,:,:) = X_ant(:,ind_d)
    comp            = comp + 1
  end if
  delete(X_ant)
  delete(X_ant_tot)
  delete(X_asmb)

  ;##############################################################################
  if INFO then
    print("### Greenland dynamics ###############################################")
  end if 

  a1_up_g           = 0.5*0.076    ; First order term (cm/y), equal to observations in 2006
  a1_lo_g           = 0.5*0.043

  if .not.CorrDYN then
    UnifDd = random_uniform(0,1,N)  ; Sample a new independent distribution
  end if

  if GRE.eq."KNMI14" then
    Delta_gre_up_2100 = 7.4
    Delta_gre_lo_2100 = 1.7
    X_gre  = proj2order(TIME2, a1_up_g, a1_lo_g, Delta_gre_up_2100, \
             Delta_gre_lo_2100, UnifDd)
    else if GRE.eq."IPCC" then
      if (SCE.eq."rcp26").or.(SCE.eq."rcp45") then
        Delta_gre_up_2100 = 6.3
        Delta_gre_lo_2100 = 1.4
        else if (SCE.eq."rcp85") then
          Delta_gre_up_2100 = 8.5
          Delta_gre_lo_2100 = 2
        end if
      end if
      X_gre  = proj2order(TIME2, a1_up_g, a1_lo_g, Delta_gre_up_2100, \
               Delta_gre_lo_2100, UnifDd)
    end if
  end if

  X_gre = X_gre + 0.15  ; Add 0.15cm, the contribution from 1995 to 2005

  ; Multiply by the fingerprint
  do t=0,nb_y2-1
    X_gre(:,t) = X_gre(:,t)*F_gdyn2(t)
  end do

  ; Compute the pdfs based on the chosen periods
  do per=0,nb_PER-1        ; Loop on the period
    X_gre_pdf(per,:)  = X_gre_pdf(per,:) + \
                        pdfx(dim_avg_n(X_gre(:,i_ysep(0,per):i_ysep(1,per)),1),nbin,opt)
  end do

  ;Update X_tot, the sum of all contributions
  if COMB.eq."DEP" then
    ; Reorder contribution in ascending order
    ip   = dim_pqsort_n(X_gre, 2, 0)
    delete(ip)
  end if

  if NoU_G then
    do y=0,nb_y2-1        ; Loop on the period
      X_tot(:,y) = X_tot(:,y) + avg(X_gre(:,y))
    end do
  else
    X_tot = X_tot + X_gre
  end if

  if Corr then
    X_all(comp,:,:) = X_gre(:,ind_d)
    comp        = comp + 1
  end if
  delete(X_gre)
  delete(UnifDd)

  ;##############################################################################
  if INFO then
    print("### Compute PDF of total SLR  ########################################")
  end if

  ; Tot = Thermal exp. and ocean dyn. + Glaciers and ice sheets + Greenland SMB +
  ;       Antractic SMB + land water + antarctic dynamic + greenland dynamics

  ;Compute the pdfs based on the chosen periods
  do per=0,nb_PER-1        ; Loop on the period
    pdfi                = pdfx(dim_avg_n(X_tot(:,i_ysep(0,per):i_ysep(1,per)),1),nbin,opt)
    bin_center_f        = pdfi@bin_center
    X_tot_pdf(per,:)  = X_tot_pdf(per,:) + pdfi
    delete(pdfi)
  end do

  X_tot_pdf@bin_center = bin_center_f  
  X_tot_pdf_i = X_tot_pdf/nb_it
  PDF_cum     = dim_cumsum_n(X_tot_pdf_i(nb_PER-1,:),0,0)
  indi        = closest_val( 99.9, PDF_cum)
  ListPush(CONV,X_tot_pdf@bin_center(indi))
  print("99.9 percentile: "+CONV[0])
  delete(indi)
  delete(PDF_cum)
  delete(X_tot_pdf_i)

  if (ListCount(CONV).ge.4) then
    dc1 = abs(CONV[0]-CONV[1])
    dc2 = abs(CONV[1]-CONV[2])
    dc3 = abs(CONV[2]-CONV[3])
    if (dc1.le.er).and.(dc2.le.er).and.(dc3.le.er).and.(MIN_IT.lt.nb_it) then
      END = True
    end if
  end if

  ;Compute the correlations (X_all: nb_comp,nb_SCE,N,nb_yd)
  ; (M_Corr_x: nb_SCE,nb_yd,nb_el)
  if Corr then
    print("Computing correlations")
    do t=0,nb_yd-1
      pl = 0
      do i=0,nb_comp-1
        do j=i+1,nb_comp-1
          M_Corr_P(t,pl) = M_Corr_P(t,pl) + escorc(X_all(i,:,t),X_all(j,:,t))
          M_Corr_S(t,pl) = M_Corr_S(t,pl) + spcorr_n(X_all(i,:,t),X_all(j,:,t),0) 
          pl = pl+1
        end do
      end do
    end do
  end if

  print("Finished iteration "+nb_it)
end do ;#### End of main loop

; Scale PDFs and correlation matrices with number of iterations:
X_O_G_pdf     = X_O_G_pdf/nb_it
X_O_A_pdf     = X_O_A_pdf/nb_it
X_B_pdf       = X_B_pdf/nb_it
X_gic_pdf     = X_gic_pdf/nb_it
X_gsmb_pdf    = X_gsmb_pdf/nb_it
X_asmb_pdf    = X_asmb_pdf/nb_it
X_landw_pdf   = X_landw_pdf/nb_it
X_ant_pdf     = X_ant_pdf/nb_it
X_gre_pdf     = X_gre_pdf/nb_it
X_ant_tot_pdf = X_ant_tot_pdf/nb_it
X_tot_pdf     = X_tot_pdf/nb_it
if Corr then
  M_Corr_P      = M_Corr_P/nb_it
  M_Corr_S      = M_Corr_S/nb_it
end if

print("### Numbers for the total distribution ###")
Perc  = (/1,5,10,17,20,50,80,83,90,95,99,99.5,99.9/)
print("### Scenario "+SCE+" ###")
p     = printPerc(X_tot_pdf(nb_PER-1,:),Perc)

if Corr then
  print("### Spearman correlations ###")
  print(M_Corr_S(nb_yd-1,:))
end if

;##############################################################################
;### Output results
nb_proc = 9+2     ; 9 different processes plus total antarctic and total

MAT_OUT           = new((/nb_PER,nb_proc,nbin/),double)
MAT_OUT(:,0,:)  = X_O_G_pdf
MAT_OUT(:,1,:)  = X_O_A_pdf
MAT_OUT(:,2,:)  = X_B_pdf
MAT_OUT(:,3,:)  = X_gic_pdf
MAT_OUT(:,4,:)  = X_gsmb_pdf
MAT_OUT(:,5,:)  = X_asmb_pdf
MAT_OUT(:,6,:)  = X_landw_pdf
MAT_OUT(:,7,:)  = X_ant_pdf
MAT_OUT(:,8,:)  = X_gre_pdf
MAT_OUT(:,9,:)  = X_ant_tot_pdf
MAT_OUT(:,10,:) = X_tot_pdf

MAT_OUT!0       = "period"
MAT_OUT!1       = "proc"
MAT_OUT!2       = "bin"
MAT_OUT@long_name = "This variable contains the pdfs of each sea level contributor and the total sea level. The list of processes included here is written in variable ProcessNames"

MAT_OUT&period  = start_periods ;Should be the start of each period, careful with
                                ; periods longer than one year.
MAT_OUT&bin     = X_tot_pdf@bin_center
delete(MAT_OUT@bin_center) ; Improve the look of ncdump

if Corr then
  MAT_OUTc1       = new((/nb_yd,nb_el/),float)
  MAT_OUTc1(:,:)  = M_Corr_P
  MAT_OUTc1!0     = "time"
  MAT_OUTc1!1     = "elem"
  MAT_OUTc1f      = MAT_OUTc1(elem|:,time|:)
  MAT_OUTc1&time  = TIME2(ind_d)
  MAT_OUTc1@long_name = "Elements of the Pearson correlation matrix. First element is the correlation between component 1 and component 2. Second element is the correlation between component 1 and 3... 8th element is the correlation between component 2 and 3..."

  MAT_OUTc2       = new((/nb_yd,nb_el/),float)
  MAT_OUTc2(:,:)  = M_Corr_S
  MAT_OUTc2!0     = "time"
  MAT_OUTc2!1     = "elem"
  MAT_OUTc2f      = MAT_OUTc2(elem|:,time|:)
  MAT_OUTc2&time  = TIME2(ind_d)
  MAT_OUTc2@long_name = "Elements of the Spearman correlation matrix. See the NameComponents variable for the list of components. First element is the correlation between component 1 and component 2. Second element is the correlation between component 1 and 3... 8th element is the correlation between component 2 and 3..."
end if

print("### Export data in NetCDF file #########################################")
system("/bin/rm -f "+NameOut)  ; remove pre-existing file
ncdf = addfile(NameOut ,"c")   ; open output netCDF file
;### Create global attributes for the file
fAtt               = True            ; assign file attributes
fAtt@title         = "Computations were done with the following options:: " + \
      "Local computations?" + LOC + \
      ", include Inverse Barometer effect: "+ IBarE + \
      ", GMST option: "+ TEMPf + \
      ", Greenland SMB and dynamics is: "+GRE + \
      ", Ocean dynamics is: " + ODYN + \
      ", Antarctic dynamics is: " + ANT_DYN + \
      ", Gamma is: " + tostring(GAM)+ \
      ", combination of processes: " + COMB + \
      ", compute correlation between processes: " + Corr + \
      ", correlation between GMST and thermal expansion is: "+ CorrGT + \
      ", measure of correlation between GMST and thermal expansion is:"+ CorrM + \
      ", correlation between surface mass balances: "+ CorrSMB + \
      ", correlation between ice sheet dynamics: "+ CorrDYN + \
      ", remove ocean dynamics uncertainty: "+ NoU_O + \
      ", remove greenland uncertainty: "+ NoU_G + \
      ", remove Antarctic uncertainty: "+ NoU_A + \
      ", remove glaciers and ice caps uncertainty: "+ NoU_Gl
fAtt@source_file   = "This NetCDF file was built from the Probabilistic Sea Level Projection code v"+VER+".ncl"
fAtt@creation_date = systemfunc ("date")
fileattdef( ncdf, fAtt )            ; copy file attributes
;### Make time an UNLIMITED dimension
;filedimdef(ncdf,"Years",-1,True)
;### Output variables directly; NCL will call appropriate functions
; to write the meta data associated with each variable
ncdf->DataPDF        = MAT_OUT
CharProcessNames     = tochar(ProcessNames)
CharProcessNames!0   = "proc"
ncdf->ProcessNames   = CharProcessNames
if Corr then
  ncdf->PearsonCor     = MAT_OUTc1
  ncdf->SpearmanCor    = MAT_OUTc2
  NameComponents_OUT   = tochar(NameComponents)
  NameComponents_OUT@long_name = "Components between which a correlation coefficient is computed"
  ncdf->NameComponents = NameComponents_OUT
end if

end 
