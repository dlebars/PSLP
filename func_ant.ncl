;###############################################################################
; func_ant.ncl: Defines functions to retun probabilistic Antarctic contribution 
;               to sea level
;###############################################################################
load "./func_misc.ncl"

;### function: ant_dyn_knmi16 ##################################################
; Compute the antarctic dynamics contribution developed for KNMI16.
; Inputs:
;     Delta_ant_up_2100: Upper bound of sea level rise in 2100 compared to 2006, as estimated by 
;                        IPCC, AR5, the period 1996-2006 is added separately.
;     Delta_ant_up_KNMI_2100: KNMI14 estimation of the SLR in 2100.
;     Delta_ant_up_new_2100[:]: New estimation from DeConto, scenario dependent
;     Delta_ant_lo_2100: Lower bound
;     a1_up_a          : Upper bound of initial trend, data from 2005-2010 (IPCC AR5) 
;     a1_lo_a          : Lower bound
;     start_date2      : Start date for integration
;     end_date2        : Last date for integration
;     TIME_loc         : Time vector
;     N                : Number of samples in the distribution
;##############################################################################
undef("ant_dyn_knmi16")
function ant_dyn_knmi16(Delta_ant_up_2100[*]:numeric, Delta_ant_up_KNMI_2100[*]:numeric, \
                       Delta_ant_up_new_2100[*]:numeric,Delta_ant_lo_2100[*]:numeric,    \
                       a1_up_a[*]:numeric, a1_lo_a[*]:numeric, start_date2[*]:numeric,   \
                       end_date2[*]:numeric, TIME_loc[*]:numeric, N[*]:numeric)
local a2_up_a, a2_lo_a, Delta_ant_up, Delta_ant_lo, a2_up_a_KNMI, Delta_ant_up_KNMI,     \
      Delta_ant_lo_KNMI, a2_up_a_new, Delta_ant_up_new, Delta_ant_lo_new, Delta_ant_cen, \
      tau_ant, mu_ant, Diff, sig_ant, mu_ant3, sig_ant3, NormD2, X_ant, nb_y_loc
begin
  nb_y_loc   = dimsizes(TIME_loc)
  alpha_95   = cdfnor_x(0.95,0,1)
  alpha_98   = cdfnor_x(0.98,0,1)

  ;Compute the second order coefficient of the equations:
  a2_up_a  = (Delta_ant_up_2100 - a1_up_a*(2100-start_date2))/(2100 - start_date2)^2
  a2_lo_a  = (Delta_ant_lo_2100 - a1_lo_a*(2100-start_date2))/(2100 - start_date2)^2
  Delta_ant_up = a1_up_a*(TIME_loc-start_date2) + a2_up_a*(TIME_loc-start_date2)^2
  Delta_ant_lo = a1_lo_a*(TIME_loc-start_date2) + a2_lo_a*(TIME_loc-start_date2)^2

  a2_up_a_KNMI  = (Delta_ant_up_KNMI_2100 - a1_up_a*(2100-start_date2))/(2100 - start_date2)^2
  Delta_ant_up_KNMI = a1_up_a*(TIME_loc-start_date2) + a2_up_a_KNMI*(TIME_loc-start_date2)^2
  Delta_ant_lo_KNMI = Delta_ant_lo

  a2_up_a_new       = (Delta_ant_up_new_2100 - a1_up_a*(2100-start_date2))/(2100 - start_date2)^2
  Delta_ant_up_new  = a1_up_a*(TIME_loc-start_date2) + a2_up_a_new*(TIME_loc-start_date2)^2
  Delta_ant_lo_new  = Delta_ant_lo

  ;Build distribution that conserves the mode of the KNMI distribution
  Delta_ant_cen = (Delta_ant_up-Delta_ant_lo)/2
  tau_ant  = Delta_ant_lo
  mu_ant   = log(Delta_ant_cen - Delta_ant_lo)
  Diff     = Delta_ant_cen-Delta_ant_lo
  Diff     = where(Diff.le.0.1,0.1,Diff)    ; Avoid 0's to divide later, not important
                                          ; because the happen during the first years
  sig_ant  = 1/alpha_98*log((Delta_ant_up_KNMI-Delta_ant_lo)/Diff)

  sig_ant3  = (-alpha_95 + sqrt(alpha_95^2 - 4*(mu_ant - sig_ant^2 - log(Delta_ant_up_new))))/2
  mu_ant3   = log(Delta_ant_up_new) - sig_ant3*alpha_95

  NormD2 = random_normal(0,1,N)
  X_ant = new((/N,nb_y_loc/),float)
  do t=0,nb_y_loc-1
    X_ant(:,t) = tau_ant(t) + exp(mu_ant3(t) + sig_ant3(t)*NormD2)
  end do
  X_ant = X_ant + 0.25 ;Add 0.25 cm for the changes between 1996 to 2006 

  return X_ant
end

;### function: ant_dyn_dc16 ####################################################
; Compute the antarctic dynamics contribution to global sea level as in IPCC AR5.
; Inputs:
;
;###############################################################################
undef("ant_dyn_dc16")
function ant_dyn_dc16(SCE:string, N[*]:numeric, TempDep[*]:numeric,\
                      ANT_SCE:numeric, a1_up_a:numeric, a1_lo_a:numeric, TIME_loc[*]:numeric)
local X_ant_m, X_ant_sig, DeltaSLR, X_ant, nb_SCE, NormD_loc, s, t
begin

  nb_y_loc   = dimsizes(TIME_loc)
  ;Dimensions of X_ant_m and X_ant_sig: 
  ; 0: Emission scenarios (RCP4.5,RCP8.5)
  ; 1: Sea level scenario from DeConto 2016, see main program
  X_ant_m    = (/(/49,58,26,41/),(/105,114,64,79/)/) ; Numbers from DeConto 2016.
  X_ant_sig  = (/(/20,28,28,30/),(/30,36,49,46/)/)
  NormD_loc  = random_normal(0,1,N)
  DeltaSLR   = X_ant_m(1,ANT_SCE) - X_ant_m(0,ANT_SCE) ; Difference between the two SLR scenarios
  X_anti     = (a1_up_a+a1_lo_a)/2 + NormD_loc*(a1_up_a-a1_lo_a)/2

  X_antf     = new(N,float)
  X_ant      = new((/N,nb_y_loc/),float)
  alp        = 5        ; Use a fifth order polynomial
  if SCE.eq."rcp45" then
    s = 0
    else if SCE.eq."rcp85"
    s = 1
    end if
  end if
  X_antf = X_ant_m(s,ANT_SCE) + NormD_loc*X_ant_sig(s,ANT_SCE) + DeltaSLR*TempDep(:)
  do t=0,nb_y_loc-1
    X_ant(:,t) =  X_anti*(TIME_loc(t)-TIME_loc(0)) + \
    ((X_antf(:) - X_anti*(2100-TIME_loc(0)))/(2100-TIME_loc(0))^alp)*(TIME_loc(t)-TIME_loc(0))^alp
  end do

  return X_ant
end

;### function: ant_dyn_lev  ####################################################
; Compute the antarctic dynamics contribution to global sea level as in 
; Levermann et al 2014, using linear response theory.
; Inputs:
;
;###############################################################################
undef("ant_dyn_lev14")
function ant_dyn_lev14(ROOTL:string, SCE:string, MOD:string, start_date2[*]:numeric, \
                      GAM:numeric, NormD:numeric, UnifDd:numeric, DIR_T:string)
local f, nb_MOD, ye, start_date, nb_y, TIME, i_ys, i_ysr_Lev, i_yer_Lev, \
      Beta_low, Beta_high,   \
      RF_AIF, RF_PS, RF_PISM, \
      RF_SICO, RF_UMISM, coeff, RF_dim, nb_bass, Ralli, Rall, TGLOB, Tref_Lev, file_n, \
      TEMP, time, dim_t, i_start, i_end, TGLOB_Lev, TGLOB_Lev_m, TGLOB_Lev_sd,  \
      Td_Lev, RMod, AlpCoeff, Beta, BMelt, modelsel, Rmat, X_ant_b, Rdist, ref
begin
  f          = addfile(ROOTL+"/RFunctions.nc","r")
  nb_MOD     = dimsizes(MOD)
  N          = dimsizes(NormD)
  ye         = 2100
  start_date = 1861   ; This is different from other runs
  Beta_low   = 7      ; Bounds to use for the basal melt rate,
  Beta_high  = 16     ; units are m.a^(-1).K^(-1)
  nb_y       = ye-start_date+1
  TIME       = ispan(start_date,ye,1)
  i_ys       = closest_val(start_date2,TIME)
  i_ysr_Lev  = closest_val(1861,TIME)  ; Reference time 1860 to 1880
  i_yer_Lev  = closest_val(1880,TIME)

  RF_AIF   = f->RF_AIF
  RF_PS    = f->RF_PS
  RF_PISM  = f->RF_PISM
  RF_SICO  = f->RF_SICO
  RF_UMISM = f->RF_UMISM
  coeff    = tofloat(f->coeff)

  RF_dim   = dimsizes(RF_AIF)
  nb_bass  = RF_dim(1)

  Ralli     = new((/5,RF_dim(0),RF_dim(1)/),float)
  Ralli(0,:,:)  = tofloat(RF_PS)
  Ralli(1,:,:)  = tofloat(RF_PISM)
  Ralli(2,:,:)  = tofloat(RF_SICO)
  Ralli(3,:,:)  = tofloat(RF_AIF)
  Ralli(4,:,:)  = tofloat(RF_UMISM)

  Ralli!0   = "model"
  Ralli!1   = "time"
  Ralli!2   = "bass"
  Rall      = Ralli(bass |:,model |:,time |:)

  ;Read Tglob, do it locally because earlier dates are necessary
  TGLOB    = new((/nb_MOD,nb_y/),float)
  Tref_Lev = new((/nb_MOD/),float)

  do m=0,nb_MOD-1
    if MOD(m).eq."BCC-CSM1-1" then
      loc_mod = "bcc-csm1-1"
      else
      loc_mod = MOD(m)
    end if
    file_n   = "global_tas_Amon_"+loc_mod+"_"+SCE+"_r1i1p1.dat"
    TEMP     = readAsciiTable(DIR_T+file_n,13,"float",3)
    time     = TEMP(:,0)
    dim_t    = dimsizes(time)
    i_start  = closest_val(start_date,time)
    i_end    = closest_val(ye,time)
    TGLOB(m,:)   = dim_avg_n(TEMP(i_start:i_end,1:12),1)   ;Data in degree Kelvin
    Tref_Lev(m)  = avg(TGLOB(m,i_ysr_Lev:i_yer_Lev))
    delete(TEMP)
    delete(time)
  end do

  ;Build the distribution of global temperature for this process
  Td_Lev = TempDist(TGLOB, Tref_Lev, GAM, NormD)

  ;Random model number: 1-19
  RMod = toint(floor(random_uniform(0,19,N))) ; Select random model indice (0 to 18)
  AlpCoeff = coeff(:,RMod,0) ; dim: bassin, N
  delete(RMod)
  ;# Use following line if Beta should have some external dependence
;  Beta     = Beta_low + UnifDd*(Beta_high - Beta_low)  ; Modify to obtain random_uniform(7,16,N)
  Beta     = random_uniform(Beta_low,Beta_high,N)

  BMelt     = new((/nb_bass,N,nb_y/),float)
  do b=0,nb_bass-1
    do t=0,nb_y-1
      BMelt(b,:,t) = AlpCoeff(b,:)*Td_Lev(:,t)*Beta
    end do
  end do
  delete(AlpCoeff)
  delete(Td_Lev)
  delete(Beta)

;  UnifDd2  = random_uniform(0,1,N)   
  Rdist    = new(N,integer)
  Rdist    = 2
  Rdist    = where(UnifDd.gt.0.33,1,Rdist)
  Rdist    = where(UnifDd.gt.0.67,0,Rdist)
  modelsel = Rdist ; Select model
;  modelsel = toint(floor(random_uniform(0,3,N))) ; Select model
  ;Rq: Select 0,1,2 for shelf models and 3,4 for other models
  Rmat     = Rall(:,modelsel,::-1) ; Reverse time to speed up loop below
  delete(modelsel)
  X_ant_b    = new((/nb_bass,N,nb_y/),float)
  do t=0,nb_y-1
    X_ant_b(:,:,t) = dim_sum_n(Rmat(:,:,RF_dim(0)-1-t:)*BMelt(:,:,:t),2)
  end do
  delete(Rmat)
  delete(BMelt)

  X_ant_b     = X_ant_b*100                   ; Convert from m to cm
  X_ant       = dim_sum_n(X_ant_b,0)          ; Sum 4 bassins

  ;Remove the uncertainty at the beginning of the projection
  ref = X_ant(:,i_ys-1)
  do t=0,nb_y-1
    X_ant(:,t) = X_ant(:,t) - ref
  end do

  return X_ant(:,i_ys:)

end

