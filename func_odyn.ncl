;###############################################################################
; func_odyn.ncl: Defines functions returning ocean dynamics probabilistic 
;                contribution to sea level
;###############################################################################

;### function: odyn_loc ########################################################
; Compute the ocean dynamics and thermal expansion contribution to local sea level.
; Inputs:
;
;###############################################################################
undef("odyn_loc")
function odyn_loc(SCE:string, MOD:string, nb_y[*]:numeric, nb_y2[*]:numeric, \
        DIR_O:string, lat_N[*]:numeric, lat_S[*]:numeric, lon_W[*]:numeric,  \
        lon_E[*]:numeric, start_date[*]:numeric, end_date2[*]:numeric,       \
        VAR:string, N[*]:numeric, i_ys[*]:numeric, Gam[*]:numeric, NormD:numeric)
local nb_MOD, MAT, MAT_G, MAT_A, fi, lon, lat, lat_Ni, lat_Si,   \
      lon_Wi, lon_Ei, TIMEt, TIMEt2, i_start, i_end, SSH, nb_y_loop, MATs,   \
      MAT_Gs, MAT_As, X_O_m, X_O_sd, X_O_G_m, X_O_G_sd, X_O_A_m, X_O_A_sd,   \
      X_O, X_O_G, X_O_A, X_O_out, s, m, t
begin
  nb_MOD = dimsizes(MOD)

  ;Initialize the SSH matrix: (scenario, model, time (years))
  MAT   = new((/nb_MOD,nb_y/),float)
  MAT_G = new((/nb_MOD,nb_y/),float) ; Global mean steric effect
  MAT_A = new((/nb_MOD,nb_y/),float) ; Local dynamics

  do m=0,nb_MOD-1
    fi      = addfile(DIR_O+MOD(m)+"_"+SCE+".nc","r")
    lon     = fi->longitude  ; Do not assume they all have the same grid, even though they DO
    lat     = fi->latitude
    lat_Ni  = closest_val(lat_N,lat)
    lat_Si  = closest_val(lat_S,lat)
    lon_Wi  = closest_val(lon_W,lon)
    lon_Ei  = closest_val(lon_E,lon)
    TIMEt   = fi->TIME
    TIMEt2  = cd_calendar(TIMEt,0)
    i_start = closest_val(start_date,TIMEt2(:,0))
    i_end   = closest_val(end_date2,TIMEt2(:,0))
    SSH     = fi->$VAR$(i_start:i_end,lat_Si:lat_Ni,lon_Wi:lon_Ei)
    nb_y_loop = i_end - i_start +1
    if nb_y_loop.eq.nb_y then
      MAT(m,:) = dim_avg_n(dim_avg_n(SSH,2),1)
             ;RQ: No scaling depending on the area, gives more weight to the southern cell
      MAT_G(m,:) = fi->$VAR$(i_start:i_end,0,0)
      else
      MAT(m,:nb_y-2)   = dim_avg_n(dim_avg_n(SSH,2),1)
      MAT(m,nb_y-1)    = MAT(m,nb_y-2)
      MAT_G(m,:nb_y-2) = fi->$VAR$(i_start:i_end,0,0)
      MAT_G(m,nb_y-1)  = MAT_G(m,nb_y-2)
    end if
    ; Remove the average SSH of the first 20 years from all models
    MAT(m,:) =  MAT(m,:) - avg(MAT(m,:20))
    ;###
    MAT_G(m,:) = MAT_G(m,:) - avg(MAT_G(m,:20))
    MAT_A(m,:) = MAT(m,:) - MAT_G(m,:)
    delete(fi)
    delete(lon)
    delete(lat)
    delete(SSH)
    delete(TIMEt)
    delete(TIMEt2)
  end do

  MATs     = MAT(:,i_ys:)*100   ; Convert from m to cm
  MAT_Gs   = MAT_G(:,i_ys:)*100 ; Convert from m to cm
  MAT_As   = MAT_A(:,i_ys:)*100 ; Convert from m to cm
  ;Build the distribution
  X_O_m    = dim_avg_n(MATs,0)    ; Compute the inter-model mean for each time
  X_O_sd   = dim_stddev_n(MATs,0) ; Compute the inter-model standard deviation
  X_O_G_m  = dim_avg_n(MAT_Gs,0)
  X_O_G_sd = dim_stddev_n(MAT_Gs,0)
  X_O_A_m  = dim_avg_n(MAT_As,0)
  X_O_A_sd = dim_stddev_n(MAT_As,0)

  X_O    = new((/N,nb_y2/),float)
  X_O_G  = new((/N,nb_y2/),float)
  X_O_A  = new((/N,nb_y2/),float)
    do t=0,nb_y2-1
      X_O(:,t)    = X_O_m(t) + Gam*NormD(:)*X_O_sd(t)
      X_O_G(:,t)  = X_O_G_m(t) + Gam*NormD(:)*X_O_G_sd(t)
      X_O_A(:,t)  = X_O_A_m(t) + Gam*NormD(:)*X_O_A_sd(t)
    end do

  X_O_out = new((/3,N,nb_y2/),float)
  X_O_out(0,:,:) = X_O
  X_O_out(1,:,:) = X_O_G
  X_O_out(2,:,:) = X_O_A

  return X_O_out
end

;### function: odyn_glob_knmi ##################################################
; Compute thermal expansion contribution to global sea level from KNMI data.
; Inputs:
;
;###############################################################################
undef("odyn_glob_knmi")
function odyn_glob_knmi(SCE:string, MOD:string, nb_y[*]:numeric, nb_y2[*]:numeric, \
      DIR_O:string, DIR_OG:string, start_date[*]:numeric, end_date2[*]:numeric, \
      VAR:string, N[*]:numeric, i_ys[*]:numeric, Gam[*]:numeric, NormD[*]:numeric)
local nb_MOD, MAT, fi, fig, TIMEt, TIMEt2, i_start, i_end, nb_y_loop, \
      MATs, X_O_m, X_O_sd, X_O, X_O_out, s, m, t
begin
  nb_MOD = dimsizes(MOD)

  ;Initialize the SSH matrix: (scenario, model, time (years))
  MAT   = new((/nb_MOD,nb_y/),float)
  print("WARNING !!!!!!! There seem to be an mistake in this script, variable fig"+ \
        " is not used, should be used instead of fi in loop?")

  do m=0,nb_MOD-1
    ;###
    fi      = addfile(DIR_O+MOD(m)+"_"+SCE+".nc","r")
    TIMEt   = fi->TIME
    TIMEt2  = cd_calendar(TIMEt,0)
    i_start = closest_val(start_date,TIMEt2(:,0))
    i_end   = closest_val(end_date2,TIMEt2(:,0))
    fig          = addfile(DIR_OG+MOD(m)+"_"+SCE+".nc","r")
    nb_y_loop = i_end - i_start +1
    if nb_y_loop.eq.nb_y then
      MAT(m,:) = fi->$VAR$(i_start:i_end,0,0)
      else
      MAT(m,:nb_y-2) = fi->$VAR$(i_start:i_end,0,0)
      MAT(m,nb_y-1)  = MAT(m,nb_y-2)
    end if
    MAT(m,:) = MAT(m,:) - avg(MAT(m,:20))
    delete(fi)
    delete(fig)
    delete(TIMEt)
    delete(TIMEt2)
  end do

  MATs     = MAT(:,i_ys:)*100   ; Convert from m to cm, and select dates after 2006
  ;Build the distribution
  X_O_m    = dim_avg_n(MATs,0)    ; Compute the inter-model mean for each time
  X_O_sd   = dim_stddev_n(MATs,0) ; Compute the inter-model standard deviation

  X_O    = new((/N,nb_y2/),float)
  do t=0,nb_y2-1
    X_O(s,:,t)    = X_O_m(s,t) + Gam*NormD(:)*X_O_sd(s,t)
  end do

  X_O_out = new((/3,N,nb_y2/),float)
  X_O_out(0,:,:) = X_O
  X_O_out(1,:,:) = X_O   ; In this case global is the same as total
  X_O_out(2,:,:) = 0     ; and anomaly is 0

  return X_O_out

end

;### function: odyn_glob_ipcc ##################################################
; Compute thermal expansion contribution to global sea level from IPCC data.
; Inputs:
;
;###############################################################################
undef("odyn_glob_ipcc")
function odyn_glob_ipcc(SCE:string ,DIR_IPCC:string, N[*]:numeric, \
                        nb_y2[*]:numeric, Gam[*]:numeric, NormD[*]:numeric)
local X_O_med, X_O_up, f_med, f_up, X_O, X_O_sd, X_O_out
begin
  X_O_med   = new(nb_y2-1,float) ; These start in 2007 instead of 2006
  X_O_up    = new(nb_y2-1,float)

  f_med     = addfile(DIR_IPCC+SCE+"_expansionmid.nc","r")
  f_up      = addfile(DIR_IPCC+SCE+"_expansionupper.nc","r")
  X_O_med   = (f_med->global_average_sea_level_change)*100
  X_O_up    = (f_up->global_average_sea_level_change)*100

  X_O    = new((/N,nb_y2/),float)
  X_O_sd = (X_O_up-X_O_med)/cdfnor_x(0.95,0,1)  ; ~1.64 

  do t=1,nb_y2-1
    X_O(:,t)    = X_O_med(t-1) + Gam*NormD(:)*X_O_sd(t-1)
  end do
  X_O(:,0)      = X_O(:,1)

  X_O_out          = new((/3,N,nb_y2/),float)
  X_O_out(0,:,:) = X_O
  X_O_out(1,:,:) = X_O   ; In this case global is the same as total
  X_O_out(2,:,:) = 0     ; and anomaly is 0
  return X_O_out

end

;### function: read_odyn_cmip5     #############################################
; Read the CMIP5 data and output time series of average and standard devation,
; to be used by odyn_cmip5 that produces the distributions and ensures continuity.
;###############################################################################
undef("read_odyn_cmip5")
function read_odyn_cmip5(DIR_OCMIP5:string, SCE:string, ys[*]:numeric, ye[*]:numeric, \
                         LOC:logical)
local fzos, fzostoga, fzossga, fzos_avg, time, iys, iye, dimt, zos, zostoga, \
      zossga, zos_avg, dimzos, OUT, zos_ModelNames, zostoga_ModelNames, zossga_ModelNames, \
      zos_avg_ModelNames, tot_sl, glob_sl, ind1, ind2, ind3
begin
    fzos     = addfile(DIR_OCMIP5+"CMIP5_SeaLevel_"+SCE+"_zos_1986-"+ye+".nc","r")
    fzostoga = addfile(DIR_OCMIP5+"CMIP5_SeaLevel_"+SCE+"_zostoga_1986-"+ye+".nc","r")
    fzossga  = addfile(DIR_OCMIP5+"CMIP5_SeaLevel_"+SCE+"_zossga_1986-"+ye+".nc","r")
    fzos_avg = addfile(DIR_OCMIP5+"CMIP5_SeaLevel_"+SCE+"_zos_avg_1986-"+ye+".nc","r")

    time     = fzos->time
    iys      = closest_val(ys,time)
    iye      = closest_val(ye,time)
    dimt     = dimsizes(time)

    zos      = fzos->LocalSeaLevel(:,iys:iye)      ;!!! Change name
    zostoga  = fzostoga->GlobalSeaLevel(:,iys:iye)
    zossga   = fzossga->GlobalSeaLevel(:,iys:iye)
    zos_avg  = fzos_avg->AverageSeaLevel(:,iys:iye)
    dimzos   = dimsizes(zos)
    OUT      = new((/4,dimzos(1)/),float)

    zos_ModelNames     = tostring(fzos->ModelNames)
    zostoga_ModelNames = tostring(fzostoga->ModelNames)
    zossga_ModelNames  = tostring(fzossga->ModelNames)
    zos_avg_ModelNames = tostring(fzos_avg->ModelNames)

    tot_sl  = zos
    glob_sl = zos
    glob_sl = 0
    if .not.LOC then
      tot_sl = 0 ; Keep the same dimensions but remove local sea level
    end if
    do m=0,dimzos(0)-1
      ind1 = ind(zostoga_ModelNames.eq.zos_ModelNames(m))
      ind2 = ind(zossga_ModelNames.eq.zos_ModelNames(m))
      ind3 = ind(zos_avg_ModelNames.eq.zos_ModelNames(m))
      if .not.ismissing(ind1)
        tot_sl(m,:)  = tot_sl(m,:) + zostoga(ind1,:)
        glob_sl(m,:) = glob_sl(m,:) + zostoga(ind1,:)
        else if .not.ismissing(ind2)
          tot_sl(m,:)  = tot_sl(m,:) + zossga(ind2,:)
          glob_sl(m,:) = glob_sl(m,:) + zossga(ind2,:)
          else
            tot_sl(m,:)  = tot_sl@_FillValue
            glob_sl(m,:) = tot_sl@_FillValue
        end if
      end if
      if ismissing(ind3)
;      print("Missing zos_avg for model "+zos_ModelNames(m))
        else
        tot_sl(m,:) = tot_sl(m,:) - zos_avg(m,:)
      end if
    end do
    OUT(0,:) = tofloat(dim_avg_n(tot_sl,0))
    OUT(1,:) = tofloat(dim_stddev_n(tot_sl,0))
    OUT(2,:) = tofloat(dim_avg_n(glob_sl,0))
    OUT(3,:) = tofloat(dim_stddev_n(glob_sl,0))

    return OUT

end

;### function: odyn_cmip5     ##################################################
; Compute thermal expansion contribution to global sea level from IPCC data.
; Output should be needed from 2006 to 2100
; Fields are already computed using PlotThermalExp_loc.ncl and PlotThermalExp.ncl
; The global average zos is removed.
; Warning!!! There is no check of the area for zos, the files should be the right 
; ones in the DIR_OCMIP5 directory.
; TS: 0: total local sea level ensemble avg
;     1: total local sea level ensemble standard deviation
;     2: global sea level ensemble avg
;     3: global sea level ensemble standard deviation
;###############################################################################
undef("odyn_cmip5")
function odyn_cmip5(SCE:string, LOC:logical, DIR_OCMIP5:string, N[*]:numeric,  \
                    ys[*]:numeric, ye[*]:numeric, Gam[*]:numeric, NormD[*]:numeric)
local nb_y2, X_O_out, tot_sl_avg, tot_sl_std, glob_sl_avg, glob_sl_std, X_O_out,
      TS2100, TS2300
begin
  nb_y2     = ye - ys +1
  X_O_out   = new((/3,N,nb_y2/),float)

  tot_sl_avg  = new(nb_y2,float)
  tot_sl_std  = new(nb_y2,float)
  glob_sl_avg = new(nb_y2,float)
  glob_sl_std = new(nb_y2,float)
  TS2100    = read_odyn_cmip5(DIR_OCMIP5, SCE, ys, 2100, LOC)
  dim2100   = dimsizes(TS2100)
  if ye.eq.2100 then
    tot_sl_avg  = TS2100(0,:)
    tot_sl_std  = TS2100(1,:)
    glob_sl_avg = TS2100(2,:)
    glob_sl_std = TS2100(3,:)
    else if ye.gt.2100 then
    TS2300           = read_odyn_cmip5(DIR_OCMIP5, SCE, ys, 2300, LOC)
    tot_sl_avg_2300  = TS2300(0,:)
    tot_sl_std_2300  = TS2300(1,:)
    glob_sl_avg_2300 = TS2300(2,:)
    glob_sl_std_2300 = TS2300(3,:)

    tot_sl_avg(:dim2100(1)-1)  = TS2100(0,:)
    tot_sl_std(:dim2100(1)-1)  = TS2100(1,:)
    glob_sl_avg(:dim2100(1)-1) = TS2100(2,:)
    glob_sl_std(:dim2100(1)-1) = TS2100(3,:)
    tot_sl_avg(dim2100(1)-1:)  = TS2300(0,dim2100(1)-1:nb_y2-1) - TS2300(0,dim2100(1)-1) + \
                                 tot_sl_avg(dim2100(1)-1)
    tot_sl_std(dim2100(1)-1:)  = TS2300(1,dim2100(1)-1:nb_y2-1) - TS2300(1,dim2100(1)-1) + \
                                 tot_sl_std(dim2100(1)-1)
    glob_sl_avg(dim2100(1)-1:) = TS2300(2,dim2100(1)-1:nb_y2-1) - TS2300(2,dim2100(1)-1) + \
                                 glob_sl_avg(dim2100(1)-1)
    glob_sl_std(dim2100(1)-1:) = TS2300(3,dim2100(1)-1:nb_y2-1) - TS2300(3,dim2100(1)-1) + \
                                 glob_sl_std(dim2100(1)-1)

    else
      print("Error: Value of "+ye+" unsupported")
      exit()
  end if
  end if

  if LOC then
    do t=0,nb_y2-1
      X_O_out(0,:,t) = tofloat(tot_sl_avg(t) + Gam*NormD(:)*tot_sl_std(t))
      X_O_out(1,:,t) = tofloat(glob_sl_avg(t) + Gam*NormD(:)*glob_sl_std(t))
    end do
    X_O_out(2,:,:)   = X_O_out(0,:,:) - X_O_out(1,:,:)
  else
    do t=0,nb_y2-1
      X_O_out(0,:,t) = tofloat(glob_sl_avg(t) + Gam*NormD(:)*glob_sl_std(t))
    end do
    X_O_out(1,:,:)   = X_O_out(0,:,:)
    X_O_out(2,:,:)   = 0
  end if
  delete(tot_sl_avg)
  delete(tot_sl_std)
  delete(glob_sl_avg)
  delete(glob_sl_std)

  X_O_out = X_O_out*100 ;Convert from meters to cm

  return X_O_out

end
