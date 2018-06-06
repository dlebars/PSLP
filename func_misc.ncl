;###############################################################################
; func_misc.ncl: Miscellaneous functions
;###############################################################################

;### function: finger1D ########################################################
; Inputs: 
;      lat_in: latitude of location of interest
;      lon_in: longitude of location of interest
;      lat[:]: 1D array of latitude of fingerprint
;      lon[:]: 1D array of longitude of fingerprint
;      fingerprint[time (or not), lat, lon]: in %, missing value 0
; Outputs:
;      output[time or not]: local fingerprint (%)
;###############################################################################
undef("finger1D")
function finger1D(lat_in[*]:numeric, lon_in[*]:numeric , lat1D[*]:numeric , lon1D[*]:numeric, \
         fingerprint:numeric)
local dim_f, mask2D, lat2D, lon2D, loc, output
begin
dim_f = dimsizes(fingerprint)
ndim  = dimsizes(dim_f)
if ndim.eq.2 then
  output = new(1,float)
  mask2D = where(fingerprint.eq.fingerprint@_FillValue,0,1)
  else if ndim.eq.3 then
    output = new(dim_f(0),float)
    mask2D = where(fingerprint(1,:,:).eq.fingerprint@_FillValue,0,1)
  end if
end if

lat2D = new((/dim_f(ndim-2),dim_f(ndim-1)/),float)
lon2D = new((/dim_f(ndim-2),dim_f(ndim-1)/),float)
do i=0,dim_f(ndim-1)-1
  lat2D(:,i) = lat1D
end do
do i=0,dim_f(ndim-2)-1
  lon2D(i,:) = lon1D
end do

lat2D = lat2D*mask2D
lat2D@_FillValue = 0
lon2D = lon2D*mask2D
lon2D@_FillValue = 0

loc = getind_latlon2d(lat2D,lon2D,lat_in,lon_in)
;print(loc)

if ndim.eq.2 then
  output = fingerprint(loc(0,0),loc(0,1))
  else if ndim.eq.3 then
    output = fingerprint(:,loc(0,0),loc(0,1))
  end if
end if

return output

end

;### function: printPerc #######################################################
; Compute percentiles from a PDF and print.
; Inputs:
;     PDF        : A pdf to compute
;     Percentiles: The percentiles to compute
;###############################################################################
undef("printPerc")
function printPerc(InPDF, Perc)
local PDF_cum, dimP, i, indi
begin 
  PDF_cum = dim_cumsum_n(InPDF,0,0)
  dimP    = dimsizes(Perc)
  do i=0, dimP-1
    print("Percentile: "+Perc(i))
    indi = closest_val( Perc(i), PDF_cum)
    print(InPDF@bin_center(indi))
  end do
  return 1
end

;### function: proj2order ####################################################
; Project future values of sea level using present day uncertanty range of the 
; contribution in cm/year and uncertainty of total contribution in 2100 in cm.
; The uncertainty is represented by a uniform distribution.
;#############################################################################
undef("proj2order")
function proj2order(TIME_loc[*]:numeric, a1_up:numeric, a1_lo:numeric,  \
   Delta_up_2100:numeric, Delta_lo_2100:numeric, Unif[*]:numeric)
local nb_y_loc, Delta_up, Delta_lo, N
begin

  nb_y_loc   = dimsizes(TIME_loc)
  N          = dimsizes(Unif)

  ;Compute the second order coefficient of the equations:
  a2_up  = (Delta_up_2100 - a1_up*(2100-TIME_loc(0)))/(2100 - TIME_loc(0))^2
  a2_lo  = (Delta_lo_2100 - a1_lo*(2100-TIME_loc(0)))/(2100 - TIME_loc(0))^2

  Delta_up = a1_up*(TIME_loc-TIME_loc(0)) + a2_up*(TIME_loc-TIME_loc(0))^2
  Delta_lo = a1_lo*(TIME_loc-TIME_loc(0)) + a2_lo*(TIME_loc-TIME_loc(0))^2

  X_out = new((/N,nb_y_loc/),float)  ; Independent of the scenario, size is to add up easily.
  do t=0,nb_y_loc-1
    X_out(:,t) = Unif*Delta_up(t) + (1-Unif)*Delta_lo(t)
  end do

  return X_out
end

;### function: TempDist ########################################################
; Build a distribution of global temperature for a contributor (reference periods 
; are different of each contributors)
;###############################################################################
undef("TempDist")
function TempDist(TGLOBs[*][*]:numeric, Tref[*]:numeric, GAM:numeric,  \
  NormD[*]:numeric)
local N, dimTGLOB, nb_MOD, nb_y2, TGLOBl, m, Td
begin

  N        = dimsizes(NormD)
  dimTGLOB = dimsizes(TGLOBs)
  nb_MOD   = dimTGLOB(0)
  nb_y2    = dimTGLOB(1)

  TGLOBl   = new((/nb_MOD,nb_y2/),float)
  do m=0,nb_MOD-1
    TGLOBl(m,:)    = TGLOBs(m,:) - Tref(m)
  end do
  TGLOB_m  = dim_avg_n(TGLOBl,0)    ; Compute the inter-model mean for each time
  TGLOB_sd = dim_stddev_n(TGLOBl,0) ; Compute the inter-model standard deviation

  Td       = new((/N,nb_y2/),float)
  do t=0,nb_y2-1
    Td(:,t)  = TGLOB_m(t) + GAM*NormD(:)*TGLOB_sd(t)
  end do

  return Td
end

