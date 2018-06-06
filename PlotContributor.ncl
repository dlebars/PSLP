;###############################################################################
; PlotContributor.ncl: Test individual contributor functions
;###############################################################################
load "./func_ant.ncl"

begin
SCE = "rcp85"
MOD = (/"ACCESS1-0","BCC-CSM1-1","CanESM2","CNRM-CM5","CSIRO-Mk3-6-0","GFDL-ESM2G", \
        "GFDL-ESM2M","GISS-E2-R","HadGEM2-CC","HadGEM2-ES","inmcm4","IPSL-CM5A-LR", \
        "IPSL-CM5A-MR","MIROC5","MIROC-ESM-CHEM","MIROC-ESM","MPI-ESM-LR","MPI-ESM-MR", \
        "MRI-CGCM3","NorESM1-ME","NorESM1-M"/)
ys     = 2006
ye     = 2100   ; End year for computation
TIME2  = ispan(ys,ye,1)
start_periods = ispan(2006,2100,1) ; start dates of periods
end_periods   = ispan(2006,2100,1) ; end dates of periods
period        = (/start_periods, \
                  end_periods /)
nb_PER = dimsizes(period(0,:))

i_ysep     = new((/2,nb_PER/),integer)      ; Start and end periods index
do i=0,nb_PER-1                             ; Loop periods
  do j=0,1                                  ; Start and end years (included)
    i_ysep(j,i)    = closest_val(period(j,i),TIME2)
  end do
end do

N      = tointeger(1e3)
GAM    = 1.64
ROOT   = "/nobackup/users/bars/SeaLevelFromHylke/CMIP5_OCEAN/"
DIR_T  = ROOT+"Tglobal/"
NormD  = random_normal(0,1,N)
UnifDd = random_uniform(0,1,N)

;### Parameters to produce PDF
;Range to produce PDF of the distributions
opt           = True
opt@bin_min   = -20.5
opt@bin_max   = 500.5     ;260.5
nbin          = toint(opt@bin_max-opt@bin_min)

X_ant  = ant_dyn_lev14(SCE, MOD, ys, GAM, NormD, UnifDd, DIR_T) ; dim: N, nb_y

;Compute the pdfs based on user defined periods of time
X_ant_pdf     = new((/nb_PER,nbin/),double)
do per=0,nb_PER-1        ; Loop on the period
  X_ant_pdf(per,:)  = pdfx(dim_avg_n(X_ant(:,i_ysep(0,per):i_ysep(1,per)),1),nbin,opt)
end do

Perc     = (/5,50,95/)
PDF_cum  = dim_cumsum_n(X_ant_pdf,0,1)
dimP     = dimsizes(Perc)
PercPDF  = new((/dimP, nb_PER/),integer)
bin1     = X_ant_pdf@bin_center

; Remove the 0s to fix an issue arising from closest_val returning the first closest
; value in the array.
PDF_cum  = where(PDF_cum.eq.0,PDF_cum@_FillValue,PDF_cum)

do t=0,nb_PER-1
  do i=0, dimP-1
    indi          = closest_val(Perc(i),PDF_cum(t,:))
    PercPDF(i, t) = toint(bin1(indi))
  end do
end do

;########### Plot #############################################################
  wks                      = gsn_open_wks("pdf" ,"PlotContributor")
  res                      = True
  res@gsnMaximize          = True
  res@xyMarkLineModes      = (/"Lines","Lines","Lines","Lines"/)
  res@pmLegendDisplayMode  = "Always"
  res@tiYAxisString        = "PDF (%)"
  res@xyLineThicknessF     = 5
  res@gsnPaperOrientation  = "portrait"
  res@xyLineColors         = (/"red","blue","green","black","purple","orange","grey"/)


  res@tiMainString   = "Plot "+ SCE + ", perc: 5-50-95, Antartica Levermann 2014"
  res@xyExplicitLegendLabels = (/"5th perc.", "50th perc.", "95th perc."/)
  plot               = gsn_csm_xy(wks,start_periods,PercPDF, res)

end
