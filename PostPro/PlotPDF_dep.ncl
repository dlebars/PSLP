;##############################################################################
; PlotPDF_dep.ncl: Read NetCDF files from the Sea Level Projection Model and
;                  make plots and tables used in Le Bars 2018.
;##############################################################################
load "./func_postpro.ncl"
begin

DIR           = "../Exp_PaperDependence/"
SCE           = "rcp45"
;Read data
;PDF format: period, process, bin
; ProcessNames = (/"Global steric", "Local ocean", "Inverse barometer", "Glaciers",    \
;                 "Greenland SMB", "Antarctic SMB", "Landwater", "Antarctic dynamics",\
;                 "Greenland dynamics", "sum anta.", "Total"/)

;files    = systemfunc("ls "+DIR+"*"+SCE+"*.nc")
files = (/"IPCC1","IPCC2","IPCC3","Prob1","Prob2","Prob3","Prob4",\
           "Prob5","Prob6","Prob7","Prob8","Prob9"/)
print("Reading file: "+files)
f1       = addfile(DIR+"/SeaLevelPDF_namelist"+files(0)+"_"+SCE+".nc","r")
PDF1     = f1->DataPDF
Pcc1     = f1->PearsonCor
Scc1     = f1->SpearmanCor
bin1     = f1->bin
Years    = PDF1&period
ELEM     = tostring(f1->NameComponents)
dimPDF   = dimsizes(PDF1)
dimPcc   = dimsizes(Pcc1)
TimeCC   = f1->time

nb_cases       = dimsizes(files)
PDF            = new((/nb_cases,dimPDF(0),dimPDF(1),dimPDF(2)/),double)
PDF(0,:,:,:) = PDF1
PCC            = new((/nb_cases,dimPcc(0),dimPcc(1)/),double) ; Pearson cross-correlation
PCC(0,:,:)   = Pcc1
SCC            = new((/nb_cases,dimPcc(0),dimPcc(1)/),double) ; Spearman cross-correlation 
SCC(0,:,:)   = Scc1

do i=1,nb_cases-1
  f       = addfile(DIR+"/SeaLevelPDF_namelist"+files(i)+"_"+SCE+".nc","r")
  PDF(i,:,:,:) = f->DataPDF
  PCC(i,:,:)   = f->PearsonCor
  SCC(i,:,:)   = f->SpearmanCor
  delete(f)
end do

PCC!0 = "cases"
PCCr = PCC(cases | :, elem | :, time | :)
SCC!0 = "cases"
SCCr = SCC(cases | :, elem | :, time | :)

printVarSummary(PCC)

RUN      = (/"IPCC AR5 partial","IPCC AR5 ind.", "IPCC AR5 dep.",                     \
             "LEV14 partial, Gam=1.64", "LEV14 ind, Gam=1.64", "LEV14 dep, Gam=1.64", \ 
             "LEV14 partial, Gam=1.64, Low", "LEV14 partial, Gam=1.64, High",         \
             "LEV14 partial, Gam=1.64, NoU_O", "LEV14 partial, Gam=1.64, NoU_G",      \
             "LEV14 partial, Gam=1.64, NoU_A", "LEV14 partial, Gam=1.64, NoU_Gl"/)

nb_proc  = dimPDF(1)

per2100  = closest_val(2100,Years)
print(Years(per2100))

Perc     = (/5,10,17,50,83,90,95,99,99.9/)
Perct    = ispan(1,99,1)
dimP     = dimsizes(Perc)
dimPt    = dimsizes(Perct)
PDF_cum  = dim_cumsum_n(PDF,0,3)
; Remove the 0s to fix an issue arising from closest_val returning the first closest
; value in the array.
PDF_cum  = where(PDF_cum.eq.0,PDF_cum@_FillValue,PDF_cum)

PercPDF  = new((/nb_cases,dimPDF(0),dimP/),integer)
PerctPDF = new((/nb_cases,dimPDF(0),dimPt/),integer)

pro  = 10 ;5 Asmb; Process to plot
do r=0,nb_cases-1
  print("### Run is "+RUN(r)+" ###")
  do y=0,dimPDF(0)-1
    do i=0, dimP-1
      indi = closest_val(Perc(i),PDF_cum(r,y,pro,:))
      PercPDF(r, y, i) = toint(bin1(indi))
      if y.eq.per2100 then
        print("Percentile: "+Perc(i))
        print((/PercPDF(r, y, i)/))
      end if
    end do
    do i=0, dimPt-1
      indi = closest_val(Perct(i),PDF_cum(r,y,pro,:))
      PerctPDF(r, y, i) = toint(bin1(indi))
    end do
  end do
end do

DiffHighLow = dim_avg_n(PerctPDF(7, per2100-2:per2100, :) - PerctPDF(6, per2100-2:per2100, :),0)
DiffO       = dim_avg_n(PerctPDF(3, per2100-2:per2100, :) - PerctPDF(8, per2100-2:per2100, :),0)
DiffG       = dim_avg_n(PerctPDF(3, per2100-2:per2100, :) - PerctPDF(9, per2100-2:per2100, :),0)
DiffA       = dim_avg_n(PerctPDF(3, per2100-2:per2100, :) - PerctPDF(10, per2100-2:per2100, :),0)
DiffGl      = dim_avg_n(PerctPDF(3, per2100-2:per2100, :) - PerctPDF(11, per2100-2:per2100, :),0)

;### Likely range
up          = 83-1
low         = 17-1
DiffLik     = diffPDFr(low,up,PerctPDF)

;### Very likely range
up           = 95-1
low          = 5-1
DiffVLik     = diffPDFr(low,up,PerctPDF)

;### 1-99
up           = 99-1
low          = 1-1
DiffVC       = diffPDFr(low,up,PerctPDF)

indPercIPCC = ind((Perc.eq.5).or.(Perc.eq.50).or.(Perc.eq.95))

;### Write a text file of Latex table #########################################
col1     = decimalPlaces(Perc(indPercIPCC) ,1,True)
col2     = decimalPlaces(Perc ,1,True)
PercPDFs = tostring(PercPDF)

TableAR5 = [/col1, PercPDFs(0,per2100,indPercIPCC), PercPDFs(1,per2100,indPercIPCC), \
                 PercPDFs(2,per2100,indPercIPCC) /]
write_table("TableAR5_"+SCE+"_PaperDependence.txt", "w", TableAR5, "%4.1f & %s & %s & %s & %s & %s & %s \\")

TableProb  = [/col2, PercPDFs(3,per2100,:), PercPDFs(6,per2100,:), PercPDFs(7,per2100,:), \
             PercPDFs(4,per2100,:), PercPDFs(5,per2100,:) /]
write_table("Table_"+SCE+"_PaperDependence.txt", "w", TableProb, "%4.1f & %s & %s & %s & %s & %s \\")

;Write the correlation matrices:
SCCrd = decimalPlaces(SCCr, 2, True )
dt    = dimPcc(0)-1            ; Time to print out
fout  = "%5.2f"
ELEM2 = (/"GMST      ", "TE        ", "GIC       ", "GSMB      ", "ASMB      ", \
          "Land Water", "AD        ", "GD        " /)
ST1   = " 1.00"
do r=0,nb_cases-1
  NameTableFile = "SpearmanCorr_namelist"+files(r)+"_"+SCE+".txt"
  write_table(NameTableFile, "w", [/ELEM2(0), ST1,SCCrd(r,0,dt), \
  SCCrd(r,1,dt), SCCrd(r,2,dt), SCCrd(r,3,dt), SCCrd(r,4,dt), SCCrd(r,5,dt),     \
  SCCrd(r,6,dt)/], "%s & %s & "+fout+" & "+fout+" & "+fout+" & "+fout+" & "+fout+" & "+fout+" & "+fout+" \\")

  write_table(NameTableFile, "a", [/ELEM2(1), "  -  ", ST1, SCCrd(r,7,dt), \
  SCCrd(r,8,dt), SCCrd(r,9,dt), SCCrd(r,10,dt), SCCrd(r,11,dt), SCCrd(r,12,dt)/], \
  "%s & %s & %s & "+fout+" & "+fout+" & "+fout+" & "+fout+" & "+fout+" & "+fout+" \\")

  write_table(NameTableFile, "a", [/ELEM2(2),"  -  ", "  -  ", ST1, SCCrd(r,13,dt), \
  SCCrd(r,14,dt), SCCrd(r,15,dt), SCCrd(r,16,dt), SCCrd(r,17,dt)/], \
  "%s & %s & %s & %s & "+fout+" & "+fout+" & "+fout+" & "+fout+" & "+fout+" \\")

  write_table(NameTableFile, "a", [/ELEM2(3), "  -  ", "  -  ", "  -  ", ST1, \
  SCCrd(r,18,dt), SCCrd(r,19,dt), SCCrd(r,20,dt), SCCrd(r,21,dt)/], \
  "%s & %s & %s & %s & %s & "+fout+" & "+fout+" & "+fout+" & "+fout+" \\")

  write_table(NameTableFile, "a", [/ELEM2(4), "  -  ", "  -  ", "  -  ", "  -  ", ST1, \
  SCCrd(r,22,dt), SCCrd(r,23,dt), SCCrd(r,24,dt)/], \
  "%s & %s & %s & %s & %s & %s & "+fout+" & "+fout+" & "+fout+" \\")

  write_table(NameTableFile, "a", [/ELEM2(5), "  -  ", "  -  ", "  -  ", "  -  ", "  -  ", \
  ST1, SCCrd(r,25,dt), SCCrd(r,26,dt)/], \
  "%s & %s & %s & %s & %s & %s & %s & "+fout+" & "+fout+" \\")

  write_table(NameTableFile, "a", [/ELEM2(6), "  -  ", "  -  ", "  -  ", "  -  ", "  -  ", \
  "  -  ", ST1, SCCrd(r,27,dt)/], \
  "%s & %s & %s & %s & %s & %s & %s & %s & "+fout+" \\")

  write_table(NameTableFile, "a", [/ELEM2(7), "  -  ", "  -  ", "  -  ", "  -  ", "  -  ", \
  "  -  ", "  -  ", ST1/], \
  "%s & %s & %s & %s & %s & %s & %s & %s & %s \\")

end do

;########### Plot #############################################################
  wks                      = gsn_open_wks("pdf" ,"PlotPDF_dep_"+SCE)
  res                      = True
  res@gsnMaximize          = True
  res@xyMarkLineModes      = (/"Lines","Lines","Lines","Lines"/)
  res@pmLegendDisplayMode  = "Always"
  res@tiYAxisString        = "PDF (%)"
  res@xyLineThicknessF     = 5
  res@gsnPaperOrientation  = "portrait"
  res@xyLineColors         = (/"red","blue","green","black","purple","orange","grey"/)

  ;### Plot comparison of total sea level: Main plot 
  ;delete(res@xyExplicitLegendLabels)
  res@trXMinF              = -10
  res@trXMaxF              = 150
  res@tiXAxisString  = "Global sea level rise (cm)"
  res@tiMainString   = "Plot IPCC "+ SCE
  res@xyExplicitLegendLabels = RUN
  plot               = gsn_csm_xy(wks,bin1,PDF(:,per2100,10,:), res)

  i = 6
  res@trXMinF         = 2006
  res@trXMaxF         = 2100
  res@tiYAxisString   = "Global sea level rise (cm)"
  res@tiXAxisString   = "Time (years)"
  res@tiMainString    = "Plot IPCC "+ SCE + ", perc: "+Perc(i)
  res@xyExplicitLegendLabels = RUN
  plot                = gsn_csm_xy(wks,Years,PercPDF(:,:,i), res)

  res@xyDashPatterns       = 1
  ;### Look into the cross-correlations
  res@tiXAxisString    = "Time (years)"
 
  ;### Look into the correlation time series
  ;### IPCC case
  res@trYMinF              = -1.1
  res@trYMaxF              = 1.1
  res@tiYAxisString    = "Spearman correlation with GMST"
  res@tiMainString   = "Plot IPCC "+ SCE
  delete(res@xyExplicitLegendLabels)
  delete(res@xyLineColors)
  delete(res@xyDashPatterns)
  res@xyLineColors         = (/"red","blue","green","black","purple"/)
  res@xyDashPatterns       = (/1,2,3,4,5/)
  res@xyExplicitLegendLabels = ELEM((/1,2,3,4,6/))
  plot               = gsn_csm_xy(wks,TimeCC,SCCr(0,(/0,1,2,3,5/),:), res)

  res@tiYAxisString    = "Spearman correlation with Glaciers"
  delete(res@xyExplicitLegendLabels)
  delete(res@xyLineColors)
  delete(res@xyDashPatterns)
  res@xyLineColors           = (/"red","green","black","purple"/)
  res@xyExplicitLegendLabels = ELEM((/1,3,4,6/))
  res@xyDashPatterns         = (/1,3,4,5/)
  plot               = gsn_csm_xy(wks,TimeCC,SCCr(0,(/7,13,14,16/),:), res)

  ;### IPCC + Levermann 2014 case
  res@tiYAxisString  = "Spearman correlation with GMST"
  res@tiMainString   = "Plot IPCC+Lev14 "+ SCE
  delete(res@xyExplicitLegendLabels)
  delete(res@xyLineColors)
  delete(res@xyDashPatterns)
  res@xyLineColors         = (/"red","blue","green","black","purple"/)
  res@xyDashPatterns       = (/1,2,3,4,5/)
  res@xyExplicitLegendLabels = ELEM((/1,2,3,4,6/))
  plot               = gsn_csm_xy(wks,TimeCC,SCCr(3,(/0,1,2,3,5/),:), res)

  res@tiYAxisString    = "Spearman correlation with Glaciers"
  delete(res@xyExplicitLegendLabels)
  delete(res@xyLineColors)
  delete(res@xyDashPatterns)
  res@xyLineColors         = (/"red","green","black","purple"/)
  res@xyExplicitLegendLabels = ELEM((/1,3,4,6/))
  res@xyDashPatterns       = (/1,3,4,5/)
  plot               = gsn_csm_xy(wks,TimeCC,SCCr(3,(/7,13,14,16/),:), res)

  ;### IPCC + Levermann 2014 high correlation case
  res@tiYAxisString    = "Spearman correlation with GMST"
  res@tiMainString     = "Plot IPCC+Lev14 "+ SCE+", high"
  delete(res@xyExplicitLegendLabels)
  delete(res@xyLineColors)
  delete(res@xyDashPatterns)
  res@xyLineColors         = (/"red","blue","green","black","purple"/)
  res@xyDashPatterns       = (/1,2,3,4,5/)
  res@xyExplicitLegendLabels = ELEM((/1,2,3,4,6/))
  plot               = gsn_csm_xy(wks,TimeCC,SCCr(7,(/0,1,2,3,5/),:), res)

  res@tiYAxisString    = "Spearman correlation with Glaciers"
  delete(res@xyExplicitLegendLabels)
  delete(res@xyLineColors)
  delete(res@xyDashPatterns)
  res@xyLineColors         = (/"red","green","black","purple"/)
  res@xyExplicitLegendLabels = ELEM((/1,3,4,6/))
  res@xyDashPatterns       = (/1,3,4,5/)
  plot               = gsn_csm_xy(wks,TimeCC,SCCr(7,(/7,13,14,16/),:), res)

  ;### Uncertainty comparison in 2100
  res1 = True
  res1@gsnMaximize          = True
;  res1@xyMarkLineModes      = (/"Lines","Lines","Lines","Lines"/)
  res1@pmLegendDisplayMode  = "Always"
  res1@xyLineThicknessF     = 5
  res1@gsnPaperOrientation  = "portrait"
  res1@xyLineColors         = (/"red","blue","green","black","purple","orange","grey"/)
  res1@tiMainString         = "Uncertainty comparison "+ SCE
  res1@tiYAxisString        = "Uncertainty contribution (cm)"
  res1@tiXAxisString        = "Percentile of total distribution"
  res1@xyExplicitLegendLabels = (/"Dependence","Antarctica","Greenland","Thermal expansion","Glaciers and Ice Caps"/)

  plot               = gsn_csm_xy(wks,Perct,(/DiffHighLow,DiffA,DiffG,DiffO,DiffGl/),res1)

  ;### Time evolution of uncertainty
  res1@tiMainString  = "Uncertainty comp. over time, likely range, "+ SCE
  res1@tiYAxisString = "Uncertainty contribution (cm)"
  res1@tiXAxisString = "Time (years)"
  plot               = gsn_csm_xy(wks, Years, (/DiffLik(0,:),DiffLik(3,:),DiffLik(2,:), \
                                              DiffLik(1,:),DiffLik(4,:)/), res1)

  res1@tiMainString  = "Uncertainty comp. over time, very likely range, "+ SCE
  plot               = gsn_csm_xy(wks, Years, (/DiffVLik(0,:),DiffVLik(3,:),DiffVLik(2,:), \
                                              DiffVLik(1,:),DiffVLik(4,:)/), res1)

  res1@tiMainString  = "Uncertainty comp. over time, 1-99 range, "+ SCE
  plot               = gsn_csm_xy(wks, Years, (/DiffVC(0,:),DiffVC(3,:),DiffVC(2,:), \
                                              DiffVC(1,:),DiffVC(4,:)/), res1)

end
