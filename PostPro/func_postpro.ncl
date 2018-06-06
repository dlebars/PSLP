;###############################################################################
; func_postpro.ncl: Functions used in post-processing
;###############################################################################

;### function: diffPDFr ########################################################
; Compute the difference of a given probability range using a PDF as an input
;###############################################################################
undef("diffPDFr")
function diffPDFr(low, up, PerctPDF)
  local dimPDF, DiffMat
  begin

  dimPDF         = dimsizes(PerctPDF)
  DiffMat        = new((/5, dimPDF(1)/),float)

  DiffMat(0,:)   = runave((PerctPDF(7, :, up) - PerctPDF(7, :, low)) \
                   - (PerctPDF(6, :, up) - PerctPDF(6, :, low)),5,0)
  DiffMat(1,:)   = runave((PerctPDF(3, :, up) - PerctPDF(3, :, low)) \
                   - (PerctPDF(8, :, up) - PerctPDF(8, :, low)),5,0)
  DiffMat(2,:)   = runave((PerctPDF(3, :, up) - PerctPDF(3, :, low)) \
                   - (PerctPDF(9, :, up) - PerctPDF(9, :, low)),5,0)
  DiffMat(3,:)   = runave((PerctPDF(3, :, up) - PerctPDF(3, :, low)) \
                   - (PerctPDF(10, :, up) - PerctPDF(10, :, low)),5,0)
  DiffMat(4,:)   = runave((PerctPDF(3, :, up) - PerctPDF(3, :, low)) \
                   - (PerctPDF(11, :, up) - PerctPDF(11, :, low)),5,0)

  return DiffMat
end
