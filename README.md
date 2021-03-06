# PSLP
The Probabilistic Sea Level Projection model makes sea level projections using the process-based method. It is based on the IPCC AR5 process-based projection with additional capabilities. 

To run the model:

1) Install the NCL language (http://www.ncl.ucar.edu/)
2) Unzip the InputData file
3) Prepare a namelist (see examples in the Exp_PaperDependence folder)
4) Modify the LoopProjections.bash script with the name of your namelist
5) Run ./LoopProjections.bash , it will read in the namelist and run main.ncl
6) Output data will be in the same folder as the namelist
7) Data can then be analyzed using the code in the PostPro folder.

Only a few options are now included mostly for reproducibility of the results discussed in Le Bars et al. 2017 and in Le Bars et al. 2018. The input data for regional projections are not yet included. Contact me if you are interested in this option (bars .at. knmi.nl)

The PSLP model was used in:

Le Bars, D., Drijfhout, S., & de Vries, H. (2017). A high-end sea level rise probabilistic projection including rapid Antarctic ice sheet mass loss. Environmental Research Letters, 12(4), 044013. http://doi.org/10.1088/1748-9326/aa6512

Le Bars, D. (2018, March 8). Uncertainty in sea level rise projections due to the dependence between contributors. http://doi.org/10.17605/OSF.IO/UVW3S
