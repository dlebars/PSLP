#!/bin/bash
### Loop over a few cases of sea level projections, read namelist and call main.ncl

#Exp_dir="ExpTest/"
Exp_dir="Exp_PaperDependence/"

#for namelist in namelistIPCC{2..3}
for namelist in namelistIPCC{1..3}  namelistProb{1..9} 
do
  echo '#######################################################################'
  echo '#############' Using $Exp_dir$namelist.txt '###########################'
  opt="$(grep -v '^#' $Exp_dir$namelist".txt")"
  for SCE in rcp45 rcp85
  do
    echo Running NCL with options: SCE='"'$SCE'"' $opt
    NameOutput="SeaLevelPDF_"$namelist"_"$SCE".nc"
    ncl main.ncl SCE='"'$SCE'"' $opt NameOut='"'$Exp_dir$NameOutput'"'
  done
done
