;******************************************************************************
; NCL script to calcualte geopotential from other variables 
;
; This script does not look for specific dates in a file - simply convert
; all dates in the input file to IM format, so it is the users responsibility
; to ensure the input data is consistent
;******************************************************************************

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/popRemap.ncl"

begin

 year = new((/11/),"string")
 year = (/"1950","1951","1952","1953","1954","1955","1956","1957","1958","1959","1960"/)

;************************************************
; file handling
;************************************************
  print("opening files")
  path_input = "../all/"
  path_input1 = "/g/data/gh02/CMIP6/N2/et_scripts/gcms/CSIRO-BOM-ACCESS-ESM1-5/preprocessed/cmorized_regridded/fx/"
  path_output ="./"
  path_input2 = "/g/data/fs38/publications/CMIP6/CMIP/CSIRO/ACCESS-ESM1-5/historical/r6i1p1f1/6hrLev/ta/gn/latest/"
do ny = 1, 9 
  yyyy = year(ny)
  yyyy1 = year(ny+1)
  yyyy2 = "1950"
  yyyy3 = "1960"
  print(yyyy)
  in_ta    = addfile(path_input+"ta_6hrLev_ACCESS-ESM1-5_historical_r6i1p1f1_gn_"+yyyy+"01010000-"+yyyy+"12311800.nc","r")      ; 6-hourly 3-d T              

  in_ta1    = addfile(path_input2+"ta_6hrLev_ACCESS-ESM1-5_historical_r6i1p1f1_gn_195001010600-195101010000.nc","r")      ; 6-hourly 3-d T              

  in_orog    = addfile(path_input1+"orog_fx_ACCESS-ESM1-5_historical_r6i1p1f1_gn.nc","r")      ; 6-hourly 3-d T              


;************************************************
; Hypsometric eqestion for calculating 3D gz
;************************************************
; p2 = p1/exp(g(z2-z1)/R/Tv)
; z2 = z1+RdTv/g*ln(P1/P2)

;************************************************
; read needed variables from file
;************************************************
  print("read variables in")

  time   = in_ta->time                            
  lat   = in_ta->lat                            
  lat_bnds   = in_ta->lat_bnds                            
  lon   = in_ta->lon                            
  lon_bnds   = in_ta->lon_bnds                            
  ta   = in_ta->ta   ; from bottom to top                         
  lev   = in_ta->lev   ; from bottom to top                         
  orog = in_orog->orog

  b = in_ta1->b ; from bottom to top

 
 dim = dimsizes(ta)
 dim0 = dim(0)
 dim1 = dim(1)
 dim2 = dim(2)
 dim3 = dim(3)

 Tmean = new((/dim0, dim1, dim2, dim3/),"float")
 Umean = new((/dim0, dim1, dim2, dim3/),"float")
 Vmean = new((/dim0, dim1, dim2, dim3/),"float")
 zg = new((/dim0, dim1, dim2, dim3/),"double")

 zg = conform(ta,lev,1) + conform(ta,b,1) * conform(ta,orog,(/2,3/)) 

  zg1 = doubletofloat(zg)
copy_VarCoords(ta,zg1)
  fout    = addfile("zg_6hrLev_ACCESS-ESM1-5_historical_r6i1p1f1_gn_"+yyyy+"01010000-"+yyyy+"12311800.nc","c")      ; 6-hourly 3-d zg              
  fatt = True
  fatt@title = "calculated 3D geopotential on pressure levels"
  fatt@date = systemfunc("date")
  fileattdef(fout, fatt)
  fout->time = time
  fout->lat = lat
  fout->lat_bnds = lat_bnds
  fout->lon = lon
  fout->lon_bnds = lon_bnds
  fout->zg = zg1 ;  bottom to top

delete(fout)

delete(lev)
delete(orog)
delete(b)
delete(zg)
delete(zg1)
delete(lat)
delete(lon)
delete(in_ta)
delete(in_ta1)
delete(in_orog)
end do


end
