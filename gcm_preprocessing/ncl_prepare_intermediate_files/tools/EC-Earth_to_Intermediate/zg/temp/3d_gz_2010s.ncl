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

 year = new((/10/),"string")
 year = (/"2010","2011","2012","2013","2014","2015","2016","2017","2018","2019"/)

;************************************************
; file handling
;************************************************
  print("opening files")
  path_input = "/g/data/gh02/CMIP6/N2/et_scripts/gcms/EC-Earth3-Veg/preprocessed/regridded2/historical/"
  path_output ="./"
do ny = 0, 9 
  yyyy = year(ny)
  print(yyyy)
  in_ta    = addfile(path_input+"ta_6hrLev_EC-Earth3-Veg_historical_r1i1p1f1_gr_"+yyyy+"01010000-"+yyyy+"12311800.nc","r")      ; 6-hourly 3-d T              
  in_hus    = addfile(path_input+"hus_6hrLev_EC-Earth3-Veg_historical_r1i1p1f1_gr_"+yyyy+"01010000-"+yyyy+"12311800.nc","r")      ; 6-hourly 3-d T              
  in_ps    = addfile(path_input+"ps_3hr_EC-Earth3-Veg_historical_r1i1p1f1_gr_"+yyyy+"01010000-"+yyyy+"12312100.nc","r")      ; 6-hourly 3-d T              
  in_tas    = addfile(path_input+"tas_3hr_EC-Earth3-Veg_historical_r1i1p1f1_gr_"+yyyy+"01010000-"+yyyy+"12312100.nc","r")      ; 6-hourly 3-d T              
  in_orog    = addfile("/g/data/gh02/CMIP6/N2/et_scripts/gcms/EC-Earth3-Veg/preprocessed/regridded2/orog_fx_EC-Earth3_historical_r1i1p1f1_gr.nc","r")      ; 6-hourly 3-d T              


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
  plev   = in_ta->plev(::-1)                            
;  lev_bnds   = in_ta->lev_bnds                            
  ta   = in_ta->ta(:,::-1,:,:)   ; from bottom to top                         
  tas   = in_tas->tas(0::2,:,:) ; start from 0UTC with an interval of 6 hr                             
;  printVarSummary(tas)
  ps   = in_ps->ps(0::2,:,:)                            
  hus   = in_hus->hus(:,::-1,:,:) ; from bottom to top                           
  orog = in_orog->orog
 
 dim = dimsizes(ta)
 dim0 = dim(0)
 dim1 = dim(1)
 dim2 = dim(2)
 dim3 = dim(3)

 Tmean = new((/dim0, dim1, dim2, dim3/),"float")
 zg = new((/dim0, dim1, dim2, dim3/),"double")
 copy_VarCoords(ta, zg)
;  Tmean = hus
;  zg = hus

   Tmean(:,0,:,:) = (tas + ta(:,0,:,:)*(1+0.63*hus(:,0,:,:)))/2.
  do n =1, 74
   Tmean(:,n,:,:) = (ta(:,n,:,:)*(1+0.63*hus(:,n,:,:))+ta(:,n-1,:,:)*(1+0.63*hus(:,n-1,:,:)))/2.
  end do

  g0 = 9.81
  R0 = 287.04; 8.314462

  d_orog = conform_dims(dimsizes(tas),orog,(/1,2/)) 
 
  p0 = ps/(/plev(0)/)
;  printVarSummary(d_orog)
;  printVarSummary(zg(:,0,:,:))
;  printVarSummary(Tmean(:,0,:,:))
;  printVarSummary(ps)
;  printVarSummary(p0)
;  printVarSummary(plev(0))
 

  zg(:,0,:,:) = d_orog + R0*Tmean(:,0,:,:)*log(p0)/g0

  do n = 1,74
  print(n)
  pn = (/plev(n-1)/plev(n)/)
  zg(:,n,:,:) = zg(:,n-1,:,:)+R0*Tmean(:,n,:,:)*log(pn)/g0
  end do


  fout    = addfile("zg_6hrLev_EC-Earth3-Veg_historical_r1i1p1f2_gr_"+yyyy+"01010000-"+yyyy+"12311800.nc","c")      ; 6-hourly 3-d zg              
  fatt = True
  fatt@title = "calculated 3D geopotential on pressure levels"
  fatt@date = systemfunc("date")
  fileattdef(fout, fatt)
  fout->time = time
  fout->lat = lat
  fout->lat_bnds = lat_bnds
  fout->lon = lon
  fout->lon_bnds = lon_bnds
  fout->lev = plev(::-1)
;  fout->lev_bnds = lev_bnds
  fout->zg = zg(:,::-1,:,:)

delete(ta)
delete(tas)
delete(ps)
delete(hus)
delete(plev)
delete(Tmean)
delete(yyyy)
delete(zg)
delete(lat)
delete(lon)
delete(p0)
delete(pn)
delete(d_orog)
delete(time)
delete(fout)
delete(in_ta)
delete(in_hus)
delete(in_ps)
delete(in_tas)
delete(in_orog)
end do


end
