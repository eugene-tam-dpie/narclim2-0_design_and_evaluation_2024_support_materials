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
 year = (/"1950","1951","1952","1953","1954","1955","1956","1957","1958","1959"/)

;************************************************
; file handling
;************************************************
  print("opening files")
  path_input = "/g/data/gh02/CMIP6/N2/scripts/tools/MPI-ESM1-2-HR_to_Intermediate/MPI_TO_INT_1950_model/data_temp/"
  path_output ="./"
do ny = 0, 0 
  yyyy = year(ny)
  print(yyyy)
  in_ta    = addfile(path_input+"ta_6hrLev_1950.nc","r")      ; 6-hourly 3-d T              
  in_hus    = addfile(path_input+"hus_6hrLev_1950.nc","r")      ; 6-hourly 3-d T              
  in_huss    = addfile(path_input+"huss_6hrPlevPt_1950.nc","r")      ; 6-hourly 3-d T              
  in_ps    = addfile(path_input+"ps_3hr_1950.nc","r")      ; 6-hourly 3-d T              
  in_tas    = addfile(path_input+"tas_6hrPlevPt_1950.nc","r")      ; 6-hourly 3-d T              
  in_orog    = addfile(path_input+"orog_fx_MPI-ESM1-2-HR_historical_r1i1p1f1_gn.nc","r")      ; 6-hourly 3-d T              


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
  tas   = in_tas->tas ; start from 0UTC with an interval of 6 hr                             
;  ps   = in_ps->ps(0::2,:,:)                            
  hus   = in_hus->hus ; from bottom to top                           
  huss   = in_huss->huss ; from bottom to top                           
  orog = in_orog->orog
 
  b = in_ta->b ; from bottom to top
  ap = in_ta->ap; from bottom to top
  ps =in_ta->ps

 plev = conform(ta,ap,1) + conform(ta,b,1) * conform(ta,ps,(/0,2,3/)) 

 
 dim = dimsizes(ta)
 dim0 = dim(0)
 dim1 = dim(1)
 dim2 = dim(2)
 dim3 = dim(3)

 Tmean = new((/dim0, dim1, dim2, dim3/),"float")
 zg = new((/dim0, dim1, dim2, dim3/),"float")
 plog = new((/dim0, dim1, dim2, dim3/),"float")
 Tmean = new((/dim0, dim1, dim2, dim3/),"float")
 copy_VarCoords(ta, zg)
 copy_VarCoords(ta, plev)
 copy_VarCoords(ta, plog)
 copy_VarCoords(ta, Tmean)
;  Tmean = hus
;  zg = hus

   Tmean(:,0,:,:) = (tas*(1+0.61*huss) + ta(:,0,:,:)*(1+0.61*hus(:,0,:,:)))/2.
  do n =1, dim1 -1
   Tmean(:,n,:,:) = (ta(:,n,:,:)*(1+0.61*hus(:,n,:,:))+ta(:,n-1,:,:)*(1+0.61*hus(:,n-1,:,:)))/2.
  end do

  g0 = 9.81
  R0 = 287.04; 8.314462

  d_orog = conform_dims(dimsizes(tas),orog,(/1,2/)) 
 
  p0 = ps/plev(:,0,:,:)
;  printVarSummary(d_orog)
;  printVarSummary(zg(:,0,:,:))
;  printVarSummary(Tmean(:,0,:,:))
;  printVarSummary(ps)
;  printVarSummary(p0)
;  printVarSummary(plev(0))
 

  zg(:,0,:,:) =doubletofloat(d_orog + R0*Tmean(:,0,:,:)*log(p0)/g0)

 plog(:,0,:,:) = doubletofloat(log(p0))

  do n = 1,dim1 -1
  print(n)
  pn = (/plev(:,n-1,:,:)/plev(:,n,:,:)/)
  plog(:,n,:,:) = doubletofloat(log(pn))
  zg(:,n,:,:) =doubletofloat(zg(:,n-1,:,:)+R0*Tmean(:,n,:,:)*log(pn)/g0)
  end do


  fout    = addfile("zg_6hrLev_1950_for_ET.nc","c")      ; 6-hourly 3-d zg              
  fatt = True
  fatt@title = "calculated 3D geopotential on model levels"
  fatt@date = systemfunc("date")
  fileattdef(fout, fatt)
  fout->time = time
  fout->lat = lat
  fout->lat_bnds = lat_bnds
  fout->lon = lon
  fout->lon_bnds = lon_bnds
  fout->zg = zg ;  bottom to top
delete(fatt)
delete(fout)

  fout    = addfile("plev_6hrLev_1950_for_ET.nc","c")      ; 6-hourly 3-d zg              
  fatt = True
  fatt@title = "calculated 3D pressure levels"
  fatt@date = systemfunc("date")
  fileattdef(fout, fatt)
  fout->time = time
  fout->lat = lat
  fout->lat_bnds = lat_bnds
  fout->lon = lon
  fout->lon_bnds = lon_bnds
  fout->plev = plev  ; from bottom to top
delete(fatt)
delete(fout)

  fout    = addfile("plog_6hrLev_1950_for_ET.nc","c")      ; 6-hourly 3-d zg              
  fatt = True
  fatt@title = "calculated 3D log pressure levels"
  fatt@date = systemfunc("date")
  fileattdef(fout, fatt)
  fout->time = time
  fout->lat = lat
  fout->lat_bnds = lat_bnds
  fout->lon = lon
  fout->lon_bnds = lon_bnds
  fout->plog = plog  ; from bottom to top
delete(fatt)
delete(fout)

  fout    = addfile("Tmean_6hrLev_1950_for_ET.nc","c")      ; 6-hourly 3-d zg              
  fatt = True
  fatt@title = "calculated 3D Tmean on model levels"
  fatt@date = systemfunc("date")
  fileattdef(fout, fatt)
  fout->time = time
  fout->lat = lat
  fout->lat_bnds = lat_bnds
  fout->lon = lon
  fout->lon_bnds = lon_bnds
  fout->Tmean = Tmean  ; from bottom to top

delete(fout)
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
