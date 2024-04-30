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
 year = (/"2080","2081","2082","2083","2084","2085","2086","2087","2088","2089"/)

;************************************************
; file handling
;************************************************
  print("opening files")
  path_input = "/g/data/gh02/CMIP6/N2/et_scripts/gcms/EC-Earth3-Veg/preprocessed/model_lev/ssp370/"
  path_input1 = "/g/data/gh02/CMIP6/N2/et_scripts/gcms/EC-Earth3-Veg/preprocessed/regridded0/"
  path_output ="./"
do ny = 0, 9 
  yyyy = year(ny)
  print(yyyy)
  in_ua    = addfile(path_input+"ua_6hrLev_EC-Earth3-Veg_ssp370_r1i1p1f1_gr_"+yyyy+"01010000-"+yyyy+"12311800.nc","r")      ; 6-hourly 3-d T              
  in_va    = addfile(path_input+"va_6hrLev_EC-Earth3-Veg_ssp370_r1i1p1f1_gr_"+yyyy+"01010000-"+yyyy+"12311800.nc","r")      ; 6-hourly 3-d T              
  in_ta    = addfile(path_input+"ta_6hrLev_EC-Earth3-Veg_ssp370_r1i1p1f1_gr_"+yyyy+"01010000-"+yyyy+"12311800.nc","r")      ; 6-hourly 3-d T              
  in_hus    = addfile(path_input+"hus_6hrLev_EC-Earth3-Veg_ssp370_r1i1p1f1_gr_"+yyyy+"01010000-"+yyyy+"12311800.nc","r")      ; 6-hourly 3-d T              

  in_uas    = addfile(path_input+"uas_E3hr_EC-Earth3-Veg_ssp370_r1i1p1f1_gr_"+yyyy+"01010130-"+yyyy+"12312230.nc","r")      ; 6-hourly 3-d T              
  in_vas    = addfile(path_input+"vas_E3hr_EC-Earth3-Veg_ssp370_r1i1p1f1_gr_"+yyyy+"01010130-"+yyyy+"12312230.nc","r")      ; 6-hourly 3-d T              
  in_tas    = addfile(path_input+"tas_3hr_EC-Earth3-Veg_ssp370_r1i1p1f1_gr_"+yyyy+"01010000-"+yyyy+"12312100.nc","r")      ; 6-hourly 3-d T              
  in_huss    = addfile(path_input+"huss_3hr_EC-Earth3-Veg_ssp370_r1i1p1f1_gr_"+yyyy+"01010000-"+yyyy+"12312100.nc","r")      ; 6-hourly 3-d T              


;  in_ua    = addfile(path_input+"ua_6hrLev_2080.nc","r")      ; 6-hourly 3-d T              
;  in_va    = addfile(path_input+"va_6hrLev_2080.nc","r")      ; 6-hourly 3-d T              
;  in_ta    = addfile(path_input+"ta_6hrLev_2080.nc","r")      ; 6-hourly 3-d T              
;  in_hus    = addfile(path_input+"hus_6hrLev_2080.nc","r")      ; 6-hourly 3-d T              
;  in_ps    = addfile(path_input+"ps_3hr_2080.nc","r")      ; 6-hourly 3-d T              
;  in_tas    = addfile(path_input+"tas_6hrPlevPt_2080.nc","r")      ; 6-hourly 3-d T              
;  in_uas    = addfile(path_input+"uas_6hrPlevPt_2080.nc","r")      ; 6-hourly 3-d T              
;  in_vas    = addfile(path_input+"vas_6hrPlevPt_2080.nc","r")      ; 6-hourly 3-d T              

  in_orog    = addfile(path_input1+"orog_fx_EC-Earth3_historical_r1i1p1f1_gr.nc","r")      ; 6-hourly 3-d T              


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
  ua   = in_ua->ua   ; from bottom to top                         
  va   = in_va->va   ; from bottom to top                         

  tas   = in_tas->tas(::2,:,:) ; start from 0UTC with an interval of 6 hr                             
  uas   = in_uas->uas(::2,:,:) ; start from 0UTC with an interval of 6 hr                             
  vas   = in_vas->vas(::2,:,:) ; start from 0UTC with an interval of 6 hr                             
  huss   = in_huss->huss(::2,:,:) ; start from 0UTC with an interval of 6 hr                             
;  ps   = in_ps->ps(0::2,:,:)                            
  hus   = in_hus->hus ; from bottom to top                           
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
 Umean = new((/dim0, dim1, dim2, dim3/),"float")
 Vmean = new((/dim0, dim1, dim2, dim3/),"float")
 zg = new((/dim0, dim1, dim2, dim3/),"double")
 copy_VarCoords(ta, zg)
 copy_VarCoords(ta, plev)
;  Tmean = hus
;  zg = hus

   Tmean(:,0,:,:) = (tas*(1+0.63*huss) + ta(:,0,:,:)*(1+0.63*hus(:,0,:,:)))/2.
   Umean(:,0,:,:) = (uas + ua(:,0,:,:))/2.
   Vmean(:,0,:,:) = (vas + va(:,0,:,:))/2.
  do n =1, dim1 -1
   Tmean(:,n,:,:) = (ta(:,n,:,:)*(1+0.63*hus(:,n,:,:))+ta(:,n-1,:,:)*(1+0.63*hus(:,n-1,:,:)))/2.
   Umean(:,n,:,:) = (ua(:,n,:,:)+ua(:,n-1,:,:))/2.
   Vmean(:,n,:,:) = (va(:,n,:,:)+va(:,n-1,:,:))/2.
  end do

  g0 = 9.81
  R0 = 287.04; 8.314462
  lat1 = lat*3.1415/180
  Omega = 7.2921159e-5
  r = 6.3781e6
  d_orog = conform_dims(dimsizes(tas),orog,(/1,2/)) 
 
  do n =1, dim1 -1
   Tmean(:,n,:,:) = (ta(:,n,:,:)*(1+0.63*hus(:,n,:,:))+ta(:,n-1,:,:)*(1+0.63*hus(:,n-1,:,:)))/2.
   Umean(:,n,:,:) = (ua(:,n,:,:)+ua(:,n-1,:,:))/2.
   Vmean(:,n,:,:) = (va(:,n,:,:)+va(:,n-1,:,:))/2.
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
   cosf = cos(lat1)
   printVarSummary(cosf)
   cosf_c = conform_dims(dimsizes(tas), cosf, 1)
   printVarSummary(cosf_c) 
   Umean02 = Umean(:,0,:,:)*Umean(:,0,:,:)
   Vmean02 = Umean(:,0,:,:)*Vmean(:,0,:,:)


  A0= -1*(2*Omega*Umean(:,0,:,:)*cosf_c + (Umean02+Vmean02)/r)/g0
 
  zg(:,0,:,:) = d_orog + R0*Tmean(:,0,:,:)*log(p0)/(g0*(1+A0))


  do n = 1,dim1 -1
  print(n)
  pn = (/plev(:,n-1,:,:)/plev(:,n,:,:)/)

   Umean2 = Umean(:,n,:,:)*Umean(:,n,:,:)
   Vmean2 = Vmean(:,n,:,:)*Vmean(:,n,:,:)

  A= -1*(2*Omega*Umean(:,n,:,:)*cosf_c+(Umean2+Vmean2)/r)/g0
  zg(:,n,:,:) = zg(:,n-1,:,:)+R0*Tmean(:,n,:,:)*log(pn)/(g0*(1+A))
  end do

  plev1 = doubletofloat(plev)
  zg1 = doubletofloat(zg)
copy_VarCoords(ta,zg1)
copy_VarCoords(ta,plev1)
  fout    = addfile("zg_6hrLev_EC-Earth3-Veg_ssp370_r1i1p1f1_gr_"+yyyy+"01010000-"+yyyy+"12311800.nc","c")      ; 6-hourly 3-d zg              
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

  fout    = addfile("plev_6hrLev_EC-Earth3-Veg_ssp370_r1i1p1f1_gr_"+yyyy+"01010000-"+yyyy+"12311800.nc","c")      ; 6-hourly 3-d zg              
  fatt = True
  fatt@title = "calculated 3D geopotential on pressure levels"
  fatt@date = systemfunc("date")
  fileattdef(fout, fatt)
  fout->time = time
  fout->lat = lat
  fout->lat_bnds = lat_bnds
  fout->lon = lon
  fout->lon_bnds = lon_bnds
  fout->plev = plev1  ; from bottom to top

delete(ta)
delete(tas)
delete(ps)
delete(hus)
delete(huss)
delete(ua)
delete(uas)
delete(va)
delete(vas)
delete(cosf_c)
delete(plev)
delete(plev1)
delete(Tmean)
delete(Umean)
delete(Vmean)
delete(yyyy)
delete(zg)
delete(zg1)
delete(lat)
delete(lon)
delete(A0)
delete(A)
delete(p0)
delete(pn)
delete(d_orog)
delete(time)
delete(fout)
delete(in_ta)
delete(in_tas)
delete(in_hus)
delete(in_huss)
delete(in_ua)
delete(in_uas)
delete(in_va)
delete(in_vas)
delete(in_orog)
delete(Umean02)
delete(Umean2)
delete(Vmean02)
delete(Vmean2)

end do


end
