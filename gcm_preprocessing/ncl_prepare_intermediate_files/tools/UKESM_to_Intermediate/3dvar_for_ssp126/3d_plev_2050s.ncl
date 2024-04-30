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
 year = (/"2050","2051","2052","2053","2054","2055","2056","2057","2058","2059"/)

;************************************************
; file handling
;************************************************
  print("opening files")
  path_input = "/g/data/gh02/CMIP6/N2/et_scripts/gcms/UKESM1-0-LL/preprocessed/model_level/ssp126/"
  path_input2 ="/g/data/gh02/CMIP6/N2/scripts/tools/UKESM_to_Intermediate/data_ssp126/temp/"
  path_input1 = "/g/data/gh02/CMIP6/N2/et_scripts/gcms/UKESM1-0-LL/preprocessed/model_level/piControl/"
  path_output ="/g/data/gh02/CMIP6/N2/scripts/tools/UKESM_to_Intermediate/data_ssp126/"
do ny = 0, 9 
  yyyy = year(ny)
  print(yyyy)
  in_ua    = addfile(path_input+"ua_6hrLev_UKESM1-0-LL_ssp126_r1i1p1f2_gn_"+yyyy+"01010000-"+yyyy+"12311800.nc","r")      ; 6-hourly 3-d T              
  in_va    = addfile(path_input+"va_6hrLev_UKESM1-0-LL_ssp126_r1i1p1f2_gn_"+yyyy+"01010000-"+yyyy+"12311800.nc","r")      ; 6-hourly 3-d T              
  in_ta    = addfile(path_input+"ta_6hrLev_UKESM1-0-LL_ssp126_r1i1p1f2_gn_"+yyyy+"01010000-"+yyyy+"12311800.nc","r")      ; 6-hourly 3-d T              
  in_hus    = addfile(path_input+"hus_6hrLev_UKESM1-0-LL_ssp126_r1i1p1f2_gn_"+yyyy+"01010000-"+yyyy+"12311800.nc","r")      ; 6-hourly 3-d T              
  in_ps    = addfile(path_input+"ps_6hrLev_UKESM1-0-LL_ssp126_r1i1p1f2_gn_"+yyyy+"01010000-"+yyyy+"12311800.nc","r")      ; 6-hourly 3-d T              

  in_ta1    = addfile(path_input2+"ta_6hrLev_UKESM1-0-LL_ssp126_r1i1p1f2_gn_"+yyyy+"01010000-"+yyyy+"12311800.nc","r")      ; 6-hourly 3-d T              
;  in_uas    = addfile(path_input+"uas_E3hr_UKESM1-0-LL_ssp126_r1i1p1f2_gn_"+yyyy+"01010130-"+yyyy+"12312230.nc","r")      ; 6-hourly 3-d T              
;  in_vas    = addfile(path_input+"vas_E3hr_UKESM1-0-LL_ssp126_r1i1p1f2_gn_"+yyyy+"01010130-"+yyyy+"12312230.nc","r")      ; 6-hourly 3-d T              
  in_tas    = addfile(path_input+"tas_3hr_UKESM1-0-LL_ssp126_r1i1p1f2_gn_"+yyyy+"01010000-"+yyyy+"12312100.nc","r")      ; 6-hourly 3-d T              
;  in_huss    = addfile(path_input+"huss_3hr_UKESM1-0-LL_ssp126_r1i1p1f2_gn_"+yyyy+"01010000-"+yyyy+"12312100.nc","r")      ; 6-hourly 3-d T              


  in_orog    = addfile(path_input1+"orog_fx_UKESM1-0-LL_piControl_r1i1p1f2_gn_fx.nc","r")      ; 6-hourly 3-d T              


;************************************************
; Hypsometric eqestion for calculating 3D gz
;************************************************
;   p2 = p1/exp(g(z2-z1)/R/Tv)
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
;  uas   = in_uas->uas(::2,:,:) ; start from 0UTC with an interval of 6 hr                             
;  vas   = in_vas->vas(::2,:,:) ; start from 0UTC with an interval of 6 hr                             
;  huss   = in_huss->huss(::2,:,:) ; start from 0UTC with an interval of 6 hr                             
  ps   = in_ps->ps                            
  hus   = in_hus->hus ; from bottom to top                           
  orog = in_orog->orog


 
  b = in_ta1->b ; from bottom to top
  lev = in_ta1->lev; from bottom to top
  ;ap = in_ta1->ap; from bottom to top
printVarSummary(ta)
printVarSummary(orog)

 zg = conform(ta,lev,1) + conform(ta,b,1) * conform(ta,orog,(/2,3/)) 
 ;plev = conform(ta,ap,1) + conform(ta,b,1) * conform(ta,ps,(/0,2,3/)) 

 
 dim = dimsizes(ta)
 dim0 = dim(0)
 dim1 = dim(1)
 dim2 = dim(2)
 dim3 = dim(3)

 Tmean = new((/dim0, dim1, dim2, dim3/),"float")
 Umean = new((/dim0, dim1, dim2, dim3/),"float")
 Vmean = new((/dim0, dim1, dim2, dim3/),"float")
 plev = new((/dim0, dim1, dim2, dim3/),"double")
 ;zg = new((/dim0, dim1, dim2, dim3/),"double")
 copy_VarCoords(ta, zg)
 copy_VarCoords(ta, plev)

   ;Tmean(:,0,:,:) = (tas*(1+0.63*huss) + ta(:,0,:,:)*(1+0.63*hus(:,0,:,:)))/2.
   ;Umean(:,0,:,:) = (uas + ua(:,0,:,:))/2.
   ;Vmean(:,0,:,:) = (vas + va(:,0,:,:))/2.
   Tmean(:,0,:,:) =  ta(:,0,:,:)*(1+0.63*hus(:,0,:,:))
   Umean(:,0,:,:) =  ua(:,0,:,:)
   Vmean(:,0,:,:) =  va(:,0,:,:)
  do n =1, dim1 -1
   Tmean(:,n,:,:) = (ta(:,n,:,:)*(1+0.63*hus(:,n,:,:))+ta(:,n-1,:,:)*(1+0.63*hus(:,n-1,:,:)))/2.
   Umean(:,n,:,:) = (ua(:,n,:,:)+ua(:,n-1,:,:))/2.
   Vmean(:,n,:,:) = (va(:,n,:,:)+va(:,n-1,:,:))/2.
  end do

  g0 = 9.81
  R0 = 287.04; 8.314462
  d_orog = conform_dims(dimsizes(tas),orog,(/1,2/)) 
 
  zg0 = zg(:,0,:,:) - d_orog
printVarSummary(ps)
printVarSummary(zg0)
printVarSummary(Tmean(:,0,:,:))
printVarSummary(plev(:,0,:,:))


  plev(:,0,:,:) = ps/(exp(g0*zg0/R0/Tmean(:,0,:,:)))

  do n = 1,dim1 -1
;  print(n)
  zgn = (/zg(:,n,:,:)-zg(:,n-1,:,:)/)
  plev(:,n,:,:) = plev(:,n-1,:,:)/(exp(g0*zgn/R0/Tmean(:,n,:,:)))
  end do

  plev1 = doubletofloat(plev)
  zg1 = doubletofloat(zg)
copy_VarCoords(ta,zg1)
copy_VarCoords(ta,plev1)
  fout    = addfile(path_output+"zg_6hrLev_UKESM1-0-LL_ssp126_r1i1p1f2_gn_"+yyyy+"01010000-"+yyyy+"12311800.nc","c")      ; 6-hourly 3-d zg              
  fatt = True
  fatt@title = "calculated 3D geopotential on pressure levels"
  fatt@date = systemfunc("date")
  fileattdef(fout, fatt)
  fout->time = time
  fout->lat = lat
  fout->lat_bnds = lat_bnds
  fout->lon = lon
  fout->lon_bnds = lon_bnds
;  fout->zg = (/zg1/) ;  bottom to top
  fout->zg = zg1 ;  bottom to top

delete(fout)

  fout    = addfile(path_output+"plev_6hrLev_UKESM1-0-LL_ssp126_r1i1p1f2_gn_"+yyyy+"01010000-"+yyyy+"12311800.nc","c")      ; 6-hourly 3-d zg              
  fatt = True
  fatt@title = "calculated 3D  pressure on hybrid model levels"
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
delete(ua)
delete(va)
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
delete(d_orog)
delete(time)
delete(fout)
delete(in_ta)
delete(in_ta1)
delete(in_tas)
delete(in_hus)
delete(in_ps)
delete(in_ua)
delete(in_va)
delete(in_orog)


end do


end
