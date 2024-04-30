;******************************************************************************
; NCL script to convert NorEMS2 (ds316.0) data to Intermediate file format
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
 year = (/"2100","2101","2102","2103","2104","2105","2106","2107","2108","2109","2020"/)

 path1 = "/g/data/gh02/CMIP6/N2/et_scripts/gcms/UKESM1-0-LL/preprocessed/model_level/ssp370/"
 path2 = "tas/"

 do nn = 0, 9
 yyyy1 = year(nn)
 yyyy2 = year(nn+1)

;************************************************
; file handling
;************************************************
  print("opening files")
  in_ta    = addfile(path1+"ta_6hrLev_UKESM1-0-LL_ssp370_r1i1p1f2_gn_"+yyyy1+"01010600-"+yyyy2+"01010000.nc","r")      ; 6-hourly 3-d T              
  in_hus    = addfile(path1+"hus_6hrLev_UKESM1-0-LL_ssp370_r1i1p1f2_gn_"+yyyy1+"01010600-"+yyyy2+"01010000.nc","r")      ; 6-hourly 3-d T              
  in_ps   = addfile(path2+"ps_6hrLev_UKESM1-0-LL_ssp370_r1i1p1f2_gn_"+yyyy1+"01010600-"+yyyy2+"01010000.nc","r")      ; 6-hourly 3-d T              
  in_tas   = addfile(path2+"tas_3hr_UKESM1-0-LL_ssp370_r1i1p1f2_gn_"+yyyy1+"01010300-"+yyyy2+"01010000.nc","r")      ; 6-hourly 3-d T              
  in_delta = addfile("/scratch/gh02/eut599/pdiff_test/delta_z.nc","r")

;************************************************
; Hypsometric eqestion for calculating 3D ps
;************************************************
; p2 = p1/exp(g(z2-z1)/R/Tv)


;************************************************
; read needed variables from file
;************************************************
  print("read variables in")

  time   = in_ta->time                            
  lat   = in_ta->lat                            
  lat_bnds   = in_ta->lat_bnds                            
  lon   = in_ta->lon                            
  lon_bnds   = in_ta->lon_bnds                            
  lev   = in_ta->lev                            
  lev_bnds   = in_ta->lev_bnds                            
  ta   = in_ta->ta                            
  tas   = in_tas->tas(1::2,:,:) ; start from 6UTC with an interval of 6 hr                             
  printVarSummary(tas)
  ps   = in_ps->ps                            
  hus   = in_hus->hus                            
  delta   = in_delta->delta_z                            
 
  depth = new((/85, 64,77/), "double") ; calculate depth between two levels
  lev1 = conform_dims(dimsizes(delta), lev, (/0/))+delta
  printVarSummary(lev1)
  
;  depth(0,:,:) =  lev(0)
  depth(0,:,:) =  20
  ;depth(0,:,:) =  lev1(0,:,:)
  do n1 = 1, 84 ; from 2 level to top
  depth(n1,:,:) = lev1(n1,:,:) - lev1(n1-1,:,:)
  end do 

  depth1 = conform_dims(dimsizes(ta),depth,(/1,2,3/))
printVarSummary(depth1)

 Tmean = new((/1440, 85, 64, 77/),"float")
 p = new((/1440, 85, 64, 77/),"double")


   Tmean(:,0,:,:) = (tas + ta(:,0,:,:))*(1+0.63*hus(:,0,:,:))/2.
  do n2 =1, 84
   Tmean(:,n2,:,:) = (ta(:,n2,:,:)*(1+0.63*hus(:,n2,:,:))+ta(:,n2-1,:,:)*(1+0.63*hus(:,n2-1,:,:)))/2.
  end do

  g0 = 9.81
  R0 = 287.04; 8.314462

  p(:,0,:,:) = ps/(exp(g0*depth1(:,0,:,:)/R0/Tmean(:,0,:,:)))
  do n3 = 1,84
  p(:,n3,:,:) = p(:,n3-1,:,:)/(exp(g0*depth1(:,n3,:,:)/R0/Tmean(:,n3,:,:)))
  end do

  p1 = doubletofloat(p)

  copy_VarCoords(ta, p1)
  fout    = addfile("pfull_6hrLev_UKESM1-0-LL_ssp370_r1i1p1f2_gn_"+yyyy1+"01010600-"+yyyy2+"01010000.nc","c")      ; 6-hourly 2-d ps                
  fatt = True
  fatt@title = "calculated 3D pressure on height levels"
  fatt@date = systemfunc("date")
  fileattdef(fout, fatt)
  fout->time = time
  fout->lat = lat
  fout->lat_bnds = lat_bnds
  fout->lon = lon
  fout->lon_bnds = lon_bnds
  fout->lev = lev
  fout->lev_bnds = lev_bnds
  fout->pfull = p1
delete(fout)
delete(fatt)
delete(in_ta)
delete(in_hus)
delete(in_delta)
delete(in_ps)
delete(in_tas)
delete(time)
delete(lev)
delete(ta)
delete(tas)
delete(ps)
delete(hus)
delete(lev1)
delete(depth)
delete(depth1)
delete(p)
delete(Tmean)
delete(lat)
delete(lon)
end do

end
