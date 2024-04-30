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



;************************************************
; file handling
;************************************************
  print("opening files")
  in_ta    = addfile("orig/1950/ta_6hrLev_UKESM1-0-LL_historical_r1i1p1f2_gn_195001010600-195101010000.nc","r")      ; 6-hourly 3-d T              
  in_ps    = addfile("orig/1950/ps_6hrLev_UKESM1-0-LL_historical_r1i1p1f2_gn_195001010600-195101010000.nc","r")      ; 6-hourly 2-d ps                
  in_tas    = addfile("orig/1950/tas_3hr_UKESM1-0-LL_historical_r1i1p1f2_gn_195001010300-195101010000.nc","r")      ; 6-hourly 2-d ps                


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
  lon   = in_ta->lon                            
  lev   = in_ta->lev                            
  lev_bnds   = in_ta->lev_bnds                            
  ta   = in_ta->ta                            
  tas   = in_tas->tas(1::2,:,:) ; start from 6UTC with an interval of 6 hr                             
  printVarSummary(tas)
  ps   = in_ps->ps                            
 
  depth = new((/85/), "double") ; calculate depth between two levels
 
  depth(0) = lev(0)
  do n = 1, 84 ; from 2 level to top
  depth(n) = lev(n) - lev (n-1)
  end do 

 Tmean = new((/1440, 85, 64, 77/),"float")
 p = new((/1440, 85, 64, 77/),"double")


   Tmean(:,0,:,:) = (tas + ta(:,0,:,:))/2.
  do n =1, 84
   Tmean(:,n,:,:) = (ta(:,n,:,:)+ta(:,n-1,:,:))/2.
  end do

  g0 = 9.81
  R0 = 287.04; 8.314462

  p(:,0,:,:) = ps/(exp(g0*depth(0)/R0/Tmean(:,0,:,:)))
  do n = 1,84
  print(n)
  p(:,n,:,:) = p(:,n-1,:,:)/(exp(g0*depth(n)/R0/Tmean(:,n,:,:)))
  end do


  fout    = addfile("new/1950/p3d_6hrLev_UKESM1-0-LL_historical_r1i1p1f2_gn_195001010600-200001010000.nc","c")      ; 6-hourly 2-d ps                
  fatt = True
  fatt@title = "calculated 3D pressure on height levels"
  fatt@date = systemfunc("date")
  fileattdef(fout, fatt)
  fout->time = time
  fout->lat = lat
  fout->lon = lon
  fout->lev = lev
  fout->lev_bnds = lev_bnds
  fout->p = p

end
