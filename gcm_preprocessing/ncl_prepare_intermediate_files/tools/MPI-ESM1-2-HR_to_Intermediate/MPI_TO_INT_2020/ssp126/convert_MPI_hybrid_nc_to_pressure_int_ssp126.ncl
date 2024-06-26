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


;******************************************************************************
; Which CASE are we doing. 
; CASE is only used for naming conventions, but make sure this is correct, 
; else you could clobber exisiting files.
;******************************************************************************

;CASE = "20THC"
;CASE = "RCP85"
;CASE = "RCP60"
;CASE = "RCP45"
if ( .not. isvar("CASE") ) then
  print("We need to know which case we are running")
  exit
end if
  
if ( .not. isvar("IM_root_name") ) then
  IM_root_name = "NorESM2-MM_"+CASE
; IM_root_name = "CESM_CMIP5_MOAR_"+CASE
else
  IM_root_name = IM_root_name+"_"+CASE
end if


if ( .not. isvar("outDIR") ) then
  outDIR = "OUTPUT"
end if


;******************************************************************************
; Don't change anything below this line
;******************************************************************************


;************************************************
; file handling
;************************************************
  print("opening files")
  in_ta    = addfile("atmos_ta_1.nc","r")      ; 6-hourly 3-d T              
  in_ua    = addfile("atmos_ua_1.nc","r")      ; 6-hourly 3-d U                
  in_va    = addfile("atmos_va_1.nc","r")      ; 6-hourly 3-d V                  
  in_hus   = addfile("atmos_hus_1.nc","r")     ; 6-hourly 3-d Q                   
  in_zg   = addfile("atmos_zg_1.nc","r")       ; 6-hourly 3-d zg                   

  in_zsfc  = addfile("atmos_zsfc.nc","r")    ; static surface geopotential        
  in_lmask = addfile("atmos_lmask.nc","r")   ; static land mask

  in_ps    = addfile("atmos_ps_1.nc","r")      ; 3-hourly surface pressure
  in_psl    = addfile("atmos_psl_1.nc","r")      ; 6-hourly surface pressure
  in_uas    = addfile("atmos_uas_1.nc","r")      ; 6-hourly surface pressure
  in_vas    = addfile("atmos_vas_1.nc","r")      ; 6-hourly surface pressure
  in_snw   = addfile("atmos_snw_1.nc","r")     ; daily SWE 
  in_mrsol = addfile("atmos_mrsol_1.nc","r")     ; 6-hourly soil moisture
  in_tas    = addfile("atmos_tas_1.nc","r")    ; 6-hourly surface temp
  in_tsl   = addfile("atmos_tsl_1.nc","r")   ; 6-hourly soil temp
  in_tos   = addfile("atmos_tos_1.nc","r")   ; 3hrly SST on gcm grids instead of pop grid (gaussian)
  in_sic   = addfile("atmos_tos_1.nc","r")   ; 3hrly SEAICE % on POP grid (guassian)   no sea ice within domain

;************************************************
; read needed variables from file
;************************************************
  print("read variables in")

  time   = in_ta->time                            
  P0 = 1000.       ; value of P0 is 1000.(units=hPa); it is used in both Pa and hPa units below depending on function
  lon = in_ta->lon
  lat = in_ta->lat

;  T = in_ta->T(:,::-1,:,:)                               
;  U = in_ua->U(:,::-1,:,:)                          
;  V = in_va->V(:,::-1,:,:)                     
;  Q = in_hus->Q(:,::-1,:,:)                    

  T = in_ta->ta     ; data are from bottom to top                         
  U = in_ua->ua
  V = in_va->va                     
  Q = in_hus->hus
  Z = in_zg->zg
 
; LMASK= in_lmask->LANDMASK*1.0
  LMASK= in_lmask->sftlf
  LMASK@missing = 1E+20
  LMASK = LMASK*0.01
  LMASK = where(LMASK .gt. 0.5, 1, 0)
  LMASK = where( ismissing(LMASK), 1, LMASK)
;  PHIS = in_zsfc->PHIS
  PHIS = in_zsfc->orog
;  ZSFC = PHIS/9.81
  ZSFC = PHIS  ; orog provide actual evaluation

;  SST_DAY_POP   = in_tos->tos(::2,:,:) ; tos is 3hrly data                  
  ICE_DAY_POP   = in_sic->tos(::2,:,:) ; convert from % to fraction                        
  ;ICE_DAY_POP   = in_sic->aice_d*0.01 ; convert from % to fraction                        
  ICE_DAY_POP   = ICE_DAY_POP *0 ; convert from % to fraction                        

;  PS   = in_ps->PS                        
  PS   = in_ps->ps  ; ps is 6hrly 
                        
  LAT  = in_ta->lat                          
;  hyam = in_ta->hyam(::-1)   ; no need to use them as ET has convert data from model levels to pressure levels                         
;  hybm = in_ta->hybm(::-1)                      
;  hyai = in_ta->hyai(::-1)                            
;  hybi = in_ta->hybi(::-1)                      

;  hyam = in_ta->ap    ;check whether it's required to upside down.                        
;  hybm = in_ta->b                      
;  hyai = in_ta->ap_bnds                            
;  hybi = in_ta->b_bnds                      
;
;  hyam = hyam/100000.
;  hyai = hyai/100000.
  ;print(hyai)
  ;print(hybi)

; get dimensions from dummy variable
  dsizes_x = dimsizes(T)
  ntim     = dsizes_x(0)
  nlev     = dsizes_x(1)
  nlat     = dsizes_x(2)
  nlon     = dsizes_x(3)

; Interpolate SST and SEA ICE fields to CESM Atmospheric domain from POP Ocean domain
;  SST_DAY = PopLatLon(SST_DAY_POP,"gx3v4","1x1d","bilin","da","020604")
;  SST_DAY = PopLatLon(SST_DAY_POP,"gx1v6","fv0.9x1.25","bilin","da","090309")  ; work on this 
  print("interploate sst")
  SST_DAY   = in_tos->tos(::2,:,:)  ; tos is remapped on the GCM grids at 3hrly.
  SST_DAY = SST_DAY+273.15 ; convert from C to K
  SST6 = new((/ntim,nlat,nlon/),float) ; create a 6 hourly variable and repeat 4x daily
;  SST6(::4,:,:)=SST_DAY  ;  SST_DAY is already at 6 hourly 
;  SST6(1::4,:,:)=SST_DAY
;  SST6(2::4,:,:)=SST_DAY
;  SST6(3::4,:,:)=SST_DAY
  SST6 = SST_DAY
  SST6 = where (ismissing(SST6), -1.E30, SST6) ; set missing value to default for METGRID.TBL

;  ICE_DAY = PopLatLon(ICE_DAY_POP,"gx1v6","fv0.9x1.25","bilin","da","090309")
  ICE_DAY = ICE_DAY_POP 
  ICE6 = new((/ntim,nlat,nlon/),float) ; create a 6 hourly variable and repeat 4x daily
;  ICE6(::4,:,:)=ICE_DAY
;  ICE6(1::4,:,:)=ICE_DAY
;  ICE6(2::4,:,:)=ICE_DAY
;  ICE6(3::4,:,:)=ICE_DAY
  ICE6 = ICE_DAY
  ;ICE6 = where (ismissing(ICE6), 1.e+20, ICE6) ; set missing value to default for METGRID.TBL
  ICE6 = where (ismissing(ICE6), -1.E30, ICE6) ; set missing value to default for METGRID.TBL

; read in monthly variables and make into 6-hourly
  print("working on monthly variables")
  TSKIN = in_tas->tas             ;6-hourly variable
;  TSKIN = in_ts->ts             ;6-hourly variable
  SNOWW = in_snw->snw           ;dailly variable
  TSOIL = in_tsl->tsl           ;6-hourly variable
  MSOIL = in_mrsol->mrsol       ;6-hourly variable
  SOILD = in_mrsol->depth  ;daily variable
   
;convert MSOIL from kg m-2 to fraction
  MSOIL = doubletofloat(MSOIL/(1000.*conform(MSOIL,SOILD,1)))

; now make into 6 hourly variables that change every ~30 days
; note that a *key* assumption here is that we are working with
; data in 3 month chunks each time we run this script!!
; note that soil M & T have 15 layers.  We only need the 4 that
; most closely correspond to the NOAH LSM 0-10,10-40,40-100,100-200 cm layers.
; these were predetermined based on inspection of the layers
; they aren't very important because they get spun up anyway...
;  TSKIN6   = PS ; dummy
  TSKIN6 = TSKIN ; 6-hourly variable

  TSOIL1_6 = new((/ntim,nlat,nlon/),float) ; create a 6 hourly variable and repeat 4x daily  soil temperature
  TSOIL2_6 = new((/ntim,nlat,nlon/),float) ; create a 6 hourly variable and repeat 4x daily
  TSOIL3_6 = new((/ntim,nlat,nlon/),float) ; create a 6 hourly variable and repeat 4x daily
  TSOIL4_6 = new((/ntim,nlat,nlon/),float) ; create a 6 hourly variable and repeat 4x daily

TSOIL1_6 = TSOIL
TSOIL2_6 = TSOIL
TSOIL3_6 = TSOIL
TSOIL4_6 = TSOIL


;  TSOIL1_6(::4,:,:) = TSOIL(:,0,:,:) ; 4cmm
;  TSOIL1_6(1::4,:,:) = TSOIL(:,0,:,:)
;  TSOIL1_6(2::4,:,:) = TSOIL(:,0,:,:)
;  TSOIL1_6(3::4,:,:) = TSOIL(:,0,:,:)

;  TSOIL2_6(::4,:,:) = TSOIL(:,1,:,:) ;26cmm
;  TSOIL2_6(1::4,:,:) = TSOIL(:,1,:,:)
;  TSOIL2_6(2::4,:,:) = TSOIL(:,1,:,:)
;  TSOIL2_6(3::4,:,:) = TSOIL(:,1,:,:)

;  TSOIL3_6(::4,:,:) = TSOIL(:,2,:,:) ;60cmm
;  TSOIL3_6(1::4,:,:) = TSOIL(:,2,:,:)
;  TSOIL3_6(2::4,:,:) = TSOIL(:,2,:,:)
;  TSOIL3_6(3::4,:,:) = TSOIL(:,2,:,:)

;  TSOIL4_6(::4,:,:) = TSOIL(:,3,:,:)  ;1.36m
;  TSOIL4_6(1::4,:,:) = TSOIL(:,3,:,:)
;  TSOIL4_6(2::4,:,:) = TSOIL(:,3,:,:)
;  TSOIL4_6(3::4,:,:) = TSOIL(:,3,:,:)

  MSOIL1_6 = new((/ntim,nlat,nlon/),float) ; create a 6 hourly variable and repeat 4x daily   soil moisture
  MSOIL2_6 = new((/ntim,nlat,nlon/),float) ; create a 6 hourly variable and repeat 4x daily
  MSOIL3_6 = new((/ntim,nlat,nlon/),float) ; create a 6 hourly variable and repeat 4x daily
  MSOIL4_6 = new((/ntim,nlat,nlon/),float) ; create a 6 hourly variable and repeat 4x daily

;  MSOIL1_6(::4,:,:) = MSOIL(:,0,:,:) ; 4cmm
;  MSOIL1_6(1::4,:,:) = MSOIL(:,0,:,:)
;  MSOIL1_6(2::4,:,:) = MSOIL(:,0,:,:)
;  MSOIL1_6(3::4,:,:) = MSOIL(:,0,:,:)

;  MSOIL2_6(::4,:,:) = MSOIL(:,1,:,:) ;26cmm
;  MSOIL2_6(1::4,:,:) = MSOIL(:,1,:,:)
;  MSOIL2_6(2::4,:,:) = MSOIL(:,1,:,:)
;  MSOIL2_6(3::4,:,:) = MSOIL(:,1,:,:)

;  MSOIL3_6(::4,:,:) = MSOIL(:,2,:,:) ;60cmm
;  MSOIL3_6(1::4,:,:) = MSOIL(:,2,:,:)
;  MSOIL3_6(2::4,:,:) = MSOIL(:,2,:,:)
;  MSOIL3_6(3::4,:,:) = MSOIL(:,2,:,:)

;  MSOIL4_6(::4,:,:) = MSOIL(:,3,:,:)  ;1.36m
;  MSOIL4_6(1::4,:,:) = MSOIL(:,3,:,:)
;  MSOIL4_6(2::4,:,:) = MSOIL(:,3,:,:)
;  MSOIL4_6(3::4,:,:) = MSOIL(:,3,:,:)

MSOIL1_6 = MSOIL(:,0,:,:)  ; four different levels
MSOIL2_6 = MSOIL(:,1,:,:)
MSOIL3_6 = MSOIL(:,2,:,:)
MSOIL4_6 = MSOIL(:,3,:,:)


  SNOWW6   = PS ; dummy
;  TSOIL1_6 = PS ; dummy 
;  TSOIL2_6 = PS ; dummy
;  TSOIL3_6 = PS ; dummy
;  TSOIL4_6 = PS ; dummy
;  MSOIL1_6 = PS ; dummy
;  MSOIL2_6 = PS ; dummy
;  MSOIL3_6 = PS ; dummy
;  MSOIL4_6 = PS ; dummy
   
;  TSKIN6(0:119,:,:)     = conform(PS(0:119,:,:),  TSKIN(0,:,:),(/1,2/))
;  TSKIN6(120:239,:,:)   = conform(PS(120:239,:,:),TSKIN(1,:,:),(/1,2/))
;  TSKIN6(240::,:,:)     = conform(PS(240::,:,:),  TSKIN(2,:,:),(/1,2/))
   
  TAVGSFC6              = TSKIN6              ; This can be used for inland lakes
  
   SNOWW6 = SNOWW  ;dail variable 
;   SNOWW6(0::4,:,:)  = SNOWW  ;dail variable 
;   SNOWW6(1::4,:,:) =  SNOWW
;   SNOWW6(2::4,:,:) = SNOWW
;   SNOWW6(3::4,:,:) = SNOWW

;  SNOWW6(0:119,:,:)     = conform(PS(0:119,:,:),  SNOWW(0,:,:),(/1,2/))
;  SNOWW6(120:239,:,:)   = conform(PS(120:239,:,:),SNOWW(1,:,:),(/1,2/))
;  SNOWW6(240::,:,:)     = conform(PS(240::,:,:),  SNOWW(2,:,:),(/1,2/))
   
;  TSOIL1_6(0:119,:,:)   = conform(PS(0:119,:,:),  TSOIL(0,2,:,:),(/1,2/))
;  TSOIL1_6(120:239,:,:) = conform(PS(120:239,:,:),TSOIL(1,2,:,:),(/1,2/))
;  TSOIL1_6(240::,:,:)   = conform(PS(240::,:,:),  TSOIL(2,2,:,:),(/1,2/))
   
;  TSOIL2_6(0:119,:,:)   = conform(PS(0:119,:,:),  TSOIL(0,4,:,:),(/1,2/))
;  TSOIL2_6(120:239,:,:) = conform(PS(120:239,:,:),TSOIL(1,4,:,:),(/1,2/))
;  TSOIL2_6(240::,:,:)   = conform(PS(240::,:,:),  TSOIL(2,4,:,:),(/1,2/))
   
;  TSOIL3_6(0:119,:,:)   = conform(PS(0:119,:,:),  TSOIL(0,6,:,:),(/1,2/))
;  TSOIL3_6(120:239,:,:) = conform(PS(120:239,:,:),TSOIL(1,6,:,:),(/1,2/))
;  TSOIL3_6(240::,:,:)   = conform(PS(240::,:,:),  TSOIL(2,6,:,:),(/1,2/))
   
;  TSOIL4_6(0:119,:,:)   = conform(PS(0:119,:,:),  TSOIL(0,8,:,:),(/1,2/))
;  TSOIL4_6(120:239,:,:) = conform(PS(120:239,:,:),TSOIL(1,8,:,:),(/1,2/))
;  TSOIL4_6(240::,:,:)   = conform(PS(240::,:,:),  TSOIL(2,8,:,:),(/1,2/))
   
;  MSOIL1_6(0:119,:,:)   = conform(PS(0:119,:,:),  MSOIL(0,2,:,:),(/1,2/))
;  MSOIL1_6(120:239,:,:) = conform(PS(120:239,:,:),MSOIL(1,2,:,:),(/1,2/))
;  MSOIL1_6(240::,:,:)   = conform(PS(240::,:,:),  MSOIL(2,2,:,:),(/1,2/))
   
;  MSOIL2_6(0:119,:,:)   = conform(PS(0:119,:,:),  MSOIL(0,4,:,:),(/1,2/))
;  MSOIL2_6(120:239,:,:) = conform(PS(120:239,:,:),MSOIL(1,4,:,:),(/1,2/))
;  MSOIL2_6(240::,:,:)   = conform(PS(240::,:,:),  MSOIL(2,4,:,:),(/1,2/))
   
;  MSOIL3_6(0:119,:,:)   = conform(PS(0:119,:,:),  MSOIL(0,6,:,:),(/1,2/))
;  MSOIL3_6(120:239,:,:) = conform(PS(120:239,:,:),MSOIL(1,6,:,:),(/1,2/))
;  MSOIL3_6(240::,:,:)   = conform(PS(240::,:,:),  MSOIL(2,6,:,:),(/1,2/))
   
;  MSOIL4_6(0:119,:,:)   = conform(PS(0:119,:,:),  MSOIL(0,8,:,:),(/1,2/))
;  MSOIL4_6(120:239,:,:) = conform(PS(120:239,:,:),MSOIL(1,8,:,:),(/1,2/))
;  MSOIL4_6(240::,:,:)   = conform(PS(240::,:,:),  MSOIL(2,8,:,:),(/1,2/))

; Set missing values to default for METGRID.TBL
  TSOIL1_6 = where (ismissing(TSOIL1_6), -1.E30, TSOIL1_6)
  TSOIL2_6 = where (ismissing(TSOIL2_6), -1.E30, TSOIL2_6)
  TSOIL3_6 = where (ismissing(TSOIL3_6), -1.E30, TSOIL3_6)
  TSOIL4_6 = where (ismissing(TSOIL4_6), -1.E30, TSOIL4_6)

  MSOIL1_6 = where (ismissing(MSOIL1_6), -1.E30, MSOIL1_6)
  MSOIL2_6 = where (ismissing(MSOIL2_6), -1.E30, MSOIL2_6)
  MSOIL3_6 = where (ismissing(MSOIL3_6), -1.E30, MSOIL3_6)
  MSOIL4_6 = where (ismissing(MSOIL4_6), -1.E30, MSOIL4_6)

; Get bottom level temperature for vertical extrapolation
  TBOT = T(:,0,:,:) 

; Calculate geopotential height
;  print ("calculate Z")
;  TV = T
;  TV = T*(1.+0.61*Q)
;  Z = T
;  Z = cz2ccm(PS,PHIS,TV(:,::-1,:,:),P0*100.,hyam,hybm,hyai,hybi)
;  Z = Z(:,::-1,:,:) ;reorder vertical so consistent with T,U,V, and Q (bottom up)

;  Z = in_zg->zg(:,::-1,:,:)  ; 6 hourly variable is available                              

; Calculate the pressures on each hybrid level (bottom up)
;  print ("calculate P")
;  P = T
;  P = pres_hybrid_ccm(PS,P0*100.,hyam,hybm) ; bottom to top
;  print(P(0,:,50,50))
; Uses stock RH function in ncl
;  print("calculate RH")
;  R = T
;  R = relhum(T,Q,P)
;  R = where(R.gt.100.0, 100.0, R)

; Calculate near surface variables (do this instead of take from files do to greater availability of 3d)
; linearly extrapolate T2 in log(P) from 2 lowest hybrid levels)
; Assume RH at 2 m (R2) is same as RH at lowest hybrid level
; calculate Q2 from R2 and T2
; assume U10 and V10 follow a power law for a neutrally stable atmosphere and
; extrapolate them from the lowest hybrid level. Alpha = 0.143
; Hsu, S.A., et al, 1994,J. Appl. Meteor., Vol. 33, pp. 757-765

; create an array of desired pressure levels:
  pnew = (/ 1000.0, 975.0, 950.0, 925.0, 900.0, 875.0, 850.0, 825.0, 800.0, 775.0, 750.0, 700.0, 650.0, 600.0, 550.0, 500.0, \
             450.0, 400.0, 350.0, 300.0, 250.0, 225.0, 200.0, 175.0, 150.0, 125.0, 100.0, 50.0 /)

  P = conform(T, pnew*100, 1)
  print("calculate near surface variables")
 ; T2 = PS
 ; T2 = T(:,1,:,:)-(T(:,1,:,:)-T(:,0,:,:))*((log(PS)-log(P(:,1,:,:)))/(log(P(:,0,:,:))-log(P(:,1,:,:))))
 ; T2 is from tas

  T2  = in_tas->tas             ;6-hourly variable
  U10 = in_uas->uas             ;6-hourly variable
  V10 = in_vas->vas             ;6-hourly variable
  SLP = in_psl->psl                       ;6-hourly variable 
  ;U10 = doubletofloat(U(:,0,:,:)*((10./(Z(:,0,:,:)-conform(PS,ZSFC,(/1,2/))))^0.145))
  ;V10 = doubletofloat(V(:,0,:,:)*((10./(Z(:,0,:,:)-conform(PS,ZSFC,(/1,2/))))^0.145))

;************************************************
; define other arguments required by vinth2p
;************************************************
; type of interpolation: 1 = linear, 2 = log, 3 = loglog
  interp = 2 

; is extrapolation desired if data is outside the range of PS
  extrap = True
   

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; WRITE OUT DATA TO WRF INTERMEDIATE FORMAT (CALL FORTRAN SUBROUTINE VIA WRAPIT)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  print("write out data to WRF intermediate files")

; Set the IM header information
; These are all passed to the write script as resources

; set the format to use for the date
  date_format  = "%Y-%N-%D_%H:00:0000000"        		; 24 chars exact for binary file
  date_format2 = "%Y-%N-%D_%H"                                  ; for filename
  yyyy_format  = "%Y" 						; for output directory


; Set up field specific header information
  FIELD_T       ="TT"
  UNITS_T       ="K"
  DESC_T        ="Temperature"

  FIELD_U       ="UU"
  UNITS_U       ="m s-1"
  DESC_U        ="Zonal Wind Speed"

  FIELD_V       ="VV"
  UNITS_V       ="m s-1"
  DESC_V        ="Meridional Wind Speed"

  FIELD_Q       ="SPECHUMD"
  UNITS_Q       ="kg kg-1"
  DESC_Q        ="Specific Humidity"

  FIELD_R       ="RH"
  UNITS_R       ="%"
  DESC_R        ="Relative Humidity"

  FIELD_Z       ="HGT"
  UNITS_Z       ="m"
  DESC_Z        ="Geopotential Height"

  FIELD_SLP     ="PMSL"
  UNITS_SLP     ="Pa"
  DESC_SLP      ="Sea level pressure"
  XLVL_SLP      = 201300.

  FIELD_PS      ="PSFC"
  UNITS_PS      ="Pa"
  DESC_PS       ="Surface pressure"

  FIELD_U10     ="UU"
  UNITS_U10     ="m s-1"
  DESC_U10      ="10 m Zonal Wind Speed"

  FIELD_V10     ="VV"
  UNITS_V10     ="m s-1"
  DESC_V10      ="10 m Meridional Wind Speed"

  FIELD_T2      ="TT"
  UNITS_T2      ="K"
  DESC_T2       ="2 m Temperature"

  FIELD_Q2      ="SPECHUMD"
  UNITS_Q2      ="kg kg-1"
  DESC_Q2       ="2 m Specific Humidity"

  FIELD_R2      ="RH"
  UNITS_R2      ="%"
  DESC_R2       ="2 m Relative Humidity"

  FIELD_SNOWW6  ="SNOW"
  UNITS_SNOWW6  ="kg m-2"
  DESC_SNOWW6   ="Water Equivalent of Accum Snow Depth"

  FIELD_TSKIN6  ="SKINTEMP"
  UNITS_TSKIN6  ="K"
  DESC_TSKIN6   ="Skin Temperature"

  FIELD_TAVGSFC6="TAVGSFC"
  UNITS_TAVGSFC6="K"
  DESC_TAVGSFC6 ="Daily mean of surface air temperature"

  FIELD_SST6    ="SST"
  UNITS_SST6    ="K"; 25 chars exact
  DESC_SST6     ="Sea Surface Temperature"

  FIELD_ICE6    ="SEAICE"
  UNITS_ICE6    ="fraction"
  DESC_ICE6     ="Sea-Ice-Fraction"

  FIELD_TSOIL1_6="ST000010"
  UNITS_TSOIL1_6="K"
  DESC_TSOIL1_6 ="Soil Temperature 0-10 cm layer"

  FIELD_TSOIL2_6="ST010040"
  UNITS_TSOIL2_6="K"
  DESC_TSOIL2_6 ="Soil Temperature 10-40 cm layer"

  FIELD_TSOIL3_6="ST040100"
  UNITS_TSOIL3_6="K"
  DESC_TSOIL3_6 ="Soil Temperature 40-100 cm layer"

  FIELD_TSOIL4_6="ST100200"
  UNITS_TSOIL4_6="K"
  DESC_TSOIL4_6 ="Soil Temperature 100-200 cm layer"

  FIELD_MSOIL1_6="SM000010"
  UNITS_MSOIL1_6="fraction"
  DESC_MSOIL1_6 ="Soil Moisture 0-10 cm layer"

  FIELD_MSOIL2_6="SM010040"
  UNITS_MSOIL2_6="fraction"
  DESC_MSOIL2_6 ="Soil Moisture 10-40 cm layer"

  FIELD_MSOIL3_6="SM040100"
  UNITS_MSOIL3_6="fraction"
  DESC_MSOIL3_6 ="Soil Moisture 40-100 cm layer"

  FIELD_MSOIL4_6="SM100200"
  UNITS_MSOIL4_6="fraction"
  DESC_MSOIL4_6 ="Soil Moisture 100-200 cm layer"

  FIELD_LMASK   ="LANDSEA"
  UNITS_LMASK   ="proprtn"
  DESC_LMASK    ="Land/Sea Flag; 0=Ocean; 1=Land"

  FIELD_ZSFC    ="SOILHGT"
  UNITS_ZSFC    ="m"
  DESC_ZSFC     ="Terrain Elevation"

   RonP= new((/ntim,nlev,nlat,nlon/),float) ; create a 6 hourly variable for R relative humidity
   R2= new((/ntim,nlat,nlon/),float) ; create a 6 hourly variable for  relative humidity
;   U10= new((/ntim,nlat,nlon/),float) ; create a 6 hourly variable for U10 
;   V10= new((/ntim,nlat,nlon/),float) ; create a 6 hourly variable for V10
;   SLP= new((/ntim,nlat,nlon/),float) ; create a 6 hourly variable for Sea Level Pressure

; Loop through each time period - do final calculations and write data to the IM format
  do TIM = 0,ntim-1 ; loop to export data for each time step 

    HDATE = cd_string(time(TIM), date_format)                          ; set dates
    OFILE = IM_root_name+":" + cd_string(time(TIM), date_format2)      ; figure the output file name so that we can clobber it
    yyyy  = cd_string(time(TIM), yyyy_format)                          ; figure out the year so that we can make directories
    print("HDATE= " + HDATE)
    system ("mkdir -p " + outDIR + "/" + yyyy)	               	       ; make the output directories
    system ("rm " + OFILE)                                             ; delete the file, so that we don't append to an existing file
    system ("rm " + outDIR + "/" +yyyy+ "/" + OFILE)                   ; delete the file, so that we don't append to an existing file

    OFILE1 = "3Dtempo"+":" + cd_string(time(TIM), date_format2)      ; figure the output file name so that we can clobber it
    OFILE2 = "SFCtempo"+":" + cd_string(time(TIM), date_format2)      ; figure the output file name so that we can clobber it
   

    ; calculate 3-d variables on pressure levels (need to do this here because ncl doesn't like multiple times)
    print("extrapolate variables to pressure levels")
    varflg = 1
;    TonP = vinth2p_ecmwf(T(TIM,::-1,:,:),hyam(::-1),hybm(::-1),pnew(::-1),PS(TIM,:,:),interp,P0,1,extrap,varflg,TBOT(TIM,:,:),PHIS)  ; top to bottom
;    varflg = -1
;    ZonP = vinth2p_ecmwf(Z(TIM,::-1,:,:),hyam(::-1),hybm(::-1),pnew,PS(TIM,:,:),interp,P0,1,extrap,varflg,TBOT(TIM,:,:),PHIS)
     ZonP = Z(TIM,:,:,:)
;    varflg = 0
   

    ; calculate 3-d variables on pressure levels (need to do this here because ncl doesn't like multiple times)
    print("extrapolate variables to pressure levels")
    varflg = 1
;    TonP = vinth2p_ecmwf(T(TIM,::-1,:,:),hyam(::-1),hybm(::-1),pnew(::-1),PS(TIM,:,:),interp,P0,1,extrap,varflg,TBOT(TIM,:,:),PHIS)  ; top to bottom
;    varflg = -1
;    ZonP = vinth2p_ecmwf(Z(TIM,::-1,:,:),hyam(::-1),hybm(::-1),pnew,PS(TIM,:,:),interp,P0,1,extrap,varflg,TBOT(TIM,:,:),PHIS)
     ZonP = Z(TIM,:,:,:)
;    varflg = 0
;    UonP = vinth2p_ecmwf(U(TIM,::-1,:,:),hyam(::-1),hybm(::-1),pnew(::-1),PS(TIM,:,:),interp,P0,1,extrap,varflg,TBOT(TIM,:,:),PHIS) ; top to bottom
;    VonP = vinth2p_ecmwf(V(TIM,::-1,:,:),hyam(::-1),hybm(::-1),pnew(::-1),PS(TIM,:,:),interp,P0,1,extrap,varflg,TBOT(TIM,:,:),PHIS) ; top to bottom
;    RonP = vinth2p_ecmwf(R(TIM,::-1,:,:),hyam(::-1),hybm(::-1),pnew,PS(TIM,:,:),interp,P0,1,extrap,varflg,TBOT(TIM,:,:),PHIS)
    
;    TonP = TonP(::-1,:,:) ; bottom to top
;    UonP = UonP(::-1,:,:) ; bottom to top
;    VonP = VonP(::-1,:,:) ; bottom to top

  print("calculate RH") ; calculate RH on pressure level
    pnew1 = conform(T(TIM,:,:,:), pnew, 0)
    pnew1 = pnew1*100
    RonP(TIM,:,:,:) = relhum(T(TIM,:,:,:),Q(TIM,:,:,:),pnew1)
    RonP = where(RonP.gt.100.0, 100.0, RonP)
    RonP = where(RonP.lt.0.0, 0.0, RonP)

     R2(TIM,:,:) = RonP(TIM,0,:,:)
;    U10 = doubletofloat(U(:,0,:,:)*((10./(Z(:,0,:,:)-conform(pnew1,ZSFC,(/1,2/))))^0.145))
;     V10 = doubletofloat(V(:,0,:,:)*((10./(Z(:,0,:,:)-conform(pnew1,ZSFC,(/1,2/))))^0.145))
;      U10(TIM,:,:) = U(TIM,0,:,:)
;      V10(TIM,:,:) = V(TIM,0,:,:)
;      HUSS = Q(TIM,0,:,:)

    ; calcualte slp using ecmwf formulation (no 6-hourly or 3-hourly slp output is otherwise available)
;    SLP = doubletofloat(pslec(TBOT(TIM,:,:),PHIS,PS(TIM,:,:),P(TIM,0,:,:)))  
     SLP(TIM,:,:) =in_psl->psl(TIM,:,:) 

  opt = True
  opt@date = HDATE
  opt@map_source             = "CESM 0.9 x 1.25"
  opt@projection             = 0
  opt@startloc               = "SWCORNER"
  opt@startlon               = doubletofloat(lon(0))
  opt@startlat               = doubletofloat(lat(0))
 ; opt@deltalon               = doubletofloat(lon(1) - lon(0))
 ; opt@deltalat               = doubletofloat(lat(1) - lat(0))
  opt@deltalon               = (doubletofloat(lon(116)) - doubletofloat(lon(0)))/116
  opt@deltalat               = (doubletofloat(lat(84)) - doubletofloat(lat(0)))/84
  opt@is_wind_earth_relative = False



    ; write any non-level variables (including treament of soil temp/moisture profiles as single levels each)
    ; Due to the fact that other programs will read these files, and there we assume the fields to in a in given order,
    ; these order of the writes below are import for this application

    do LEV=0,dimsizes(pnew)-1
      opt@level = pnew(LEV)*100.
      wrf_wps_write_int(IM_root_name,FIELD_Z,UNITS_Z,DESC_Z, Z(TIM,LEV,:,:),opt)
      wrf_wps_write_int(IM_root_name,FIELD_T,UNITS_T,DESC_T, T(TIM,LEV,:,:),opt)
      wrf_wps_write_int(IM_root_name,FIELD_U,UNITS_U,DESC_U, U(TIM,LEV,:,:),opt)
      wrf_wps_write_int(IM_root_name,FIELD_V,UNITS_V,DESC_V, V(TIM,LEV,:,:),opt)
      wrf_wps_write_int(IM_root_name,FIELD_R,UNITS_R,DESC_R, RonP(TIM,LEV,:,:),opt)
    end do 

      print(Z(TIM, 2, 50, 50))

    opt@level = 200100.
    wrf_wps_write_int(IM_root_name,FIELD_ZSFC,UNITS_ZSFC,DESC_ZSFC, doubletofloat(ZSFC(:,:)),opt)
;    wrf_wps_write_int(IM_root_name,FIELD_LMASK,UNITS_LMASK,DESC_LMASK, LMASK(:,:),opt)
    wrf_wps_write_int(IM_root_name,FIELD_T2,UNITS_T2,DESC_T2, T2(TIM,:,:),opt) 	         
    wrf_wps_write_int(IM_root_name,FIELD_R2,UNITS_R2,DESC_R2, R2(TIM,:,:),opt) 	       
    wrf_wps_write_int(IM_root_name,FIELD_U10,UNITS_U10,DESC_U10, U10(TIM,:,:),opt) 	           
    wrf_wps_write_int(IM_root_name,FIELD_V10,UNITS_V10,DESC_V10, V10(TIM,:,:),opt) 	          
    wrf_wps_write_int(IM_root_name,FIELD_PS,UNITS_PS,DESC_PS, PS(TIM,:,:),opt) 	            

    opt@level = XLVL_SLP	; Need a different level for SLP
    wrf_wps_write_int(IM_root_name,FIELD_SLP,UNITS_SLP,DESC_SLP, SLP(TIM,:,:),opt)                      

    opt@level = 200100.
    wrf_wps_write_int(IM_root_name,FIELD_SST6,UNITS_SST6,DESC_SST6, SST6(TIM,:,:),opt)    
    wrf_wps_write_int(IM_root_name,FIELD_TSKIN6,UNITS_TSKIN6,DESC_TSKIN6, TSKIN6(TIM,:,:),opt)        
    wrf_wps_write_int(IM_root_name,FIELD_SNOWW6,UNITS_SNOWW6,DESC_SNOWW6, SNOWW6(TIM,:,:),opt)  

;    wrf_wps_write_int(IM_root_name,FIELD_TAVGSFC6,UNITS_TAVGSFC6,DESC_TAVGSFC6, TAVGSFC6(TIM,:,:),opt) 
;    wrf_wps_write_int(IM_root_name,FIELD_ICE6,UNITS_ICE6,DESC_ICE6, ICE6(TIM,:,:),opt)   
    wrf_wps_write_int(IM_root_name,FIELD_TSOIL1_6,UNITS_TSOIL1_6,DESC_TSOIL1_6, TSOIL1_6(TIM,:,:),opt) 
    wrf_wps_write_int(IM_root_name,FIELD_TSOIL2_6,UNITS_TSOIL2_6,DESC_TSOIL2_6, TSOIL2_6(TIM,:,:),opt) 
    wrf_wps_write_int(IM_root_name,FIELD_TSOIL3_6,UNITS_TSOIL3_6,DESC_TSOIL3_6, TSOIL3_6(TIM,:,:),opt) 
    wrf_wps_write_int(IM_root_name,FIELD_TSOIL4_6,UNITS_TSOIL4_6,DESC_TSOIL4_6, TSOIL4_6(TIM,:,:),opt) 
    wrf_wps_write_int(IM_root_name,FIELD_MSOIL1_6,UNITS_MSOIL1_6,DESC_MSOIL1_6, MSOIL1_6(TIM,:,:),opt) 
    wrf_wps_write_int(IM_root_name,FIELD_MSOIL2_6,UNITS_MSOIL2_6,DESC_MSOIL2_6, MSOIL2_6(TIM,:,:),opt) 
    wrf_wps_write_int(IM_root_name,FIELD_MSOIL3_6,UNITS_MSOIL3_6,DESC_MSOIL3_6, MSOIL3_6(TIM,:,:),opt) 
    wrf_wps_write_int(IM_root_name,FIELD_MSOIL4_6,UNITS_MSOIL4_6,DESC_MSOIL4_6, MSOIL4_6(TIM,:,:),opt) 

; fout = addfile( "input.nc","c")
;fatt = True
;fatt@title = "outputs for intermediate files"
;fatt@date = systemfunc("date")
;fileattdef(fout, fatt)
;fout->u = U(TIM,:,:,:)
;fout->v = V(TIM,:,:,:)
;fout->t = T(TIM,:,:,:)
;fout->rh = RonP(TIM,:,:,:)
;fout->z = Z(TIM,:,:,:)
;fout->u10 = U10(TIM,:,:)
;fout->v10 = V10(TIM,:,:)
;fout->t2 = T2(TIM,:,:)
;fout->rh2 = R2(TIM,:,:)
;fout->tsoil1_6 = TSOIL1_6(TIM,:,:)
;fout->msoil1_6 = MSOIL1_6(TIM,:,:)


    ; move the file to the outDIR directory
      print ("mv " + OFILE + "   " + outDIR + "/" +yyyy+ "/.")
      system ("mv " + OFILE + "   " + outDIR + "/" +yyyy+ "/.")
      print("  ")

     delete(HDATE)
     delete(OFILE)
     delete(opt)
  end do 

end

