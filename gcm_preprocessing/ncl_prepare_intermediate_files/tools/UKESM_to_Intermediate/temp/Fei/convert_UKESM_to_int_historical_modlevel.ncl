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
  IM_root_name = "UKESM_"+CASE
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

  in_zsfc  = addfile("/g/data/gh02/cr3059/UKESM_test/data_tmp/" +year +"/atmos_zsfc.nc","r")    ; static surface geopotential        
  in_lmask = addfile("/g/data/gh02/cr3059/UKESM_test/data_tmp/" +year +"/atmos_lmask.nc","r")   ; static land mask

  in_tas   = addfile("/g/data/gh02/cr3059/UKESM_test/data_tmp/" +year +"/atmos_tas.nc","r")     ; 3-hourly surface temp
  in_mrsos = addfile("/g/data/gh02/cr3059/UKESM_test/data_tmp/" +year +"/atmos_mrsos.nc","r")   ; 3-hourly soil moisture

  in_ta    = addfile("/g/data/gh02/cr3059/UKESM_test/data_tmp/" +year +"/atmos_ta.nc","r")      ; 6-hourly 3-d T              
  in_ua    = addfile("/g/data/gh02/cr3059/UKESM_test/data_tmp/" +year +"/atmos_ua.nc","r")      ; 6-hourly 3-d U                
  in_va    = addfile("/g/data/gh02/cr3059/UKESM_test/data_tmp/" +year +"/atmos_va.nc","r")      ; 6-hourly 3-d V                  
  in_hus   = addfile("/g/data/gh02/cr3059/UKESM_test/data_tmp/" +year +"/atmos_hus.nc","r")     ; 6-hourly 3-d Q                   

  in_ps     = addfile("/g/data/gh02/cr3059/UKESM_test/data_tmp/" +year +"/atmos_ps.nc","r")     ; 6-hourly surface pressure
  in_pfull  = addfile("/g/data/gh02/cr3059/UKESM_test/data_tmp/" +year +"/atmos_pfull.nc","r")  ; 6-hourly 3-d pressure

  in_tos   = addfile("/g/data/gh02/cr3059/UKESM_test/data_tmp/" +year +"/atmos_tos.nc","r")     ; daily SST on gcm grids instead of pop grid (gaussian)
  in_sic   = addfile("/g/data/gh02/cr3059/UKESM_test/data_tmp/" +year +"/atmos_tos.nc","r")     ; no sea ice within domain, SST used instead

  ; "b" needs to be read from uncropped files (it doesn't exist after cropping)
  ; the following file is used to extract a and b in a+b*orog
  in_ta1    = addfile("/g/data/gh02/CMIP6/N2/scripts/UKESM1-0-LL/historical/ta_6hrLev_UKESM1-0-LL_historical_r1i1p1f2_gn_195001010600-195101010000.nc","r")
  ; and this one is used to extract orog
  in_orog   = addfile("/g/data/gh02/CMIP6/N2/et_scripts/gcms/UKESM1-0-LL/preprocessed/model_level/fx/orog_fx_UKESM1-0-LL_piControl_r1i1p1f2_gn.nc","r")



;************************************************
; read needed variables from file
;************************************************
  print("read variables in")

  time   = in_ta->time                            
  lon = in_ta->lon
  lat = in_ta->lat

  T = in_ta->ta                               
  U = in_ua->ua                          
  V = in_va->va                     
  Q = in_hus->hus                    
  P = in_pfull->pfull
  plev = P 
  pnew1 = P

  printVarSummary(T)
  printVarSummary(P)
 
  LMASK= in_lmask->sftlf
  LMASK@missing = 1E+20
  LMASK = LMASK*0.01
  LMASK = where(LMASK .gt. 0.5, 1, 0)
  LMASK = where( ismissing(LMASK), 1, LMASK)
  PHIS = in_zsfc->orog
  ZSFC = PHIS  ; orog provide actual evaluation

  SST_DAY_POP   = in_tos->tos ; tos is daily data        
  ICE_DAY_POP   = in_sic->tos ; use SST instead, all gets turned to 0 below                      
  ICE_DAY_POP   = ICE_DAY_POP *0 ;

  PS   = in_ps->ps                        
  LAT  = in_ta->lat                          

; get dimensions from dummy variable
  dsizes_x = dimsizes(T)
  ntim     = dsizes_x(0)
  nlev     = dsizes_x(1)
  nlat     = dsizes_x(2)
  nlon     = dsizes_x(3)

  print("UKESM SST: interpolated to (nlat,nlon) in preprocess, converted to K and repeated to ntim (daily->6-hourly) here")
  SST_DAY = in_tos->tos
  SST_DAY = SST_DAY+273.15 ; convert from C to K
  SST6 = new((/ntim,nlat,nlon/),float) ; create a 6 hourly variable and repeat 4x daily
  SST6(::4,:,:)=SST_DAY  
  SST6(1::4,:,:)=SST_DAY
  SST6(2::4,:,:)=SST_DAY
  SST6(3::4,:,:)=SST_DAY
  SST6 = where (ismissing(SST6), -1.E30, SST6) ; set missing value to default for METGRID.TBL

  ICE_DAY = ICE_DAY_POP
  ICE6 = new((/ntim,nlat,nlon/),float) ; create a 6 hourly variable and repeat 4x daily
  ICE6(::4,:,:)=ICE_DAY
  ICE6(1::4,:,:)=ICE_DAY
  ICE6(2::4,:,:)=ICE_DAY
  ICE6(3::4,:,:)=ICE_DAY
  ICE6 = where (ismissing(ICE6), -1.E30, ICE6) ; set missing value to default for METGRID.TBL

  TSKIN = in_tas->tas(0::2,:,:)             ;3-hourly variable
  SNOWW = in_tas->tas(0::2,:,:)             ;since no snw, only getting dims here from tas and setting all 0 below 
  SNOWW = 0
  MSOIL = in_mrsos->mrsos(0::2,:,:)         ;3-hourly variable
  TSOIL = TSKIN                             ;missing - use TSKIN

; now make into 6 hourly variables that change every ~30 days
; note that a *key* assumption here is that we are working with
; data in 3 month chunks each time we run this script!!
; note that soil M & T have 15 layers.  We only need the 4 that
; most closely correspond to the NOAH LSM 0-10,10-40,40-100,100-200 cm layers.
; these were predetermined based on inspection of the layers
; they aren't very important because they get spun up anyway...

  TSKIN6 = TSKIN ; 6-hourly variable

  TSOIL1_6 = new((/ntim,nlat,nlon/),float) ; create a 6 hourly variable and repeat 4x daily  soil temperature
  TSOIL2_6 = new((/ntim,nlat,nlon/),float) ; create a 6 hourly variable and repeat 4x daily
  TSOIL3_6 = new((/ntim,nlat,nlon/),float) ; create a 6 hourly variable and repeat 4x daily
  TSOIL4_6 = new((/ntim,nlat,nlon/),float) ; create a 6 hourly variable and repeat 4x daily

  TSOIL1_6 = TSOIL
  TSOIL2_6 = TSOIL
  TSOIL3_6 = TSOIL
  TSOIL4_6 = TSOIL

  MSOIL1_6 = new((/ntim,nlat,nlon/),float) ; create a 6 hourly variable and repeat 4x daily   soil moisture
  MSOIL2_6 = new((/ntim,nlat,nlon/),float) ; create a 6 hourly variable and repeat 4x daily
  MSOIL3_6 = new((/ntim,nlat,nlon/),float) ; create a 6 hourly variable and repeat 4x daily
  MSOIL4_6 = new((/ntim,nlat,nlon/),float) ; create a 6 hourly variable and repeat 4x daily

  MSOIL1_6 = MSOIL  ; four different levels
  MSOIL2_6 = MSOIL
  MSOIL3_6 = MSOIL
  MSOIL4_6 = MSOIL

  SNOWW6 = SNOWW  

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

; calculate zg
  lev   = in_ta->lev   ; from bottom to top
  orog = in_orog->orog
  b = in_ta1->b ; from bottom to top
  zg = conform(T,lev,1) + conform(T,b,1) * conform(T,orog,(/2,3/))
  GHT = doubletofloat(zg)
  copy_VarCoords(T, GHT)

  T2  = in_tas->tas(::2,:,:)           ;3-hourly variable
  U10 = U(:,0,:,:)             ;already 6-hourly, here we assume that index 0 is == 10m, we could do it more accurately as below:
;  U10 = doubletofloat(U(:,0,:,:)*((10./(Z(:,0,:,:)-conform(pnew1,ZSFC,(/1,2/))))^0.145))
  V10 = V(:,0,:,:)             ;

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

  FIELD_Z       ="GHT"
  UNITS_Z       ="m"
  DESC_Z        ="Geopotential Height"

  FIELD_PRES       ="PRES"
  UNITS_PRES       ="Pa"
  DESC_PRES        ="3D Pressure"

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
   SLP= new((/ntim,nlat,nlon/),float) ; create a 6 hourly variable for Sea Level Pressure

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
   
  print("calculate RH") ; calculate RH on pressure level
    RonP(TIM,:,:,:) = relhum(T(TIM,:,:,:),Q(TIM,:,:,:),pnew1(TIM,:,:,:))
    RonP = where(RonP.gt.100.0, 100.0, RonP)
    RonP = where(RonP.lt.0.0, 0.0, RonP)
    R2(TIM,:,:) = RonP(TIM,0,:,:)

    ; calculate slp using ecmwf formulation (no 6-hourly or 3-hourly slp output is otherwise available)
    SLP(TIM,:,:) = doubletofloat(pslec(TBOT(TIM,:,:),PHIS,PS(TIM,:,:),P(TIM,0,:,:)))  

  opt = True
  opt@date = HDATE
  opt@map_source             = "CESM 0.9 x 1.25"
  opt@projection             = 0
  opt@startloc               = "SWCORNER"
  opt@startlon               = doubletofloat(lon(0))
  opt@startlat               = doubletofloat(lat(0))
  opt@deltalon               = doubletofloat(lon(1) - lon(0))
  opt@deltalat               = doubletofloat(lat(1) - lat(0))
  opt@is_wind_earth_relative = False

    ; write any non-level variables (including treament of soil temp/moisture profiles as single levels each)
    ; Due to the fact that other programs will read these files, and there we assume the fields to in a in given order,
    ; these order of the writes below are import for this application

     do LEV=0, nlev -1
      opt@level = nlev - LEV 
      wrf_wps_write_int(IM_root_name,FIELD_Z,UNITS_Z,DESC_Z, GHT(TIM,LEV,:,:),opt)
      wrf_wps_write_int(IM_root_name,FIELD_T,UNITS_T,DESC_T, T(TIM,LEV,:,:),opt)
      wrf_wps_write_int(IM_root_name,FIELD_U,UNITS_U,DESC_U, U(TIM,LEV,:,:),opt)
      wrf_wps_write_int(IM_root_name,FIELD_V,UNITS_V,DESC_V, V(TIM,LEV,:,:),opt)
      wrf_wps_write_int(IM_root_name,FIELD_R,UNITS_R,DESC_R, RonP(TIM,LEV,:,:),opt)
      wrf_wps_write_int(IM_root_name,FIELD_PRES,UNITS_PRES,DESC_PRES, plev(TIM,LEV,:,:),opt) ; model level pressure
    end do 

    delete(opt@level)
    opt@level = 200100.
    wrf_wps_write_int(IM_root_name,FIELD_ZSFC,UNITS_ZSFC,DESC_ZSFC, doubletofloat(ZSFC(:,:)),opt)
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

