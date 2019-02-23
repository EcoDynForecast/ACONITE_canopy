module aconite_function_mod


	public :: acm, mass_resp, daylength, GDD_calc,area_resp,atkin_resp,ryan_resp, clm_resp

contains



!-------PHENOLOGY GROWING DEGREE DAY CALCULATOR
double precision function GDD_calc(tmax,tmin,GDD,doy)
	double precision, intent(in):: tmax, tmin, GDD
	integer, intent(in):: doy
	double precision :: tav
	
    tav = (tmax + tmin)/2.
    if(doy == 1) then
      GDD_calc = 0
    elseif(tav - 8 > 0) then
      GDD_calc = GDD + tav - 8
    endif
  End Function
!-------------------------------------------------
  
  
  
!--------DAY LENGTH CALCULATOR
double precision function  daylength(Lat,doy)
	double precision, parameter :: pi = 3.14159265359
	double precision, intent(in) :: Lat
	integer, intent(in):: doy
	double precision :: LatRad, r, z, decl,z2,h,TA,AC
	
    LatRad = Lat * (2.0 * pi) / 360.0
    r = 1 - (0.0167 * COS(0.0172 * (doy - 3)))
    z = 0.39785 * sin(4.868961 + 0.017203 * doy + 0.033446 *SIN(6.224111 + 0.017202 * doy))
    if (ABS(z) < 0.7) then
      decl = ATAN(z / (sqrt(1.0 - (z**2))))
    else
      decl = pi / 2.0 - atan(sqrt(1 - z**2) / z)
    endif
    
    if (abs(LatRad) >= (pi/2.0)) then
      if (Lat < 0) then
        LatRad = (-1.0) * (pi/2.0 - 0.01);
      else
        LatRad = (1.0) * (pi/2.0 - 0.01);
      endif
    endif
    
    z2 = -TAN(decl) * TAN(LatRad);
    if (z2 >= 1.0) then
      h = 0;
    elseif (z2 <= -1.0) then
      h = pi;
    else
      TA = abs(z2);
      if (TA < 0.7)then
        AC = 1.570796 - atan(TA / sqrt(1.0 - (TA**2)));
      else
        AC = atan(sqrt(1 - (TA**2)) / TA);
      endif
      if (z2 < 0) then
        h = pi-AC;
      else
        h = AC;
      endif
    endif
    daylength = 2.0 * (h*24.0) / (2.0*pi);
End function
  
!-------  ACM PHOTOSYNTHESIS MODEL
  double precision function acm(N,lai,tmax,tmin,DL,rad,co2,param,psid,rtot,leafgrow_doy,leaf_drop_occured,doy)
  	double precision, intent(in) :: N, lai, tmax,tmin,DL,rad,co2, param(16), psid,rtot, &
  									leafgrow_doy,leaf_drop_occured
  	double precision ::  tav,gs,pp,trange,qq,ci,e0,cps 
  		integer, intent(in):: doy
  
    tav = (tmax + tmin)/2.
    if(tmin > param(15) .AND. doy >= leafgrow_doy .AND. doy < leaf_drop_occured) then  
      trange = tmax - tmin ! daily temperature range (C)
      gs = abs(psid)**param(10) / ((param(6)* rtot + 0.5*trange))
      pp = (N * param(1) * exp(param(8) * tmax))/gs
      qq = param(3) - param(4)
      ci = 0.5 * (co2 + qq - pp + sqrt((co2 + qq - pp)**2 - 4 * (co2 * qq - pp * param(3))))
      e0 = param(7) * lai**2/(lai**2 + param(9))
      cps = (e0* rad * gs * (co2 - ci))/(e0 * rad + gs * (co2-ci))
      acm = cps * (param(2) * DL + param(5))
    else
      acm = 0.
    endif
    


end function
  
!-------NON-LINEAR RESPIRATION MODEL (REICH ET AL. 2008): MASS
  double precision function  mass_resp(tissueC,tissueN,resp_parm1,resp_parm2,resp_parm3,twq) 
  	double precision, intent(in) ::   tissueC,tissueN,resp_parm1,resp_parm2,twq,resp_parm3 
  	double precision :: total_biomass,mmolesN_per_biomass,nmolesC_per_biomass_per_sec, &
  						gC_per_biomass_per_sec,gC_per_sec
  						
    if(tissueC > 0.0 .AND. tissueN > 0.0) then
      total_biomass = tissueC/0.5
      mmolesN_per_biomass = ((tissueN/14.)*1000)/total_biomass
      
      !nmolesC_per_biomass_per_sec =exp(mass_param1+mass_param2*log(mmolesN_per_biomass))
      nmolesC_per_biomass_per_sec = resp_parm1 + log10(mmolesN_per_biomass)*resp_parm2
      gC_per_biomass_per_sec =((10**nmolesC_per_biomass_per_sec)/10**9)*12.
      gC_per_sec = gC_per_biomass_per_sec * total_biomass
      mass_resp = gC_per_sec * (60.*60.*24.) 
      !print *, mmolesN_per_biomass,10**nmolesC_per_biomass_per_sec
      !nmolesC_per_biomass_per_sec = resp_parm1 * (mmolesN_per_biomass ** resp_parm2)
      !gC_per_biomass_per_sec =(nmolesC_per_biomass_per_sec/10**9)*12.
      !gC_per_sec = gC_per_biomass_per_sec * total_biomass
      !mass_resp = gC_per_sec * (60.*60.*24.) 
      !print *, mmolesN_per_biomass,nmolesC_per_biomass_per_sec
    else
      mass_resp = 0.0
    endif
end function

!-------NON-LINEAR RESPIRATION MODEL (REICH ET AL. 2008) :AREA
  double precision function  area_resp(LAI,tissueN,resp_parm1,resp_parm2,resp_parm3,twq) 
  	double precision, intent(in) ::   LAI,tissueN,resp_parm1,resp_parm2,twq,resp_parm3
  	double precision :: total_biomass,N_per_area,gC_per_area_per_sec, &
  						gC_per_sec
  						
    if(LAI > 0.0 .AND. tissueN > 0.0) then
       N_per_area = (tissueN/LAI)
       umolC_per_area_per_sec = ((resp_parm1-resp_parm3*twq)* (N_per_area ** resp_parm2))
       gC_per_area_per_sec = (umolC_per_area_per_sec/10**6)*12.
       gC_per_sec = gC_per_area_per_sec*LAI
       area_resp = gC_per_sec * (60.*60.*24.)
    else
      area_resp = 0.0
    endif
end function

!-------AKTIN EQUATION 
  double precision function  atkin_resp(LAI,tissueN,twq,resp_parm1,resp_parm2,resp_parm3) 
  	double precision, intent(in) ::   LAI,tissueN,twq,resp_parm1,resp_parm2,resp_parm3
  	double precision :: total_biomass,N_per_area,gC_per_area_per_sec, &
  						gC_per_sec
  						
    if(tissueC > 0.0 .AND. tissueN > 0.0) then
       N_per_area = (tissueN/LAI)
       umolC_per_area_per_sec  = resp_parm1 + (N_per_area * resp_parm2) - resp_parm3*twq 
       gC_per_area_per_sec = (umolC_per_area_per_sec/10**6)*12.
       gC_per_sec = gC_per_area_per_sec*LAI
       atkin_resp = gC_per_sec * (60.*60.*24.)
       atkin_resp = max(atkin_resp,0.0)
    else
      atkin_resp = 0.0
    endif
end function

!-------RYAN EQUATION
  double precision function  ryan_resp(tissueN,tissueC,resp_parm1) 
  	double precision, intent(in) ::   tissueN,tissueC, resp_parm1 
  	double precision ::  gC_per_sec,mol_C,mmol_N,mmolN_per_molC,mmol_per_mol_per_hr, &
  		mmol_per_mol_per_second, mmol_per_second
  						
    if(tissueC > 0.0 .AND. tissueN > 0.0) then
       mol_C = tissueC / 12.0
       mmol_N = (tissueN/14.0)*1000
       mmolN_per_molC = mmol_N/mol_C
 
       mmol_per_mol_per_hr = resp_parm1 * mmolN_per_molC
       mmol_per_mol_per_second = mmol_per_mol_per_hr/(60.*60.)
       mmol_per_second = mmol_per_mol_per_second*mol_C
       gC_per_sec = (mmol_per_second*12.0)/1000
       ryan_resp = gC_per_sec * (60.*60.*24.)
    else
      ryan_resp = 0.0
    endif

end function

!-------CLM EQUATION
  double precision function  clm_resp(tissueN,tissueC,resp_parm1) 
  	double precision, intent(in) ::   tissueN,tissueC, resp_parm1 
  	double precision ::  gC_per_sec,mol_C,mmol_N,mmolN_per_molC,mmol_per_mol_per_hr, &
  		mmol_per_mol_per_second, mmol_per_second

    if(tissueC > 0.0 .AND. tissueN > 0.0) then
       gC_per_sec = TissueN * resp_parm1
       clm_resp = gC_per_sec * (60.*60.*24.)
    else
       clm_resp = 0.0
    endif
           
end function




end module




