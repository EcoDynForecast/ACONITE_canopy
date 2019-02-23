subroutine LAI_N_Rinterface(LAI_N_dataframe,par,LAI_N,climdata,constants, &
							LAI_N_data_frame_dims,LAI_N_dims,climdata_dims)
							
	use aconite_function_mod
	
 	implicit none
							
	integer, intent(in) :: LAI_N_data_frame_dims(2),LAI_N_dims(2),climdata_dims(2)
	double precision, intent(inout) ::LAI_N_dataframe(LAI_N_data_frame_dims(1),LAI_N_data_frame_dims(2))
	double precision, intent(inout) ::LAI_N(LAI_N_dims(1),LAI_N_dims(2))
	double precision, intent(inout) ::climdata(climdata_dims(1),climdata_dims(2))
		   
	double precision, intent(in) :: constants(13), par(20)
		
	double precision :: leafC_return(365),leafN_return(365),gpp_day(365), &
  					  gpp_leafc_marg(365),gpp_leafn_marg(365),RA(365), &
  					  RA_leafn_marg(365),RA_leafc_marg(365),growth_respiration(365), &
  					  leaf_allocation(365)
  
  	double precision :: leafC_cohort(2),leafN_cohort(2),cohort(365,2),add_c(365),add_n(365)
  	
  	double precision :: Lat,psid,rtot,leafgrow_doy,leafdrop_doy,base_resp_temp,lca, &
  						LL,ncohorts,nyears,resp_parm1,resp_parm2,Q10,grow_resp,psn_temp_thres, &
  						m,temp_max(365),temp_min(365),radiation(365),co2(365),doy(365), &
  						sla_top,maxleafC,curr_lca,leafC,leafN,leafCN,GDD,leaf_out_occured, &
  						z,tav,DL,curr_lai,add_lai,CN,fluxC,fluxN,leaf_horizon, &
  						t_leaf,leaflife_span,T_response,leaf_C_ramp_up,leaf_N_ramp_up, &
  						leaf_C_ramp_down,leaf_N_ramp_down,start_leafC, &
  						start_leafN,leaf_drop_occured,ramp_up,max_add_c,max_add_n, &
  						ramp_down,ramp_down_doy,use_area_resp,t10_vector(10),t10, resp_parm3,twq
  						
  	integer :: h,j,i
		   	
  	Lat = constants(1)
  	psid = constants(2)
  	rtot = constants(3)
  	leafgrow_doy = constants(4)
  	leafdrop_doy = constants(5)
  	max_add_c = constants(6)
  	base_resp_temp = constants(7)
  	lca = constants(8)
  	LL= constants(9)
  	ncohorts= constants(10)
  	nyears= constants(11)
  	use_area_resp = constants(12)
  	twq = constants(13)
  
  	resp_parm1 = par(11)
  	resp_parm2= par(12)
  	Q10 = par(13)
  	grow_resp = par(14)
  	psn_temp_thres = par(15)
  	m = par(16)
  	ramp_up = par(17)
  	ramp_down = par(18)
  	ramp_down_doy=par(19)
  	resp_parm3 = par(20)
  	  	
  	temp_max(:) = climdata(:,6)
  	temp_min(:) = climdata(:,7)
  	radiation(:) = climdata(:,4)
  	co2(:) = climdata(:,10)
  	doy(:) = climdata(:,3)
  	

  	max_add_c = 1.0
  	max_add_n = 0.1
  
  	do h = 1,LAI_N_data_frame_dims(1)
  	
  		!LAI_N(h,2) = LAI_N(h,2)*LAI_N(h,1)
  	
  		maxleafC=0.0
  		
  		sla_top = 1/lca
    	maxleafC = ((log(LAI_N(h,1)*m + sla_top) - log(sla_top))/m)
  		curr_lca = maxleafC/LAI_N(h,1)
  	
    	if(Lat > 30) then ! Seasonal Phenology
    		if(LL <= 1) then
    			start_leafC =  0.0
    			start_leafN = 0.0
    		else
    	    	start_leafC = max(maxleafC * (1-1/LL),0.0)
    			start_leafN = max(LAI_N(h,2) *(1-1/LL),0.0)
    		endif
    		leafC = start_leafC
    		leafN = start_leafN
    	else ! Tropical Penology
    	    start_leafC = maxleafC
    		start_leafN = LAI_N(h,2)
    		leafC = start_leafC
    		leafN = start_leafN
    	endif    	
    	
    	leafCN = maxleafC/LAI_N(h,2)
    	GDD = 0.0
    	leaf_out_occured = 0.0
    	leaf_drop_occured = 365
      	leaf_C_ramp_up = (maxleafC - leafC)/ramp_up
      	leaf_N_ramp_up = (LAI_N(h,2) - leafN)/ramp_up
    	leaf_C_ramp_down = (maxleafC - start_leafC)/ramp_down
      	leaf_N_ramp_down = (LAI_N(h,2) - start_leafN)/ramp_down
      	t10_vector(:) = (temp_max(1) + temp_min(1))/2. 
      	t10 = sum(t10_vector)/10
      	
    		do i = 1,365
      			tav = (temp_max(i) + temp_min(i))/2.
      			t10_vector(1:9) = t10_vector(2:10)
      			t10_vector(10) =  (temp_max(i) + temp_min(i))/2.
      			t10 = sum(t10_vector)/10
      			GDD = GDD_calc(temp_max(i),temp_min(i),GDD,i)
      			DL = daylength(Lat,i)
      			sla_top = 1/lca
      			curr_lai = (leafC)/curr_lca
      			add_c(i) = max_add_c*(leafC/(maxleafC))
      			add_lai = (leafC+add_c(i))/curr_lca
      			!add_lai = (sla_top*(exp(m*(leafC+add_c))-1))/m
      			add_n(i) = max_add_n*(leafC/(maxleafC))
      			
      			!print *, h, leafC, leafN
      			
      			if(leafC > 0) then
        				CN = leafC/leafN
        				!T_response = exp(0.1012*(MAX(tav,0.0)-base_resp_temp) - 0.0005*(MAX(tav,0.0)**2-base_resp_temp**2))
        				T_response = (Q10**((tav-base_resp_temp)/10)) !*(0.67+0.33*(1-DL/24))
        				gpp_day(i) = acm(leafN,curr_lai,temp_max(i),temp_min(i),DL,radiation(i),co2(i),&
        					par,psid,rtot,leafgrow_doy,leafdrop_doy,i)
        				if(use_area_resp == 0) then
        					RA(i) = mass_resp(leafC,leafN,resp_parm1,resp_parm2,resp_parm3,twq)*T_response
        				elseif(use_area_resp == 1) then
        					RA(i) = area_resp(curr_lai,leafN,resp_parm1,resp_parm2,resp_parm3,twq)*T_response			     				
        				elseif(use_area_resp == 2) then
        					RA(i) = atkin_resp(curr_lai,leafN,twq,resp_parm1,resp_parm2,resp_parm3)*T_response
        					if(i == 200) then
        						!print *, RA(i)/T_response,curr_lai, &
        						!	leafN,twq,resp_parm1,resp_parm2,resp_parm3,T_response
        					endif        				
        				elseif(use_area_resp == 3) then
        					RA(i) = ryan_resp(leafN,leafC,resp_parm1)*T_response
        				elseif(use_area_resp == 4) then
        					RA(i) = clm_resp(leafN,leafC,resp_parm1)*T_response

        				endif

        				
        				gpp_leafc_marg(i) = acm(leafN,add_lai,temp_max(i),temp_min(i),DL,&
        					radiation(i),co2(i),par,psid,rtot,leafgrow_doy,leafdrop_doy,i)
        				gpp_leafn_marg(i) = acm((leafN+add_n(i)),curr_lai,temp_max(i),temp_min(i),&
        					DL,radiation(i),co2(i),par,psid,rtot,leafgrow_doy,leafdrop_doy,i)
        				if(use_area_resp == 0) then
        					RA_leafc_marg(i)= mass_resp((leafC+add_c(i)),leafN,&
        						resp_parm1,resp_parm2,resp_parm3,twq)*T_response  
        					RA_leafn_marg(i)= mass_resp(leafC,(leafN+add_n(i)), &
        						resp_parm1,resp_parm2,resp_parm3,twq)*T_response	
        				elseif(use_area_resp == 1) then
        					RA_leafc_marg(i)= area_resp((add_lai),leafN,&
        						resp_parm1,resp_parm2,resp_parm3,twq)*T_response  
        					RA_leafn_marg(i)= area_resp(curr_lai,(leafN+add_n(i)), &
        						resp_parm1,resp_parm2,resp_parm3,twq)*T_response	
        				elseif(use_area_resp == 2) then
        					RA_leafc_marg(i)= atkin_resp((add_lai),leafN,twq,resp_parm1, &
        						resp_parm2,resp_parm3)*T_response  
        					RA_leafn_marg(i)= atkin_resp(curr_lai,(leafN+add_n(i)),twq,resp_parm1, &
        						resp_parm2,resp_parm3)*T_response	
        				elseif(use_area_resp == 3) then
        					RA_leafc_marg(i)= ryan_resp(leafN,leafC+add_c(i),resp_parm1)*T_response  
        					RA_leafn_marg(i)= ryan_resp((leafN+add_n(i)),leafC,resp_parm1)*T_response	
      					elseif(use_area_resp == 4) then
        					RA_leafc_marg(i)= clm_resp(leafN,leafC+add_c(i),resp_parm1)*T_response  
        					RA_leafn_marg(i)= clm_resp((leafN+add_n(i)),leafC,resp_parm1)*T_response	
      					endif
      			else	
        				gpp_day(i) = 0.0
        				gpp_leafc_marg(i) = 0.0
        				gpp_leafn_marg(i) = 0.0
        				RA(i) = 0.0
        				RA_leafc_marg(i)= 0.0
        				RA_leafn_marg(i)= 0.0
      			endif
      				
      				growth_respiration(i) = 0.0
      				leaf_allocation(i) = 0.0 
      				      				
      				if(Lat > 30) then  ! Seasonal Phenology
       					growth_respiration(i) = 0.0
      					leaf_allocation(i) = 0.0      					
      					if (i > leafgrow_doy .AND. leaf_out_occured == 0)  then
      						if(maxleafC-0.1 <= leafC) then
        						leaf_out_occured = i
        					else
        					leafC = leafC + leaf_C_ramp_up
        					leaf_allocation(i) = leaf_C_ramp_up
        					growth_respiration(i) = leaf_allocation(i)*grow_resp 
        					leafN = leafN + leaf_N_ramp_up
							endif
      					endif

      					if(i > (ramp_down_doy) .AND. leafC > start_leafC+0.1) then
      						leafC = leafC - leaf_C_ramp_down
      						leafN = leafN - leaf_N_ramp_down
      						if(leafC <= start_leafC+0.1) then 
      							leafC = start_leafC
      							leaf_drop_occured = i
      						endif
      						if(leafN < 0.001) then 
      							leafN  = 0.0
      						endif   					
      					endif
      				else  ! Tropical Phenology
      					leaf_allocation(i) = (maxleafC * (1/LL))/365
      					growth_respiration(i) = leaf_allocation(i)*grow_resp 
      					leaf_out_occured = 1
      					leaf_drop_occured = 365
      				endif
      							
      				if(leafN <0) then
      					print *,h, leafN,LAI_N(h,1),LAI_N(h,2)
      				endif
      				
    		end do
  
  		if(LL > 1 .OR. Lat < 30) then
    		leaf_horizon = LL
  		else
    		leaf_horizon = 1.0
  		endif
  		
  		t_leaf = 1/leaf_horizon
  		leaflife_span = ((leafdrop_doy-leaf_out_occured + 1/LL))
  		LAI_N_dataframe(h,1) = sum(gpp_day)
  		LAI_N_dataframe(h,2) = sum(RA)
  		LAI_N_dataframe(h,3)= (LAI_N_dataframe(h,1)- LAI_N_dataframe(h,2))  &
  			- sum(growth_respiration) - sum(leaf_allocation)
  	  	LAI_N_dataframe(h,4)= ((sum(gpp_leafc_marg)-sum(gpp_day)-(sum(RA_leafc_marg) &
  			- sum(RA)))) - (max_add_c*grow_resp + max_add_c)/(leaf_horizon)	  		
  		LAI_N_dataframe(h,5) = ((sum(gpp_leafn_marg)-sum(gpp_day) - (sum(RA_leafn_marg) &
  			- sum(RA))))
  		LAI_N_dataframe(h,6) = maxleafC
  		LAI_N_dataframe(h,7)=  maxleafC/LAI_N(h,2)
  		LAI_N_dataframe(h,8) = t_leaf
  		LAI_N_dataframe(h,11)=  leaflife_span
  		LAI_N_dataframe(h,12) = sum(leaf_allocation)
  		LAI_N_dataframe(h,13) = leaf_horizon
  		LAI_N_dataframe(h,14) = sum(growth_respiration)
  		
  		!print *, LAI_N(h,1),LAI_N(h,2),LAI_N_dataframe(h,2)
  		
 	end do
 		
end subroutine
