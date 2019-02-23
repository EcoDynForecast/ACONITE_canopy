aconite_canopy_model <- function(fortran_code,resp_type,
                                 tropic_LAI_data,
                                 tropic_TFN_data,
                                 evergreen_LAI_data,
                                 evergreen_TFN_data,
                                 decid_LAI_data,
                                 decid_TFN_data,
                                 output_file,
                                 plot_name,
                                 resp_parm1_sens,
                                 resp_parm2_sens,  
                                 resp_parm3_sens){
  
  grow_resp <- 0.28
  psn_temp_thres <- 0
  biome_acm <- TRUE
  biome_Q10 <- FALSE
  save_contours <- TRUE
  use_pdf <- FALSE
  create_plots <- FALSE
  biome <- c(1,2,3)
  LL_sens <- 1
  LCA_sens <- 1
  ACM_1_sens <- 1
  Q10_sens <- 1
  temp_increase <- 0
  grow_season_increase <- 0
  co2_increase <- 0
  
  decid_fit <- lm(decid_TFN_data~decid_LAI_data)
  evergreen_fit <- lm(evergreen_TFN_data~evergreen_LAI_data)
  tropic_fit <- lm(tropic_TFN_data~tropic_LAI_data)
  
  decid_CI <- confint(decid_fit)
  evergreen_CI <- confint(evergreen_fit)
  tropic_CI <- confint(tropic_fit)
  
  if(resp_type==0){
    base_resp_temp = 20 #Reich et al. 2008 normalized to 20C
  }else if(resp_type == 1){
    base_resp_temp = 20 #Reich et al. 2008 normalized to 20C
  }else if(resp_type == 2){
    base_resp_temp = 25  #Aktin et al 2016 normalized to 25C
  }else if(resp_type == 3){
    base_resp_temp = 20 #Ryan 1991 normalized to 20C
  }
  
  #Parameters that arent important
  ncohorts = 2    # number of tracked leaf cohorts (maximum foliage age)  #not currently used
  max_add_c = 1.0
  
  # run duration
  nyears = 1    # length of simulation in years
  
  #Create arrays to storage output from the simulations
  LAI_TFN_mean_export = array(NA,c(length(biome),length(resp_parm1_sens),length(resp_parm2_sens),length(LL_sens),length(LCA_sens),length(ACM_1_sens),length(Q10_sens)))
  LAI_TFN_mean_GPP = array(NA,c(length(biome),length(resp_parm1_sens),length(resp_parm2_sens),length(LL_sens),length(LCA_sens),length(ACM_1_sens),length(Q10_sens)))
  LAI_TFN_mean_ra_mass =array(NA,c(length(biome),length(resp_parm1_sens),length(resp_parm2_sens),length(LL_sens),length(LCA_sens),length(ACM_1_sens),length(Q10_sens)))
  LAI_TFN_mean_la_rg = array(NA,c(length(biome),length(resp_parm1_sens),length(resp_parm2_sens),length(LL_sens),length(LCA_sens),length(ACM_1_sens),length(Q10_sens)))
  LAI_TFN_mean_la = array(NA,c(length(biome),length(resp_parm1_sens),length(resp_parm2_sens),length(LL_sens),length(LCA_sens),length(ACM_1_sens),length(Q10_sens)))
  LAI_TFN_mean_rg = array(NA,c(length(biome),length(resp_parm1_sens),length(resp_parm2_sens),length(LL_sens),length(LCA_sens),length(ACM_1_sens),length(Q10_sens)))
  max_export_sens = array(NA,c(length(biome),length(resp_parm1_sens),length(resp_parm2_sens),length(LL_sens),length(LCA_sens),length(ACM_1_sens),length(Q10_sens)))
  max_gpp_sens = array(NA,c(length(biome),length(resp_parm1_sens),length(resp_parm2_sens),length(LL_sens),length(LCA_sens),length(ACM_1_sens),length(Q10_sens)))
  max_ra_mass_sens = array(NA,c(length(biome),length(resp_parm1_sens),length(resp_parm2_sens),length(LL_sens),length(LCA_sens),length(ACM_1_sens),length(Q10_sens)))
  max_la_rg_sens = array(NA,c(length(biome),length(resp_parm1_sens),length(resp_parm2_sens),length(LL_sens),length(LCA_sens),length(ACM_1_sens),length(Q10_sens)))
  max_la_sens = array(NA,c(length(biome),length(resp_parm1_sens),length(resp_parm2_sens),length(LL_sens),length(LCA_sens),length(ACM_1_sens),length(Q10_sens)))
  max_rg_sens = array(NA,c(length(biome),length(resp_parm1_sens),length(resp_parm2_sens),length(LL_sens),length(LCA_sens),length(ACM_1_sens),length(Q10_sens)))
  max_LAI_sens = array(NA,c(length(biome),length(resp_parm1_sens),length(resp_parm2_sens),length(LL_sens),length(LCA_sens),length(ACM_1_sens),length(Q10_sens)))
  max_N_sens = array(NA,c(length(biome),length(resp_parm1_sens),length(resp_parm2_sens),length(LL_sens),length(LCA_sens),length(ACM_1_sens),length(Q10_sens)))
  data_proportion = array(NA,c(length(biome),length(resp_parm1_sens),length(resp_parm2_sens),length(LL_sens),length(LCA_sens),length(ACM_1_sens),length(Q10_sens)))
  data_proportion_1LAI = array(NA,c(length(biome),length(resp_parm1_sens),length(resp_parm2_sens),length(LL_sens),length(LCA_sens),length(ACM_1_sens),length(Q10_sens)))
  model_slope = array(NA,c(length(biome),length(resp_parm1_sens),length(resp_parm2_sens),length(LL_sens),length(LCA_sens),length(ACM_1_sens),length(Q10_sens)))
  slope_in = array(NA,c(length(biome),length(resp_parm1_sens),length(resp_parm2_sens),length(LL_sens),length(LCA_sens),length(ACM_1_sens),length(Q10_sens)))
  
  #loop through biome/PFT combinations
  for(kk in 1:length(biome)){
    #print(kk)
    if(biome[kk] == 1){
      #Toolik decidious shrubs
      climdata = read.csv("input/climate_TL_2007_2008_noleap.csv")
      twq = 9.324752  #Temperature of warmest quarter
      mean_LAI = mean(decid_LAI_data) #Get LAI from dataset
      mean_TFN = mean(decid_TFN_data) #Get TCN from dataset
      CI = decid_CI
      obs_LAI = decid_LAI_data
      obs_TFN =  decid_TFN_data
      Lat = 68 #Latitude
      leafgrow_doy = 160 - grow_season_increase/2 #Day of year that spring phenology starts
      leafdrop_doy = 260 + grow_season_increase/2 #Day of year that fall phenology starts
      lca= 41 #Leaf carbon per area
      LL = 1 #Leaf lifespan (years)
      ramp_up = 20  #Number of days used to leaf out
      ramp_down = 20   #Number of days used to drop leaves
      ramp_down_doy = leafdrop_doy-ramp_down
      LAI = c(0.05,seq(0.1,1,0.1),seq(1.25,8,0.1))  #Values of LAI and N to simulation
      N = c(0.05,seq(0.1,1,0.1),seq(1.25,20,0.1))
      m = 0.0000000000001  #Small number means that M is effectively zero.  This allows for future studies that have mean SLA change with total LAI
      Q10 = 2.0
      if(biome_acm == FALSE){
        load("input/acm_recal_with_spa_200pixels_continuous_parameters.RData")
      }else{
        load("input/acm_recal_with_spa_continuous_arctic_optimum_for_arctic_parameters.RData")        
      }
      a1 = median(parameters[1,,])
      a2 = median(parameters[2,,])
      a3 = median(parameters[3,,])
      a4 = median(parameters[4,,])
      a5 = median(parameters[5,,])
      a6 = median(parameters[6,,])
      a7 = median(parameters[7,,])
      a8 = median(parameters[8,,])
      a9 = median(parameters[9,,])
      a10 = median(parameters[10,,])
      psid = median(parameters[11,,])
      rtot = median(parameters[12,,])
      acm_par=c(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10)
    }else if(biome[kk] == 2){
      #Toolik evergreen shrubs
      climdata = read.csv("input/climate_TL_2007_2008_noleap.csv")
      twq = 9.324752
      mean_LAI = mean(evergreen_LAI_data)
      mean_TFN = mean(evergreen_TFN_data)
      CI = evergreen_CI
      obs_LAI = evergreen_LAI_data
      obs_TFN =  evergreen_TFN_data
      Lat = 68
      leafgrow_doy = 160 - grow_season_increase/2
      leafdrop_doy = 260 + grow_season_increase/2
      lca= 85
      LL = 2.191781
      ramp_up = 60
      ramp_down = 120
      ramp_down_doy = leafgrow_doy+ramp_up
      LAI = seq(1,4,0.01)
      N = seq(2,8,0.01)
      LAI = c(0.05,seq(0.1,1,0.05),seq(1.25,8,0.1))
      N = c(0.05,seq(0.1,1,0.05),seq(1.25,20,0.1))
      m = 0.0000000000001
      Q10 = 2.0
      if(biome_acm == FALSE){
        load("input/acm_recal_with_spa_200pixels_continuous_parameters.RData")
      }else{
        load("input/acm_recal_with_spa_continuous_arctic_optimum_for_arctic_parameters.RData")        
      }
      a1 = median(parameters[1,,])
      a2 = median(parameters[2,,])
      a3 = median(parameters[3,,])
      a4 = median(parameters[4,,])
      a5 = median(parameters[5,,])
      a6 = median(parameters[6,,])
      a7 = median(parameters[7,,])
      a8 = median(parameters[8,,])
      a9 = median(parameters[9,,])
      a10 = median(parameters[10,,])
      psid = median(parameters[11,,])
      rtot = median(parameters[12,,])
      acm_par=c(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10)
    }else if(biome[kk] == 3){
      #Tropical broadleaf evergreen La Selva
      climdata = read.csv("input/climate_LaSelva_2007.csv")
      twq = 27.15946
      mean_LAI = mean(tropic_LAI_data)
      mean_TFN = mean(tropic_TFN_data)
      CI = tropic_CI
      obs_LAI = tropic_LAI_data
      obs_TFN =  tropic_TFN_data
      Lat = 10.42
      leafgrow_doy = 1
      leafdrop_doy = 365
      lca= 44    
      LL = 1.328767
      ramp_up = 20
      ramp_down = 20
      ramp_down_doy = leafdrop_doy-ramp_down
      LAI = c(0.05,seq(0.1,1,0.1),seq(1.25,13,0.1))
      N = c(0.05,seq(0.1,1,0.1),seq(1.25,26,0.1))
      m = 0.00000000001
      Q10 = 2.0
      load("input/acm_recal_with_spa_200pixels_continuous_parameters.RData")
      a1 = median(parameters[1,,])
      a2 = median(parameters[2,,])
      a3 = median(parameters[3,,])
      a4 = median(parameters[4,,])
      a5 = median(parameters[5,,])
      a6 = median(parameters[6,,])
      a7 = median(parameters[7,,])
      a8 = median(parameters[8,,])
      a9 = median(parameters[9,,])
      a10 = median(parameters[10,,])
      psid = median(parameters[11,,])
      rtot = median(parameters[12,,])
      acm_par=c(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10)
    }
    
    climdata$TEMP_MAX = climdata$TEMP_MAX + temp_increase
    climdata$TEMP_MIN = climdata$TEMP_MIN + temp_increase
    climdata$CO2 = climdata$CO2 + co2_increase
    
    if(length(ACM_1_sens) == 1) {
      ACM_1_sens = a1
    }else{
      ACM_1_sens = ACM_1_sens*a1
    }
    if(length(Q10_sens) == 1) {
      Q10_sens = Q10
    }  
    if(length(LCA_sens)==1){
      LCA_sens = lca
    }
    if(length(LL_sens)==1){
      LL_sens = LL
    }
    
    for(i1 in 1:length(resp_parm1_sens)){
      for(i2 in 1:length(resp_parm2_sens)){ 
        for(i3 in 1:length(LL_sens)){ 
          for(i4 in 1:length(LCA_sens)){ 
            for(i5 in 1:length(ACM_1_sens)){ 
              for(i6 in 1:length(Q10_sens)){ 
                
                resp_parm1 = resp_parm1_sens[i1]
                resp_parm2 = resp_parm2_sens[i2]    
                
                LL = LL_sens[i3]
                lca = LCA_sens[i4]
                acm_par[1] = ACM_1_sens[i5]
                Q10 = Q10_sens[i6]
                
                constants = c(Lat,psid,rtot,leafgrow_doy,leafdrop_doy,max_add_c,base_resp_temp,lca,LL,ncohorts,nyears,resp_type,twq)
                par=c(acm_par,resp_parm1,resp_parm2,Q10,grow_resp,psn_temp_thres,m,ramp_up,ramp_down,ramp_down_doy,resp_parm3_sens)
                
                LAI_N = matrix(-99,length(N)*length(LAI),2)
                index = 0
                for(j in 1:length(LAI)){
                  for(k in 1:length(N)){
                    index = index + 1
                    LAI_N[index,]=c(LAI[j],N[k])
                  }
                }
                
                gpp= array(-99, dim=c(length(LAI_N[,1]),1))
                ra_mass= array(-99, dim=c(length(LAI_N[,1]),1))
                NetPsn_mass= array(-99, dim=c(length(LAI_N[,1]),1))
                leafC_return_mass= array(-99, dim=c(length(LAI_N[,1]),1))
                leafN_return_mass= array(-99, dim=c(length(LAI_N[,1]),1))
                leafC_used_in_return = array(-99, dim=c(length(LAI_N[,1]),1))
                leafCN_used_in_return = array(-99, dim=c(length(LAI_N[,1]),1))
                optimal = array(0, dim=c(length(LAI_N[,1]),1))
                optimal_in_data = array(0, dim=c(length(LAI_N[,1]),1))
                t_leaf_used = array(0, dim=c(length(LAI_N[,1]),1))
                leaflife_span= array(0, dim=c(length(LAI_N[,1]),1))
                leaf_allocation = array(-99, dim=c(length(LAI_N[,1]),1))
                leaf_horizon = array(-99, dim=c(length(LAI_N[,1]),1))
                leaf_growth_respiration = array(-99, dim=c(length(LAI_N[,1]),1))
                
                LAI_N_dataframe <- data.frame(gpp,ra_mass,NetPsn_mass,leafC_return_mass,leafN_return_mass,
                                              leafC_used_in_return,leafCN_used_in_return,t_leaf_used,optimal,optimal_in_data,leaflife_span,
                                              leaf_allocation,leaf_horizon,leaf_growth_respiration)
                
                dyn.load(fortran_code)
                
                fortran_output=.Fortran("LAI_N_Rinterface",
                                        LAI_N_dataframe = as.matrix(LAI_N_dataframe),
                                        par = as.matrix(par),
                                        LAI_N = as.matrix(LAI_N),
                                        climdata = as.matrix(climdata),
                                        constants = as.matrix(constants),
                                        LAI_N_data_frame_dims = as.integer(c(dim(LAI_N_dataframe)[1],dim(LAI_N_dataframe)[2])),
                                        LAI_N_dims = as.integer(c(dim(LAI_N)[1],dim(LAI_N)[2])),
                                        climdata_dims = as.integer(c(dim(climdata)[1],dim(climdata)[2])))                            
                
                par_out = array(fortran_output$LAI_N_dataframe, dim=c(dim(LAI_N_dataframe)[1],dim(LAI_N_dataframe)[2]))
                
                gpp= array(-99, dim=c(length(LAI),length(N)))
                ra_mass= array(-99, dim=c(length(LAI),length(N)))
                NetPsn_mass= array(-99, dim=c(length(LAI),length(N)))
                leafC_return_mass= array(-99, dim=c(length(LAI),length(N)))
                leafN_return_mass= array(-99, dim=c(length(LAI),length(N)))
                leafC_used_in_return = array(-99, dim=c(length(LAI),length(N)))
                leafCN_used_in_return = array(-99, dim=c(length(LAI),length(N)))
                t_leaf_used = array(-99, dim=c(length(LAI),length(N)))
                leaflife_span = array(-99, dim=c(length(LAI),length(N)))
                leaf_allocation = array(-99, dim=c(length(LAI),length(N)))
                leaf_horizon = array(-99, dim=c(length(LAI),length(N)))
                leaf_growth_respiration = array(-99, dim=c(length(LAI),length(N)))
                
                index=0
                for(j in 1:length(LAI)){
                  for(k in 1:length(N)){
                    index = index + 1
                    gpp[j,k]= par_out[index,1]
                    ra_mass[j,k]= par_out[index,2]
                    NetPsn_mass[j,k]= par_out[index,3]
                    leafC_return_mass[j,k]= par_out[index,4]
                    leafN_return_mass[j,k]= par_out[index,5]
                    leafC_used_in_return[j,k] = par_out[index,6]
                    leafCN_used_in_return[j,k] = par_out[index,7]
                    t_leaf_used[j,k] = par_out[index,8]
                    leaflife_span[j,k]=par_out[index,11]
                    leaf_allocation[j,k]=par_out[index,12]
                    leaf_horizon[j,k]=par_out[index,13]
                    leaf_growth_respiration[j,k]=par_out[index,14]
                  }
                }
                
                max_export=-9999
                max_gpp=-9999
                max_ra_mass=-9999
                max_la_rg=-9999
                max_la=-9999
                max_rg=-9999
                max_LAI=-9999
                max_N=-9999
                
                for(j in 1:length(LAI)){
                  for(k in 1:length(N)){
                    if(NetPsn_mass[j,k] > max_export){
                      max_export = NetPsn_mass[j,k]
                      max_gpp = gpp[j,k]
                      max_ra_mass = ra_mass[j,k]
                      max_la = leaf_allocation[j,k]
                      max_rg = leaf_growth_respiration[j,k]
                      max_LAI = LAI[j]
                      max_N = N[k]
                    }
                  }
                }
                
                j = which.min(abs(LAI-mean_LAI))
                k = which.min(abs(N - mean_TFN))
                LAI_TFN_mean_export[kk,i1,i2,i3,i4,i5,i6] = NetPsn_mass[j,k]
                LAI_TFN_mean_GPP[kk,i1,i2,i3,i4,i5,i6] = gpp[j,k]
                LAI_TFN_mean_ra_mass[kk,i1,i2,i3,i4,i5,i6] =ra_mass[j,k]
                LAI_TFN_mean_la_rg[kk,i1,i2,i3,i4,i5,i6] = leaf_allocation[j,k] + leaf_growth_respiration[j,k]
                LAI_TFN_mean_la[kk,i1,i2,i3,i4,i5,i6] = leaf_allocation[j,k]
                LAI_TFN_mean_rg[kk,i1,i2,i3,i4,i5,i6] = leaf_growth_respiration[j,k]
                
                max_export_sens[kk,i1,i2,i3,i4,i5,i6] = max_export
                max_gpp_sens[kk,i1,i2,i3,i4,i5,i6] = max_gpp
                max_ra_mass_sens[kk,i1,i2,i3,i4,i5,i6] = max_ra_mass
                max_la_rg_sens[kk,i1,i2,i3,i4,i5,i6] = max_la_rg
                max_la_sens[kk,i1,i2,i3,i4,i5,i6] = max_la
                max_rg_sens[kk,i1,i2,i3,i4,i5,i6] = max_rg
                max_LAI_sens[kk,i1,i2,i3,i4,i5,i6] = max_LAI
                max_N_sens[kk,i1,i2,i3,i4,i5,i6] = max_N
                leafC_return_mass= array(-99, dim=c(length(LAI),length(N)))
                leafN_return_mass= array(-99, dim=c(length(LAI),length(N)))
                
                index=0
                for(j in 1:length(LAI)){
                  for(k in 1:length(N)){
                    index = index + 1
                    leafC_return_mass[j,k]= par_out[index,4]
                    leafN_return_mass[j,k]= par_out[index,5]
                    
                  }
                }
                zero_iso_leafC = -99
                zero_iso_leafN = -99
                leafC.linelist=contourLines(LAI,N,leafC_return_mass)
                leafN.linelist=contourLines(LAI,N,leafN_return_mass)
                for(i in 1:length(leafC.linelist)){
                  if(leafC.linelist[[i]]$level==0){
                    zero_iso_leafC = i
                  }
                }
                for(i in 1:length(leafN.linelist)){
                  if(leafN.linelist[[i]]$level==0){
                    zero_iso_leafN = i
                  }
                }
                
                if(zero_iso_leafC != -99 && zero_iso_leafN !=-99){
                  x_isolineC = leafC.linelist[[zero_iso_leafC]]$x[which(leafC.linelist[[zero_iso_leafC]]$x>0.02)]
                  y_isolineC = leafC.linelist[[zero_iso_leafC]]$y[which(leafC.linelist[[zero_iso_leafC]]$x>0.02)]
                  x_isolineN = leafN.linelist[[zero_iso_leafN]]$x
                  y_isolineN = leafN.linelist[[zero_iso_leafN]]$y
                }
                
                data_in_range = 0
                tot_data_points = 0
                tot_data_points_1LAI = 0
                data_in_range_1LAI = 0
                
                if(zero_iso_leafC != -99 && zero_iso_leafN !=-99){
                  for(d in 1:length(obs_TFN)){
                    tot_data_points = tot_data_points + 1
                    if(obs_LAI[d] > 0.5) {tot_data_points_1LAI = tot_data_points_1LAI + 1}
                    if(obs_TFN[d]>=y_isolineC[which.min(abs(x_isolineC-obs_LAI[d]))] && 
                       obs_TFN[d]<=y_isolineN[which.min(abs(x_isolineN-obs_LAI[d]))] &&
                       obs_LAI[d]<=x_isolineC[which.min(abs(y_isolineC-obs_TFN[d]))] && 
                       obs_LAI[d]>=x_isolineN[which.min(abs(y_isolineN-obs_TFN[d]))]){
                      data_in_range = data_in_range + 1
                      if(obs_LAI[d] > 0.5) {data_in_range_1LAI = data_in_range_1LAI + 1}
                    }
                  }
                }
                
                data_proportion[kk,i1,i2,i3,i4,i5,i6]= data_in_range/tot_data_points
                data_proportion_1LAI[kk,i1,i2,i3,i4,i5,i6]= data_in_range_1LAI/tot_data_points_1LAI
                
                if(max_LAI_sens[kk,i1,i2,i3,i4,i5,i6] < max(LAI) && max_N_sens[kk,i1,i2,i3,i4,i5,i6] < max(N)){
                  model_slope[kk,i1,i2,i3,i4,i5,i6] = max_N_sens[kk,i1,i2,i3,i4,i5,i6]/max_LAI_sens[kk,i1,i2,i3,i4,i5,i6]
                }else{
                  model_slope[kk,i1,i2,i3,i4,i5,i6] = -9999
                } 
                
                ###############################
                
                ylabel = expression(paste('Total canopy nitrogen (g N ',m^-2,')',sep=''))
                xlabel = expression(paste('Leaf area index (',m^2,' ',m^-2,')',sep=''))
                
                
                leafC.linelist=contourLines(LAI,N,leafC_return_mass)
                leafN.linelist=contourLines(LAI,N,leafN_return_mass)
                zero_iso_leafC = -99
                zero_iso_leafN = -99
                for(i in 1:length(leafC.linelist)){
                  if(leafC.linelist[[i]]$level==0){
                    zero_iso_leafC = i
                  }
                }
                for(i in 1:length(leafN.linelist)){
                  if(leafN.linelist[[i]]$level==0){
                    zero_iso_leafN = i
                  }
                }
                
                if(length(obs_LAI) > 1){
                  slope_in[kk,i1,i2,i3,i4,i5,i6] = 0
                  predicted_N_low = CI[1,1] + max_LAI_sens[kk,i1,i2,i3,i4,i5,i6]*CI[2,1]
                  predicted_N_high = CI[1,2] + max_LAI_sens[kk,i1,i2,i3,i4,i5,i6]*CI[2,2]
                  if(max_N_sens[kk,i1,i2,i3,i4,i5,i6] >= predicted_N_low && max_N_sens[kk,i1,i2,i3,i4,i5,i6] <= predicted_N_high){
                    slope_in[kk,i1,i2,i3,i4,i5,i6] = 1
                  }
                }
                if(biome[kk]==1)canopy_name = 'arctic_decid'
                if(biome[kk]==2)canopy_name = 'arctic_evergreen'
                if(biome[kk]==3)canopy_name = 'tropical'
                
                if(save_contours){
                  save(LAI,N,gpp,ra_mass,leaf_growth_respiration,leaf_allocation,leafC.linelist,zero_iso_leafC,leafN.linelist,zero_iso_leafN,mean_LAI,mean_TFN,CI,max_LAI_sens,max_N_sens,obs_LAI,obs_TFN,file= paste(figure_directory,'/',plot_name,'_',canopy_name,'_',as.character(resp_parm1),'_',as.character(resp_parm2),'_contours.Rdata',sep=''))
                }
                
                if(create_plots == TRUE){
                  
                  
                  
                  file_figures = paste(figure_directory,'/',plot_name,'_',canopy_name,'_',as.character(resp_parm1),'_',as.character(resp_parm2),'.pdf',sep='')
                  if(use_pdf) pdf(file_figures,height = 7.007874*0.7,width = 7.007874*2.5)
                  par(mfrow=c(1,4), tcl=-0.4) 
                  interval = 100
                  par(mar = c(5,5,2,1),oma=c(0,0,0,0))
                  contour_cex = 1.1
                  if(biome[kk] ==3){
                    xlim=c(0,14)
                    ylim=c(0,28)
                    interval = 300
                    upper_gpp = 4000
                    upper_ra = 3000
                    upper_alloc = 3000
                  }else{
                    xlim=c(0,8)
                    ylim=c(0,16)  
                    upper_gpp = 1200
                    upper_ra = 1200
                    upper_alloc = 5000
                  }
                  
                  contour(LAI,N,((gpp-ra_mass)-leaf_growth_respiration-leaf_allocation),xlab='',ylab='',xlim=xlim,ylim=ylim,levels=seq(0,5000,interval),lwd=2.0,labcex=contour_cex,cex.lab=2.0,cex.axis=2.0)
                  title(main = 'Carbon export',cex.main = 1.5)
                  #contour(LAI,N,((gpp-ra_mass)-leaf_growth_respiration-leaf_allocation),xlab='',ylab='',levels=seq(-600,-100,interval),add=TRUE,col='gray',lwd=2.0,labcex=contour_cex,cex.lab=2.0,cex.axis=2.0)
                  points(leafC.linelist[[zero_iso_leafC]]$x,leafC.linelist[[zero_iso_leafC]]$y,type='l',col='red',lwd=4,lty='dashed')
                  points(leafN.linelist[[zero_iso_leafN]]$x,leafN.linelist[[zero_iso_leafN ]]$y,type='l',col='blue',lwd=4,lty='longdash')
                  
                  text(0.3,14.5,'(a)',cex = 1.5)
                  contour(LAI,N,gpp,xlab='',ylab='',levels=seq(0,upper_gpp,interval),lwd=1.5,labcex=contour_cex,xlim=xlim,ylim=ylim,cex.lab=2.0,cex.axis=2.0)
                  rect(xleft=-0.1, ybottom = 14.8 , xright = 0.7, ytop = 16.2, col = 'white', border = 'white')
                  text(0.3,14.5,'(b)',cex = 1.5)
                  title(main = 'GPP',cex.main = 1.5)
                  contour(LAI,N,ra_mass,xlab='',ylab='',levels=seq(-600,upper_ra,interval),lwd=2.0,labcex=contour_cex,xlim=xlim,ylim=ylim,cex.lab=2.0,cex.axis=2.0)
                  text(0.3,14.5,'(c)',cex = 1.5)
                  title(main = 'Rm',cex.main = 1.5)
                  contour(LAI,N,leaf_allocation + leaf_growth_respiration,xlab='',ylab='',levels=seq(-600,upper_alloc,interval),lwd=2.0,labcex=contour_cex,xlim=xlim,ylim=ylim,cex.lab=2.0,cex.axis=2.0)
                  
                  text(0.3,14.5,'(d)',cex = 1.5)
                  title(main = 'AL + Rg',cex.main = 1.5)
                  
                  
                  interval = 100
                  contour_cex = 1.1
                  contour(LAI,N,((gpp-ra_mass)-leaf_growth_respiration-leaf_allocation),xlab=xlabel,ylab=ylabel,levels=seq(0,5000,interval),lwd=2.0,labcex=contour_cex,cex.lab=2.0,cex.axis=2.0)
                  contour(LAI,N,((gpp-ra_mass)-leaf_growth_respiration-leaf_allocation),xlab=xlabel,ylab=ylabel,levels=seq(-600,-100,interval),add=TRUE,col='gray',lwd=2.0,labcex=contour_cex,cex.lab=2.0,cex.axis=2.0)
                  title(main = 'Net Canopy Carbon Export')
                  points(leafC.linelist[[zero_iso_leafC]]$x,leafC.linelist[[zero_iso_leafC]]$y,type='l',col='red',lwd=2,lty='dashed')
                  points(leafN.linelist[[zero_iso_leafN]]$x,leafN.linelist[[zero_iso_leafN ]]$y,type='l',col='blue',lwd=2,lty='longdash')
                  
                  contour(LAI,N,((gpp-ra_mass)-leaf_growth_respiration-leaf_allocation),xlab=xlabel,ylab=ylabel,levels=seq(0,5000,interval),lwd=2.0,labcex=contour_cex,cex.lab=2.0,cex.axis=2.0)
                  contour(LAI,N,((gpp-ra_mass)-leaf_growth_respiration-leaf_allocation),xlab=xlabel,ylab=ylabel,levels=seq(-600,-100,interval),add=TRUE,col='gray',lwd=2.0,labcex=contour_cex,cex.lab=2.0,cex.axis=2.0)
                  title(main = 'Net Canopy Carbon Export')
                  points(leafC.linelist[[zero_iso_leafC]]$x,leafC.linelist[[zero_iso_leafC]]$y,type='l',col='red',lwd=2,lty='dashed')
                  points(leafN.linelist[[zero_iso_leafN]]$x,leafN.linelist[[zero_iso_leafN ]]$y,type='l',col='blue',lwd=2,lty='longdash')
                  if(length(obs_LAI) > 1)points(mean_LAI,mean_TFN,pch=19)
                  if(length(obs_LAI) > 1){
                    abline(CI[1,1],CI[2,1],lty='dotted',col='brown',lwd=2)
                    abline(CI[1,2],CI[2,2],lty='dotted',col='brown',lwd=2)   
                  }
                  plot.new()
                  legend1=expression(paste('contour (postive values: g C ',m^-2,' ',yr^-1,')',sep=''))
                  legend2=expression(paste('contour (negative values: g C ',m^-2,' ',yr^-1,')',sep=''))
                  legend('top',c(legend1,legend2,'C allocation marginal return = 0','N allocation marginal return = 0','Observations: 95% C.I.'),col=c('black','gray','red','blue','brown'),lty=c('solid','solid','dashed','longdash','dotted'),bty='n',lwd=c(2,2,2,2),cex=0.90)
                  legend('top',c(' ',' ',' ','','','Site Mean                                       '),col=c('transparent','transparent','transparent','transparent','transparent','black'),pch = c(1,1,1,1,1,19),bty='n',box.col='transparent',bg = 'transparent',text.col='black',cex=0.90)
                  if(use_pdf) dev.off()
                }
              }
            }
          }
        }
      }
    }
  }
  
  save(LAI_TFN_mean_export,LAI_TFN_mean_GPP,LAI_TFN_mean_ra_mass,LAI_TFN_mean_la_rg,
       LAI_TFN_mean_rg,max_export_sens,max_gpp_sens,max_ra_mass_sens,max_la_rg_sens,
       max_la_sens,max_rg_sens,NetPsn_mass,max_LAI_sens,max_N_sens,data_proportion,
       data_proportion_1LAI,model_slope,slope_in,resp_parm1_sens,resp_parm2_sens,
       LL_sens,LCA_sens,ACM_1_sens,Q10_sens,biome,file = output_file)
}
