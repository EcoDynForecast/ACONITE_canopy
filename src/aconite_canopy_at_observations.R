aconite_canopy_at_observations <- function(){
  
  grow_resp <- 0.28
  psn_temp_thres <- 0
  biome_acm <- TRUE
  biome_Q10 <- FALSE
  save_contours <- TRUE
  biome <- c(1,2,3)
  temp_increase <- 0
  grow_season_increase <- 0
  co2_increase <- 0
  
  ncohorts <- 2    # number of tracked leaf cohorts (maximum foliage age)  #not currently used
  max_add_c <- 1.0
  nyears <- 1    # length of simulation in years
  
  #loop through biome/PFT combinations
  for(kk in 1:length(biome)){
    #print(kk)
    if(biome[kk]  ==  1){
      #Toolik decidious shrubs
      climdata <- read.csv("input/climate_TL_2007_2008_noleap.csv")
      twq <- 9.324752  #Temperature of warmest quarter
      mean_LAI <- mean(decid_LAI_data) #Get LAI from dataset
      mean_TFN <- mean(decid_TFN_data) #Get TCN from dataset
      obs_LAI <- decid_LAI_data
      obs_TFN <-  decid_TFN_data
      Lat <- 68 #Latitude
      leafgrow_doy <- 160 - grow_season_increase/2 #Day of year that spring phenology starts
      leafdrop_doy <- 260 + grow_season_increase/2 #Day of year that fall phenology starts
      lca <- 41 #Leaf carbon per area
      LL <- 1 #Leaf lifespan (years)
      ramp_up <- 20  #Number of days used to leaf out
      ramp_down <- 20   #Number of days used to drop leaves
      ramp_down_doy <- leafdrop_doy-ramp_down
      LAI <- c(0.05,seq(0.1,1,0.1),seq(1.25,8,0.1))  #Values of LAI and N to simulation
      N <- c(0.05,seq(0.1,1,0.1),seq(1.25,16,0.1))
      m <- 0.0000000000001  #Small number means that M is effectively zero.  This allows for future studies that have mean SLA change with total LAI
      Q10 <- 2.0
      if(biome_acm  ==  FALSE){
        load("input/acm_recal_with_spa_200pixels_continuous_parameters.RData")
      }else{
        load("input/acm_recal_with_spa_continuous_arctic_optimum_for_arctic_parameters.RData")        
      }
    }else if(biome[kk]  ==  2){
      #Toolik evergreen shrubs
      climdata <- read.csv("input/climate_TL_2007_2008_noleap.csv",header  =  TRUE)
      twq <- 9.324752
      mean_LAI <- mean(evergreen_LAI_data)
      mean_TFN <- mean(evergreen_TFN_data)
      obs_LAI <- evergreen_LAI_data
      obs_TFN <-  evergreen_TFN_data
      Lat <- 68
      leafgrow_doy <- 160 - grow_season_increase/2
      leafdrop_doy <- 260 + grow_season_increase/2
      lca <- 85
      LL <- 2.191781
      ramp_up <- 60
      ramp_down <- 120
      ramp_down_doy <- leafgrow_doy+ramp_up
      LAI <- seq(1,4,0.01)
      N <- seq(2,8,0.01)
      LAI <- c(0.05,seq(0.1,1,0.05),seq(1.25,8,0.1))
      N <- c(0.05,seq(0.1,1,0.05),seq(1.25,16,0.1))
      m <- 0.0000000000001
      Q10 <- 2.0
      if(biome_acm  ==  FALSE){
        load("input/acm_recal_with_spa_200pixels_continuous_parameters.RData")
      }else{
        load("input/acm_recal_with_spa_continuous_arctic_optimum_for_arctic_parameters.RData")        
      }
    }else if(biome[kk]  ==  3){
      #Tropical broadleaf evergreen La Selva
      climdata <- read.csv("input/climate_LaSelva_2007.csv")
      twq <- 27.15946
      mean_LAI <- mean(tropic_LAI_data)
      mean_TFN <- mean(tropic_TFN_data)
      obs_LAI <- tropic_LAI_data
      obs_TFN <-  tropic_TFN_data
      Lat <- 10.42
      leafgrow_doy <- 1
      leafdrop_doy <- 365
      lca <- 44    
      LL <- 1.328767
      ramp_up <- 20
      ramp_down <- 20
      ramp_down_doy <- leafdrop_doy-ramp_down
      LAI <- c(0.05,seq(0.1,1,0.1),seq(1.25,13,0.1))
      N <- c(0.05,seq(0.1,1,0.1),seq(1.25,26,0.1))
      m <- 0.00000000001
      Q10 <- 2.0
      load("input/acm_recal_with_spa_200pixels_continuous_parameters.RData")
    }
    
    a1 <- median(parameters[1,,])
    a2 <- median(parameters[2,,])
    a3 <- median(parameters[3,,])
    a4 <- median(parameters[4,,])
    a5 <- median(parameters[5,,])
    a6 <- median(parameters[6,,])
    a7 <- median(parameters[7,,])
    a8 <- median(parameters[8,,])
    a9 <- median(parameters[9,,])
    a10 <- median(parameters[10,,])
    psid <- median(parameters[11,,])
    rtot <- median(parameters[12,,])
    acm_par <- c(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10)
    
    climdata$TEMP_MAX <- climdata$TEMP_MAX + temp_increase
    climdata$TEMP_MIN <- climdata$TEMP_MIN + temp_increase
    climdata$CO2 <- climdata$CO2 + co2_increase
    
    LAI_N <- cbind(obs_LAI,obs_TFN)
    
    # RYAN RESPIRATION 
    resp_type <- 3
    resp_parm1 <- 0.0106 #2.525e-6
    resp_parm2 <- 1.0  
    resp_parm3_sens <- 0.0
    base_resp_temp <- 20 #Reich et al. 2008 normalized to 20C
    
    gpp <- array(-99, dim  =  c(length(LAI_N[, 1]), 1))
    ra_mass <- array(-99, dim  =  c(length(LAI_N[, 1]), 1))
    NetPsn_mass <- array(-99, dim  =  c(length(LAI_N[, 1]), 1))
    leafC_return_mass <- array(-99, dim  =  c(length(LAI_N[,1]), 1))
    leafN_return_mass <- array(-99, dim  =  c(length(LAI_N[,1]), 1))
    leafC_used_in_return <- array(-99, dim  =  c(length(LAI_N[, 1]), 1))
    leafCN_used_in_return <- array(-99, dim  =  c(length(LAI_N[, 1]), 1))
    optimal <- array(0, dim  =  c(length(LAI_N[, 1]), 1))
    optimal_in_data <- array(0, dim  =  c(length(LAI_N[, 1]), 1))
    t_leaf_used <- array(0, dim  =  c(length(LAI_N[, 1]), 1))
    leaflife_span <- array(0, dim  =  c(length(LAI_N[, 1]), 1))
    leaf_allocation <- array(-99, dim  =  c(length(LAI_N[, 1]), 1))
    leaf_horizon <- array(-99, dim = c(length(LAI_N[, 1]), 1))
    leaf_growth_respiration <- array(-99, dim  =  c(length(LAI_N[, 1]), 1))
    
    LAI_N_dataframe <- data.frame(gpp, ra_mass, NetPsn_mass, leafC_return_mass, 
                                  leafN_return_mass, leafC_used_in_return,
                                  leafCN_used_in_return, t_leaf_used,optimal,
                                  optimal_in_data, leaflife_span,
                                  leaf_allocation, leaf_horizon, leaf_growth_respiration)
    
    constants <- c(Lat, psid, rtot, leafgrow_doy, leafdrop_doy, max_add_c, base_resp_temp, 
                   lca, LL, ncohorts, nyears, resp_type, twq)
    par <- c(acm_par, resp_parm1, resp_parm2, Q10, grow_resp, psn_temp_thres, m, ramp_up,
             ramp_down, ramp_down_doy, resp_parm3_sens)
    
    
    dyn.load(fortran_code)
    fortran_output <- .Fortran("LAI_N_Rinterface",
                               LAI_N_dataframe  =  as.matrix(LAI_N_dataframe),
                               par  =  as.matrix(par),
                               LAI_N  =  as.matrix(LAI_N),
                               climdata  =  as.matrix(climdata),
                               constants  =  as.matrix(constants),
                               LAI_N_data_frame_dims  =  as.integer(c(dim(LAI_N_dataframe)[1],
                                                                      dim(LAI_N_dataframe)[2])),
                               LAI_N_dims  =  as.integer(c(dim(LAI_N)[1],
                                                           dim(LAI_N)[2])),
                               climdata_dims  =  as.integer(c(dim(climdata)[1],
                                                              dim(climdata)[2])))                            
    
    par_out_ryan <- array(fortran_output$LAI_N_dataframe, dim = c(dim(LAI_N_dataframe)[1],
                                                                  dim(LAI_N_dataframe)[2]))
    
    # AKTIN RESPIRATION 
    resp_type <- 2
    resp_parm1 <- 1.7560
    resp_parm2 <- 0.2061 
    resp_parm3_sens <- 0.0402
    base_resp_temp <- 25 
    
    constants <- c(Lat, psid, rtot, leafgrow_doy, leafdrop_doy, max_add_c, base_resp_temp, 
                   lca, LL, ncohorts, nyears, resp_type, twq)
    par <- c(acm_par, resp_parm1, resp_parm2, Q10, grow_resp, psn_temp_thres, m, ramp_up,
             ramp_down, ramp_down_doy, resp_parm3_sens)
    
    gpp <- array(-99, dim  =  c(length(LAI_N[, 1]), 1))
    ra_mass <- array(-99, dim  =  c(length(LAI_N[, 1]), 1))
    NetPsn_mass <- array(-99, dim  =  c(length(LAI_N[, 1]), 1))
    leafC_return_mass <- array(-99, dim  =  c(length(LAI_N[,1]), 1))
    leafN_return_mass <- array(-99, dim  =  c(length(LAI_N[,1]), 1))
    leafC_used_in_return <- array(-99, dim  =  c(length(LAI_N[, 1]), 1))
    leafCN_used_in_return <- array(-99, dim  =  c(length(LAI_N[, 1]), 1))
    optimal <- array(0, dim  =  c(length(LAI_N[, 1]), 1))
    optimal_in_data <- array(0, dim  =  c(length(LAI_N[, 1]), 1))
    t_leaf_used <- array(0, dim  =  c(length(LAI_N[, 1]), 1))
    leaflife_span <- array(0, dim  =  c(length(LAI_N[, 1]), 1))
    leaf_allocation <- array(-99, dim  =  c(length(LAI_N[, 1]), 1))
    leaf_horizon <- array(-99, dim = c(length(LAI_N[, 1]), 1))
    leaf_growth_respiration <- array(-99, dim  =  c(length(LAI_N[, 1]), 1))
    
    LAI_N_dataframe <- data.frame(gpp, ra_mass, NetPsn_mass, leafC_return_mass,
                                  leafN_return_mass, leafC_used_in_return,
                                  leafCN_used_in_return, t_leaf_used, optimal,
                                  optimal_in_data, leaflife_span, leaf_allocation,
                                  leaf_horizon, leaf_growth_respiration)
    
    fortran_output <- .Fortran("LAI_N_Rinterface",
                               LAI_N_dataframe  =  as.matrix(LAI_N_dataframe),
                               par  =  as.matrix(par),
                               LAI_N  =  as.matrix(LAI_N),
                               climdata  =  as.matrix(climdata),
                               constants  =  as.matrix(constants),
                               LAI_N_data_frame_dims  =  as.integer(c(dim(LAI_N_dataframe)[1],
                                                                      dim(LAI_N_dataframe)[2])),
                               LAI_N_dims  =  as.integer(c(dim(LAI_N)[1],
                                                           dim(LAI_N)[2])),
                               climdata_dims  =  as.integer(c(dim(climdata)[1],
                                                              dim(climdata)[2])))    
    
    par_out_atkin  =  array(fortran_output$LAI_N_dataframe, dim = c(dim(LAI_N_dataframe)[1],dim(LAI_N_dataframe)[2]))
    
    # REICH RESPIRATION
    resp_type  =  0  
    resp_parm1  =  0.691  #All leaves from Reich Table
    resp_parm2  =  1.639 #All leaves from Reich Table
    resp_parm3_sens  =  0.0
    base_resp_temp  =  20 #Reich et al. 2008 normalized to 20C
    
    constants <- c(Lat, psid, rtot, leafgrow_doy, leafdrop_doy, max_add_c, base_resp_temp, 
                   lca, LL, ncohorts, nyears, resp_type, twq)
    par <- c(acm_par, resp_parm1, resp_parm2, Q10, grow_resp, psn_temp_thres, m, ramp_up,
             ramp_down, ramp_down_doy, resp_parm3_sens)
    
    gpp <- array(-99, dim  =  c(length(LAI_N[, 1]), 1))
    ra_mass <- array(-99, dim  =  c(length(LAI_N[, 1]), 1))
    NetPsn_mass <- array(-99, dim  =  c(length(LAI_N[, 1]), 1))
    leafC_return_mass <- array(-99, dim  =  c(length(LAI_N[,1]), 1))
    leafN_return_mass <- array(-99, dim  =  c(length(LAI_N[,1]), 1))
    leafC_used_in_return <- array(-99, dim  =  c(length(LAI_N[, 1]), 1))
    leafCN_used_in_return <- array(-99, dim  =  c(length(LAI_N[, 1]), 1))
    optimal <- array(0, dim  =  c(length(LAI_N[, 1]), 1))
    optimal_in_data <- array(0, dim  =  c(length(LAI_N[, 1]), 1))
    t_leaf_used <- array(0, dim  =  c(length(LAI_N[, 1]), 1))
    leaflife_span <- array(0, dim  =  c(length(LAI_N[, 1]), 1))
    leaf_allocation <- array(-99, dim  =  c(length(LAI_N[, 1]), 1))
    leaf_horizon <- array(-99, dim = c(length(LAI_N[, 1]), 1))
    leaf_growth_respiration <- array(-99, dim  =  c(length(LAI_N[, 1]), 1))
    
    LAI_N_dataframe <- data.frame(gpp, ra_mass, NetPsn_mass, leafC_return_mass,
                                  leafN_return_mass, leafC_used_in_return,
                                  leafCN_used_in_return, t_leaf_used, optimal,
                                  optimal_in_data, leaflife_span, leaf_allocation,
                                  leaf_horizon, leaf_growth_respiration)
    
    fortran_output <- .Fortran("LAI_N_Rinterface",
                               LAI_N_dataframe  =  as.matrix(LAI_N_dataframe),
                               par  =  as.matrix(par),
                               LAI_N  =  as.matrix(LAI_N),
                               climdata  =  as.matrix(climdata),
                               constants  =  as.matrix(constants),
                               LAI_N_data_frame_dims  =  as.integer(c(dim(LAI_N_dataframe)[1],
                                                                      dim(LAI_N_dataframe)[2])),
                               LAI_N_dims  =  as.integer(c(dim(LAI_N)[1],
                                                           dim(LAI_N)[2])),
                               climdata_dims  =  as.integer(c(dim(climdata)[1],
                                                              dim(climdata)[2])))    
    
    par_out_reich  <-  array(fortran_output$LAI_N_dataframe, dim = c(dim(LAI_N_dataframe)[1],dim(LAI_N_dataframe)[2]))
    
    save(par_out_atkin,par_out_ryan, par_out_reich, 
         file  =  paste0("output/Rm_at_obs_LAI_TCN_biome", biome[kk], ".Rdata"))
  }
}