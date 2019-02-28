# Code to produce the analysis in
#  Thomas et al. 
#
library(plotrix)

working_directory <- "/Users/quinn/Dropbox (VTFRS)/Shared_folders/ACONITE/Manuscript 2/2019 ms/analysis/final_rerun_for_reproducability/ACONITE_canopy/"
figure_directory <- paste(working_directory, "/output/", sep  =  "")
fortran_code <- "src/aconite_canopy_Rinterface.so"
################
# LOAD SCRIPTS
################

setwd(working_directory)
source("src/aconite_canopy_model.R")
source("src/aconite_canopy_at_observations.R")

#######################################
# IMPORT AND ANALYZE LAI AND TFN 
#OBSERVATIONS ACROSS THE 3 CANOPY TYPES
#######################################
data <- read.csv("input/LAI_TFN_combined.csv", header = TRUE)
tropic_LAI_data <- data$LAI[which(data$species_class_id  ==  3)]
tropic_TFN_data <- data$TFN[which(data$species_class_id  ==  3)]

evergreen_LAI_data <-  data$LAI[data$species_class_id  ==  2 &
                                  data$site_code  ==  1]
evergreen_TFN_data <- data$TFN[data$species_class_id  ==  2 & 
                                 data$site_code  ==  1]
decid_LAI_data <- data$LAI[data$species_class_id  ==  1 & 
                             data$site_code  ==  1]
decid_TFN_data <- data$TFN[data$species_class_id  ==  1 & 
                             data$site_code  ==  1]

all_LAI <- c(decid_LAI_data, evergreen_LAI_data, tropic_LAI_data)
all_TFN <- c(decid_TFN_data, evergreen_TFN_data, tropic_TFN_data)

###########################################
#RUN CANOPY MODEL FOR DIFFERENT RESPIRATION MODELS
#resp_type variable code:
#0 = exponential using mass units (Reich 2008)
#1  =  NOT USED
#2  =  using Atkin et al. 2016 linear model,
#3  =  Linear with Ryan 1991
###########################################

#RYAN 1991
aconite_canopy_model(fortran_code,
                     resp_type = 3,
                     tropic_LAI_data,
                     tropic_TFN_data,
                     evergreen_LAI_data,
                     evergreen_TFN_data,
                     decid_LAI_data,
                     decid_TFN_data,
                     output_file = "output/baseline_ryan.Rdata",
                     plot_name = "baseline_ryan",
                     resp_parm1_sens = 0.0106,
                     resp_parm2_sens = 1.0,  
                     resp_parm3_sens = 0.0) #Not used in model


#AKTIN 2016
aconite_canopy_model(fortran_code,
                     resp_type = 2,
                     tropic_LAI_data,
                     tropic_TFN_data,
                     evergreen_LAI_data,
                     evergreen_TFN_data,
                     decid_LAI_data,
                     decid_TFN_data,
                     output_file = "output/baseline_atkin.Rdata",
                     plot_name = "baseline_atkin",
                     resp_parm1_sens = 1.7560,
                     resp_parm2_sens = 0.2061,  
                     resp_parm3_sens = 0.0402) #Not used in model

#REICH 2008
aconite_canopy_model(fortran_code,
                     resp_type = 0,
                     tropic_LAI_data,
                     tropic_TFN_data,
                     evergreen_LAI_data,
                     evergreen_TFN_data,
                     decid_LAI_data,
                     decid_TFN_data,
                     output_file = "output/baseline_reich.Rdata",
                     plot_name = "baseline_reich",
                     resp_parm1_sens = 0.691,  #All leaves from Reich Table X
                     resp_parm2_sens = 1.639, #All leaves from Reich Table X 
                     resp_parm3_sens = 0.0) #Not used in model

##############################################
#RUN CANOPY MODEL FOR ALL RESPIRATION MODELS 
# AT OBSERVED LAI AND TCN
#############################################

aconite_canopy_at_observations()

##############
#Figure 1
###############

#Diagram in Lucid chart

###############
#FIGURE 2
##############
ylabel <- expression(paste("Total canopy nitrogen (g N ", m^-2, ")", sep = ""))
xlabel <- expression(paste("Leaf area index (",m^2," ",m^-2,")",sep = ""))
pdf("output/Figure2_Rm_flux.pdf",height  =  7, width  =  7)
par(mfrow = c(3, 3), mar  =  c(1.1, 1.1, 1.1, 1.1), oma  =  c(5, 6, 2, 0))
#############
load("output/baseline_ryan_tropical_0.0106_1_contours.Rdata")
contour(LAI,N, ra_mass, levels = seq(0, 7000, 500))
mtext(ylabel, side = 2, adj  =  0.6, line  =  2.5, cex = 0.7)
mtext("Tropical Evergreen", side = 2, adj  =  0.6, line  =  4.5)
mtext("Ryan 1991", side = 3, adj  =  0.6, line  =  1.0)
text(0.5, 24.8, "(a)")
points(obs_LAI, obs_TFN, pch = 20, cex = 2)
#############
load("output/baseline_atkin_tropical_1.756_0.2061_contours.Rdata")
contour(LAI, N, ra_mass, levels = seq(0, 7000, 500))
points(obs_LAI, obs_TFN, pch = 20, cex = 2)
mtext("Aktin et al. 2016", side = 3, adj  =  0.6, line  =  1.0)
text(0.5, 24.8, "(b)")
#############
load("output/baseline_reich_tropical_0.691_1.639_contours.Rdata")
contour(LAI,N, ra_mass, levels = seq(0, 7000, 500))
points(obs_LAI, obs_TFN, pch = 20, cex = 2)
mtext("Reich et al. 2008", side = 3, adj  =  0.6, line  =  1.0)
text(0.5, 24.8, "(c)")
#############
load("output/baseline_ryan_arctic_evergreen_0.0106_1_contours.Rdata")
contour(LAI[which(LAI <=  1.5)], N[which(N <=  3)], ra_mass[which(LAI <=  1.5), which(N <=  3)], levels = seq(0, 100, 10), ylab  =  "Arctic Evergreen")
mtext(ylabel, side = 2, adj  =  0.6, line  =  2.5, cex = 0.7)
mtext("Arctic Evergreen", side = 2, adj  =  0.6, line  =  4.5)
points(obs_LAI, obs_TFN, pch = 20, cex = 2)
text(0.1, 2.8, "(d)")
#############
load("output/baseline_atkin_arctic_evergreen_1.756_0.2061_contours.Rdata")
contour(LAI[which(LAI <=  1.5)], N[which(N <=  3)], ra_mass[which(LAI <=  1.5), which(N <=  3)], levels = seq(0, 100, 10), ylab  =  "Arctic Evergreen")
points(obs_LAI, obs_TFN, pch = 20, cex = 2)
text(0.1, 2.8, "(e)")
#############
load("output/baseline_reich_arctic_evergreen_0.691_1.639_contours.Rdata")
contour(LAI[which(LAI <=  1.5)], N[which(N <=  3)], ra_mass[which(LAI <=  1.5), which(N <=  3)],levels = seq(0, 100, 10), ylab  =  "Arctic Evergreen")
points(obs_LAI, obs_TFN, pch = 20, cex = 2)
text(0.1, 2.8, "(f)")
#############
load("output/baseline_ryan_arctic_decid_0.0106_1_contours.Rdata")
contour(LAI[which(LAI <=  5)], N[which(N <=  10)], ra_mass[which(LAI <=  5), which(N <=  10)], levels = seq(0, 200, 20), ylab  =  "Arctic Decidious", xlab = "LAI")
mtext(ylabel, side = 2, adj  =  0.6, line  =  2.5, cex = 0.7)
mtext("Arctic Deciduous", side = 2, adj  =  0.6, line  =  4.5)
points(obs_LAI, obs_TFN, pch = 20, cex = 2)
mtext(xlabel, side = 1, adj  =  0.6, line  =  2.5, cex = 0.7)
text(0.3, 9.5, "(g)")
#############
load("output/baseline_atkin_arctic_decid_1.756_0.2061_contours.Rdata")
contour(LAI[which(LAI <=  5)], N[which(N <=  10)], ra_mass[which(LAI <=  5), which(N <=  10)],levels = seq(0, 300, 20), ylab  =  "Arctic Decidious")
points(obs_LAI, obs_TFN, pch = 20, cex = 2)
mtext(xlabel, side = 1, adj  =  0.6, line  =  2.5, cex = 0.7)
text(0.3, 9.5, "(h)")
#############
load("output/baseline_reich_arctic_decid_0.691_1.639_contours.Rdata")
contour(LAI[which(LAI <=  5)], N[which(N <=  10)], ra_mass[which(LAI <=  5), which(N <=  10)], levels = seq(0, 200, 20), ylab  =  "Arctic Decidious")
points(obs_LAI, obs_TFN, pch = 20, cex = 2)
mtext(xlabel, side = 1, adj  =  0.6, line  =  2.5, cex = 0.7)
text(0.3, 9.5, "(i)")
#############
dev.off()


############
#FIGURE 3
############

pdf("output/Figure3_Rm_histogram.pdf",height  =  7.8, width  =  7.8)
par(mfrow = c(3, 3),mar  =  c(1.1, 1.1, 1.1, 1.1), oma  =  c(5, 6, 2, 0))
xlab  =  expression(paste("annual respiration (g C ", m^2, " ", yr^-1, ")", sep = ""))
load("output/Rm_at_obs_LAI_TCN_biome3.Rdata")
ylim  =  range(c(density(par_out_reich[, 2])$y, density(par_out_ryan[,  2])$y),  density(par_out_atkin[,  2])$y)
xlim  =  c(0,  max(c(density(par_out_reich[, 2])$x, density(par_out_ryan[, 2])$x), density(par_out_atkin[, 2])$x))
plot(density(par_out_reich[, 2], cut = c(0, 10000)), xlim = xlim, ylim = ylim,  main = "", lty = "dotted")
points(density(par_out_atkin[, 2], cut = c(0, 10000)), type = "l", col = "blue", lty = "dashed")
points(density(par_out_ryan[, 2], cut = c(0, 10000)), type  =  "l", col = "red", lty = "solid")     
mtext("density", side = 2, adj  =  0.6, line  =  2.5, cex = 0.7)
mtext("Tropical Evergeen", side = 2, adj  =  0.6, line  =  4.5)
v  =  expression("Canopy R"[leaf][",m"])
mtext(v, side = 3, adj  =  0.6, line  =  1.0)
text(xlim[1], ylim[2], "\n   (a)")
legend("topright", c("Ryan 1991", "Atkin 2016", "Reich 2008"), col = c("red", "blue", "black"), bty = "n", lty = c("solid", "dashed", "dotted"))

#Blank
plot.new()
plot.new()

load("output/Rm_at_obs_LAI_TCN_biome2.Rdata")
ylim  =  range(c(density(par_out_reich[, 2])$y, density(par_out_ryan[, 2])$y), density(par_out_atkin[, 2])$y)
xlim  =  c(0, max(c(density(par_out_reich[, 2])$x, density(par_out_ryan[, 2])$x), density(par_out_atkin[, 2])$x))
plot(density(par_out_reich[, 2], cut = c(0, 10000)), xlim = xlim, ylim = ylim, xlab  =  "annual respiration (gC m-2 yr-1)", ylab  =  "Arctic Evergreen",  main = "", lty = "dotted")
points(density(par_out_atkin[, 2], cut = c(0, 10000)), type = "l", col = "blue", lty = "dashed")
points(density(par_out_ryan[, 2], cut = c(0, 10000)), type  =  "l", col = "red", lty = "solid")        
mtext("density", side = 2, adj  =  0.6, line  =  2.5, cex = 0.7)
mtext("Arctic Evergreen", side = 2, adj  =  0.6, line  =  4.5)
text(xlim[1], ylim[2], "\n   (b)")

plot.new()
plot.new()

load("output/Rm_at_obs_LAI_TCN_biome1.Rdata")
ylim  =  range(c(density(par_out_reich[, 2])$y, density(par_out_ryan[, 2])$y), density(par_out_atkin[, 2])$y)
xlim  =  c(0, max(c(density(par_out_reich[, 2])$x, density(par_out_ryan[, 2])$x), density(par_out_atkin[, 2])$x))
plot(density(par_out_reich[, 2], cut = c(0, 10000)), xlim = xlim, ylim = ylim, xlab  =  "annual respiration (gC m-2 yr-1)", ylab  =  "Arctic Deciduous",  main = "", lty = "dotted")
points(density(par_out_atkin[, 2], cut = c(0, 10000)), type = "l", col = "blue", lty = "dashed")
points(density(par_out_ryan[, 2], cut = c(0, 10000)), type  =  "l", col = "red", lty = "solid")  
mtext("density", side = 2, adj  =  0.6, line  =  2.5, cex = 0.7)
mtext("Arctic Deciduous", side = 2, adj  =  0.6, line  =  4.5)
mtext(xlab, side = 1, adj  =  0.6, line  =  2.5, cex = 0.7)
text(xlim[1], ylim[2], "\n   (c)")

#Blank
plot.new()
plot.new()

dev.off()

###################
#FIGURE 4
###################

pdf("output/Figure4_other_fluxes.pdf", height  =  7,  width  =  7)
ylabel  =  expression(paste("Total canopy nitrogen (g N ", m^-2, ")", sep = ""))
xlabel  =  expression(paste("Leaf area index (", m^2, " ", m^-2, ")", sep = ""))
par(mfrow = c(3, 3), mar  =  c(1.1, 1.1, 1.1, 1.1), oma  =  c(5, 6, 2, 0))
#####
load("output/baseline_ryan_tropical_0.0106_1_contours.Rdata")
contour(LAI, N, gpp, levels = seq(0, 7000, 500))
mtext(ylabel, side = 2, adj  =  0.6, line  =  2.5, cex = 0.7)
mtext("Tropical Evergeen", side = 2, adj  =  0.6, line  =  4.5)
mtext("GPP", side = 3, adj  =  0.6, line  =  1.0)
text(0.5, 24.8, "(a)")
points(obs_LAI, obs_TFN, pch = 20, cex = 2)
#######
#abline(lm(obs_TFN~obs_LAI))
load("output/baseline_atkin_tropical_1.756_0.2061_contours.Rdata")
contour(LAI, N, leaf_allocation, levels = seq(0, 7000, 50))
mtext(expression("A"[l]), side = 3, adj  =  0.6, line  =  1.0)
text(0.5, 24.8, "(b)")
points(obs_LAI, obs_TFN, pch = 20, cex = 2)
#######
load("output/baseline_reich_tropical_0.691_1.639_contours.Rdata")
contour(LAI, N, leaf_growth_respiration, levels = seq(0, 7000, 50))
mtext(expression("R"[g]), side = 3, adj  =  0.6, line  =  1.0)
text(0.5, 24.8, "(c)")
points(obs_LAI, obs_TFN, pch = 20, cex = 2)
#######
load("output/baseline_ryan_arctic_evergreen_0.0106_1_contours.Rdata")
contour(LAI[which(LAI <=  1.5)], N[which(N <=  3)], gpp[which(LAI <=  1.5), which(N <=  3)], levels = seq(0, 2000, 50))
mtext(ylabel, side = 2, adj  =  0.6, line  =  2.5, cex = 0.7)
mtext("Arctic Evergreen", side = 2, adj  =  0.6, line  =  4.5)
points(obs_LAI, obs_TFN, pch = 20, cex = 2)
text(0.1, 2.8, "(d)")
#######
load("output/baseline_atkin_arctic_evergreen_1.756_0.2061_contours.Rdata")
contour(LAI[which(LAI <=  1.5)], N[which(N <=  3)], leaf_allocation[which(LAI <=  1.5), which(N <=  3)], levels = seq(0, 100, 10))
text(0.1, 2.8, "(e)")
points(obs_LAI, obs_TFN, pch = 20, cex = 2)
#######
load("output/baseline_reich_arctic_evergreen_0.691_1.639_contours.Rdata")
contour(LAI[which(LAI <=  1.5)], N[which(N <=  3)], leaf_growth_respiration[which(LAI <=  1.5), which(N <=  3)], levels = seq(0, 100, 5))
points(obs_LAI, obs_TFN, pch = 20, cex = 2)
text(0.1, 2.8, "(f)")
#######
load("output/baseline_ryan_arctic_decid_0.0106_1_contours.Rdata")
contour(LAI[which(LAI <=  5)], N[which(N <=  10)], gpp[which(LAI <=  5), which(N <=  10)], levels = seq(0, 2000, 100))
mtext(ylabel, side = 2, adj  =  0.6, line  =  2.5, cex = 0.7)
mtext("Arctic Deciduous", side = 2, adj  =  0.6, line  =  4.5)
points(obs_LAI, obs_TFN, pch = 20, cex = 2)
mtext(xlabel, side = 1, adj  =  0.6, line  =  2.5, cex = 0.7)
text(0.3, 9.5, "(g)")
#######
load("output/baseline_atkin_arctic_decid_1.756_0.2061_contours.Rdata")
contour(LAI[which(LAI <=  5)], N[which(N <=  10)], leaf_allocation[which(LAI <=  5), which(N <=  10)], levels = seq(0, 200, 50))
points(obs_LAI, obs_TFN, pch = 20, cex = 2)
mtext(xlabel, side = 1, adj  =  0.6, line  =  2.5, cex = 0.7)
text(0.3, 9.5, "(h)")
#######
load("output/baseline_reich_arctic_decid_0.691_1.639_contours.Rdata")

contour(LAI[which(LAI <=  5)], N[which(N <=  10)], leaf_growth_respiration[which(LAI <=  5), which(N <=  10)], levels = seq(0, 200, 10))
points(obs_LAI, obs_TFN, pch = 20, cex = 2)
mtext(xlabel, side = 1, adj  =  0.6, line  =  2.5, cex = 0.7)
text(0.3, 9.5, "(i)")
#######
dev.off()

#####################
###Figure 5
#####################

pdf("output/Figure5_Rm_GPP_histograms.pdf", height  =  7.8,  width  =  7.8)
par(mfrow = c(3, 3), mar  =  c(1.1, 1.1, 1.1, 1.1), oma  =  c(5, 6, 2, 0))
load("output/Rm_at_obs_LAI_TCN_biome3.Rdata")
ylim  =  range(c(density(par_out_reich[, 2]/par_out_reich[, 1])$y, density(par_out_ryan[, 2]/par_out_ryan[, 1])$y), density(par_out_atkin[, 2]/par_out_atkin[, 1])$y)
xlim  =  range(c(density(par_out_reich[, 2]/par_out_reich[, 1])$x, density(par_out_ryan[, 2]/par_out_ryan[, 1])$x), density(par_out_atkin[, 2]/par_out_atkin[, 1])$x)
plot(density(par_out_reich[, 2]/par_out_reich[, 1]), xlim = xlim, ylim = ylim, main  =  "", lty = "dotted")
points(density(par_out_atkin[, 2]/par_out_atkin[, 1]), type = "l", col = "blue", lty = "dashed")
points(density(par_out_ryan[, 2]/par_out_ryan[, 1]), type  =  "l", col = "red", lty = "solid")  
v  =  expression(paste("Canopy ", R[leaf][",m"],"/GPP", sep = " "))
mtext(v, side = 3, adj  =  0.6, line  =  1.0)
text(xlim[1], ylim[2], "\n   (a)")
mtext("density", side = 2, adj  =  0.6, line  =  2.5, cex = 0.7)
mtext("Tropical Evergeen", side = 2, adj  =  0.6, line  =  4.5)
legend("topright", c("Ryan 1991", "Atkin 2016", "Reich 2008"), col = c("red", "blue", "black"), bty = "n", lty = c("solid", "dashed", "dotted"))


#Blank
plot.new()
plot.new()

load("output/Rm_at_obs_LAI_TCN_biome2.Rdata")
ylim  =  range(c(density(par_out_reich[, 2]/par_out_reich[, 1])$y, density(par_out_ryan[, 2]/par_out_ryan[, 1])$y), density(par_out_atkin[, 2]/par_out_atkin[, 1])$y)
xlim  =  range(c(density(par_out_reich[, 2]/par_out_reich[, 1])$x, density(par_out_ryan[, 2]/par_out_ryan[, 1])$x), density(par_out_atkin[, 2]/par_out_atkin[, 1])$x)
plot(density(par_out_reich[, 2]/par_out_reich[, 1]), xlim = xlim, ylim = ylim, xlab  =  "Rm/GPP Ratio",  main = "", lty = "dotted")
points(density(par_out_atkin[, 2]/par_out_atkin[, 1]), type = "l", col = "blue", lty = "dashed")
points(density(par_out_ryan[, 2]/par_out_ryan[, 1]), type  =  "l", col = "red", lty = "solid")    
mtext("density", side = 2, adj  =  0.6, line  =  2.5, cex = 0.7)
mtext("Arctic Evergreen",  side = 2, adj  =  0.6, line  =  4.5)
text(xlim[1],  ylim[2],  "\n   (b)")

plot.new()
plot.new()

load("output/Rm_at_obs_LAI_TCN_biome1.Rdata")
ylim = range(c(density(par_out_reich[, 2] / par_out_reich[,  1])$y, density(par_out_ryan[, 2] / par_out_ryan[, 1])$y), density(par_out_atkin[, 2] / par_out_atkin[, 1])$y)
xlim = range(c(density(par_out_reich[, 2] / par_out_reich[, 1])$x, density(par_out_ryan[, 2] / par_out_ryan[, 1])$x), density(par_out_atkin[, 2] / par_out_atkin[, 1])$x)
plot(density(par_out_reich[, 2] / par_out_reich[, 1]), xlim = xlim, ylim = ylim, xlab  =  "Rm/GPP Ratio", main = "", lty = "dotted")
points(density(par_out_atkin[, 2] / par_out_atkin[, 1]), type = "l", col = "blue", lty = "dashed")
points(density(par_out_ryan[, 2] / par_out_ryan[, 1]), type  =  "l", col = "red", lty = "solid") 
text(xlim[1], ylim[2], "\n   (c)")
mtext("density", side = 2, adj  =  0.6, line  =  2.5, cex = 0.7)
mtext("Arctic Deciduous", side = 2, adj  =  0.6, line  =  4.5)
mtext("ratio", side = 1, adj  =  0.6, line  =  2.5, cex = 0.7)
#Blank
plot.new()
plot.new()

dev.off()

##############
# Figure 6
##############

ylabel = expression(paste("Total canopy nitrogen (g N " ,m^-2, ")" ,sep = ""))
xlabel = expression(paste("Leaf area index (", m^2, " " ,m^-2, ")" ,sep = ""))
pdf("output/Figure6_NCCE_flux.pdf", height = 7, width  =  7)
par(mfrow = c(3, 3), mar = c(1.1, 1.1, 1.1, 1.1), oma  =  c(5, 6, 2, 0))
#########
load("output/baseline_ryan_tropical_0.0106_1_contours.Rdata")
contour(LAI, N, ((gpp - ra_mass) - leaf_growth_respiration - leaf_allocation), levels = seq(0, 7000, 200))
points(leafC.linelist[[zero_iso_leafC]]$x, leafC.linelist[[zero_iso_leafC]]$y, type = "l", col = "red", lwd = 2, lty = "dotted")
points(leafN.linelist[[zero_iso_leafN]]$x, leafN.linelist[[zero_iso_leafN ]]$y, type = "l", col = "blue", lwd = 2, lty = "longdash")
mtext(ylabel, side = 2, adj = 0.6, line  =  2.5, cex = 0.7)
mtext("Tropical Evergeen", side = 2, adj  =  0.6, line  =  4.5)
mtext("Ryan 1991", side = 3, adj  =  0.6, line  =  1.0)
text(0.5, 24.8, "(a)")
points(obs_LAI, obs_TFN, pch = 20, cex = 2)
abline(lm(obs_TFN ~ obs_LAI))
#########
load("output/baseline_atkin_tropical_1.756_0.2061_contours.Rdata")
contour(LAI, N, ((gpp - ra_mass) - leaf_growth_respiration - leaf_allocation), levels = seq(0, 7000, 200))
points(leafC.linelist[[zero_iso_leafC]]$x, leafC.linelist[[zero_iso_leafC]]$y, type = "l", col = "red", lwd = 2, lty = "dotted")
points(leafN.linelist[[zero_iso_leafN]]$x, leafN.linelist[[zero_iso_leafN ]]$y, type = "l", col = "blue", lwd = 2, lty = "longdash")
points(obs_LAI, obs_TFN, pch = 20, cex = 2)
mtext("Aktin et al. 2016", side = 3, adj  =  0.6, line  =  1.0)
text(0.5 ,24.8, "(b)")
abline(lm(obs_TFN ~ obs_LAI))
#########
load("output/baseline_reich_tropical_0.691_1.639_contours.Rdata")
contour(LAI, N, ((gpp - ra_mass) - leaf_growth_respiration - leaf_allocation), levels = seq(0, 7000, 200))
points(leafC.linelist[[zero_iso_leafC]]$x, leafC.linelist[[zero_iso_leafC]]$y, type = "l", col = "red", lwd = 2, lty = "dotted")
points(leafN.linelist[[zero_iso_leafN]]$x, leafN.linelist[[zero_iso_leafN ]]$y, type = "l", col = "blue", lwd = 2, lty = "longdash")
points(obs_LAI, obs_TFN, pch = 20, cex = 2)
mtext("Reich et al. 2008", side = 3, adj  =  0.6, line  =  1.0)
text(0.5, 24.8, "(c)")
abline(lm(obs_TFN ~ obs_LAI))
#########
load("output/baseline_ryan_arctic_evergreen_0.0106_1_contours.Rdata")
contour(LAI, N, ((gpp - ra_mass) - leaf_growth_respiration - leaf_allocation), levels = seq(0, 7000, 100), ylab  =  "Arctic Evergreen")
points(leafC.linelist[[zero_iso_leafC]]$x, leafC.linelist[[zero_iso_leafC]]$y, type = "l", col = "red", lwd = 2, lty = "dotted")
points(leafN.linelist[[zero_iso_leafN]]$x, leafN.linelist[[zero_iso_leafN ]]$y, type = "l", col = "blue", lwd = 2, lty = "longdash")
mtext(ylabel, side = 2, adj  =  0.6, line  =  2.5, cex = 0.7)
mtext("Arctic Evergreen", side = 2, adj  =  0.6, line  =  4.5)
points(obs_LAI, obs_TFN, pch = 20, cex = 2)
text(0.4, 19, "(d)")
abline(lm(obs_TFN ~ obs_LAI))
#########
load("output/baseline_atkin_arctic_evergreen_1.756_0.2061_contours.Rdata")
contour(LAI, N, ((gpp - ra_mass) - leaf_growth_respiration - leaf_allocation), levels = seq(0, 7000, 100))
points(leafC.linelist[[zero_iso_leafC]]$x, leafC.linelist[[zero_iso_leafC]]$y, type = "l", col = "red", lwd = 2, lty = "dotted")
points(leafN.linelist[[zero_iso_leafN]]$x, leafN.linelist[[zero_iso_leafN ]]$y, type = "l", col = "blue", lwd = 2, lty = "longdash")
text(0.4, 19, "(e)")
points(obs_LAI, obs_TFN, pch = 20, cex = 2)
abline(lm(obs_TFN ~ obs_LAI))
#########
load("output/baseline_reich_arctic_evergreen_0.691_1.639_contours.Rdata")
contour(LAI, N, ((gpp - ra_mass) - leaf_growth_respiration - leaf_allocation), levels = seq(0, 7000, 100))
points(leafC.linelist[[zero_iso_leafC]]$x, leafC.linelist[[zero_iso_leafC]]$y, type = "l", col = "red", lwd = 2, lty = "dotted")
points(leafN.linelist[[zero_iso_leafN]]$x, leafN.linelist[[zero_iso_leafN ]]$y, type = "l", col = "blue", lwd = 2, lty = "longdash")
points(obs_LAI, obs_TFN, pch = 20, cex = 2)
text(0.4, 19, "(f)")
abline(lm(obs_TFN ~ obs_LAI))
#########
load("output/baseline_ryan_arctic_decid_0.0106_1_contours.Rdata")
contour(LAI, N, ((gpp-ra_mass) - leaf_growth_respiration - leaf_allocation), levels = seq(0, 7000, 100),ylab  =  "Arctic Decidious")
points(leafC.linelist[[zero_iso_leafC]]$x, leafC.linelist[[zero_iso_leafC]]$y, type = "l", col = "red", lwd = 2, lty = "dotted")
points(leafN.linelist[[zero_iso_leafN]]$x, leafN.linelist[[zero_iso_leafN ]]$y, type = "l", col = "blue", lwd = 2, lty = "longdash")
mtext(ylabel, side = 2, adj  =  0.6, line  =  2.5, cex = 0.7)
mtext("Arctic Deciduous", side = 2, adj  =  0.6, line  =  4.5)
points(obs_LAI, obs_TFN, pch = 20, cex = 2)
mtext(xlabel, side = 1, adj  =  0.6, line  =  2.5, cex = 0.7)
text(0.4, 19, "(g)")
abline(lm(obs_TFN ~ obs_LAI))
#########
load("output/baseline_atkin_arctic_decid_1.756_0.2061_contours.Rdata")
contour(LAI, N, ((gpp-ra_mass) - leaf_growth_respiration - leaf_allocation), levels = seq(0, 7000, 100))
points(leafC.linelist[[zero_iso_leafC]]$x, leafC.linelist[[zero_iso_leafC]]$y, type = "l", col = "red", lwd = 2, lty = "dotted")
points(leafN.linelist[[zero_iso_leafN]]$x, leafN.linelist[[zero_iso_leafN ]]$y, type = "l", col = "blue", lwd = 2, lty = "longdash")
points(obs_LAI, obs_TFN, pch = 20, cex = 2)
mtext(xlabel, side = 1, adj  =  0.6, line  =  2.5, cex = 0.7)
text(0.4, 19, "(h)")
abline(lm(obs_TFN ~ obs_LAI))
#########
load("output/baseline_reich_arctic_decid_0.691_1.639_contours.Rdata")
contour(LAI, N, ((gpp-ra_mass) - leaf_growth_respiration - leaf_allocation), levels = seq(0, 7000, 100))
points(leafC.linelist[[zero_iso_leafC]]$x,leafC.linelist[[zero_iso_leafC]]$y, type = "l", col = "red", lwd = 2, lty = "dotted")
points(leafN.linelist[[zero_iso_leafN]]$x, leafN.linelist[[zero_iso_leafN ]]$y, type = "l", col = "blue", lwd = 2, lty = "longdash")
points(obs_LAI, obs_TFN, pch = 20, cex = 2)
mtext(xlabel, side = 1, adj  =  0.6, line  =  2.5, cex = 0.7)
text(0.4, 19, "(i)")
abline(lm(obs_TFN ~ obs_LAI))
#########
dev.off()

##########
#Rm/GPP Ratios in text
##########

out_table <- array(NA, dim  =  c(6, 3))

load("output/Rm_at_obs_LAI_TCN_biome3.Rdata")

tropical_aktin_Rmleaf <- mean(par_out_atkin[, 2]/par_out_atkin[, 1])
tropical_aktin_Rleaf <- mean((par_out_atkin[,2]+par_out_atkin[,14])/par_out_atkin[, 1])

tropical_ryan_Rmleaf <- mean(par_out_ryan[, 2] / par_out_ryan[, 1])
tropical_ryan_Rleaf <- mean((par_out_ryan[, 2] + par_out_ryan[, 14]) / par_out_ryan[, 1])

tropical_reich_Rmleaf <- mean(par_out_reich[, 2] / par_out_reich[, 1])
tropical_reich_Rleaf <- mean((par_out_reich[, 2] + par_out_reich[, 14])/ par_out_reich[, 1])

out_table[1:2, 1] <- c(tropical_aktin_Rmleaf, tropical_aktin_Rleaf)
out_table[1:2, 2] <- c(tropical_ryan_Rmleaf, tropical_ryan_Rleaf)
out_table[1:2, 3] <- c(tropical_reich_Rmleaf, tropical_reich_Rleaf)


tropical_aktin_Rmleaf_mean <- mean(par_out_atkin[, 2])
tropical_ryan_Rmleaf_mean <-mean(par_out_ryan[, 2])
tropical_reich_Rmleaf_mean <-mean(par_out_reich[, 2])

load("output/Rm_at_obs_LAI_TCN_biome2.Rdata")

ever_aktin_Rmleaf <- mean(par_out_atkin[, 2] / par_out_atkin[, 1])
ever_aktin_Rleaf <- mean((par_out_atkin[, 2] + par_out_atkin[, 14]) / par_out_atkin[, 1])

ever_ryan_Rmleaf <- mean(par_out_ryan[, 2] / par_out_ryan[, 1])
ever_ryan_Rleaf <- mean((par_out_ryan[, 2] + par_out_ryan[, 14]) / par_out_ryan[, 1])

ever_reich_Rmleaf <- mean(par_out_reich[, 2] / par_out_reich[, 1])
ever_reich_Rleaf <- mean((par_out_reich[, 2] + par_out_reich[, 14]) / par_out_reich[, 1])

out_table[3:4, 1] <- c(ever_aktin_Rmleaf, ever_aktin_Rleaf)
out_table[3:4, 2] <- c(ever_ryan_Rmleaf, ever_ryan_Rleaf)
out_table[3:4, 3] <- c(ever_reich_Rmleaf, ever_reich_Rleaf)

ever_aktin_Rmleaf_mean <- mean(par_out_atkin[, 2])
ever_ryan_Rmleaf_mean <-mean(par_out_ryan[, 2])
ever_reich_Rmleaf_mean <-mean(par_out_reich[, 2])

load("output/Rm_at_obs_LAI_TCN_biome1.Rdata")

decid_aktin_Rmleaf <- mean(par_out_atkin[, 2] / par_out_atkin[, 1])
decid_aktin_Rleaf <- mean((par_out_atkin[, 2] + par_out_atkin[, 14]) / par_out_atkin[, 1])

decid_ryan_Rmleaf <- mean(par_out_ryan[, 2] / par_out_ryan[, 1])
decid_ryan_Rleaf <- mean((par_out_ryan[, 2] + par_out_ryan[, 14]) / par_out_ryan[, 1])

decid_reich_Rmleaf <- mean(par_out_reich[, 2] / par_out_reich[, 1])
decid_reich_Rleaf <- mean((par_out_reich[, 2] + par_out_reich[, 14]) / par_out_reich[, 1])

out_table[5:6, 1] <- c(decid_aktin_Rmleaf, decid_aktin_Rleaf)
out_table[5:6, 2] <- c(decid_ryan_Rmleaf, decid_ryan_Rleaf)
out_table[5:6, 3] <- c(decid_reich_Rmleaf, decid_reich_Rleaf)

out_table<-round(out_table, 2)

decid_aktin_Rmleaf_mean <- mean(par_out_atkin[, 2])
decid_ryan_Rmleaf_mean <-mean(par_out_ryan[, 2])
decid_reich_Rmleaf_mean <-mean(par_out_reich[, 2])

colnames(out_table) <- c("Atkin","Ryan","Reich")
row.names(out_table) <- c("Tropical (Rm,leaf)/GPP", "Tropical (Rleaf)/GPP", 
                          "Evergreen (Rm,leaf)/GPP","Evergreen (Rleaf)/GPP",
                          "Decid (Rm,leaf)/GPP", "decid (Rleaf)/GPP")


(tropical_ryan_Rmleaf_mean-tropical_reich_Rmleaf_mean) / tropical_reich_Rmleaf_mean
(tropical_ryan_Rmleaf_mean-tropical_aktin_Rmleaf_mean) / tropical_aktin_Rmleaf_mean 

(ever_ryan_Rmleaf_mean-ever_reich_Rmleaf_mean) / ever_reich_Rmleaf_mean
(ever_ryan_Rmleaf_mean-ever_aktin_Rmleaf_mean) / ever_aktin_Rmleaf_mean 

(decid_ryan_Rmleaf_mean-decid_reich_Rmleaf_mean) / decid_reich_Rmleaf_mean
(decid_ryan_Rmleaf_mean-decid_aktin_Rmleaf_mean) / decid_aktin_Rmleaf_mean 

############################  
# EXTRACT OPTIMAL LAI AND TFN
############################

canopy_simulations <- c("output/baseline_ryan.Rdata",
                        "output/baseline_atkin.Rdata",
                        "output/baseline_reich.Rdata")

max_LAI <- array(-9999, dim = c(3,3))
max_N <- array(-9999, dim = c(3,3))
for(i in 1:3){
  load(canopy_simulations[i])
  for(j in 1:3){
    max_LAI[i, j] <- max_LAI_sens[j]
    max_N[i, j] <- max_N_sens[j]
  }
}

row.names(max_LAI) <- c("Ryan", "Atkin", "Reich")
row.names(max_N) <- c("Ryan", "Atkin", "Reich")
colnames(max_LAI) <- c("Arctic Evergreen", "Arctic Deciduous", "Tropical")
colnames(max_N) <- c("Arctic Deciduous", "Arctic Evergreen", "Tropical")

max_LAI
max_N

