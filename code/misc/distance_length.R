setwd("/Volumes/Work_drive/prj/BAV_TAV/data/raw_internal/Interaction_march_corrected/")

BAVTAV.Proximities.Probe_Distal <- read.delim("BAVTAV.Proximities.Probe_Distal_SP4_p001_filtered.txt")

head (BAVTAV.Proximities.Probe_Distal)

BAVTAV.Proximities.Probe_Distal$breadth = (BAVTAV.Proximities.Probe_Distal$Interactor_End - 
                                             BAVTAV.Proximities.Probe_Distal$Interactor_Start)

BAVTAV.Proximities.Probe_Distal_SP4_p05_filtered <- read.delim("BAVTAV.Proximities.Probe_Distal_SP4_p05_filtered.txt")

BAVTAV.Proximities.Probe_Distal_SP4_p05_filtered$breadth = (BAVTAV.Proximities.Probe_Distal_SP4_p05_filtered$Interactor_End-
                                                              BAVTAV.Proximities.Probe_Distal_SP4_p05_filtered$Interactor_Start)



BAVTAV.Proximities.Probe_Distal_SP4_p001_distance_length_filtered <- read.delim("BAVTAV.Proximities.Probe_Distal_SP4_p001_distance_length_filtered.txt",
                                                                                header= FALSE)

BAVTAV.Proximities.Probe_Distal_SP4_p001_distance_length_filtered$breadth = BAVTAV.Proximities.Probe_Distal_SP4_p001_distance_length_filtered$V11 -
                                                                              BAVTAV.Proximities.Probe_Distal_SP4_p001_distance_length_filtered$V10


plot(
     BAVTAV.Proximities.Probe_Distal_SP4_p05_filtered$breadth,
     abs(BAVTAV.Proximities.Probe_Distal_SP4_p05_filtered$distance),
     pch = 19,  cex= 0.75,
     col ="red")
     
plot(
  BAVTAV.Proximities.Probe_Distal$breadth,
  abs(BAVTAV.Proximities.Probe_Distal$distance),
  pch = 19,  cex= 0.75,
  col ="red")

plot(
  BAVTAV.Proximities.Probe_Distal_SP4_p001_distance_length_filtered$breadth,
  abs(BAVTAV.Proximities.Probe_Distal_SP4_p001_distance_length_filtered$V12),
  pch = 19,  cex= 0.75,
  col ="red")


plot(BAVTAV.Proximities.Probe_probe_SP4_p001_filtered$abs.Distance,
     pch =20, cex=0.2,
     col ="darkgreen")

plot(  abs(BAVTAV.Proximities.Probe_Distal$distance),
       pch = 19,  cex= 0.75,
       col ="red")
)


plot(BAVTAV.Proximities.Probe_probe_SP4_p05_filtered$abs.Distance,
     pch =20, cex=0.2,
     col = "darkgreen")

plot(BAVTAV.Proximities.Probe_probe_SP4_p001_filtered$abs.Distance,
     pch =20, cex=0.2,
     col ="darkgreen")