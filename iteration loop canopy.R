# call LeafEB_LAI function
# Including canopy depth
# Note: this formulation follows Leuning et al. (1995) PC&E 18:1183-1200

# Checking by comparing to alternative functions:
# PhotosynEB from plantecophys
# Lian paper formulation
# tealeaves

library(plantecophys)
library(tealeaves)

# Environmental settings
Tair <- 35
PPFD <- 1500
RH <- 50
VPD <- RHtoVPD(RH,Tair)
Wind <- 2 
Wleaf <- 0.2
LeafAbsPAR <- 0.85
LeafAbsNIR <- 0.3
LAI <- 0.1
kd <- 0.5
kw = 0.5

# Iteration loop
Tleaf <- gs <- ET <- Rnetiso <- c()
Tleaf[1] <- Tair

for (i in 1:20) {
  
  # Call Photosyn to get gs at given Tleaf
  AT <- Photosyn(Tleaf = Tleaf[i], PPFD = LeafAbsPAR*PPFD*kd*exp(-kd*LAI), 
                 VPD = VPD)
  gs[i] <- AT$GS
  ET[i] <- AT$ELEAF
  # Call LeafEB with Tleaf and gs to get new Tleaf
  EB <- LeafEB_LAI(Tleaf = Tleaf[i], Tair = Tair, gs = AT$GS,
               PPFD = PPFD, VPD = VPD, Wind = Wind, Wleaf = Wleaf,
               LeafAbsPAR = LeafAbsPAR, LeafAbsNIR = LeafAbsNIR,
               LAI = LAI,kd = kd)
  Rnetiso[i] <- EB$Rnetiso
  Tleaf[i+1] <- EB$Tleaf2
  
}

# Check outputs - with lower radiation, Tleaf < Tair
Tleaf
gs
ET
Rnetiso

# Compare with other formulations

# tealeaves - this appears consistent with PhotosynEB
# Set parameters to above inputs
leaf_par <- make_leafpar(
  replace = list(
    leafsize = set_units(Wleaf, "m"),
    abs_s = set_units((LeafAbsPAR + LeafAbsNIR)/2, "1"),
    g_sw = set_units(gs[10]*10, "umol/m^2/s/Pa")
  )
)
UMOLPERJ <- 4.57 # umol J-1 
RH <- VPDtoRH(VPD,Tair)
enviro_par <- make_enviropar(
  replace = list(
    T_air = set_units(Tk(Tair), "K"),
    wind = set_units(Wind*exp(-kw*LAI), "m/s"),
    S_sw = set_units(PPFD/UMOLPERJ*2*kd*exp(-kd*LAI), "W/m2"),
    RH = set_units(RH/100, "1")
  )
)
constants <- make_constants()
# Solve for T_leaf over a range of T_air
T_leaves <- tleaves(leaf_par, enviro_par, constants,
                    quiet = TRUE)
drop_units(T_leaves$T_air) - 273.15
drop_units(T_leaves$T_leaf) - 273.15
T_leaves$g_sw
