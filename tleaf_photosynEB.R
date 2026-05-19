
# Functions from LeafEnergyBalance in plantecophys package
# Based on Leuning et al. (1995) PC&E 18:1186-1200
# Some key differences from PhotosynEB, but consistent with MAESPA

# saturation vapour pressure (Pa)
satur <- function (Tair) {
  613.75*exp(17.502*Tair/(240.97+Tair))
} 

# temperature in Kelvin (Degrees)
Tk <- function (Tair) (Tair + 273.15)

## Leaf T formulation from PhotosynEB
LeafEB_2 <- function (Tleaf = 21.5, Tair = 20, gs = 0.15, PPFD = 1500, VPD = 2, 
          Patm = 101, Wind = 2.5, Wleaf = 0.02, 
          StomatalRatio = 1, LeafAbsPAR = 0.85, LeafAbsNIR = 0.3) {

  # PPFD = photosynthetic photon flux density = incident PAR (umol m-2 s-1)  
  # LeafAbsPAR = leaf absorptance of PAR (fraction)
  # LeafAbsNIR = leaf absorptance of NIR (fraction) ### Changed from RAD
  # gs = stomatal conductance to water vapour (mol m-2 s-1)
  # VPD = vapour pressure deficit (kPa)
  # Wind = windspeed (m s-1) - default 2.5
  # Wleaf = leaf width (m) - PFT specific
  # Patm = atmospheric pressure (kPa) - default 101.5
  
# constants
  # Stefan-Boltzmann
Boltz <- 5.67 * 10^-8 #  W m-2 K-4
# Leaf emissivity
# Wang & Leuning (1998) AFM has Emleaf = 0.96. RAD 0.95. Leuning simplifies to 1
Emissivity <- 0.96 
# Heat capacity of air
CPAIR <- 1010 # J kg-1 K-1
# Latent heat of vapourisation at T = 0 deg C
H2OLV0 <- 2501000 # J kg-1
# Molecular weight of water
H2OMW <- 0.018 # kg mol-1
# Molecular weight of air
AIRMA <- 0.029 # kg mol-1
# Conversion of PPFD from umol to J
UMOLPERJ <- 4.57 # umol J-1 
# Molecular diffusivity for heat
DHEAT <- 2.15e-05
# Universal gas constant
RGAS <- 8.314 # J mol-1 K-1

# short calculations
# Density of air
AIRDENS <- Patm * 1000/(287.058 * Tk(Tair)) # kg m-3 
# Latent heat of vapourisation
LHV <- (H2OLV0 - 2365 * Tair) * H2OMW # J kg-1 * kg mol-1 = J mol-1
# Term s, rate of change of saturated vapour pressure with temperature
SLOPE <- (satur(Tair + 0.1) - satur(Tair))/0.1 # Pa K-1

# Conversion for heat conductance from m s-1 to mol m-2 s-1
# Units: 1 J/m3 = 1 Pa. So Pa / (J mol-1 K-1) / K = Pa mol / J = J m-3 mol /J = mol m-3
CMOLAR <- Patm * 1000/(RGAS * Tk(Tair)) # mol m-3

# conductances - all in mol m-2 s-1
# Radiation conductance
# NB in canopy needs to be multiplied by absorptance kd (1-kd*LAI)
# Leuning eqn D7
Gradiation <- 4 * Boltz * Tk(Tair)^3 * Emissivity/(CPAIR * AIRMA)
# Boundary layer conductance for heat from forced convection
# NB in canopy needs to be corrected for decline in wind speed with depth
Gbhforced <- 0.003 * sqrt(Wind/Wleaf) * CMOLAR 
# Grashof number
GRASHOF <- 1.6e+08 * abs(Tleaf - Tair) * (Wleaf^3)
# Boundary layer conductance for heat from free convection
Gbhfree <- 0.5 * DHEAT * (GRASHOF^0.25)/Wleaf * CMOLAR
# Total boundary layer conductance for heat for both sides of the leaf
Gbh <- 2 * (Gbhfree + Gbhforced)
# Boundary layer conductance for water vapour
# StomatalRatio = 1 for amphistomatous leaves, 0.5 for hypostomatous
Gbw <- StomatalRatio * 1.075 * Gbh

# removing below eqn and following Leuning exactly
# Gbhr <- Gbh + 2 * Gradiation
# Total conductance to water vapour
gw <- gs * Gbw/(gs + Gbw)

# Net radiation
# Below calculation not used? 
# Rlongup <- Emissivity * Boltz * Tk(Tleaf)^4
# Rnet <- Rsol - Rlongup

# Absorbed SW radiation. Assume incident NIR = incident PAR, but use different absorptances
Rsol <- (LeafAbsPAR + LeafAbsNIR) * PPFD/UMOLPERJ # J m-2 s-1 - absorbed SW radiation = 2 * APAR
# Sky emissivity - use the Brutsaert formulation
ea <- satur(Tair) - 1000 * VPD
ema <- 0.642 * (ea/Tk(Tair))^(1/7)
# Isothermal net radiation
# Note - longwave needs to be corrected for transmission through the canopy
# See Leuning eqn D1
Rnetiso <- Rsol - (Emissivity - ema) * Boltz * Tk(Tair)^4

# Penman-Monteith equation (Leuning eqn 10)
# GAMMA is the psychrometer constant
# Jones p111: gamma = P cp / (0.622 lambda). 
# Here *AIRMA and LHV is *H20MW so 1/0.622 accounted for as AIRMA/H2OMW = 29/18
GAMMA <- CPAIR * AIRMA * Patm * 1000/LHV
# Transpiration rate (mol H2O m-2 s-1) : in place of gs * VPD
Y <- 1/(1 + Gradiation/Gbh) # Leuning eqn D6
ET <- (1/LHV) * (SLOPE * Y * Rnetiso + 1000 * VPD * Gbh * CPAIR * 
                   AIRMA)/(SLOPE * Y + GAMMA * Gbh/gw) # Leuning eqn 10

# calculate leaf temperature
lambdaET <- LHV * ET
H2 <- Y * (Rnetiso - lambdaET)
# H <- -CPAIR * AIRDENS * (Gbh/CMOLAR) * (Tair - Tleaf)
# Tleaf2 <- Tair + H2/(CPAIR * AIRDENS * (Gbh/CMOLAR))
Tleaf2 <- Tair + H2/(CPAIR * AIRMA * Gbh) # Units I have been using throughout
EnergyBal <- Tleaf - Tleaf2

# return info
l <- data.frame(ELEAFeb = 1000 * ET, Gradiation = Gradiation, 
                Rsol = Rsol, Rnetiso = Rnetiso,
                lambdaET = lambdaET, gw = gw, Gbh = Gbh, 
                H2 = H2, Tleaf2 = Tleaf2)
return(l)
}


