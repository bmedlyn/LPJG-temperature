
# Damage toymodel

# Arnold formulation: damage increases linearly over time
temps <- seq(20,60)
dams <- seq(0,100,by=5)

# temperature causing 50% damage after mins minutes
T50 <- function(T50par,z,mins=60) {
   T50par - z * log10(mins)
}

# temperature causing x% damage after mins minutes
Tdam <- function(dam,T50par,z,mins=60) {
  T50par - z * log10(mins*50/dam)
}

### THIS FUNCTION DREW
# damage caused by sitting at temperature Tleaf for mins minutes
damT <- function(Tleaf,T50par,z,mins=60) {
  dam <- 50*mins*10^(-1/z*(T50par - Tleaf))
  dam <- ifelse(dam>100,100,dam)
  return(dam)
}

# visualise relationships
T50par <- 55
z <- 4
mins <- 60

# temperature giving dams damage after an hour
plot(dams,Tdam(dams,T50par = T50par, z = z, mins=300))
# temperature giving dams damage after 30 mins
points(dams,Tdam(dams,T50par = T50par, z = z, mins=30),col="red")

# damage after being at temperature Tleaf for an hour
plot(temps,damT(Tleaf=temps,T50par = T50par, z = z, mins=300),ylim=c(0,100))
# damage after being at temperature Tleaf for 30 mins
points(temps,damT(Tleaf=temps,T50par = T50par, z = z, mins=30),col="red")



# Camille suggested an alternative function for damage 
# Based on shape of Tcrit relationship
Tcrit <- T50par - 5
damage <- function(Tleaf,Tcrit,T50) {
  r <- 2/(T50 - Tcrit)
  100 / (1 + exp(-r*(Tleaf - T50)))
}  

# But the problem with this is that it is not additive
# Twice the time spent at the same temperature does not cause twice the damage
# 1 minute
plot(temps,damage(temps,Tcrit,T50par),type="l",ylim=c(0,100))
# 1 hour
points(temps,damage(temps,T50(Tcrit,z,mins=60),T50(T50par,z,mins=60)),
       type="l",col="red")
# 2 hours
points(temps,damage(temps,T50(Tcrit,z,mins=120),T50(T50par,z,mins=120)),
       type="l",col="orange")
# 3 hours
points(temps,damage(temps,T50(Tcrit,z,mins=180),T50(T50par,z,mins=180)),
       type="l",col="blue")

# compare with damage from Arnold formulation
# When z = 7 they are similar up to 50% damage but the Arnold formulation goes to 100% damage in two steps basically
points(temps,damT(Tleaf=temps,T50par = T50par, z = z, mins=60),col="red")
points(temps,damT(Tleaf=temps,T50par = T50par, z = z, mins=120),col="orange")

# When z is lower the Arnold formulation gives a much stronger threshold
z<-2
plot(temps,damage(temps,T50(Tcrit,z,mins=60),T50(T50par,z,mins=60)),
       type="l",col="red")
points(temps,damT(Tleaf=temps,T50par = T50par, z = z, mins=60),col="red")

# But I think we have to go with the Arnold formulation 
# because of the need for damage to be additive
