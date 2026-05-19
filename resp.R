
# respiration toymodel

temps <- seq(10,40,by=1)

resp <- function(BC4,Temp,Tgrow) {
  exp(0.067*(Temp-25)) * exp (BC4*(Tgrow - 25))
}

resph <- function(BC4,Temp,Tgrow) {
  exp(0.1012*(Temp-25) - 0.0005*(Temp^2 - 25^2)) * exp (BC4*(Tgrow - 25))
}

plot(temps,resph(BC4 = -0.05,temps, Tgrow = 25),type="l")
points(temps,resp(BC4 = -0.05,temps,Tgrow = 25),type="l",lty=2)
points(temps,resph(BC4 = -0.05,temps, Tgrow = 15),type="l",col="blue")
points(temps,resph(BC4 = -0.05,temps, Tgrow = 35),type="l",col="red")

abline(v=15)
abline(v=25)

resph(BC4 = -0.05,15,Tgrow = 15)
resph(BC4 = -0.05,25,Tgrow = 25)
resph(BC4 = -0.05,35,Tgrow = 35)

