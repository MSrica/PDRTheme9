rm(list = ls())
setwd('.')

library(stringr)

# ------------------------------------------------------------------------------

# constants
fiPole <- 78.3 * pi / 180
lambdaPole <- 291 * pi / 180

ippHeight <- 350
earthRadius <- 6378
lightSpeed <- 2.99792458e+08

# ------------------------------------------------------------------------------

# receiver (user) cooridinates (Pevex Kukuljanovo)
userFiLambda <- data.frame(
  fiUser = 45.331360488178596*(pi/180), # longitude
  lambdaUser = 14.509391791906435*(pi/180) # latitude
)
print(userFiLambda)

# real Cartesian receiver (user) coordinates
userCoords <- data.frame(
  X = earthRadius * cos(userFiLambda$lambdaUser) * cos(userFiLambda$fiUser),
  Y = earthRadius * cos(userFiLambda$lambdaUser) * sin(userFiLambda$fiUser),
  Z = earthRadius * sin(userFiLambda$lambdaUser)
)
print(userCoords)

# ------------------------------------------------------------------------------

rinexData <- read.delim('data/rinex.txt', header = FALSE)[[1]]

alphaIonRow  <- str_trim(gsub('\\s+',' ', rinexData[4]))
alphaIonData <- str_split(alphaIonRow, ' ')[[1]][1:4]
alphaIon <- c(
  as.double(sub('D', 'e', alphaIonData[1])),
  as.double(sub('D', 'e', alphaIonData[2])),
  as.double(sub('D', 'e', alphaIonData[3])),
  as.double(sub('D', 'e', alphaIonData[4]))
)
print(alphaIon)

betaIonRow  <- str_trim(gsub('\\s+',' ', rinexData[5]))
betaIonData <- str_split(betaIonRow, ' ')[[1]][1:4]
betaIon <- c(
  as.double(sub('D', 'e', betaIonData[1])),
  as.double(sub('D', 'e', betaIonData[2])),
  as.double(sub('D', 'e', betaIonData[3])),
  as.double(sub('D', 'e', betaIonData[4]))
)
print(betaIon)

# ------------------------------------------------------------------------------

sp3Data <- read.delim('data/sp3.txt', header=FALSE)

sp3DataFormatted <- list()
rowIndex <- 23

# instead of indexing by time, index by gps id !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
while(TRUE) {
  row = str_trim(gsub('\\s+', ' ', sp3Data[rowIndex, ]))
  rowIndex = rowIndex+1
  
  rowDataStr = str_split(row, ' ')[[1]]
  # maybe useless
  rowData = str_trim(gsub('\\s+', ' ', rowDataStr))
  
  # EOF
  if (length(rowData) == 1) break
  
  # * year month day hour minutes seconds
  if (length(rowData) == 7) { 
    key = row
    sp3DataFormatted[key] = c()
    next
  }
  
  # id x y z clock deviations
  sp3DataFormatted[[key]] = append(sp3DataFormatted[[key]], row)
}

# ------------------------------------------------------------------------------

ionVector <- c()
selectedSatelite <- '30'

exportEmptySpaces <- strrep(' ', 32)

timestampHeader  <- substr(str_interp('TIME [s]${exportEmptySpaces}'), 0, 10)
sateliteIdHeader <- substr(str_interp('SATELITE_ID${exportEmptySpaces}'), 0, 10)
distanceHeader   <- substr(str_interp('DISTANCE [m]${exportEmptySpaces}'), 0, 20)

header <- str_interp('${timestampHeader} ${sateliteIdHeader} ${distanceHeader}')
lines  <- c()
lines  <- append(lines, header)

for(sp3Time in names(sp3DataFormatted)) {
  sp3TimeFormatted <- str_trim(sub('\\*', '', sp3Time))
  sp3TimeValues    <- str_split(sp3TimeFormatted, ' ')[[1]]
  
  timeHours <- as.integer(sp3TimeValues[4])
  timeMinutes <- as.integer(sp3TimeValues[5])
  timeSeconds <- as.integer(sp3TimeValues[6])
  timeGps <- (timeHours * 60 + timeMinutes) * 60 + timeSeconds
  
  sp3TimeFormattedLength <- substr(str_interp('${timeGps}${exportEmptySpaces}'), 0, 10)
  
  for(sp3SateliteValueForTime in sp3DataFormatted[sp3Time][[1]]) {
    values <- str_split(sp3SateliteValueForTime, ' ')[[1]]
    
    sateliteId <- values[1]
    sateliteX <- as.double(values[2])
    sateliteY <- as.double(values[3])
    sateliteZ <- as.double(values[4])
    
    sateliteIdFormatted <- str_trim(sub('PG', '', sateliteId))
    sateliteIdFormattedLength <- substr(str_interp('${sateliteIdFormatted}${exportEmptySpaces}'), 0, 10)
    
    elevation <- sateliteZ - userCoords$Z
    azimuth <- atan((sateliteX - userCoords$X) - (sateliteY - userCoords$Y))
    
    psi <- pi / 2 - elevation - asin(earthRadius / (earthRadius + ippHeight) * cos(elevation))
    fiIon <- asin(sin(userFiLambda$fiUser) * cos(psi) + cos(userFiLambda$fiUser) * sin(psi) * cos(azimuth))
    lambdaIon <- userFiLambda$lambdaUser + psi * sin(azimuth) / cos(fiIon)
    fiMag <- asin(sin(fiIon) * sin(fiPole) + cos(fiIon) * cos(fiPole) * cos(lambdaIon - lambdaPole))
    
    time <- 43200 * lambdaIon / pi + timeGps
    if (time >= 86400) time <- time - 86400
    if (time < 0) time <- time + 86400
    
    azimuthIon <- 0
    for (i in seq(1, 4, 1)) azimuthIon <- azimuthIon + alphaIon[i] * `^`(fiMag / pi, i - 1)
    if (azimuthIon < 0) azimuthIon <- 0
    
    psiIon <- 0
    for (i in seq(1, 4, 1)) psiIon <- psiIon + as.double(betaIon[i]) * `^`(fiMag / pi, i - 1)
    if (psiIon > 72000) psiIon <- 72000
    
    xIon <- 2 * pi * (time - 50400) / psiIon
    Fun <- `^`(1 - `^`(earthRadius / (earthRadius + ippHeight) * cos(elevation), 2), -1 / 2)
    
    if (abs(xIon) < pi / 2) dIon <- (5e-9 + azimuthIon * cos(xIon)) * Fun
    else dIon <- 5e-9 * Fun
    dIon <- dIon * lightSpeed
    
    if (sateliteIdFormatted == selectedSatelite) ionVector <- append(ionVector, dIon)
    
    diffXSquared <- (userCoords$X - sateliteX) ** 2
    diffYSquared <- (userCoords$Y - sateliteY) ** 2
    diffZSquared <- elevation ** 2
    
    distance <- sqrt(diffXSquared + diffYSquared + diffZSquared) * 1000 + dIon
    distanceFormattedLength <- substr(str_interp('${distance}${exportEmptySpaces}'), 0, 20)
    
    line <- str_interp('${sp3TimeFormattedLength} ${sateliteIdFormattedLength} ${distanceFormattedLength}')
    lines <- append(lines, line)
  }
}

# ------------------------------------------------------------------------------

fileConn <- file('output.txt')
writeLines(lines, fileConn)
close(fileConn)

print(ionVector)
time <- seq(0, 24, len = 24 * 4)
plot(time, ionVector, type = 'l', main=str_interp('Klobucharov model za ${selectedSatelite}'), xlab='time[h]', ylab='dIon[m]')
