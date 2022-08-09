library(ggplot2)
library(data.table)
library(dplyr)
library(tidyr)

base <- read.csv("C:/Users/Shixt/OneDrive/Desktop/INTERNFILES2/rnaseq_donor10021/10021Base.csv")
base <- base[,colSums(is.na(base))<nrow(base)]
base <- na.omit(base)

#ANG 2 10
angComp <- data.frame(base)
angComp <- angComp %>% relocate(contains("AnG"), .after = X)
tests1 <- c()
logDiff1 <- rnorm(20748)
for(i in 1:nrow(angComp)){
  if(i%%1000 == 0){
    cat(i, "\n")
  }
  tests1 <- c(tests1, t.test(angComp[i,2:10],angComp[i,11:ncol(angComp)])$p.value)
}
for(i in 1:nrow(angComp)){
  if(i%%100 == 0){
    cat(i, "\n")
  }
  if(rowMeans(angComp[i,2:10])==0){
    logDiff1[i] <- log(1e-20/rowMeans(angComp[11:ncol(angComp)]))
  }else{
    logDiff1[i] <- log(rowMeans(angComp[i,2:10])/rowMeans(angComp[11:ncol(angComp)]))
  }
}
data1 <- data.frame(base$X, logDiff1, tests1)
write.csv(data1, "C:/Users/Shixt/OneDrive/Desktop/INTERNFILES2/rnaseq_donor10021/angdata")

#CgG 2 5
cggComp <- data.frame(base)
cggComp <- cggComp %>% relocate(contains("CgG"), .after = X)
tests3 <- c()
logDiff3 <- rnorm(20748)
for(i in 1:nrow(angComp)){
  if(i%%1000 == 0){
    cat(i, "\n")
  }
  tests3 <- c(tests3, t.test(cggComp[i,2:5],cggComp[i,6:ncol(angComp)])$p.value)
}
for(i in 1:nrow(angComp)){
  if(i%%100 == 0){
    cat(i, "\n")
  }
  if(rowMeans(cggComp[i,2:5])==0){
    logDiff3[i] <- log(1e-20/rowMeans(angComp[6:ncol(angComp)]))
  }else{
    logDiff3[i] <- log(rowMeans(cggComp[i,2:5])/rowMeans(cggComp[6:ncol(angComp)]))
  }
}
data3 <- data.frame(base$X, logDiff3, tests3)
write.csv(data3, "C:/Users/Shixt/OneDrive/Desktop/INTERNFILES2/cggData")

#Cun 2 5
cunComp <- data.frame(base)
cunComp <- cunComp %>% relocate(contains("Cun"), .after = X)
tests4 <- c()
logDiff4 <- rnorm(20748)
for(i in 1:nrow(angComp)){
  if(i%%1000 == 0){
    cat(i, "\n")
  }
  tests4 <- c(tests4, t.test(cunComp[i,2:5],cunComp[i,6:ncol(angComp)])$p.value)
}
for(i in 1:nrow(angComp)){
  if(i%%100 == 0){
    cat(i, "\n")
  }
  if(rowMeans(cunComp[i,2:5])==0){
    logDiff4[i] <- log(1e-20/rowMeans(angComp[6:ncol(angComp)]))
  }else{
    logDiff4[i] <- log(rowMeans(cunComp[i,2:5])/rowMeans(cunComp[6:ncol(angComp)]))
  }
}
data4 <- data.frame(base$X, logDiff4, tests4)
write.csv(data4, "C:/Users/Shixt/OneDrive/Desktop/INTERNFILES2/cunData")

#FuG 2 5
fugComp <- data.frame(base)
fugComp <- fugComp %>% relocate(contains("FuG"), .after = X)
tests5 <- c()
logDiff5 <- rnorm(20748)
for(i in 1:nrow(angComp)){
  if(i%%1000 == 0){
    cat(i, "\n")
  }
  tests5 <- c(tests5, t.test(fugComp[i,2:5],fugComp[i,6:ncol(angComp)])$p.value)
}
for(i in 1:nrow(angComp)){
  if(i%%100 == 0){
    cat(i, "\n")
  }
  if(rowMeans(fugComp[i,2:5])==0){
    logDiff5[i] <- log(1e-20/rowMeans(fugComp[6:ncol(fugComp)]))
  }else{
    logDiff5[i] <- log(rowMeans(fugComp[i,2:5])/rowMeans(fugComp[6:ncol(angComp)]))
  }
}
data5 <- data.frame(base$X, logDiff5, tests5)
write.csv(data5, "C:/Users/Shixt/OneDrive/Desktop/INTERNFILES2/fugData")

#GP(e+i) 2 5
gpComp <- data.frame(base)
gpComp <- gpComp %>% relocate(contains("GPe") | contains("GPi"), .after = X)
tests6 <- c()
logDiff6 <- rnorm(20748)
for(i in 1:nrow(angComp)){
  if(i%%1000 == 0){
    cat(i, "\n")
  }
  tests6 <- c(tests6, t.test(gpComp[i,2:5],gpComp[i,6:ncol(angComp)])$p.value)
}
for(i in 1:nrow(angComp)){
  if(i%%100 == 0){
    cat(i, "\n")
  }
  if(rowMeans(gpComp[i,2:5])==0){
    logDiff6[i] <- log(1e-20/rowMeans(gpComp[6:ncol(angComp)]))
  }else{
    logDiff6[i] <- log(rowMeans(gpComp[i,2:5])/rowMeans(gpComp[6:ncol(angComp)]))
  }
}
data6 <- data.frame(base$X, logDiff6, tests6)
write.csv(data6, "C:/Users/Shixt/OneDrive/Desktop/INTERNFILES2/gpData")

#GRe 2 5
greComp <- data.frame(base)
greComp <- greComp %>% relocate(contains("GRe"), .after = X)
tests7 <- c()
logDiff7 <- rnorm(20748)
for(i in 1:nrow(angComp)){
  if(i%%1000 == 0){
    cat(i, "\n")
  }
  tests7 <- c(tests7, t.test(greComp[i,2:5],greComp[i,6:ncol(angComp)])$p.value)
}
for(i in 1:nrow(angComp)){
  if(i%%100 == 0){
    cat(i, "\n")
  }
  if(rowMeans(greComp[i,2:5])==0){
    logDiff7[i] <- log(1e-20/rowMeans(greComp[6:ncol(angComp)]))
  }else{
    logDiff7[i] <- log(rowMeans(greComp[i,2:5])/rowMeans(greComp[6:ncol(angComp)]))
  }
}
data7 <- data.frame(base$X, logDiff7, tests7)
write.csv(data7, "C:/Users/Shixt/OneDrive/Desktop/INTERNFILES2/greData")

#HCd 2 5
hcdComp <- data.frame(base)
hcdComp <- hcdComp %>% relocate(contains("HCd"), .after = X)
tests8 <- c()
logDiff8 <- rnorm(20748)
for(i in 1:nrow(angComp)){
  if(i%%1000 == 0){
    cat(i, "\n")
  }
  tests8 <- c(tests8, t.test(hcdComp[i,2:5],hcdComp[i,6:ncol(angComp)])$p.value)
}
for(i in 1:nrow(angComp)){
  if(i%%100 == 0){
    cat(i, "\n")
  }
  if(rowMeans(hcdComp[i,2:5])==0){
    logDiff8[i] <- log(1e-20/rowMeans(hcdComp[6:ncol(angComp)]))
  }else{
    logDiff8[i] <- log(rowMeans(hcdComp[i,2:5])/rowMeans(hcdComp[6:ncol(angComp)]))
  }
}
data8 <- data.frame(base$X, logDiff8, tests8)
write.csv(data8, "C:/Users/Shixt/OneDrive/Desktop/INTERNFILES2/hcdData")

#Hem 2 3
hemComp <- data.frame(base)
hemComp <- hemComp %>% relocate(contains("He"), .after = X)
tests9 <- c()
logDiff9 <- rnorm(20748)
for(i in 1:nrow(angComp)){
  if(i%%1000 == 0){
    cat(i, "\n")
  }
  tests9 <- c(tests9, t.test(hemComp[i,2:3],hemComp[i,4:ncol(angComp)])$p.value)
}
for(i in 1:nrow(angComp)){
  if(i%%100 == 0){
    cat(i, "\n")
  }
  if(rowMeans(greComp[i,2:3])==0){
    logDiff9[i] <- log(1e-20/rowMeans(hemComp[4:ncol(angComp)]))
  }else{
    logDiff9[i] <- log(rowMeans(hemComp[i,2:3])/rowMeans(hemComp[4:ncol(angComp)]))
  }
}
data9 <- data.frame(base$X, logDiff9, tests9)
write.csv(data9, "C:/Users/Shixt/OneDrive/Desktop/INTERNFILES2/hemData")

#ITG 2 5
itgComp <- data.frame(base)
itgComp <- itgComp %>% relocate(contains("ITG"), .after = X)
tests10 <- c()
logDiff10 <- rnorm(20748)
for(i in 1:nrow(angComp)){
  if(i%%1000 == 0){
    cat(i, "\n")
  }
  tests10 <- c(tests10, t.test(itgComp[i,2:5],itgComp[i,6:ncol(angComp)])$p.value)
}
for(i in 1:nrow(angComp)){
  if(i%%100 == 0){
    cat(i, "\n")
  }
  if(rowMeans(itgComp[i,2:5])==0){
    logDiff10[i] <- log(1e-20/rowMeans(itgComp[6:ncol(angComp)]))
  }else{
    logDiff10[i] <- log(rowMeans(itgComp[i,2:5])/rowMeans(itgComp[6:ncol(angComp)]))
  }
}
data10 <- data.frame(base$X, logDiff10, tests10)
write.csv(data10, "C:/Users/Shixt/OneDrive/Desktop/INTERNFILES2/itgData")

#LIG 2 3
LIGComp <- data.frame(base)
LIGComp <- LIGComp %>% relocate(contains("LIG"), .after = X)
tests11 <- c()
logDiff11 <- rnorm(20748)
for(i in 1:nrow(angComp)){
  if(i%%1000 == 0){
    cat(i, "\n")
  }
  tests11 <- c(tests11, t.test(ligComp[i,2:3],ligComp[i,4:ncol(angComp)])$p.value)
}
for(i in 1:nrow(angComp)){
  if(i%%100 == 0){
    cat(i, "\n")
  }
  if(rowMeans(ligComp[i,2:3])==0){
    logDiff11[i] <- log(1e-20/rowMeans(ligComp[4:ncol(angComp)]))
  }else{
    logDiff11[i] <- log(rowMeans(ligComp[i,2:3])/rowMeans(ligComp[4:ncol(angComp)]))
  }
}
data11 <- data.frame(base$X, logDiff11, tests11)
write.csv(data11, "C:/Users/Shixt/OneDrive/Desktop/INTERNFILES2/LIGData")

#LiG 2 5
ligComp <- data.frame(base)
ligComp <- ligComp %>% relocate(contains("LiG"), .after = X)
tests12 <- c()
logDiff12 <- rnorm(20748)
for(i in 1:nrow(angComp)){
  if(i%%1000 == 0){
    cat(i, "\n")
  }
  tests11 <- c(tests12, t.test(ligComp[i,2:5],ligComp[i,6:ncol(angComp)])$p.value)
}
for(i in 1:nrow(angComp)){
  if(i%%100 == 0){
    cat(i, "\n")
  }
  if(rowMeans(ligComp[i,2:5])==0){
    logDiff12[i] <- log(1e-20/rowMeans(ligComp[6:ncol(angComp)]))
  }else{
    logDiff12[i] <- log(rowMeans(ligComp[i,2:5])/rowMeans(ligComp[6:ncol(angComp)]))
  }
}
data11 <- data.frame(base$X, logDiff12, tests12)
write.csv(data11, "C:/Users/Shixt/OneDrive/Desktop/INTERNFILES2/ligData")

#MFG 2 5
mfgComp <- data.frame(base)
mfgComp <- mfgComp %>% relocate(contains("MFG"), .after = X)
tests14 <- c()
logDiff14 <- rnorm(20748)
for(i in 1:nrow(angComp)){
  if(i%%1000 == 0){
    cat(i, "\n")
  }
  tests14 <- c(tests14, t.test(mfgComp[i,2:5],mfgComp[i,6:ncol(angComp)])$p.value)
}
for(i in 1:nrow(angComp)){
  if(i%%100 == 0){
    cat(i, "\n")
  }
  if(rowMeans(mfgComp[i,2:5])==0){
    logDiff14[i] <- log(1e-20/rowMeans(mfgComp[6:ncol(angComp)]))
  }else{
    logDiff14[i] <- log(rowMeans(mfgComp[i,2:5])/rowMeans(mfgComp[6:ncol(angComp)]))
  }
}
data14 <- data.frame(base$X, logDiff14, tests14)
write.csv(data14, "C:/Users/Shixt/OneDrive/Desktop/INTERNFILES2/mfgData")

#MOrG 2 4
mogComp <- data.frame(base)
mogComp <- mogComp %>% relocate(contains("MOrG"), .after = X)
tests15 <- c()
logDiff15 <- rnorm(20748)
for(i in 1:nrow(angComp)){
  if(i%%1000 == 0){
    cat(i, "\n")
  }
  tests15 <- c(tests15, t.test(mogComp[i,2:4],mogComp[i,5:ncol(angComp)])$p.value)
}
for(i in 1:nrow(angComp)){
  if(i%%100 == 0){
    cat(i, "\n")
  }
  if(rowMeans(mogComp[i,2:4])==0){
    logDiff15[i] <- log(1e-20/rowMeans(mogComp[5:ncol(angComp)]))
  }else{
    logDiff15[i] <- log(rowMeans(mogComp[i,2:4])/rowMeans(mogComp[5:ncol(angComp)]))
  }
}
data15 <- data.frame(base$X, logDiff15, tests15)
write.csv(data15, "C:/Users/Shixt/OneDrive/Desktop/INTERNFILES2/mogData")

#MTG 2 5
mtgComp <- data.frame(base)
mtgComp <- mtgComp %>% relocate(contains("MTG"), .after = X)
tests16 <- c()
logDiff16 <- rnorm(20748)
for(i in 1:nrow(angComp)){
  if(i%%1000 == 0){
    cat(i, "\n")
  }
  tests16 <- c(tests16, t.test(mtgComp[i,2:5],mtgComp[i,6:ncol(angComp)])$p.value)
}
for(i in 1:nrow(angComp)){
  if(i%%100 == 0){
    cat(i, "\n")
  }
  if(rowMeans(mtgComp[i,2:5])==0){
    logDiff16[i] <- log(1e-20/rowMeans(mtgComp[6:ncol(angComp)]))
  }else{
    logDiff16[i] <- log(rowMeans(mtgComp[i,2:5])/rowMeans(mtgComp[6:ncol(angComp)]))
  }
}
data16 <- data.frame(base$X, logDiff16, tests16)
write.csv(data16, "C:/Users/Shixt/OneDrive/Desktop/INTERNFILES2/mtgData")

#orIFG 2 5
ifgComp <- data.frame(base)
ifgComp <- ifgComp %>% relocate(contains("orIFG"), .after = X)
tests17 <- c()
logDiff17 <- rnorm(20748)
for(i in 1:nrow(angComp)){
  if(i%%1000 == 0){
    cat(i, "\n")
  }
  tests17 <- c(tests17, t.test(ifgComp[i,2:5],ifgComp[i,6:ncol(angComp)])$p.value)
}
for(i in 1:nrow(angComp)){
  if(i%%100 == 0){
    cat(i, "\n")
  }
  if(rowMeans(ifgComp[i,2:5])==0){
    logDiff17[i] <- log(1e-20/rowMeans(ifgComp[6:ncol(angComp)]))
  }else{
    logDiff17[i] <- log(rowMeans(ifgComp[i,2:5])/rowMeans(ifgComp[6:ncol(angComp)]))
  }
}
data17 <- data.frame(base$X, logDiff17, tests17)
write.csv(data17, "C:/Users/Shixt/OneDrive/Desktop/INTERNFILES2/ifgData")

#PCLa 2 10 
pclComp <- data.frame(base)
pclComp <- pclComp %>% relocate(contains("PCL"), .after = X)
tests18 <- c()
logDiff18 <- rnorm(20748)
for(i in 1:nrow(angComp)){
  if(i%%1000 == 0){
    cat(i, "\n")
  }
  tests18 <- c(tests18, t.test(pclComp[i,2:10],pclComp[i,11:ncol(angComp)])$p.value)
}
for(i in 1:nrow(angComp)){
  if(i%%100 == 0){
    cat(i, "\n")
  }
  if(rowMeans(pclComp[i,2:10])==0){
    logDiff18[i] <- log(1e-20/rowMeans(pclComp[11:ncol(angComp)]))
  }else{
    logDiff18[i] <- log(rowMeans(pclComp[i,2:10])/rowMeans(pclComp[11:ncol(angComp)]))
  }
}
data18 <- data.frame(base$X, logDiff18, tests18)
write.csv(data18, "C:/Users/Shixt/OneDrive/Desktop/INTERNFILES2/pclData")

#Pcu 2 5 
pcuComp <- data.frame(base)
pcuComp <- pcuComp %>% relocate(contains("Pcu"), .after = X)
tests19 <- c()
logDiff19 <- rnorm(20748)
for(i in 1:nrow(angComp)){
  if(i%%1000 == 0){
    cat(i, "\n")
  }
  tests19 <- c(tests19, t.test(pcuComp[i,2:5],pcuComp[i,6:ncol(angComp)])$p.value)
}
for(i in 1:nrow(angComp)){
  if(i%%100 == 0){
    cat(i, "\n")
  }
  if(rowMeans(pcuComp[i,2:5])==0){
    logDiff19[i] <- log(1e-20/rowMeans(pcuComp[6:ncol(angComp)]))
  }else{
    logDiff19[i] <- log(rowMeans(pcuComp[i,2:5])/rowMeans(pcuComp[6:ncol(angComp)]))
  }
}
data19 <- data.frame(base$X, logDiff19, tests19)
write.csv(data19, "C:/Users/Shixt/OneDrive/Desktop/INTERNFILES2/pcuData")

#PHG 2 5
phgComp <- data.frame(base)
phgComp <- phgComp %>% relocate(contains("PHG"), .after = X)
tests20 <- c()
logDiff20 <- rnorm(20748)
for(i in 1:nrow(angComp)){
  if(i%%1000 == 0){
    cat(i, "\n")
  }
  tests20 <- c(tests20, t.test(phgComp[i,2:5],phgComp[i,6:ncol(angComp)])$p.value)
}
for(i in 1:nrow(angComp)){
  if(i%%100 == 0){
    cat(i, "\n")
  }
  if(rowMeans(phgComp[i,2:5])==0){
    logDiff20[i] <- log(1e-20/rowMeans(phgComp[6:ncol(angComp)]))
  }else{
    logDiff20[i] <- log(rowMeans(phgComp[i,2:5])/rowMeans(phgComp[6:ncol(angComp)]))
  }
}
data20 <- data.frame(base$X, logDiff20, tests20)
write.csv(data20, "C:/Users/Shixt/OneDrive/Desktop/INTERNFILES2/phgData")

#PoG 2 9
pogComp <- data.frame(base)
pogComp <- pogComp %>% relocate(contains("PoG"), .after = X)
tests21 <- c()
logDiff21 <- rnorm(20748)
for(i in 1:nrow(angComp)){
  if(i%%1000 == 0){
    cat(i, "\n")
  }
  tests21 <- c(tests21, t.test(pogComp[i,2:9],pogComp[i,10:ncol(angComp)])$p.value)
}
for(i in 1:nrow(angComp)){
  if(i%%100 == 0){
    cat(i, "\n")
  }
  if(rowMeans(pogComp[i,2:9])==0){
    logDiff21[i] <- log(1e-20/rowMeans(pogComp[10:ncol(angComp)]))
  }else{
    logDiff21[i] <- log(rowMeans(pogComp[i,2:5])/rowMeans(pogComp[10:ncol(angComp)]))
  }
}
data21 <- data.frame(base$X, logDiff21, tests21)
write.csv(data21, "C:/Users/Shixt/OneDrive/Desktop/INTERNFILES2/pogData")

#PrG 2 5
prgComp <- data.frame(base)
prgComp <- prgComp %>% relocate(contains("PrG"), .after = X)
tests22 <- c()
logDiff22 <- rnorm(20748)
for(i in 1:nrow(angComp)){
  if(i%%1000 == 0){
    cat(i, "\n")
  }
  tests22 <- c(tests22, t.test(prgComp[i,2:5],prgComp[i,6:ncol(angComp)])$p.value)
}
for(i in 1:nrow(angComp)){
  if(i%%100 == 0){
    cat(i, "\n")
  }
  if(rowMeans(prgComp[i,2:5])==0){
    logDiff22[i] <- log(1e-20/rowMeans(prgComp[6:ncol(angComp)]))
  }else{
    logDiff22[i] <- log(rowMeans(prgComp[i,2:5])/rowMeans(prgComp[6:ncol(angComp)]))
  }
}
data22 <- data.frame(base$X, logDiff22, tests22)
write.csv(data22, "C:/Users/Shixt/OneDrive/Desktop/INTERNFILES2/prgData")

#Pu 2 5
puComp <- data.frame(base)
puComp <- puComp %>% relocate(contains("Pu"), .after = X)
tests23 <- c()
logDiff23 <- rnorm(20748)
for(i in 1:nrow(angComp)){
  if(i%%1000 == 0){
    cat(i, "\n")
  }
  tests23 <- c(tests23, t.test(puComp[i,2:5],puComp[i,6:ncol(angComp)])$p.value)
}
for(i in 1:nrow(angComp)){
  if(i%%100 == 0){
    cat(i, "\n")
  }
  if(rowMeans(puComp[i,2:5])==0){
    logDiff20[i] <- log(1e-20/rowMeans(puComp[6:ncol(angComp)]))
  }else{
    logDiff20[i] <- log(rowMeans(puComp[i,2:5])/rowMeans(puComp[6:ncol(angComp)]))
  }
}
data23 <- data.frame(base$X, logDiff23, tests23)
write.csv(data23, "C:/Users/Shixt/OneDrive/Desktop/INTERNFILES2/puData")

#PV 2 4
pvComp <- data.frame(base)
pvComp <- pvComp %>% relocate(contains("Pv"), .after = X)
tests235 <- c()
logDiff235 <- rnorm(20748)
for(i in 1:nrow(angComp)){
  if(i%%1000 == 0){
    cat(i, "\n")
  }
  tests235 <- c(tests235, t.test(puComp[i,2:4],puComp[i,5:ncol(angComp)])$p.value)
}
for(i in 1:nrow(angComp)){
  if(i%%100 == 0){
    cat(i, "\n")
  }
  if(rowMeans(puComp[i,2:4])==0){
    logDiff235[i] <- log(1e-20/rowMeans(pvComp[5:ncol(angComp)]))
  }else{
    logDiff235[i] <- log(rowMeans(pvComp[i,2:4])/rowMeans(pvComp[5:ncol(angComp)]))
  }
}
data235 <- data.frame(base$X, logDiff235, tests235)
write.csv(data235, "C:/Users/Shixt/OneDrive/Desktop/INTERNFILES2/puData")

#SFG 2 9
sfgComp <- data.frame(base)
sfgComp <- sfgComp %>% relocate(contains("SFG"), .after = X)
tests25 <- c()
logDiff25 <- rnorm(20748)
for(i in 1:nrow(angComp)){
  if(i%%1000 == 0){
    cat(i, "\n")
  }
  tests25 <- c(tests25, t.test(sfgComp[i,2:9],sfgComp[i,10:ncol(angComp)])$p.value)
}
for(i in 1:nrow(angComp)){
  if(i%%100 == 0){
    cat(i, "\n")
  }
  if(rowMeans(sfgComp[i,2:9])==0){
    logDiff25[i] <- log(1e-20/rowMeans(phgComp[10:ncol(angComp)]))
  }else{
    logDiff25[i] <- log(rowMeans(sfgComp[i,2:9])/rowMeans(sfgComp[10:ncol(angComp)]))
  }
}
data25 <- data.frame(base$X, logDiff25, tests25)
write.csv(data25, "C:/Users/Shixt/OneDrive/Desktop/INTERNFILES2/sfgData")

#SIG 2 3
sigComp <- data.frame(base)
sigComp <- sigComp %>% relocate(contains("SiG"), .after = X)
tests26 <- c()
logDiff26 <- rnorm(20748)
for(i in 1:nrow(angComp)){
  if(i%%1000 == 0){
    cat(i, "\n")
  }
  tests26 <- c(tests20, t.test(sigComp[i,2:3],sigComp[i,4:ncol(angComp)])$p.value)
}
for(i in 1:nrow(angComp)){
  if(i%%100 == 0){
    cat(i, "\n")
  }
  if(rowMeans(sigComp[i,2:3])==0){
    logDiff26[i] <- log(1e-20/rowMeans(sigComp[4:ncol(angComp)]))
  }else{
    logDiff26[i] <- log(rowMeans(sigComp[i,2:3])/rowMeans(sigComp[4:ncol(angComp)]))
  }
}
data26 <- data.frame(base$X, logDiff26, tests26)
write.csv(data26, "C:/Users/Shixt/OneDrive/Desktop/INTERNFILES2/sigData")

#SMG 2 5
smgComp <- data.frame(base)
smgComp <- smgComp %>% relocate(contains("SMG"), .after = X)
tests27 <- c()
logDiff27 <- rnorm(20748)
for(i in 1:nrow(angComp)){
  if(i%%1000 == 0){
    cat(i, "\n")
  }
  tests27 <- c(tests27, t.test(smgComp[i,2:5],smgComp[i,6:ncol(angComp)])$p.value)
}
for(i in 1:nrow(angComp)){
  if(i%%100 == 0){
    cat(i, "\n")
  }
  if(rowMeans(smgComp[i,2:5])==0){
    logDiff27[i] <- log(1e-20/rowMeans(smgComp[6:ncol(angComp)]))
  }else{
    logDiff27[i] <- log(rowMeans(smgComp[i,2:5])/rowMeans(smgComp[6:ncol(angComp)]))
  }
}
data27 <- data.frame(base$X, logDiff27, tests27)
write.csv(data27, "C:/Users/Shixt/OneDrive/Desktop/INTERNFILES2/smgData")

#SPL 2 6
splComp <- data.frame(base)
splComp <- splComp %>% relocate(contains("SPL"), .after = X)
tests28 <- c()
logDiff28 <- rnorm(20748)
for(i in 1:nrow(angComp)){
  if(i%%1000 == 0){
    cat(i, "\n")
  }
  tests28 <- c(tests28, t.test(splComp[i,2:6],splComp[i,7:ncol(angComp)])$p.value)
}
for(i in 1:nrow(angComp)){
  if(i%%100 == 0){
    cat(i, "\n")
  }
  if(rowMeans(splComp[i,2:6])==0){
    logDiff28[i] <- log(1e-20/rowMeans(splComp[7:ncol(angComp)]))
  }else{
    logDiff28[i] <- log(rowMeans(splComp[i,2:6])/rowMeans(splComp[7:ncol(angComp)]))
  }
}
data28 <- data.frame(base$X, logDiff28, tests28)
write.csv(data28, "C:/Users/Shixt/OneDrive/Desktop/INTERNFILES2/splData")

#STG 2 5
stgComp <- data.frame(base)
stgComp <- stgComp %>% relocate(contains("STG"), .after = X)
tests29 <- c()
logDiff29 <- rnorm(20748)
for(i in 1:nrow(angComp)){
  if(i%%1000 == 0){
    cat(i, "\n")
  }
  tests29 <- c(tests29, t.test(stgComp[i,2:5],stgComp[i,6:ncol(angComp)])$p.value)
}
for(i in 1:nrow(angComp)){
  if(i%%100 == 0){
    cat(i, "\n")
  }
  if(rowMeans(stgComp[i,2:5])==0){
    logDiff29[i] <- log(1e-20/rowMeans(stgComp[6:ncol(angComp)]))
  }else{
    logDiff29[i] <- log(rowMeans(stgComp[i,2:5])/rowMeans(stgComp[6:ncol(angComp)]))
  }
}
data29 <- data.frame(base$X, logDiff29, tests29)
write.csv(data29, "C:/Users/Shixt/OneDrive/Desktop/INTERNFILES2/stgData")








