###Modelos de Proje??o Integarl (IPM)###

#Clean memory
rm(list=ls(all=TRUE))

##Pacotes necess?rios para an?lise###

install.packages(IPMpack)
install.packages(gamm4)
install.packages("popbio")
install.packages(fields)
library(IPMpack)
library (gamm4)
library(popbio)
require(fields)

###Carregando tabela de dados#######
###A esp?cie Mora paraensis (Ducke) Ducke pertence a subfam?lia Caesalpinioideae DC.da fam?lia Fabaceae Lindl. Conhecida popularmente como pracu?ba, ? nativa da Amaz?nia e considerada end?mica da regi?o estuarina (Wittmann et al. 2013). Trata-se de uma ?rvore de grande porte, que possui tronco cil?ndrico volumoso, sustentado por enormes sapopemas, copa extensa, ritidoma esverdeado, rugoso, com presen?a de lenticelas e casca interna amarelada com aproximadamente 2 cm de espessura (Loureiro et al. 1979). ####
#d<-read.csv("C:\\R\\Mora paraensis.csv",h=T)
d <- IPM_SH_0506
#View(d)
head(d)

###Calculando a fecundidade (F) e crescimento e sobrevivencia (P), que entra na formula do kernel)###

###Fecundidade###

fec0=d$fec0
d$fec0[which(d$sizeDiam>0)]=0
d$fec0[which(fec0==1)]=1

###Crescimento - growth###

d$alive=d$aliveNext=NA
d$alive[which(d$sizelogDiam>0)]=1
d$aliveNext[which(d$sizeNextlogDiam>0)]=1
d$aliveNext[which(is.na(d$sizeNextlogDiam))]=0

###Sobreviv?ncia###

d$surv=NA
d$surv[which(d$alive==1 & d$aliveNext==1)]=1
d$surv[which(d$alive==1 & d$aliveNext==0)]=0


###Distribui??o de mortalidade e sobreviv?ncia por tamanho###
###Pode se usar apenas o Diam sem logaritimizar. LEMBRAR que para calcular o lambda ? melhor logaritimizar###

par(mfrow=c(1,1))

plot(d$sizeDiam,jitter(d$surv,.2),xlab="Di?metro (cm)", ylab="Sobreviv?ncia",xlim=c(0,150))
plot(d$sizelogDiam,jitter(d$surv,.1))

###Crescimento####
plot(d$sizeDiam,d$sizeNextDiam, xlab="Di?metro t (mm)", ylab="Di?metro t+1 (mm)",
xlim=c(0,150))
abline(0,1,lty=2)

plot(d$sizelogDiam,d$sizeNextlogDiam)

###Individuos - em DIAM - que reproduziram (0 ? sem reprodu??o e 1 reprodutivo)###
plot(d$sizeDiam,d$fec0, xlab="Di?metro (mm)", ylab="Reprodu??o",
     xlim=c(0,150))
plot(d$sizelogDiam,d$fec0)

###Fecundidade por diametro do reprodutor, ou seja, quantas pl?ntulas cada individuo adulto produz por ano###
plot(d$sizeDiam,jitter(d$fec1,1),ylim=c(0,25), 
     xlab="Di?metro (mm)", ylab="Fecundidade", xlim=c(0,150))
plot(d$sizelogDiam,jitter(d$fec1,1),ylim=c(0,25))

###Determining min size, max size and mesh size (tamanho do kernel, que varia do min ao max) of juvs and adults###
hist(d$sizeDiam)

##################################################################
## histograma sem log###################
myhist1 <- hist(d$sizeNextDiam, breaks=20)
multiplier1 <- myhist1$counts / myhist1$density
mydensity1 <- density(d$sizeDiam, na.rm = T)
mydensity1$y <- mydensity1$y * multiplier1[1]
plot(myhist1, axes = F, col = "gray", xlab = "Diameter at ground level (mm)", ylab = "Number of individuals", breaks = 2, main = "", xlim = c(0, 160), ylim = c(0, 50))
axis(1, at=seq(0, 160, by=20))
axis(2)
lines(mydensity1, col = "red", lwd = 2)
#hist(d$sizeDiam, prob = T)
#############################################################################
## histograma com log
#hist(d$sizelogDiam, # histogram
     #col="gray", # column color
     #border="black",
     #prob = T, # show densities instead of frequencies
     #xlab = "Log-Diameter at ground level (mm)",
     #ylab= "Variation in the number of individuals",
     #breaks = 10, main = "", xlim = c(0, 2.5), ylim = c(0, 3))
#hist(d$sizeDiam, prob = T, col = "gray", xlab = "Diameter at ground level (mm)", ylab = "Variation in the number of individuals", breaks = 10, main = "", xlim = c(0, 160), ylim = c(0, 0.05))
#lines(density(d$sizeDiam, na.rm = T), # density plot
      #lwd = 1, # thickness of line
      #col = "red" )
#lines(frequency(d$sizeDiam, na.rm = T), # density plot
     # lwd = 1, # thickness of line
     # col = "red" )


##############################################################################
range(d$sizelogDiam,na.rm=T)
range(d$sizeDiam,na.rm=T)
minDiam=min(d$sizelogDiam,na.rm=T)
maxDiam=max(d$sizelogDiam,na.rm=T)

#IPM for all based on basal diameter
d$size=d$sizelogDiam
d$sizeNext=d$sizeNextlogDiam
minSize=min(d$size,na.rm=T)
maxSize=max(d$size,na.rm=T)
meshSize=110

###Escolher o modelo com menor valor de AIC para sobrevivencia, ou o mais simples com o valor semelhante ao do menor###
soComp=survModelComp(d,expVars = c(surv~1, surv~size, surv~size+size2, surv~size+size2+size3), testType="AIC", makePlot=T, legendPos="bottomright")

###Jogar o modelo (sobrevivencia segundo o melhor modelo)
#so=makeSurvObj(d,surv~1)
so=makeSurvObj(d,surv~size)
#so=makeSurvObj(d,surv~size+size2)
#so=makeSurvObj(d,surv~size+size2+size3)

surv.reg=lm(surv~1, data=d)
summary(surv.reg)

###Escolher o modelo com menor valor de AIC o crescimento, ou o mais ismples com o valor semelhante ao do menor###
goComp=growthModelComp(d,expVars = c(sizeNext~1, sizeNext~size, sizeNext~size+size2, sizeNext~size+size2+size3), testType="AIC", makePlot=T, legendPos="bottomright")

#go=makeGrowthObj(d,sizeNext~1) #,regType="changingVar"
go=makeGrowthObj(d,sizeNext~size)
#go=makeGrowthObj(d,sizeNext~size+size2)
#go=makeGrowthObj(d,sizeNext~size+size2+size3)

growth.reg=lm(sizeNext~size,data=d)
summary(growth.reg)


###Inserindo a fecundidade no kernel###
fo=makeFecObj(d, Formula=c(fec0~size, fec1~size),
              Family=c("binomial", "poisson"),
              Transform=c("none", "none"),
              meanOffspringSize=mean(d[is.na(d$size)==TRUE & is.na(d$sizeNext)==FALSE,"sizeNext"]),
              sdOffspringSize=sd(d[is.na(d$size)==TRUE & is.na(d$sizeNext)==FALSE,"sizeNext"]))

###Juntar sobrevivencia e crescimento####
Pmatrix=makeIPMPmatrix(survObj=so,growObj=go,minSize=minSize,maxSize=maxSize,nBigMatrix=meshSize,correction="discretizeExtremes")


###Diagnoze do modelo de sobrevivencia e crescimento####
###Vai gerar graficos, observar linhas, se estiverem todas "umas em cima das outras" ta ok, se n?o mexer no mesh###

diagnosticsPmatrix(Pmatrix, growObj = go, survObj = so, correction="constant")


#Fazendo a fecundidade para entrar no Kernel

Fmatrix=makeIPMFmatrix(fecObj=fo,minSize=minSize,maxSize=maxSize,nBigMatrix=meshSize,correction="discretizeExtremes")


###Calculo do Kernel (com todo mundo incluido, sobrevivencia, crescimento e fecundidade)
#IPM_bigfrag0506=Pmatrix+Fmatrix
IPM_bigfrag=Pmatrix+Fmatrix



###Diagnoze do kernel###
diagnosticsPmatrix(Pmatrix, growObj = go, survObj = so, correction="constant")


###Calculo do IPM, que ? o equivalente da matriz de transi??o###
par(mfrow=c(1,1))

image.plot(Pmatrix@meshpoints, 
           Pmatrix@meshpoints, 
           t(IPM_bigfrag), 
           zlim=c(0,.1),
           main = "",
           xlab = "Size (mm) at time t", 
           ylab = "Size (mm) at time t + 1")

abline(0,1,lwd=1, lty=2)

###Calculo do lambda e lan?ar no grafico###
lambda1=round(lambda(IPM_bigfrag),2)
mtext(expression(paste(lambda,"       ")),3,-2,col="white")
mtext(paste("       ","=",lambda1),3,-2,col="white")


###Calculo da elasticidade, ou seja, impacto de cada individuo (nesse caso diametro), no lambda
elasbigfrag=elasticity(IPM_bigfrag)
par(mfrow=c(1,1))
image.plot(Pmatrix@meshpoints, 
           Pmatrix@meshpoints, 
           t(elasbigfrag), 
           zlim=c(0,0.0005),
           main = "",
           xlab = "Size (mm) at time t", 
           ylab = "Size (mm) at time t + 1")
abline(0,1,lwd=1, lty=2)
#locator()    ### adicionei o ocator para sabermos as coordenadas dos pontos

#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#
#
# Intervalo de Confiança - Deterministico (utiliza apenas um intervalo de tempo)
# Atenção ao mesh size e ao modelo ideal!!

corte<-function(X){
  dezp<-length(X)/10
  dezp<-round(dezp)
  return(length(X)-dezp)
}

#sample(IPM_SH_0506$id,corte(IPM_SH_0506$id))
#IPM_SH_0506[sample(1:dim(IPM_SH_0506)[1],corte(IPM_SH_0506$id)),]

d<-read.csv("C:\\R\\Mora paraensis.csv",h=T)
d2<-d <- IPM_AJ_1819


Multi.lambda<-NULL

for (i in 1:500){
  d<-d2
  d<-d[sample(1:dim(d)[1],corte(d$id)),]
  
  ###Calculando a fecundidade (F) e crescimento e sobrevivencia (P), que entra na formula do kernel)###
  
  ###Fecundidade###
  
  fec0=d$fec0
  d$fec0[which(d$sizeDiam>0)]=0
  d$fec0[which(fec0==1)]=1
  
  ###Crescimento - growth###
  
  d$alive=d$aliveNext=NA
  d$alive[which(d$sizelogDiam>0)]=1
  d$aliveNext[which(d$sizeNextlogDiam>0)]=1
  d$aliveNext[which(is.na(d$sizeNextlogDiam))]=0
  
  ###Sobreviv?ncia###
  
  d$surv=NA
  d$surv[which(d$alive==1 & d$aliveNext==1)]=1
  d$surv[which(d$alive==1 & d$aliveNext==0)]=0
  
  
  ###Distribui??o de mortalidade e sobreviv?ncia por tamanho###
  ###Pode se usar apenas o Diam sem logaritimizar. LEMBRAR que para calcular o lambda ? melhor logaritimizar###
  
  par(mfrow=c(1,1))
  
  plot(d$sizeDiam,jitter(d$surv,.2),xlab="Di?metro (cm)", ylab="Sobreviv?ncia",xlim=c(0,150))
  plot(d$sizelogDiam,jitter(d$surv,.1))
  
  ###Crescimento####
  plot(d$sizeDiam,d$sizeNextDiam, xlab="Di?metro t (mm)", ylab="Di?metro t+1 (mm)",
       xlim=c(0,150))
  abline(0,1,lty=2)
  
  plot(d$sizelogDiam,d$sizeNextlogDiam)
  
  ###Individuos - em DIAM - que reproduziram (0 ? sem reprodu??o e 1 reprodutivo)###
  plot(d$sizeDiam,d$fec0, xlab="Di?metro (mm)", ylab="Reprodu??o",
       xlim=c(0,150))
  plot(d$sizelogDiam,d$fec0)
  
  ###Fecundidade por diametro do reprodutor, ou seja, quantas pl?ntulas cada individuo adulto produz por ano###
  plot(d$sizeDiam,jitter(d$fec1,1),ylim=c(0,25), 
       xlab="Di?metro (mm)", ylab="Fecundidade", xlim=c(0,150))
  plot(d$sizelogDiam,jitter(d$fec1,1),ylim=c(0,25))
  
  ###Determining min size, max size and mesh size (tamanho do kernel, que varia do min ao max) of juvs and adults###
  hist(d$sizeDiam)
  range(d$sizelogDiam,na.rm=T)
  range(d$sizeDiam,na.rm=T)
  minDiam=min(d$sizelogDiam,na.rm=T)
  maxDiam=max(d$sizelogDiam,na.rm=T)
  
  #IPM for all based on basal diameter
  d$size=d$sizelogDiam
  d$sizeNext=d$sizeNextlogDiam
  minSize=min(d$size,na.rm=T)
  maxSize=max(d$size,na.rm=T)
  meshSize=160
  
  ###Escolher o modelo com menor valor de AIC para sobrevivencia, ou o mais simples com o valor semelhante ao do menor###
  soComp=survModelComp(d,expVars = c(surv~1, surv~size, surv~size+size2, surv~size+size2+size3), testType="AIC", makePlot=T, legendPos="bottomright")
  
  ###Jogar o modelo (sobrevivencia segundo o melhor modelo)
  #so=makeSurvObj(d,surv~1)
  so=makeSurvObj(d,surv~size)
  #so=makeSurvObj(d,surv~size+size2)
  #so=makeSurvObj(d,surv~size+size2+size3)
  
  ###Escolher o modelo com menor valor de AIC o crescimento, ou o mais ismples com o valor semelhante ao do menor###
  goComp=growthModelComp(d,expVars = c(sizeNext~1, sizeNext~size, sizeNext~size+size2, sizeNext~size+size2+size3), testType="AIC", makePlot=T, legendPos="bottomright")
  
  #go=makeGrowthObj(d,sizeNext~1) #,regType="changingVar"
  go=makeGrowthObj(d,sizeNext~size)
  #go=makeGrowthObj(d,sizeNext~size+size2)
  #go=makeGrowthObj(d,sizeNext~size+size2+size3)
  
  ###Inserindo a fecundidade no kernel###
  fo=makeFecObj(d, Formula=c(fec0~size, fec1~size),
                Family=c("binomial", "poisson"),
                Transform=c("none", "none"),
                meanOffspringSize=mean(d[is.na(d$size)==TRUE & is.na(d$sizeNext)==FALSE,"sizeNext"]),
                sdOffspringSize=sd(d[is.na(d$size)==TRUE & is.na(d$sizeNext)==FALSE,"sizeNext"]))
  
  ###Juntar sobrevivencia e crescimento####
  Pmatrix=makeIPMPmatrix(survObj=so,growObj=go,minSize=minSize,maxSize=maxSize,nBigMatrix=meshSize,correction="discretizeExtremes")
  
  
  ###Diagnoze do modelo de sobrevivencia e crescimento####
  ###Vai gerar graficos, observar linhas, se estiverem todas "umas em cima das outras" ta ok, se n?o mexer no mesh###
  #diagnosticsPmatrix(Pmatrix, growObj = go, survObj = so, correction="constant")
  
  
  #Fazendo a fecundidade para entrar no Kernel
  
  Fmatrix=makeIPMFmatrix(fecObj=fo,minSize=minSize,maxSize=maxSize,nBigMatrix=meshSize,correction="discretizeExtremes")
  
  
  ###Calculo do Kernel (com todo mundo incluido, sobrevivencia, crescimento e fecundidade)
  IPM_bigfrag=Pmatrix+Fmatrix
  
  Multi.lambda[i]<-round(lambda(IPM_bigfrag),2)
  
}

hist(Multi.lambda)

Conf.interval<-function(X){
  dse <- 1.96 * sqrt(var(X)/length(X))
  CI <- c(mean(X)- dse, mean(X) + dse)
  return(list(lambda=mean(X), confidence=CI))
}

Conf.interval(Multi.lambda)



########################################### Intervalo de confian;a, inspirado no pacote popgrowthrate do popbio
######## Intervalo de COnfiança - Estocástico (utiliza diversos intervalos de tempo)

# utiliza mais de um intervalo de ano, ai gera uma media desses intervalos
# ir acrscentando os IPMs de diferentes intervalos
## primeiro rodar até o IPM_bigfrag e acrscentar o que falta
args(stochGrowthRateSampleList)

exp(stochGrowthRateSampleList(list(IPM_bigfrag0506,IPM_bigfrag0607),nRunIn = 100, tMax = 5000))

Replicas<-replicate(100,exp(stochGrowthRateSampleList(list(IPM_bigfrag0506,IPM_bigfrag0607),nRunIn = 100, tMax = 5000)))

Conf.interval<-function(X){
  dse <- 1.96 * sqrt(var(X)/length(X))
  CI <- c(mean(X)- dse, mean(X) + dse)
  return(list(lambda_S=mean(X), confidence=CI))
}

Conf.interval(Replicas)


edit(popbio::stoch.growth.rate)
################################################


