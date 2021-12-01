ArchTest<-function(x, lags=20, demean=FALSE)
{
  # Capture name of x for documentation in the output  
  xName <- deparse(substitute(x))
  # 
  x <- as.vector(x)
  if(demean) x <- scale(x, center = TRUE, scale = FALSE)
  #  
  lags <- lags + 1
  mat <- embed(x^2, lags)
  arch.lm <- summary(lm(mat[, 1] ~ mat[, -1]))
  STATISTIC <- arch.lm$r.squared * length(resid(arch.lm))
  names(STATISTIC) <- "Chi-squared"
  PARAMETER <- lags - 1
  names(PARAMETER) <- "df"
  PVAL <- 1 - pchisq(STATISTIC, df = PARAMETER)
  METHOD <- "ARCH LM-test;  Null hypothesis:  no ARCH effects"
  result <- list(statistic = STATISTIC, parameter = PARAMETER, 
                 p.value = PVAL, method = METHOD, data.name =
                   xName)
  class(result) <- "htest"
  return(result)
}
library(PerformanceAnalytics)
library(quantmod)
library(fGarch)
library(stats)
library(lubridate)
library(forecast)
library(tseries)
library(normtest)
library(sjPlot)

clave<-c("ARA.MX","AUTLANB.MX","AXTELCPO.MX")

datos <- new.env()
precio <- c()

for (i in clave){
  getSymbols(i,from="2015-01-01", to=today(), warnings ="FALSE", env = datos)
  activo<-datos[[i]]
  chartSeries(activo, theme="white",name=i)
  precio[[i]]<-activo[,6]
}


# Estimar rendimiento y volatilidad por periodos
periodos<-c(2015,2016,2017,2018,2019)
mu<-matrix(0,length(periodos),length(clave))
sigma<-matrix(0,length(periodos),length(clave))
precio_usar<-c()
rtn<-c()
rend<-c()
precios<-c()
lnprecio<-c()
lna<-c()

for (j in 1:length(periodos)){
  for (i in 1:length(precio)){
    precio_usar[[i]]<-precio[[i]][year(precio[[i]])==periodos[j] & month(precio[[i]])>=1 & month(precio[[i]])<=12]
    rtn[[i]]<-na.omit(diff(log(precio_usar[[i]])))
    lna[[i]]<- log(precio_usar[[i]])
    mu[j,i]<-mean(rtn[[i]])
    sigma[j,i]<-as.numeric(sqrt(var(rtn[[i]])))
  }
  lnprecio[[j]]<-do.call(cbind, lna)
  rend[[j]]<-do.call(cbind, rtn)
  precios[[j]]<-do.call(cbind,precio_usar)
}
preciotodos<-do.call(cbind,precio)
rtn<-na.omit(diff(log(preciotodos)))
lna<- log(preciotodos)

for (i in 1:length(clave)){
  chartSeries(rtn[,i], theme="white",name=clave[i])
  chartSeries(lna[,i], theme="white",name=clave[i])
}


lnprecio[[j+1]]<-lna
rend[[j+1]]<-rtn
precios[[j+1]]<-preciotodos

mu<-rbind(mu,colMeans(rtn))
sigma<-rbind(sigma,sqrt(sapply(rtn,var)))

periodos[length(periodos)+1]<-"todo"
rownames(sigma)<-periodos
rownames(mu)<-periodos

muanual = mu*252
sigmaanual = sigma*(sqrt(252))

# Pruebas de Normalidad
test<-function(i,j,x) ifelse((jb.norm.test(x[[i]][,j])[[1]]>9.2134), "no normalidad", "normalidad")
rendm <- matrix(0,6,3)
preciom <- matrix(0,6,3)
lnpreciom <- matrix(0,6,3)
for (i in 1:length(precios)){
  for (j in 1:ncol(precios[[1]])){
    
    preciom[i,j]<-test(i,j,precios)
    rendm[i,j]<-test(i,j,rend)
    lnpreciom[i,j]<-test(i,j,lnprecio)
}
}

colnames(preciom)<-clave
colnames(rendm)<-clave
colnames(lnpreciom)<-clave
rownames(preciom)<-periodos
rownames(rendm)<-periodos
rownames(lnpreciom)<-periodos

#Prueba de rendimiento estadisticamente igual a cero (T student)

TS<-matrix(0,length(periodos),length(clave))
rownames(TS)<-periodos
colnames(TS)<-clave

for (i in 1:length(rend)){
  for (j in 1:length(clave)){
    t<-t.test(as.vector(rend[[i]][,j]))[["statistic"]]
    n<-nrow(rend[[i]][,j])
    tr<-qt(.99,n-1)
  
    if(t>tr){
      b<-"media dif 0"
    }else {b<-"media igual a cero"}
    TS[i,j]<-b
  }
}


#Parte 2 
#montecarlo

preciomont<-c()
rtnmont<-c()
mumont<-c()
sigmamont<-c()

for (i in 1:length(precio)){
    preciomont[[i]]<-precio[[i]][year(precio[[i]])==2019 & month(precio[[i]])>=4 & month(precio[[i]])<=6]
    rtnmont[[i]]<-na.omit(diff(log(preciomont[[i]])))
    mumont[i]<-mean(rtnmont[[i]])
    sigmamont[i]<-as.numeric(sqrt(var(rtnmont[[i]])))
}
preciomont <- do.call(cbind, preciomont)

#Usando siempre el mismo valor inicial
So<-as.numeric(tail(preciomont,1))
ne<-1000
dias<-20
Fechas_est<-seq(date(tail(preciomont,1))+1,as.Date("2019-11-20"),by="days")
Fechas_Selec<-subset(Fechas_est, wday(Fechas_est)>=2 & wday(Fechas_est)<=6)
Ssel<-na.omit(preciotodos[Fechas_Selec])
dias_est<-nrow(Ssel)
##Aleatorios
epsilon<-rnorm(ne)
plot(epsilon, type="l")
St1 <- c()

for (i in 1:length(So)){
  St1[[i]] <- matrix(0,ne,dias_est)
}
dt<-seq(1:dias_est) 

for (x in 1:length(So)){
  for (i in 1:dias_est) {
    for (j in 1:ne) {
      St1[[x]][j,i]<-So[x]*exp((mumont[x]-.5*sigmamont[x]^2)*dt[i]+sigmamont[x]*sqrt(dt[i])*epsilon[j])
    }}}


St2 <- c()
for (i in 1:length(So)){
  St2[[i]] <- matrix(So[i],ne,dias_est)
}
dt<-1
for (x in 1:length(clave)){
  for (i in 2:dias_est) {
    for (j in 1:ne) {
      St2[[x]][j,i]<-Ssel[i,x]*exp((mumont[x]-.5*sigmamont[x]^2)*dt+sigmamont[x]*sqrt(dt)*epsilon[j])
    }}}

PRO<-c()
Promedio<-c()
PR<-c()
Pr<-c()
St<-list(St1,St2)
boxplot(St1[i])
for (j in 1:length(St)){
  for (i in 1:length(St1)){
    P_Est<-xts(t(St[[j]][[i]]), order.by=date(Ssel))
    inicio<-date(P_Est[1,])
    fin<-date(P_Est[dias_est,])
    PR[[i]]<-subset(precio[[i]], date(precio[[i]])>=inicio & date(precio[[i]])<=fin)
    Pr[[i]]<-xts(rowMeans(P_Est), order.by=date(Ssel))
    print(plot(P_Est,main=clave[i]))
    
  }
  PRO[[j]]<-do.call(cbind,PR)
  Promedio[[j]]<-do.call(cbind,Pr)
  colnames(Promedio[[j]])<-colnames(PRO)
  for (i in 1:length(clave)){
    Quantiles<-apply(xts(t(St[[j]][[i]]), order.by=date(Ssel)), 1, FUN=quantile) 
    Quantiles<-xts(t(Quantiles), order.by=date(Ssel))
    plot(Quantiles,main=clave[i])
    lines(PRO[[j]][,i],col="black",lwd=5)
    print(lines(Promedio[[j]][,i],col="blue",lwd=5))
  }
}
#Promedio movil
mav <- function(x,n){filter(x,rep(1/n,n), sides=1)} 
mavback <- function(x,n){ filter(x, c(0, rep(1/n,n)), sides=1) }

Var_Real<-rtn^2
Vol_Real<-sqrt(Var_Real)
n<-nrow(Var_Real)
mopt<-c()
Var_Est<-c()

for (i in 1:ncol(Var_Real)){
  Var_Est5<-xts(mavback(Var_Real[,i], n=5),order.by=date(Var_Real))
  Var_Est10<-xts(mavback(Var_Real[,i], n=10),order.by=date(Var_Real))
  Var_Est20<-xts(mavback(Var_Real[,i], n=20),order.by=date(Var_Real))
  Var_Est40<-xts(mavback(Var_Real[,i], n=40),order.by=date(Var_Real))
  Var_Est[[i]]<-cbind(Var_Est5,Var_Est10,Var_Est20,Var_Est40)
  names(Var_Est[[i]])<-c("M5","M10","M20","M40")
  
  Dif_Est<-cbind((Var_Real[,i]-Var_Est5)^2, (Var_Real[,i]-Var_Est10)^2, (Var_Real[,i]-Var_Est20)^2, (Var_Real[,i]-Var_Est40)^2)
  names(Dif_Est)<-c("M5","M10","M20","M40")
  
  n<-nrow(rtn)
  H<-n-colSums(is.na(Dif_Est))
  RSME<-sqrt(colMeans(na.omit(Dif_Est)))
  mopt[i]<-which.min(RSME)
}

tit<-c("Estimación M40","Estimación M20","Estimación M40")

for (i in 1:length(clave)){
  par(mfrow=c(1,2))
  print(plot(sqrt(Var_Real[,i]),main=clave[i]))
  print(plot(sqrt(Var_Est[[i]][,mopt[i]]), main=tit[i], col = "blue"))
}

#EWMA
#Varianza y Volatilidad real observada

Var_Est<-matrix(0,n-1,ncol(rtn))
Func<-matrix(0,n-1,ncol(rtn))
Var_Est[1,]<-Var_Real[2,]

a<-seq(.1,.99,by=.01)

L<-.99
FMaxAnt<-0
Loptimo<-c()
num<-c()
for (j in 1:length(precio)){
  FMaxAnt<-0
  for (L in a){
    for (i in 2:(n-1)) {
      Var_Est[i,j] <- (1-L)*Var_Real[i-1,j]+L*Var_Est[i-1,j]
      Func[i,j]<--log(Var_Est[i,j])-Var_Real[i,j]/Var_Est[i,j]
    }
    FMaxAct<-colSums(Func)[j]
    if(FMaxAct>FMaxAnt){
      FMaxAnt<-FMaxAct
      FMax<-FMaxAnt
      Loptimo[[j]]<-L
    }}}

for (j in 1:length(precio)){
  for (i in 2:(n-1)) {
    Var_Est[i,j] <- (1-Loptimo[j])*Var_Real[i-1,j]+Loptimo[j]*Var_Est[i-1,j]
  }
}

for (i in 1:length(clave)){
  par(mfrow=c(2,1))
  print(plot(sqrt(Var_Real[,i]),main=clave[i]))
  print(plot(sqrt(Var_Est)[,i],type="l"))
}

#GARCH


##1 CHECAR INDEPENDENCIA DE LOS RENDIMIENTOS

#Ho:Existe independencia entre los rezagos
#Ha: No existe independencia entre los rezagos
## se debe rechazar Ho, para que exista dependencia de los rezagos

# SI Pvalue es menor a 0.05 se rechaza Ho al 5%

#coredata(datos) #para limpiar la info

pvalue<-sapply(coredata(preciotodos),Box.test,type="Ljung",lag=10)[3,]
pvalue<-lapply(coredata(preciotodos),Box.test,type="Ljung",lag=10)

Box.test(preciotodos[,1])
lapply(preciotodos,Box.test)

#X-squared = 20194, df = 10, p-value < 2.2e-16
qchisq(0.99,10) ##se distribuye como una Chi con 10 grados de libertad
#  23.20925

###2 PROBAR LOS EFECTOS ARCH 
#Ver si es adecuado el uso de un modelo de volatilidad para la serie

#Ho: No hay efectos ARCH/GARCH 
#Ha: Si hay efectos ARCH/GARCH

pvalues<-c()
for (i in 1:length(clave)){
  pvalues[i]<-ArchTest(precio[[i]],lag=10)$p.value
}
qchisq(0.99,10)
#Chi-squared = 2074.6, df = 10, p-value < 2.2e-16

##MODELO ARIMA PARA LA MEDIA

for (i in 1:length(clave)){
  par(mfrow=c(2,1))
  acf(precio[[i]],main=clave[i])
  pacf(precio[[i]],main=clave[i])
}

#Parece que no es estacionaria por lo que se hace prueba

pvalues<-c()
for (i in 1:length(clave)){
  pvalues[i]<-adf.test(precio[[i]])$p.value
}

#H0: Serie no estacionaria
#Ha: Serie estacioanria

#Entonces hay que hacer la primera diferencia de la serie

Dif_Precio<-na.omit(diff(preciotodos))

for (i in 1:length(clave)){
  pvalues[i]<-adf.test(Dif_Precio[,i])$p.value
}

for (i in 1:length(clave)){
  par(mfrow=c(2,1))
  acf(Dif_Precio[,i],main=clave[i])
  pacf(Dif_Precio[,i],main=clave[i])
}

##d es la derivada
##p se ve en la FACP
##q se ven en la FAC
##Probar con q=0, 1 y p=1

#ARA p=1 q=1 o 3
#AUTLANB p=1 q=1
#AXTELCPO p=1 q=2 probamos con 1 tambn
modelo<-c()
for (i in 1:length(clave)){
  arima311 <- arima(precio[[i]],order=c(3,1,1))
  arima211 <- arima(precio[[i]],order=c(2,1,1))
  arima111 <- arima(precio[[i]],order=c(1,1,1))
  arima310 <- arima(precio[[i]],order=c(3,1,0))
  arima210 <- arima(precio[[i]],order=c(2,1,0))
  arima110 <- arima(precio[[i]],order=c(1,2,0))
  
  aic311 <- arima311$aic
  aic211 <- arima211$aic
  aic111 <- arima111$aic
  aic310 <- arima310$aic
  aic210 <- arima210$aic
  aic110 <- arima110$aic
  
  modelo[i]<-which.min(c(aic311,aic211,aic111,aic310,aic210,aic110))
}


#Parte 3 Valor en riesgo VaR

#Portafolio
varcov<-cov(rtn)
mup<-do.call(cbind,lapply(rtn,mean))
sim<-nrow(rtn)
w1<-matrix(0,sim,3)
z<--qnorm(.99)
for (i in 1:sim){
  w1[i,]<-runif(n=3)
}
w1<-w1/rowSums(w1)
rp<-as.numeric(tail(mu,1))
sp<-as.numeric(tail(sigma,1))

res<-matrix(0,nrow(w1),ncol(w1)+1)
res[,1:3]<-w1/rowSums(w1)
colnames(res)<-c("w1","w2","w3","VaR")



for(j in 1:nrow(w1)){
  VaR<-z*sqrt(w1[j,]%*%varcov%*%w1[j,])+w1[j,]%*%t(mup)
  res[j,4]<-VaR
}

res<-as.data.frame(res)
VaR_ordenado<- res[order(-res$VaR),]
woptimo<-VaR_ordenado[1,1:3]
colnames(woptimo)<-clave
kable(woptimo,format="html",caption="Pesos optimos") %>% kable_styling(font_size = 7)

#Parametrico
dt<-c(1,5,10,20)
varpar<-matrix(0,3,4)
VaRpar<-c()
rownames(varpar)<-c("99%","95%","90%")
colnames(varpar)<-paste(dt,"dia/s")
s<-as.numeric(tail(sigma,1))
m<-as.numeric(tail(mu,1))

for (j in 1:length(clave)){
  for (i in 1:length(dt)){
    varpar[1,i]<-(qnorm(.99)*s[j] + m[j])*sqrt(dt[i])*-1
    varpar[2,i]<-(qnorm(.95)*s[j] + m[j])*sqrt(dt[i])*-1
    varpar[3,i]<-(qnorm(.90)*s[j] + m[j])*sqrt(dt[i])*-1
  }
  VaRpar[[j]]<-varpar
}
for (i in 1:length(dt)){
  varpar[1,i]<-(qnorm(.99)*sqrt(as.numeric(woptimo)%*%varcov%*%t(woptimo))+as.numeric(woptimo)%*%t(mup))*-sqrt(dt[i])
  varpar[2,i]<-(qnorm(.95)*sqrt(as.numeric(woptimo)%*%varcov%*%t(woptimo))+as.numeric(woptimo)%*%t(mup))*-sqrt(dt[i])
  varpar[3,i]<-(qnorm(.90)*sqrt(as.numeric(woptimo)%*%varcov%*%t(woptimo))+as.numeric(woptimo)%*%t(mup))*-sqrt(dt[i])
}
VaRpar[[j+1]]<-varpar

for(i in 1:length(nombres)){
  print(kable(VaRpar[[i]],format="html",caption=nombres[i]) %>% kable_styling(font_size = 7))
}

#Precios historicos

S<-(tail(preciotodos,1))
S<-cbind(S,sum(as.numeric(S)*woptimo))
dia_inicial<-date(S)
dia_est <-dia_inicial+1
dt<-c(1,5,10,20)
rtnp<-rtn*w1
nombres<-c(clave,"Portafolio")

VaRhist<-c()
varhist<-matrix(0,3,4)
rownames(varhist)<-c("99%","95%","90%")
colnames(varhist)<-paste(dt,"dia/s")

vep<-rtn%*%t(woptimo)
vep<-cbind(rtn,vep)

for (i in 1:length(nombres)){
  for (j in 1:length(dt)){
    
    varhist[1,j]<-quantile(vep[,i],.01)*sqrt(dt[j]) #al 99% de confianza la perdida maxima es VaR
    varhist[2,j]<-quantile(vep[,i],.05)*sqrt(dt[j])
    varhist[3,j]<-quantile(vep[,i],.1)*sqrt(dt[j])
  }
  VaRhist[[i]]<-as.data.frame(varhist)
}

for(i in 1:length(nombres)){
  print(kable(VaRhist[[i]],format="html",caption=nombres[i]) %>% kable_styling(font_size = 7))
}

for (i in 1:length(nombres)){
  for (j in 1:length(dt)){
    par(mai=c(1,1,1,1))
    h<-hist(vep[,i]*sqrt(dt[j]),nclass=30,main=c(nombres[i],paste(dt[j],"dia/s")),xlab="Pesos")
    y<-max(h$counts)
    h<-y/6
    x<-max(vep[,i]*sqrt(dt[j]))+min(vep[,i]*sqrt(dt[j]))/2
    abline(v=VaRhist[[i]][1,j],col="red",lwd=2)
    text(x,y,paste("Valor inicial",S[[i]]),col="black",cex=.7)
    text(x,y-h,paste("Fecha inicial",dia_inicial),col="black",cex=.7)
    text(x,y-2*h,paste("Fecha Est",dia_est),col="black",cex=.7)
    text(x,y-3*h,paste("VaR99",formatC(VaRhist[[i]][1,j],format="f",2)),col="red",cex=.7)
    abline(v=VaRhist[[i]][2,j],col="blue",lwd=2)
    text(x,y-4*h,paste("VaR95",formatC(VaRhist[[i]][2,j],format="f",2)),col="blue",cex=.7)
    abline(v=VaRhist[[i]][3,j],col="purple",lwd=2)
    text(x,y-5*h,paste("VaR90",formatC(VaRhist[[i]][3,j],format="f",2)),col="purple",cex=.7)
  }
}

#Montecarlo
n<-nrow(preciotodos)
##Aleatorios
epsilon<-rnorm(n)
pe_a1 <- matrix(0,n,length(dt))
pe<-c()
mup<-c(mu[6,],as.numeric(mu[6,]%*%t(woptimo)))
sigmap<-c(sigma[6,],sqrt(as.numeric(woptimo)%*%varcov%*%t(woptimo)))

for (x in 1:length(S)){
  for (i in 1:length(dt)){
    for (j in 1:n) {
      pe_a1[j,i]<-S[[x]]*exp((mup[x]-.5*sigmap[x]^2)*dt[i]+sigmap[x]*sqrt(dt[i])*epsilon[j])
    }
  }
  pe[[x]]<-pe_a1
}


pe[[x+1]]<-vep
var<-matrix(0,3,4)
VaRmont<-c()
rownames(var)<-c("99%","95%","90%")
colnames(var)<-paste(dt,"dia/s")

for (i in 1:length(nombres)){
  for (j in 1:length(dt)){
    var[1,j]<-quantile(pe[[i]][,j],.01) #al 99% de confianza la perdida maxima es VaR
    var[2,j]<-quantile(pe[[i]][,j],.05)
    var[3,j]<-quantile(pe[[i]][,j],.1)
  }
  VaRmont[[i]]<-var/as.numeric(S[1,i])-1
}


for(i in 1:length(nombres)){
  print(kable(VaRmont[[i]],format="html",caption=nombres[i]) %>% kable_styling(font_size = 7))
}

for (i in 1:length(nombres)){
  for (j in 1:length(dt)){
    par(mai=c(1,1,1,1))
    h<-hist((pe[[i]]/as.numeric(S[1,i])-1)*sqrt(dt[j]),nclass=30,main=c(nombres[i],paste(dt[i],"dia/s")),xlab="Pesos")
    y<-max(h$counts)
    h<-y/6
    x<-max((pe[[i]]/as.numeric(S[1,i])-1)*sqrt(dt[j]))+min((pe[[i]]/as.numeric(S[1,i])-1)*sqrt(dt[j]))/2
    abline(v=VaRmont[[i]][1,j],col="red",lwd=2)
    text(x,y,paste("Valor inicial",S[1,i]),col="black",cex=.7)
    text(x,y-h,paste("Fecha inicial",dia_inicial),col="black",cex=.7)
    text(x,y-2*h,paste("Fecha Est",dia_est),col="black",cex=.7)
    text(x,y-3*h,paste("VaR99",formatC(VaRmont[[i]][1,j],format="f",2)),col="red",cex=.7)
    abline(v=VaRmont[[i]][2,j],col="blue",lwd=2)
    text(x,y-4*h,paste("VaR95",formatC(VaRmont[[i]][2,j],format="f",2)),col="blue",cex=.7)
    abline(v=VaRmont[[i]][3,j],col="purple",lwd=2)
    text(x,y-5*h,paste("VaR90",formatC(VaRmont[[i]][3,j],format="f",2)),col="purple",cex=.7)
  }
}

#Backtest
backtest <- function(var,p,w) {
  
  S<-(tail(p,1))
  S<-cbind(S,sum(as.numeric(S)*w))
  bt<-matrix(0,3,4)
  rownames(bt)<-c("99%","95%","90%")
  colnames(bt)<-dt
  b<-c()
  for (i in 1:length(var)){
    for (j in 1:4){
      vp_real<-p[,i]
      per_gan_obs<-na.omit(diff(vp_real))
      
      btest99<-ifelse(per_gan_obs<var[[i]][1,j]*S[1,i],1,0)
      btest95<-ifelse(per_gan_obs<var[[i]][2,j]*S[1,i],1,0)
      btest90<-ifelse(per_gan_obs<var[[i]][3,j]*S[1,i],1,0)
      
      n1<-nrow(p)
      eV99<-(sum(btest99)/n1)
      eV95<-(sum(btest95)/n1)
      eV90<-(sum(btest90)/n1)
      
      bt[1,j]<-ifelse(eV99<=.01,"Adecuado","Excede")
      bt[2,j]<-ifelse(eV95<=.05,"Adecuado","Excede")
      bt[3,j]<-ifelse(eV90<=.1,"Adecuado","Excede")
    }
    b[[i]]<-bt
  }
  return(b)
}

portafolio<-preciotodos%*%t(woptimo)
preciosport<-cbind(preciotodos,portafolio)

backtestpar<-backtest(VaRpar,preciosport,woptimo)
backtesthist<-backtest(VaRhist,preciosport,woptimo)
backtestmont<-backtest(VaRmont,preciosport,woptimo)

for(i in 1:length(nombres)){
  print(kable(backtestpar[[i]],format="html",caption=paste("Parametrico",nombres[i])) %>% kable_styling(font_size = 7))
}

for(i in 1:length(nombres)){
  print(kable(backtesthist[[i]],format="html",caption=paste("Precios Historicos",nombres[i])) %>% kable_styling(font_size = 7))
}

for(i in 1:length(nombres)){
  print(kable(backtestmont[[i]],format="html",caption=paste("Montecarlo",nombres[i])) %>% kable_styling(font_size = 7))
}

# Rendimiento de portafolio y comparacion
rendimiento<-mu[6,]%*%t(woptimo)

kable(paste(round(rendimiento*100,4),"%"),format="html",col.names = "Rendimiento") %>% kable_styling()

comparacion<-matrix(0,3,4)
sumapar<-VaRpar[[1]]+VaRpar[[2]]+VaRpar[[3]]
kable(sumapar,format="html",caption= "Suma de VaR de activo por separado") %>% kable_styling()
kable(VaRpar[[4]],format="html",caption = "VaR portafolio") %>% kable_styling()


# Criterio Sharpe

v0<-100000
varcov<-cov(rtn)
sim<-nrow(rtn)
w1<-matrix(0,sim,3)
for (i in 1:sim){
  w1[i,]<-runif(n=3)
}
w1<-w1/rowSums(w1)
rp<-as.numeric(tail(mu,1))
sp<-as.numeric(tail(sigma,1))
rtnp<-rowSums(t(t(w1)*rp))
sigmap<-c()
for (i in 1:nrow(w1)){
  sigmap[i]<-sqrt(w1[i,]%*%varcov%*%w1[i,])
}
sharpe<-rtnp/sigmap
wsharpe<-do.call(cbind,list(w1,sharpe))

woptimo<-t(matrix(wsharpe[which.max(wsharpe[,4]),1:3]))
colnames(woptimo)<-clave
kable(woptimo,format="html",caption="Pesos optimos") %>% kable_styling(font_size = 7)

