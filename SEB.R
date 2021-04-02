library(raster)
library(rgdal)

albedo <- brick("C:/Users/17252302/Downloads/Sentinel2/UCD-Lyons/Albedo_all.tif")
LAI <- brick("C:/Users/17252302/Downloads/Sentinel2/UCD-Lyons/LAI_all.tif")
NDVI <- brick("C:/Users/17252302/Downloads/Sentinel2/UCD-Lyons/NDVI_all.tif")
df <- read.table(file.choose(), header=T, sep=',',fill = T)

rastercon=function(condition, trueValue,falseValue){
  return (condition*trueValue+(!condition)*falseValue)
}
##########################################

#################################################

#########  Surfex function  ########################

surfex <- function(S,Ta,pp,RH,u,sT,sm,albedo,LAI,NDVI,DEM, iter.max=10,t1=1, model="SEBAL")#UseMethod ("sebal") 
{
  
  #if(!class(S) == "RasterLayer") stop( "S is not a raster object")
 # if(!class(Ta) == "RasterLayer") stop( "Ta is not a raster object")
 # if(!class(pp) == "RasterLayer") stop( "pp is not a raster object")
 # if(!class(Td) == "RasterLayer") stop( "Td is not a raster object")
 # if(!class(u) == "RasterLayer") stop( "u is not a raster object")
  if(!class(albedo) == "RasterLayer") stop( "albedo is not a raster object")
  if(!class(LAI) == "RasterLayer") stop( "LAI is not a raster object")
  if(!class(NDVI) == "RasterLayer") stop( "NDVI is not a raster object")
  #if(!class(DEM) == "RasterLayer") stop( "DEM is not a raster object")
  #if(!class(sm) == "RasterLayer") stop( "Ts is not a raster object")
 # if(!class(sT) == "RasterLayer") stop( "Ts is not a raster object")
 # if(!class(sdrain) == "SpatialPolygonsDataFrame") stop( "sdrain is not a spatial polygon object")
  
  z=2
  k=0.41
  Ag=9
  adiab=0.01
  #pp <- (pp/1000)*(1-((0.01*DEM)/(Ta+(0.01*DEM))))
  den = (1000*pp)/(1.01*((Ta)*287))
  #den=1.225
  cp=1005
  den_cp=den*cp
  Lv=2450000
  rd=287
  rv=462
  rd_rv=rd/rv
  stef=5.67E-08
  
  
  Tac=Ta-273.15
  es=0.6108*exp(((17.27*(Ta-273.15))/((Ta-273.15)+237.3)))
  ea=RH/100*es
  #ea=0.6108*exp(((17.27*(Td-273.15))/((Td-273.15)+237.3)))
  dq=(621.9907*es/(pp-es)) - (621.9907*ea/(pp-ea))
  psyc=(1005*pp*462)/(287*2450000)
  eslope=(4098*(0.6108*exp((17.27*Tac)/(Tac+237.3))))/(Tac+237.3)^2
  
  
  ## surface roughness parameters
  print(paste("Computing surface momentum and thermal roughness lengths"))
  zom = 0.123 * (LAI/24)
  zoh = 0.1 * zom
  #d = 2/3*(LAI/24)
  ##### surface and atmospheric emissivity
  print(paste("Computing surface and atmospheric emissivities"))
  er=1.2*((ea*10)/Ta)^0.143
  e0=rastercon(NDVI<0 & albedo<0.47, 0.985, rastercon(LAI>=3, 0.98, 0.95+(LAI*0.01)))
  
  ##soil moisture stress function parameters
  print(paste("Computing environmental stress functions"))
  #y=1/(1+0.16*(dq-3))
  #Fdq=ifelse(y>1,1,y)
  Fdq = rastercon((1/(1+0.16*(dq-3)))>1,1,1/(1+0.16*(dq-3)))
  #Fdq = ifelse((1/(1+0.16*(dq-3)))>1,1,1/(1+0.16*(dq-3)))
  Fs=(770*S)/(1000*S+230*(1000-2*S))
  #Fsm=ifelse(sm<0.3,1+6.3*(sm-0.3),ifelse(1+6.3*(sm-0.3)<0,0.01,1))
  #Fm <- rastercon(levels(sdrain@data$Drainage_C) == "Poor",rastercon(sm<0.3,1+6.3*(sm-0.3), 
  #               rastercon(1+6.3*(sm-0.3)<0,0.01,1)),rastercon(sm<0.3,1+4.3*(sm-0.3), 
  #                rastercon(1+4.3*(sm-0.3)<0,0.01,1)))
  Fm = rastercon(sm<0.3,1+4.3*(sm-0.3), rastercon(1+4.3*(sm-0.3)<0,0.01,1))
  #Fm = ifelse(sm<0.3,1+4.3*(sm-0.3),1)
  
  print(paste("Computing u_star, ra, H and L for neutral condition "))
  u200=u*(log(200/2)/log(10/2))
  ustar=(u200*k)/(log(200/zom))
  #ustar=(u*k)/(log(10/zom))
  ra=(log(z/zoh)/(k*ustar))
  ra= rastercon(ra >1000,1000,ra)
  rc=(0.47 * 110/LAI)/(Fdq*Fs*Fm)
  K=(1-albedo)*S
  Lin = (er*stef*Ta^4) #+ 60*N
  R= (eslope+(psyc*(1+(rc/ra))))
  A = ((eslope+psyc)*R)-(eslope*(eslope+psyc))
  B = -1*(eslope+psyc)
  C = (eslope+psyc)*R
  D = K+Lin+(3*e0*stef*Ta^4)+(Ag*sT)-((4*e0*stef*Ta^3)+Ag)*(Ta+adiab*2)
  E = ((4*e0*stef*Ta^3)+Ag)*(ra/den_cp)
  H = (A*D+B)/(C+A*E)
  tvs=(-1*H)/(den_cp*ustar)
  L=(Ta*ustar^2)/(k*9.8*tvs)
  
  print(paste("Computing u_star, ra, H and L with stability correction"))
  i=1
  while (i<=iter.max) 
  {
    
    print(paste("Monin-Obukhov length iteration",i,"of",iter.max))
    
    x_bus200=(1-(16*(200/L)))^0.25
    x_bus01=(1-(16*(zoh/L)))^0.25
    x_bus=(1-(16*(2/L)))^0.25
    #x_bus=0.5^0.25
    
    stab_u200=rastercon(L<0,(2*log((1+x_bus200)/2))+(log((1+x_bus200^2)/2))-(2*atan(x_bus200))+(pi/2),rastercon(L>0,-5*(2/L),0))
    #stab_u=rastercon(L<0,(2*log((1+x_bus)/2))+(log((1+x_bus^2)/2))-(2*atan(x_bus))+(pi/2),rastercon(L>0,-5*(2/L),0))
    
    stab_t01=rastercon(L<0,2*log((1+x_bus01^2)/2),rastercon(L>0,-5*(zoh/L),0))
    stab_t=rastercon(L<0,2*log((1+x_bus^2)/2),rastercon(L>0,-5*(2/L),0))
    
    ustar_new = (u200*k)/(log(200/zom)-stab_u200)
    #ustar_new =(u*k)/(log(10/zom)-(stab_u*(10/L))+(stab_u*(zom/L)))
    ra=(1/(ustar_new*k))*(log(z/zoh) - stab_t + stab_t01)
    #ra=(1/(ustar_new*k))*(log(z/zoh)-(stab_t*(z/L))+(stab_t01*(zoh/L)))
    ra= rastercon(ra >1000,1000,ra)
    rc=(0.47 * 110/LAI)/(Fdq*Fs*Fm)
    K=(1-albedo)*S
    er=1.2*((ea*10)/Ta)^0.143
    Lin = (er*stef*Ta^4) #+ 60*N
    R= (eslope+(psyc*(1+(rc/ra))))
    A = ((eslope+psyc)*R)-(eslope*(eslope+psyc))
    B = -1*(eslope+psyc)
    C = (eslope+psyc)*R
    D = K+Lin+(3*e0*stef*Ta^4)+(Ag*sT)-((4*e0*stef*Ta^3)+Ag)*(Ta+adiab*2)
    E = ((4*e0*stef*Ta^3)+Ag)*(ra/den_cp)
    H = (A*D+B)/(C+A*E)
    tvs=(-1*H)/(den_cp*ustar_new)
    L=(Ta*ustar_new^2)/(k*9.8*tvs)
    
    i=i+1
  }
  
  print(paste("Computing the final outputs"))
  Ts=(Ta+((H*ra)/den_cp)+(DEM*0.01)) 
  Lou=(e0*stef*Ts^4) + (1-e0)*Lin
  Go = Ag*(Ts-sT)
  Rn = (K + (er-1)*(er*stef*Ta^4)) - (e0*stef*4*Ta^3*(Ts-Ta))
  #G= Rn*((Ts-273.15)*(0.0038+0.0074*albedo)*(1-0.98*NDVI^4))
  #G = rastercon(LAI < 0.5, (Rn*(1.8*(Ts-273.15)/Rn) + 0.084), Rn*(0.05 +(0.18*exp(-0.521*LAI))))
  Le = (eslope*(Rn - Go) + ((den_cp*(es-ea))/ra)) / (eslope + psyc*(1 + rc/ra))
  #Lee = (eslope*(Rn - G) + ((den_cp*(es-ea))/ra)) / (eslope + psyc*(1 + rc/ra))
  #Leb = Rn-H-Go
  #Lebb = Rn-H-G
  
  ETint=(t1*3600*Le)/Lv
  
  factor=list(ra=ra,rc=rc, Rn=Rn, H=H, Go=Go,Le=Le,ETins=ETint,
              Ts=Ts, Lou=Lou, Lin=Lin)
  
  factor$call<-match.call()
  
  class(factor)<-model
  factor
  
}

for (i in 1:12){
  mod<-surfex(S=df[i,2],Ta=df[i,4],pp=df[i,5],RH=df[i,6],u=df[i,7],sT=df[i,8],sm=df[i,9],
              albedo[[i]],LAI[[i]],NDVI[[i]],DEM=df[1,12], iter.max=7,t1=1, model="SEBAL")
  writeRaster(mod$Rn,filename = paste0("C:/Users/17252302/Downloads/Sentinel2/UCD-Lyons/output/Rn/", 
                                       "Rn_",stringr::str_pad(string = i,width = 4,side = "left",pad=0)), 
              format ="GTiff", overwrite=T)
  writeRaster(mod$H,filename = paste0("C:/Users/17252302/Downloads/Sentinel2/UCD-Lyons/output/H/", 
                                      "H_",stringr::str_pad(string = i,width = 4,side = "left",pad=0)), 
              format ="GTiff", overwrite=T)
  writeRaster(mod$Le,filename = paste0("C:/Users/17252302/Downloads/Sentinel2/UCD-Lyons/output/Le/", 
                                       "LE_",stringr::str_pad(string = i,width = 4,side = "left",pad=0)), 
              format ="GTiff", overwrite=T)
  writeRaster(mod$Go,filename = paste0("C:/Users/17252302/Downloads/Sentinel2/UCD-Lyons/output/G/", 
                                       "G_",stringr::str_pad(string = i,width = 4,side = "left",pad=0)), 
              format ="GTiff", overwrite=T)
  writeRaster(mod$ETins,filename = paste0("C:/Users/17252302/Downloads/Sentinel2/UCD-Lyons/output/ETins/", 
                                        "ETins_",stringr::str_pad(string = i,width = 4,side = "left",pad=0)), 
              format ="GTiff", overwrite=T)
  writeRaster(mod$Lou,filename = paste0("C:/Users/17252302/Downloads/Sentinel2/UCD-Lyons/output/Lou/", 
                                        "Lou_",stringr::str_pad(string = i,width = 4,side = "left",pad=0)), 
              format ="GTiff", overwrite=T)
  writeRaster(mod$Ts,filename = paste0("C:/Users/17252302/Downloads/Sentinel2/UCD-Lyons/output/Ts/", 
                                       "Ts_",stringr::str_pad(string = i,width = 4,side = "left",pad=0)), 
              format ="GTiff", overwrite=T)
  writeRaster(mod$rc,filename = paste0("C:/Users/17252302/Downloads/Sentinel2/UCD-Lyons/output/rc/", 
                                       "rc_",stringr::str_pad(string = i,width = 4,side = "left",pad=0)), 
              format ="GTiff", overwrite=T)
  writeRaster(mod$ra,filename = paste0("C:/Users/17252302/Downloads/Sentinel2/UCD-Lyons/output/ra/", 
                                       "ra_",stringr::str_pad(string = i,width = 4,side = "left",pad=0)), 
              format ="GTiff", overwrite=T)
}

EF <- Le/(Rn-G)

Rn_24 <- list()
ET_24 <- list()

for (i in 1:12) {
  Rn_24[[i]] <- ((1-albedo[[i]])*(df[i,11]) - (110*(0.75+0.00002*df[1,12])))
  ET_24[[i]] <- (86400000*1*EF[[i]]*Rn_24[[i]])/(1000*2450000)
}
ET_24 <- stack(ET_24)
Rn_24 <- stack(Rn_24)