triangle <-
function(NDVI,LST,Tact,Rad,
	dem = 100, ti.method = "linear", delta.method = "batra", pressure.method = "barometric",gamma.method = "wang",outname = "ndvi-ts-plot",
	rn = c(),g.flux = c(),fao.par1 = 0.0864, fao.par2 = 0.408, rn.perc = 0.55, g.perc  = 0.07, out.unit = "mm",
	write.ras= F, write.plot =F ,outpath.ras = getwd(), outname.ras = "ET_Act.if",ta.conv = "F", int.method = "ngb",classify = T,...){


Tact.Text = mean(Tact[],na.rm=T)
Rad.Text = mean(Rad[],na.rm=T)
Dem.Text = mean(dem,na.rm=T)
if(class(NDVI)[[1]]!="RasterLayer"){
NDVI = raster(NDVI)
}

if( max(NDVI[],na.rm=T) > 1){ stop("NDVI can not be > 1!")}
if( min(NDVI[],na.rm=T) < -1){ stop("NDVI can not be < -1!")}

if(class(LST)[[1]]!="RasterLayer"){
LST = raster(LST)
}
LST = resample(LST,NDVI,method = int.method)

if(!is.numeric(Tact)){
Tact.Text = names(Tact)
if(class(Tact)[[1]]!="RasterLayer"){
Tact = raster(Tact) 
}
Tact = resample(Tact,NDVI,method = int.method)
}else{
Temp = NDVI
Temp[] = Tact
Tact = Temp
}

if(!is.numeric(dem)){
Dem.Text = names(dem)
if(class(dem)[[1]]!="RasterLayer"){
dem = raster(dem) 
}
dem = resample(dem,NDVI,method = int.method)
z = dem
}else{
Temp = NDVI
Temp[] = dem
z = Temp
}

if(!is.numeric(Rad)){
Rad.Text =names(Rad)
if(class(Rad)[[1]]!="RasterLayer"){
Rad = raster(Rad) 
}
Rad =  resample(Rad, NDVI,method = int.method)
}else{
Temp = NDVI
Temp[] = Rad
Rad = Temp
}
text.ta  = paste("obs. temperature:",Tact.Text)
text.rad  = paste("obs. radiation:",Rad.Text)
text.dem = paste("Altitude:",Dem.Text)
text.ndvi = paste("NDVI-file:",names(NDVI))
text.lst = paste("LST-file:",names(LST))
InfoFile = paste(outname.ras,"_LOg.txt",sep="")
file.create(InfoFile)
if( ta.conv==T){ Tact = Tact - 273.5 }

if(classify==T){
a = seq((-1),0.95,by=0.05)
b = seq((-0.95),1,by=0.05)
 c = seq((-0.95),1,by=0.05)
 abc = cbind(a,b,c)
 class.NDVI = reclassify(NDVI,abc)
}else{
class.NDVI = NDVI
}

stats.out = ndvilst.plot(NDVI,LST,classify = classify,write = write.plot,outname=outname,...)

t = stats.out$t.max
y = stats.out$t.min
y1 = stats.out$y1
y2 = stats.out$y2
x1 = stats.out$x1
x2 = stats.out$x2
m  = stats.out$m

print(m)

text.tmin = paste("t.min:",round(y,2))
text.tmax = paste("t.max:",round(t,2))
if(write.plot==T){
text.plot = paste("Corresponding NDVI-LST plot: '",outname,".png'",sep="")
}else{
text.plot = "Corresponding NDVI-LST plot: No .png file written (write.plot = F)"
}

setwd(outpath.ras)


#==== Calculation of PhiMin and Phi following Stisen et al. 2007: ==============
#PhiMin=max(min((1.26*( (  class.NDVI - min(class.NDVI[],na.rm=T) ) / 
#( max(class.NDVI[],na.rm=T)  - min(class.NDVI[],na.rm=T) )
#)^2)
#, 1.26),0)

PhiMin=1.26*( (  class.NDVI - min(class.NDVI[],na.rm=T) ) / 
( max(class.NDVI[],na.rm=T)  - min(class.NDVI[],na.rm=T) )
)^2


if(ti.method == "tang"){
Ti.max= y1 + PhiMin/1.26 * (y2 - y1)
}else{
Ti.max = m*class.NDVI + t
}

PhiList=min(((Ti.max - LST) /
(Ti.max - y))
* (1.26 - PhiMin)+ PhiMin
, 1.26)

  
#==== Calculation of evaporative fraction (EF):=================================================
##############
# Delta:= the slope of the saturated vapor pressure at given temperature, as defined by Murray (1967), see Batra et al. 2006 for details.
##############
P0= 1013.15# atm. pressure at sea level

### After BATRA et al. 2006, Delta in [hPa]: 
if(delta.method=="batra"){
text.delta = "Delta calculated acc. to Batra et al. 2006"
print(text.delta)
delta = (26297.76/((Tact + 273.15) - 29.65)^2) * exp((17.67 * Tact)/((Tact + 273.15) - 29.65))
}
if(delta.method == "fao"){
text.delta = "Delta calculated acc. to FAO (Murray 1966)"
print(text.delta)
delta = 4098 * (0.6108 * exp((17.27*Tact)/(Tact + 237.3)))/(Tact + 237.3)^2
delta = delta * 10
}
if(delta.method == "wang"){
text.delta = "Delta calculated acc. to Wang et al. 2006"
print(text.delta)
Tr = 1 - 373.15/(Tact+273.15)
e_star = P0 * exp(13.3185*Tr - 1.976*Tr^2 - 0.6445*Tr^3 - 0.1299*Tr^4)
delta = (373.15*e_star/(Tact+273.15)^2)*(13.3185 - 3.952*Tr - 1.9335*Tr^2 - 0.5196*Tr^3)
}


####################
#  P (pressure at given temp in average catchment height z):
####################
### After WANg et al. 2006:
if(pressure.method == "wang"){
text.pressure = "Pressure calculated acc. to Wang et al. 2006"
 print(text.pressure)
P = P0*10^(((-1)*z/18400)*((Tact+273.15)/273.15)) 
}
### After ROEDEL (1994) in [hPa]:
if (pressure.method == "barometric"){
text.pressure = "Pressure calculated acc. barometric approach (Roedel 1994)"
print(text.pressure)
g = 9.81 # m/s^2
R = 8.31451 #J/(mol*K)
M = 0.02897 # kg/mol mitllere Molmasse der Luft
P = P0 * exp((-1)*(z*M*g)/(R*(Tact + 273.15)))
}


###################
# gamma (psychrometric constant):
###################
### After WANg et al . 2006, gamma in [hPa/degC]:
if(gamma.method=="wang"){
text.gamma = "gamma calculated acc. to Wang et al. 2006"
print(text.gamma)
Cp  = 1013  # [J/kgdegC], specific heat capacity of air at constant P
lamda = 4.2*(597-0.6*(Tact))*1000
 gamma =(Cp * P) / (0.622 * lamda)

}
### After FAO
if(gamma.method == "fao"){
text.gamma = "gamma calculated acc. to method of FAO"
print(text.gamma)
epsilon = 0.622 # ratio molecular weight of water vapour to dry air
lamda.fao = 2.45 # MJ/kg, latent heigh vaporisation
Cp = 1.1013 * 10^(-3) # MJ/(kg*degC), specific heat capacity of air at constant P
gamma = (Cp * P) / (epsilon*lamda.fao)
}

###################
# Final calculation of EF:
###################

print("EF")
EF = PhiList* (delta/(delta + gamma))
print("EF-done")

if(is.null(rn)){
rn = Rad * rn.perc  # MJ/(m^2*d) ### R: net radiation (~ 55% of glob. radiation) (Davies [1967])
}else{
rn = Rad *rn
}
print("rn-done")
if(is.null(g.flux)){
print("g=0")
g.flux  = Rad * g.perc  # MJ/(m^2*d) ### g: soil heat flux (~ 7% of R) (Crago and Brutsaert  [1996])
}else{
print(g.flux)
g = Rad * .flux 
print(g.flux)

}
print("g-done")
if(out.unit == "mm"){
text.format = "Unit of estimated ETact is: mm/d"
ETact  =max(((EF *(rn - g))*fao.par1)*fao.par2,0)
}
print("ETact-done1")
if(out.unit == "watt"){
text.format = "Unit of estimated ETact is: W/m^2"
ETact  =max(((EF *(rn - g))),0)
}
print("ETact-done2")

unit=c()
if(out.unit == "mm"){unit=1}
if(out.unit == "watt"){unit=1}
if(is.null(unit)){ 
print("Unknown output units - out.unit must either be 'mm' or 'watt'. Unit for ETact will now be set to default (mm per day)!")
text.format = "Unit of estimated ETact is: mm/d"
ETact  =max(((EF *(rn - g))*fao.par1)*fao.par2, 0)
}

print("ETact-done3")

if(write.ras == T){
if(!is.null(outpath.ras)){setwd(outpath.ras)}
    writeRaster(ETact , filename=outname.ras , "gTiff", full.names = F,overwrite=T, NAflag = -9999) 
	text.outname = paste("Final ETact-file: '",outname.ras,"'",sep="")
}else{
 text.outname = "Final ETact-file: No .tif file written (write.ras = F)"
}


text.abstr1 = "Summary for the calculation of actual evapotranspiration"
text.abstr2 = "using the triangle method approach by Jiang and Islam 1999"
text.empty = "======================================================"
text.time = as.character(Sys.time())
text.minet = paste("min. ET act:", round(min(ETact[],na.rm=T),2))
text.maxet = paste("max. ET act:", round(max(ETact[],na.rm=T),2))
text.meanet = paste("mean ET act:", round(mean(ETact[],na.rm=T),2))
InfoFile = file(InfoFile)
writeLines(c(text.abstr1,text.abstr2,text.empty,text.ta,text.rad,text.dem,text.ndvi,text.lst,text.empty,text.tmin,text.tmax,text.delta,text.pressure,text.gamma,text.empty,text.plot,text.outname,text.format,text.empty,text.meanet,text.maxet,text.minet,text.empty,text.time),InfoFile,sep="\n")
close(InfoFile)
print("write text file --> done")

return(ETact)

invisible(y)
}
