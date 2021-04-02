
########################################################
# This python code estimates day time surface radiation and energy fluxes over grass using a land surface scheme
# that incorporates limited routine weather data. The code also evaluate the sensitivity of the scheme
# to different parameterizations of surface resistance for contrasting soil moisture regimes.

######################################################

##########################################################################
# Please cite the following reference if you have used this code in any publication
# Citation: K. A. Ishola, G. Mills, R. M., Fealy, Ó. Ní Choncubhair, R. Fealy (2020). Improving a land
#                   surface scheme for estimating sensible and latent heat fluxes above grasslands with
#                   contrasting soil moisture zones. Agric. Fores Meteor. 294, 108151,
#                   https://doi.org/10.1016/j.agrformet.2020.108151
########################################################################

########   Define the required libraries ################
import math
from math import e
import matplotlib.pyplot as plt
import scipy.interpolate as sci
from scipy.stats import *
from matplotlib.pyplot import *
import numpy as np
from sklearn.metrics import mean_squared_error
from math import sqrt
import pandas as pd
import datetime
import matplotlib.dates as dates
####################################################

#########   Function to convert air temperature from celsius to kelvin  ################
def Tk(tc):
    x = tc + 273.15
    return(x)
##############################################################################

########## Function to compute water vapour fractional conductance ######################
def Fdq(dq):
    y = 1/(1 + 0.16*((dq) - 3))
    y = [1 if i > 1 else i for i in y]
    return(y)
######################################################################################

##########   Function to compute radiation fraction  ####################################
def Fs(S):
    z = (770*S)/ (1000*S + 230*(1000 - 2*S))
    return(z)
#########################################################################################

#######################  Function to compute soil moisture response    #############
def Fm(sm):
    u = [0.1 if i < 0.1 else i for i in sm]
    u = np.array([1 + 4.3*(i - 0.3) if  i < 0.3 else 1 for i in u])
    return(u)
############################################################################

#omit error values ############
np.seterr(divide='ignore', invalid='ignore')
##################################

#######  Function for the first loop i.e for neutral condition. ############
def it_1(ws,t24h,tc,P,S,rh):
 psyc = (1005*P*462)/(287*2450000) #(0.001005*P)/(0.622*2.45)
 eslope = 4098*((0.6108*np.exp((17.27*tc)/(tc+237.3)))/(tc+237.3)**2)
 ustar=(ws*k)/(math.log(10/zom,math.e))
 ra=(math.log(z/zoh,math.e)/(k*ustar))
 #rc=(10**4)*((es-ea)/P)*rd_rv       ###dRH99
 #rc=70                              ###FAO
 rc=25.9/(Fdq(dq)*Fs(S)*Fm(sm))      ###BB97
 K=(1-a)*S
 er=1.2*((ea*10)/Tk(tc))**0.143
 Lin = (er*stef*Tk(tc)**4 ) #+ 60*N
 R= (eslope+(psyc*(1+(rc/ra))))
 A = ((eslope+psyc)*R)-(eslope*(eslope+psyc))
 B = -1*(eslope+psyc)
 C = (eslope+psyc)*R
 D = K+Lin+(3*0.94*stef*Tk(tc)**4)+(Ag*t24h)-((4*0.94*stef*Tk(tc)**3)+Ag)*(Tk(tc)+adiab*2)
 E = ((4*0.94*stef*Tk(tc)**3)+Ag)*(ra/den_cp)
 H = (A*D+B)/(C+A*E)
 tvs=(-1*H)/(den_cp*ustar)
 L=(Tk(tc)*ustar**2)/(k*9.8*tvs)
 return(ra,H,L)
###########################################################################

####   Function for the second loop i.e for stablity correction.   ##############
def it_2(ws,t24h,tc,P,S,rh,H,L):
 psyc = (1005*P*462)/(287*2450000)
 eslope = 4098*((0.6108*np.exp((17.27*tc)/(tc+237.3)))/(tc+237.3)**2)
 stab_u=[(2*np.log((1+((1-(16*(z/i)))**0.25))/2))+(np.log((1+((1-(16*(z/i)))**0.25)**2)/2))-(2*np.arctan((1-(16*(z/i)))**0.25))+(np.pi/2) if i < 0 else -5*(z/i) for i in L]
 stab_t=[2*np.log((1+((1-(16*(zoh/i)))**0.25)**2)/2) if i < 0 else -5*(zoh/i) for i in L]
 ustar_new=(ws*k)/(math.log(10/zom,math.e)-(stab_u*(10/L))+(stab_u*(zom/L)))
 ra=(1/(ustar_new*k))*(math.log(z/zoh,math.e)-(stab_t*(z/L))+(stab_t*(zoh/L)))
 ra = np.array([1000 if i > 1000 else i for i in ra])
 #rc=(10**4)*((es-ea)/P)*rd_rv                   ##dRH99
 #rc=70                                          ###FAO
 rc=25.9/(Fdq(dq)*Fs(S)*Fm(sm))                  ###BB97
 K=(1-a)*S
 er=1.2*((ea*10)/Tk(tc))**0.143
 Lin = (er*stef*Tk(tc)**4) #+ 60*N
 R= (eslope+(psyc*(1+(rc/ra))))
 A = ((eslope+psyc)*R)-(eslope*(eslope+psyc))
 B = -1*(eslope+psyc)
 C = (eslope+psyc)*R
 D = K+Lin+(3*0.94*stef*Tk(tc)**4)+(Ag*t24h)-((4*0.94*stef*Tk(tc)**3)+Ag)*(Tk(tc)+adiab*2)
 E = ((4*0.94*stef*Tk(tc)**3)+Ag)*(ra/den_cp)
 H = (A*D+B)/(C+A*E)
 tvs=(-1*H)/(den_cp*ustar_new)
 L=(Tk(tc)*ustar_new**2)/(k*9.8*tvs)
 return(ra,H,L,ustar_new,Lin)
###################################################################################

####### import input data in csv format  (.csv) from the local directory ##################
data = pd.read_csv("C:/Users/17252302/Desktop/PhD_MU/data/point_estimated fluxes/Johnstown/test.csv")
data = data.convert_objects(convert_numeric=True)
#########################################

##### define the input and validation variables and in the data #####################################
date=data.iloc[:, 0]            #date
S=data.iloc[:, 1]               #global solar radiation (W m-2)
tc=data.iloc[:, 2]              #near-surface temperature at observation height (2 m) (oC)
P=data.iloc[:, 3]               # msl pressure (kPa)
rh=data.iloc[:, 4]              # Relative humdity (%)
ws=data.iloc[:, 5]              # Wind speed (m s-1)
t24h=data.iloc[:, 6]            # mean air temperature in the last 24hr (24hr moving average) (oK)
sm=data.iloc[:, 7]             # Measured soil moisture content at top 20 cm (m3 m-3)
Rn_obs=data.iloc[:, 8]          # Measured total radiative flux (W m-2)
H_obs=data.iloc[:, 9]           # Measured Sensible heat flux (W m-2)
Le_obs=data.iloc[:, 10]         # Measure Latent heat flux (W m-2)
##########################################################################################

####### define initial coefficients  ###########
a=0.23                   ##surface albedo for grass
z=2                     ##observation height (m)
zom=0.01                ## surface roughness length for momentum (m)
zoh=0.1*zom             ## surface roughness length for heat (m)
k=0.41                  ## von Karma constant
Ag=9                   ## soil heat flux coefficient (W m-2 K-1)
adiab=0.01              ##dry adiabatic lapse rate (K m-1)
den=1.225               ## density of dry air (kg m-3)
cp=1005                 ##specific heat capacity of dry air (J kg-1 K-1)
den_cp=den*cp
Lv=2450000              ## latent heat of vaporization (J kg-1)
rd=287                  ## specific gas constant for dry air (J kg-1 K-1)
rv=462                  ## specific gas constant for vapour (J kg-1 K-1)
rd_rv=rd/rv
stef=5.67*10**-8        ## stefan boltzmann's constant (W m-2 K-1)
###############################

###Compute moisture deficit ####################
es = 0.6108*(np.exp((17.27*tc)/(tc+237.3))) ###saturated vapour pressure (kPa)
ea = rh/100 * es              ####actual vapour pressure (kPa)
dq = (621.9907*es/(P-es)) - (621.9907*ea/(P-ea))  ###moisture deficit (kPa)
#VPD = es-ea
#########################################################################

######### Businger-Dyer calculations #####################
#x_bus=0.5**0.25 #  (1 - (19*z/L))**0.25
#stab_u=[(2*math.log((1+x_bus)/2))+(math.log((1+x_bus**2)/2))-(2*math.atan(x_bus))+(math.pi/2) if i < 0 else -5*(z/i) for i in L]
#stab_t=[2*math.log((1+x_bus**2)/2) if i < 0 else -5*(z/i) for i in L]
##########################################################

#### iteration 1 loop setup ####################
loop_setup=it_1(ws,t24h,tc,P,S,rh)
#######################################################

#### iteration2: loop action   #############################
loop_action=it_2(ws,t24h,tc,P,S,rh,loop_setup[2],loop_setup[1])
loop_action_2=it_2(ws,t24h,tc,P,S,rh,loop_action[2],loop_action[1])
loop_action_3=it_2(ws,t24h,tc,P,S,rh,loop_action_2[2],loop_action_2[1])
loop_action_4=it_2(ws,t24h,tc,P,S,rh,loop_action_3[2],loop_action_3[1])
loop_action_5=it_2(ws,t24h,tc,P,S,rh,loop_action_4[2],loop_action_4[1])
loop_action_6=it_2(ws,t24h,tc,P,S,rh,loop_action_5[2],loop_action_5[1])
##############################################################################################

####  define output objects  ############################
L=loop_action_6[2]              # Obukhov length (m)
H=loop_action_6[1]              # estimated Sensible heat flux (W m-2)
ra=loop_action_6[0]             # estimated aerodynamic resistance (s m-1)
ustar=loop_action_6[3]          # estimated friction velocity (m s-1)
Lin=loop_action_6[4]            # estimated longwave radiation downward (W m-2)
##############################################

####### Compute surface temperature, outgoing longwave, ground heat flux, and Net radiation
eslope = 4098*((0.6108*np.exp((17.27*tc)/(tc+237.3)))/(tc+237.3)**2)
er=1.2*((ea*10)/Tk(tc))**0.143          # apparent atmospheric emissivity
psyc =(1005*P*462)/(287*2450000)        # pychrometic constant (kPa K-1)
K=(1-a)*S                               # net shortwave radiation (W m-2)
ts=(Tk(tc)+((H*ra)/den_cp)+(z*0.01))    # Land surface temperature (oK)
Lou=(0.94*stef*ts**4) + (1-0.94)*Lin    # longwave radiation upward (W m-2)
Go = Ag*(ts-t24h)                       # Soil heat flux (W m-2)
Rn = (K + (er-1)*(er*stef*Tk(tc)**4)) - (0.94*stef*4*Tk(tc)**3*(ts-Tk(tc))) # estimated net radiative flux (W m-2)
##########################################################################

############### Compute latent heat flux with different parameterizations of rc ######################
#rc= (10**4)*((es-ea)/P)*rd_rv        ####### dRH99 (de Rooy and Holtslag, 1999)
#rc=70                                ###### FAO    (Allen et al., 1998)
rc=25.9/(Fdq(dq)*Fs(S)*Fm(sm))       ###### BB97    (Beljaars and Bosveld, 1997)
Le = (eslope*(Rn - Go) + ((den_cp*(es-ea))/ra)) / (eslope + psyc*(1 + rc/ra))  ## P-M approach
#Le = Rn - Go - H         #### Balance method
##############################################################

##### export data in csv format to local directory #######################
my_data=np.vstack((date,Rn, H, Le, Go, rc, ra, Lou, Lin,S, ts))
my_data=my_data.T
df= pd.DataFrame(my_data)
colnames=['date','Rn','H','Le','G','rc','ts']
df.to_csv('C:/Users/17252302/Desktop/PhD_MU/data/point_estimated fluxes/Johnstown/flux_test.csv',
         index=False, header=colnames) # write output to csv
#######################################################################

#######################################################
#####     Compute RMSE and bias with nan values for Le   ########################
mask = ~np.isnan(Le_obs) & ~np.isnan(Le)
rms = sqrt(mean_squared_error(Le_obs[mask], Le[mask]))
bias = sum(Le[mask] - Le_obs[mask])/8760
###########################################################

#######################################################
#####     Compute RMSE and bias with nan values for H   ########################
maskk = ~np.isnan(H_obs) & ~np.isnan(H)
rms_h = sqrt(mean_squared_error(H_obs[maskk], H[maskk]))
bias_h = sum(H[maskk] - H_obs[maskk])/8760
###########################################################

#######################################################
#####     Compute RMSE and bias with nan values for Rn   ########################
masck = ~np.isnan(Rn_obs) & ~np.isnan(Rn)
rms_Rn = sqrt(mean_squared_error(Rn_obs[masck], Rn[masck]))
bias_Rn = sum(Rn[masck] - Rn_obs[masck])/8760
###########################################################


############# scatterplots to visualize the relationship between measured and estimated Le########
slope, intercept, r_value, p_value, std_err = stats.linregress(Le_obs[mask], Le[mask])
print(r_value,slope,intercept,p_value)
fig = plt.figure()
ax = fig.add_subplot(111)
plt.plot(Le_obs,Le, '.',c="black",label="BB97")
plt.plot(Le_obs, slope*Le_obs + intercept, '-',c="r",label="r = 0.76")
#label axes
xlabel("Measured QE (W m^-2)")
ylabel("Estimated QE (W m^-2)")
ax.text(0.1, 0.9,'', horizontalalignment='center',
     verticalalignment='center',
     transform=ax.transAxes)
plt.legend()
plt.show()
#fig.savefig('C:\\Users\\admin\\Desktop\\Oak\\Le06_Schm4.png',dpi=600,transparent=True)
#############################################################################
