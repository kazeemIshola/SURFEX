# SURFLEX python code for simulating surface radiation and energy flux exchanges


This python script estimates daytime surface radiation and energy fluxes over grass using a land surface scheme
that incorporates limited routine weather data. The script also evaluates the sensitivity of the scheme
to different parameterizations of surface resistance for contrasting soil moisture regimes.


The scheme is adapted from de Rooy and Holtslag (1999), and the tested surface resistance
models are based on FAO method (Allen et al., 1998) and modified Jarvis approach (Jarvis, 1976), implemented by
Beljaars and Bosveld (1997). 


The python code is based on the study carried out over temperate grassland with contrasting soil moisture zones (wet and free-draining)
in Ireland (Ishola et al., 2020).


# Acknowledgment
The project is funded by Teagasc under Walsh Fellowship programme. We acknowledge the providers of data used as forcing and evaluation of the scheme. The meteorological data are available from MET Eiréann, Irish Meteorological Service (https://www.met.ie/) and eddy covariance flux measurements are available from European Fluxes Database Cluster (http://www.europe-fluxdata.eu/). 


# References

Allen, R. G., Pereira, L. S., Raes, D. and Smith, M., 1998. Crop evapotranspiration. Guidelines for computing crop water requirements. 
      Irrigation and Drainage Paper No. 56. Rome: FAO.

Beljaars, A. C. M. and Bosveld, F. C., 1997. Cabauw Data for the Validation of Land Surface Parameterization Schemes J. Climate 10, 1172 – 1193.

de Rooy, W. C.  and Holtslag, A. A. M., 1999. Estimation of surface radiation and energy ﬂux densities from single-level weather data.
      J. Appl. Meteor., 38, 526– 540.
      
Ishola, K. A., G. Mills, R. M., Fealy, Ó. Ní Choncubhair, R. Fealy (2020). Improving a land surface scheme for estimating sensible and latent heat fluxes 
      above grasslands with contrasting soil moisture zones. Agric. Fores Meteor. 294, 108151, https://doi.org/10.1016/j.agrformet.2020.108151
      
Jarvis, P., 1976. The interpretation of leaf water potential and stomatal conductance found in canopies in the ﬁeld. 
      Philos. Trans. Roy. Soc. London B, 273, 593–610.
