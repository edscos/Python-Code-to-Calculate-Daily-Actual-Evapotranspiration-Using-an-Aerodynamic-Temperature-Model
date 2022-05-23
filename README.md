# Python Code to Calculate Daily Actual Evapotranspiration Using an Aerodynamic Temperature Model

Python code to calculate actual daily evapotranspiration (ETa) from an one-source surface energy balance approach (OSEB) using point-based data inputs.

# Author: Edson Costa-Filho
# Civil and Environmental Engineering Department
# Colorado State University, Fort Collins, CO

# REFERENCES:

Allen, R. G., Pereira, L. S., Raes, D., & Smith, M. (1998). Crop evapotranspiration-Guidelines for computing crop water requirements-FAO Irrigation and drainage paper 56. Fao, Rome, 300(9), D05109.

Allen, R. G., & Pereira, L. S. (2009). Estimating crop coefficients from fraction of ground cover and height. Irrigation Science, 28(1), 17-34.

ASCE-EWRI (2005) The ASCE standardized reference evapotranspiration equation. In: Allen RG, Walter IA, Elliott R, Howell T, Itenfisu D, Jensen M (eds.) ASCE-EWRI Task committee report. American Society of Civil Engineers, Reston, VA, 171 pp.

Brunsell, N. A., & Gillies, R. R. (2002). Incorporating surface emissivity into a thermal atmospheric correction. Photogrammetric engineering and remote sensing, 68(12), 1263-1270.

Brutsaert, W. (1975). On a derivable formula for longwave radiation from clear skies. Water resources research, 11(5), 742-744.

Brutsaert, W. (1982). Evaporation into the atmosphere: Theory. History, and Applications, 1.

Buck, A. L. (1981). New equations for computing vapor pressure and enhancement factor. Journal of Applied Meteorology and Climatology, 20(12), 1527-1532.

Businger, J. A., Wyngaard, J. C., Izumi, Y., & Bradley, E. F. (1971). Flux-profile relationships in the atmospheric surface layer. Journal of Atmospheric Sciences, 28(2), 181-189.

Choudhury, B. J., & Monteith, J. L. (1988). A four?layer model for the heat budget of homogeneous land surfaces. Quarterly Journal of the Royal Meteorological Society, 114(480), 373-398.

Costa-Filho, E., Chávez, J. L., & Comas, L. (2020). Determining maize water stress through a remote sensing-based surface energy balance approach. Irrigation Science, 1-18.

Costa-Filho, E., Chávez, J. L., Zhang, H., & Andales, A. A. (2021). An optimized surface aerodynamic temperature approach to estimate maize sensible heat flux and evapotranspiration. Agricultural and Forest Meteorology, 311, 108683.

Ham, J. M. (2005). Useful equations and tables in micrometeorology. Micrometeorology in agricultural systems, 47, 533-560.

Johnson, L. F., & Trout, T. J. (2012). Satellite NDVI assisted monitoring of vegetable crop evapotranspiration in California’s San Joaquin Valley. Remote Sensing, 4(2), 439-455.

Kustas, W. P., & Daughtry, C. S. (1990). Estimation of the soil heat flux/net radiation ratio from spectral data. Agricultural and Forest Meteorology, 49(3), 205-223.

Monteith, J. (1973). Principles of environmental physics edward arnold. London, 214p.

Su, Z. (2002). The Surface Energy Balance System (SEBS) for estimation of turbulent heat fluxes. Hydrology and earth system sciences, 6(1), 85-100.

Thom, A. 1. (1972). Momentum, mass and heat exchange of vegetation. Quarterly Journal of the Royal Meteorological Society, 98(415), 124-134.

Verma, S.B. (1989) Aerodynamic resistances to transfers of heat, mass and momentum. T.A. Black, D.L. Splittlehouse, M.D. Novak, D.T. Price (Eds.), Estimation of Areal Evapotranspiration, IAHS Publ. No. 177, International Association for Hydrological Science, Wallingford, UK (1989), pp. 13-20

Webb, E. K. (1970). Profile relationships: The log?linear range, and extension to strong stability. Quarterly Journal of the Royal Meteorological Society, 96(407), 67-90.

Zhang, Y., Liu, C., Lei, Y., Tang, Y., Yu, Q., Shen, Y., & Sun, H. (2006). An integrated algorithm for estimating regional latent heat flux and daily evapotranspiration. International Journal of Remote Sensing, 27(1), 129-152.

# CODE DESCRIPTION:

The Python script labeled "ETa_OSEB_RS_Code.py" performs the calculation of daily actual ETa using an OSEB approach for latent heat flux based on an optimized aerodynamic surface temperature (see Costa-Filho et al., 2021). The remote sensing input data for running this code is point-based measurements of the following variables:
    
   Nadir-looking surface temperature or Ts (Celsius)
   Normalized Difference Vegetation Index or NDVI

The ancillary data input for running this python script are the following variables:
  
   Air temperature above the canopy or Ta (Celsius)
   Incoming short-wave solar radiation or Rs (W/m^2)
   Relative humidity above the canopy or RH (%)
   Wind speed above the canopy (m/s)
   Wind direction above the canopy (Decimal Degrees)
   Hourly alfalfa-based reference ET (mm/h)
   Hourly grass-based reference ET (mm/h)
   Daily alfalfa-based reference ET (mm/d)
   Daily grass-based reference ET (mm/d)
   
There is an spreadsheet labeled "Input_Data_v3.xlsx" that provides a file to stored the input variables needed to calculate daily maize ETa.

The alfalfa and grass-based reference ET has to be calculated or obtained before running the Python script. Agricultural weather station data often provide reference ET values for a given location. Users might find the ETo calculator from the Food and Agriculture Organization (FAO) useful to calculate daily values of reference ET using agricultural weather station data (see https://www.fao.org/land-water/databases-and-software/eto-calculator/en/).

The code is designed to provide an exported .csv containing the following output variables:

    Day of the year or DOY (Dimensionless)
    Hourly instantaneous ETa (mm/h)
    Crop coefficient based on alfalfa reference ET (Dimensionless)
    Crop coefficient based on grass reference ET (Dimensionless)
    Daily ETa based on alfalfa reference ET (mm/d)
    Daily ETa based on grass reference ET (mm/d)
    Estimated net radiation flux (W/m^2)
    Estimated soil heat flux (W/m^2)
    Estimated sensible heat flux (W/m^2)
    Estimated latent heat flux (W/m^2)
    Nadir-looking surface temperature (Celsius)
    Air temperature above the canopy (Celsius)
    Estimated surface aerodynamic temperature (Celsius)
    Aerodynamic resistance to heat transfer (s/m)
    Surface shear velocity (m/s)
    Fractional vegetation cover (Dimensionless)
    Leaf area index (m^2/m^2)
    Surface albedo (Dimensionless)

# NOTE: Make sure to change the file directory in which the spreadsheet input data file is stored.


















