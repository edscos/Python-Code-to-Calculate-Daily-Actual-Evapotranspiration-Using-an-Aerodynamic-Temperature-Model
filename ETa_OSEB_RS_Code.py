# Python Code to run OSEB ETa using and optimized aerodynamic temperature approach
# Author: Edson Costa-Filho, M.Sc.
# Date: 5/21/2022

# References:

# Allen, R. G., Pereira, L. S., Raes, D., & Smith, M. (1998).
# Crop evapotranspiration-Guidelines for computing crop
# water requirements-FAO Irrigation and drainage paper 56.
# Fao, Rome, 300(9), D05109.

# Allen, R. G., & Pereira, L. S. (2009). Estimating crop
# coefficients from fraction of ground cover and height.
# Irrigation Science, 28(1), 17-34.

# ASCE-EWRI (2005) The ASCE standardized reference evapotranspiration
# equation. In: Allen RG, Walter IA, Elliott R, Howell T, Itenfisu
# D, Jensen M (eds.) ASCE-EWRI Task committee report. American
# Society of Civil Engineers, Reston, VA, 171 pp.

# Brunsell, N. A., & Gillies, R. R. (2002).
# Incorporating surface emissivity into a thermal atmospheric correction.
# Photogrammetric engineering and remote sensing, 68(12), 1263-1270.

# Brutsaert, W. (1975). On a derivable formula for longwave radiation
# from clear skies.
# Water resources research, 11(5), 742-744.

# Brutsaert, W. (1982). Evaporation into the atmosphere: Theory.
# History, and Applications, 1.

# Buck, A. L. (1981). New equations for computing vapor pressure
# and enhancement factor.
# Journal of Applied Meteorology and Climatology, 20(12), 1527-1532.

# Businger, J. A., Wyngaard, J. C., Izumi, Y., & Bradley, E. F. (1971).
# Flux-profile relationships in the atmospheric surface layer.
# Journal of Atmospheric Sciences, 28(2), 181-189.

# Choudhury, B. J., & Monteith, J. L. (1988). A four?layer model for
# the heat budget of homogeneous land surfaces.
# Quarterly Journal of the Royal Meteorological Society,
# 114(480), 373-398.

# Costa-Filho, E., Chávez, J. L., & Comas, L. (2020). Determining maize
# water stress through a remote sensing-based
# surface energy balance approach.
# Irrigation Science, 1-18.

# Costa-Filho, E., Chávez, J. L., Zhang, H., & Andales, A. A. (2021).
# An optimized surface aerodynamic temperature approach to
# estimate maize sensible heat flux and evapotranspiration.
# Agricultural and Forest Meteorology, 311, 108683.

# Ham, J. M. (2005). Useful equations and tables in micrometeorology.
# Micrometeorology in agricultural systems, 47, 533-560.

# Johnson, L. F., & Trout, T. J. (2012).
# Satellite NDVI assisted monitoring of vegetable
# crop evapotranspiration in California’s San Joaquin Valley.
# Remote Sensing, 4(2), 439-455.

# Kustas, W. P., & Daughtry, C. S. (1990). Estimation of the soil heat flux/net radiation
# ratio from spectral data. Agricultural and Forest Meteorology, 49(3), 205-223.

# Monteith, J. (1973). Principles of environmental physics edward arnold. London, 214p.

# Su, Z. (2002). The Surface Energy Balance System (SEBS) for
# estimation of turbulent heat fluxes.
# Hydrology and earth system sciences, 6(1), 85-100.

# Thom, A. 1. (1972). Momentum, mass and heat exchange of vegetation.
# Quarterly Journal of the Royal Meteorological Society, 98(415),
# 124-134.

# Verma, S.B. (1989) Aerodynamic resistances to transfers of heat,
# mass and momentum. T.A. Black, D.L. Splittlehouse, M.D. Novak, D.T.
# Price (Eds.), Estimation of Areal Evapotranspiration,
# IAHS Publ. No. 177, International Association for
# Hydrological Science, Wallingford, UK (1989), pp. 13-20

# Webb, E. K. (1970). Profile relationships: The log?linear range, and
# extension to strong stability.
# Quarterly Journal of the Royal Meteorological Society,
# 96(407), 67-90.

# Zhang, Y., Liu, C., Lei, Y., Tang, Y., Yu, Q., Shen, Y., & Sun, H. (2006).
# An integrated algorithm for estimating regional latent heat flux and daily evapotranspiration.
# International Journal of Remote Sensing, 27(1), 129-152.

from math import *
# Packages for data analysis
import pandas as pd
import numpy as np

# Package for writing files as csv
import csv

# Packages for visualizing data
import imageio
from PIL import Image

# Defining Constant Variables Regarding the Date of the Image

L = 0.50  # Constant for Calculating SAVI
# Low vegetation: L = 1.00
# Intermediate vegetation: L = 0.50
# High vegetation: L = 0.25

Elev = 1270  # Site Ground elevation in m
Z_u = 3.5  # ON-SITE Wind reference height [m]
Z_T = 3.5  # ON-SITE Air Temp reference height [m]

# NOTE 1: These parameters do not need to be changed.

g = 9.81  # Gravitational Acceleration [m/s^2]
k = 0.41  # Von Karman constant
df = 3600  # Data record frequency [sec]
min_u = 0.50  # Lowest allowed wind speed value [m./s]

# Importing ancillary data from excel (local weather and reference ET data)
# NOTE: Add the directory where the input data spreadsheet is located
data_excel = pd.read_excel(r'file_directory\Input_Data_v3.xlsx')
all_columns = list(data_excel.columns)
print(all_columns)
Date = data_excel['Date']  # Date of the Data (MM/DD/YYYY)
Time = data_excel['Time']  # Time of the Data (HH:MM:SS)
DOY = data_excel['DOY']  # Day of the Year (DOY)
Ta = data_excel['Ta']  # Air Temperature data (C)
Rs_in = data_excel['Rs']  # Incoming Short-wave Solar radiation (W/m^2)
RH = data_excel['RH']  # Relative Humidity (%)
U = data_excel['U']  # Horizontal Wind Speed (m/s)
WD = data_excel['WD']  # Wind Direction (Degrees)
ETr_h = data_excel['ETrh']  # Hourly Alfalfa-based Reference ET (mm/hr)
ETo_h = data_excel['EToh']  # Hourly Grass-based Reference ET (mm/hr)
ETr_d = data_excel['ETrd']  # Daily Alfalfa-based Reference ET (mm/d)
ETo_d = data_excel['ETod']  # Daily Grass-based Reference ET (mm/d)
Ts = data_excel['Ts']  # Nadir-looking Surface Temperature (C)
NDVI = data_excel['NDVI']  # Normalized Difference Vegetation Index (NDVI)

# Obtaining the number of rows and columns of the dataset
size = range(len(Ta))  # Number of rows in the spreadsheet dataset

# Assigning initial values for each variable to be calculated in the code
LAI = 0
ALBEDO = 0
Rs_out = 0
Rl_in = 0
Rl_out = 0
Rn = 0
G = 0
X = 0
To = 0
U_star = 0
rah = 0
H = 0
L_mo = 0
iteration = 0
state = 0
x1 = 0
phi_h1 = 0
phi_m1 = 0
LE = 0
ET_inst = 0
Kc_alfalfa = 0
Kc_grass = 0
ETa_Kc_alfalfa = 0
ETa_Kc_grass = 0
Es = 0
d = 0
Zom = 0
Zoh = 0
hc = 0
rp = 0
fc_norman = 0
LAI_L = 0
fs_kustas = 0
fc = 0
CF = 0

# for loop to calculate SEB fluxes, instantaneous and daily actual ETa

for row in size:

    # Complementary Micrometeorological Variables

    # ON-SITE latent heat of vaporization (J./Kg)
    lambda_v = (2.501 - 0.00236 * Ta) * 10 ** 6

    # Atmospheric pressure (Pa) (As used in Ham, 2005)
    Pa = (101.3 * np.exp((-0.0342 * Elev) / (Ta + 273.15))) * 1000

    # Saturated vapor pressure (kPa) (Buck, 1981)
    es = 0.61121 * np.exp(17.502 * Ta / (Ta + 240.97))

    # ON-SITE actual vapor pressure (kPa)
    ea = es * RH / 100

    # Converting the Actual Vapor Pressure from kPa to mb
    ea_mb = 10 * ea

    # ON-SITE virtual temperature (K) (As used in Ham, 2005)
    Tkv = (Ta + 273.15) / (1 - 0.378 * ea / Pa / 1000)

    # ON-SITE reference height (m) for wind speed measurement
    Zm = Z_u

    # Gas constant for dry air (J/Kg/K)
    R_d = 287.04

    # Air density (kg/m^3) (As used in Ham, 2005)
    rho_a = (Pa / (R_d * (Ta + 273.15))) * (1 - ((0.378 * (1000 * ea)) / Pa))

    # Specific Heat of Air (J/kg/K) (As used in Ham, 2005)
    Cpa = 1004.7 * (1 + (0.522 * (1000 * ea) / Pa))

    # Converting the air temperature from Celsius degree to Kelvin
    T_aK = Ta + 273.15

    # Air emissivity Calculation (Brutsaert, 1975)
    Ea = 1.24 * ((ea_mb / T_aK) ** (1 / 7))

    # Converting surface temperature from Celsius to Kelvin
    T_sK = Ts + 273.15

    # Calculation of the turbulent mixing-row resistance (s/m): Costa-Filho et al. (2021)
    if [WD <= 90]:
        rp = ((90 - WD) / (90 + WD)) * (1 / U)
    elif [90 < WD <= 180]:
        rp = ((WD - 90) / (270 - WD)) * (1 / U)
    elif [180 < WD <= 270]:
        rp = ((270 - WD) / (WD - 90)) * (1 / U)
    else:
        rp = ((WD - 270) / (450 - WD)) * (1 / U)

    # Leaf Area Index Calculation (Zhang et al., 2006)
    LAI = 0.60 * np.exp(3.00 ** NDVI)

    # Fractional Vegetation Cover Calculation (Norman et al., 1995)
    fc_norman = 1 - np.exp(-0.50 * LAI)

    # Local LAI calculation (Kustas and Norman, 2000)
    LAI_L = LAI / fc_norman

    # Fractional soil cover (Kustas and Norman, 2000)
    fs_kustas = fc_norman * np.exp(-0.50 * LAI_L) + (1 - fc_norman)

    # Clumping factor (Kustas and Normam, 2000)
    CF = -np.log(fs_kustas) / (0.50 * LAI)

    # Updated Fractional Vegetation Cover (Kustas and Norman, 2000)
    fc = 1 - np.exp(-0.50 * CF * LAI)

    # Calculating the Surface Emissivity from Fractional Vegetation Cover
    e_v = 0.98  # Surface emissivity for a fully vegetated canopy (Campbell and Norman, 1998)
    e_s = 0.93  # Surface emissivity for a bare soil (Campbell and Norman, 1998)

    # Surface Emissivity equation (Brunsell and Gillies, 2002)
    Es = fc * e_v + (1 - fc) * e_s

    # Albedo Calculation (Costa-Filho et al., 2020)
    ALBEDO = 0.0179 * LAI ** (-3.27) + 0.1929

    # Net Radiation Flux (Rn) Calculation
    Rs_out = ALBEDO * Rs_in  # Rs_out: Outgoing Short-wave Radiation in W/m^2
    Rl_in = Ea * (5.67 * 10 ** (-8)) * T_aK ** 4  # Rl_in: Incoming Long-wave Radiation in W/m^2
    Rl_out = Es * (5.67 * 10 ** (-8)) * T_sK ** 4  # Rl_out: Outgoing Long-wave Radiation in W/m^2
    Rn = Rs_in - Rs_out + Rl_in - Rl_out  # Rn: Estimated Net Radiation Flux in W/m^2

    # Soil Heat Flux (G) Calculation (Su, 2002)
    Gamma_veg = 0.05  # G/Rn ratio for bare soil surfaces (Kustas and Daughtry, 1990)
    Gamma_soil = 0.315  # G/Rn ratio for fully vegetated surfaces (Monteith, 1973)
    G = Rn * (Gamma_veg + (1 - fc) * (Gamma_soil - Gamma_veg))

    # Canopy Height Calculation (Costa-Filho et al., 2020)
    hc = 0.697 * np.exp(0.236 * LAI) - 3.42 * np.exp(-3.177 * LAI)

    # Calculating the Canopy Roughness Length Terms
    X = 0.20 * LAI  # Fraction of green LAI

    # Zero-displacement Height Calculation (Monteith and Choudhury, 1988)
    d = hc * (np.log(1 + (X ** (1 / 6))) + 0.03 * np.log(1 + (X ** 6)))

    # Roughness Length for Momentum Transfer (Monteith and Choudhury, 1988)
    if [X < 0.20]:
        Zom = 0.01 + 0.28 * hc * np.sqrt(X)
    else:
        Zom = 0.3 * hc * (1 - d / hc)

    # Roughness length for heat transfer (Brutsaert, 1982)
    Zoh = Zom * 0.10

    # Calculating the Surface Aerodynamic Temperature (Costa-Filho et al., 2021)
    if [LAI <= 1.50]:
        To = -8.742 * fc + 0.571 * Ta + 0.529 * Ts + 0.806 * rp + 3.295
    elif [1.50 < LAI <= 2.50]:
        To = -9.168 * fc + 0.485 * Ta + 0.575 * Ts - 0.160 * rp + 6.491
    elif [2.50 < LAI <= 3.50]:
        To = 4.708 * fc + 0.350 * Ta + 0.580 * Ts + 0.086 * rp
    else:
        To = -1.912 * fc + 0.443 * Ta + 0.509 * Ts + 0.115 * rp + 5.014

    # Monin-Obhukov Stability Theory (MOST) for Atmospheric Stability Corrections
    # First Iteration, Assuming Neutral Atmospheric Conditions:
    U_star = (U * k) / (np.log((Zm - d) / Zom))
    rah = np.log((Zm - d) / Zom) * np.log((Zm - d) / Zoh) / (U * k ** 2)
    H = (rho_a * Cpa) * (To - Ta) / rah
    L_mo = -(U_star ** 3 * (Ta + 273.15) * rho_a * Cpa) / (g * k * H)
    x1 = 0
    phi_h1 = 0
    phi_m1 = 0

    # Stability according to the Monin-Obukhov Stability Length (L_mo)
    a = 0
    b = 0

    if [[L_mo < 0] and [L_mo > -500]]:  # Unstable Atmospheric Conditions (Ts > Ta)
        for i in range(50):
            a = a + 1

            # Equations from Paulson (1970) and Businger et al. (1971)
            x1 = (1 - 16 * (Zm - d) / L_mo) ** 0.25
            phi_h1 = 2 * np.log((1 + x1 ** 2) / 2)
            phi_m1 = 2 * np.log((1 + x1) / 2) + np.log((1 + x1 ** 2) / 2) - 2 * (
                    np.pi / 4 * x1 + 0.285 * x1 * (1 - np.abs(x1))) + pi / 2

            # Thom (1972); Brutsaert (1982); Verma (1989)
            U_star = (U * k) / (np.log((Zm - d) / Zom - phi_m1))  # Surface Shear Velocity (m/s)
            rah = (np.log((Zm - d) / Zoh) - phi_h1) / (U_star * k)  # Aerodynamic Resistance (s/m)
            H = (rho_a * Cpa) * (To - Ta) / rah
            L_mo = -(U_star ** 3 * (Ta + 273.15) * rho_a * Cpa) / (g * k * H)

    if [[L_mo > 0] and [L_mo < 500]]:  # Stable Atmospheric Conditions (Ts < Ta)
        for i in range(50):
            b = b + 1

            # Equations from Paulson (1970) and Businger et al. (1971)
            phi_h1 = -5 * ((Zm - d) / L_mo)
            phi_m1 = phi_h1

            # Thom (1972); Brutsaert (1982); Verma (1989)
            U_star = (U * k) / (np.log((Zm - d) / Zom - phi_m1))  # Surface Shear Velocity (m/s)
            rah = (np.log((Zm - d) / Zoh) - phi_h1) / (U_star * k)  # Aerodynamic Resistance (s/m)
            H = (rho_a * Cpa) * (To - Ta) / rah
            L_mo = -(U_star ** 3 * (Ta + 273.15) * rho_a * Cpa) / (g * k * H)

    if [[L_mo > 500] or [L_mo < -500]]:
        # Thom (1972); Brutsaert (1982); Verma (1989)
        U_star = (U * k) / (np.log((Zm - d) / Zom))
        rah = np.log((Zm - d) / Zom) * np.log((Zm - d) / Zoh) / (U * k ** 2)
        H = (rho_a * Cpa) * (To - Ta) / rah
        L_mo = -(U_star ** 3 * (Ta + 273.15) * rho_a * Cpa) / (g * k * H)

    # Calculation of Instantaneous Latent Heat Flux (LE)
    LE = Rn - G - H

    # Converting Instantaneous Latent Heat Flux to Instantaneous ET in water depth (mm/h)
    ET_inst = df * LE / lambda_v / 1

    # Calculating the Instantaneous Evaporative Fraction considering alfalfa ref ET
    Kc_alfalfa = ET_inst / ETr_h

    # Calculating the Instantaneous Evaporative Fraction considering grass ref ET
    Kc_grass = ET_inst / ETr_h

# Calculation of Daily ETa (Alfalfa and Grass Extrapolation Methods)
ETa_Kc_alfalfa = Kc_alfalfa * ETr_d
ETa_Kc_grass = Kc_grass * ETo_d

# Exporting a List Variable into the .csv File
data = [DOY, ET_inst, Kc_alfalfa, Kc_grass, ETa_Kc_alfalfa, ETa_Kc_grass,
        Rn, G, H, LE, Ts, Ta, To, rah, U_star, fc, LAI, ALBEDO]

# Converting list to array format using the numpy package
out_arr = np.array(data)
# Transposing the data from row to vector columns
OUTPUT = np.transpose(out_arr)

# Creating the .csv File gets created in the current working directory
# Adding headers to each column in the future .csv file
headerList = ['DOY', 'ET_inst', 'Kc_alfalfa', 'Kc_grass', 'ETa_Kc_alfalfa', 'ETa_Kc_grass',
              'Rn', 'G', 'H', 'LE', 'Ts', 'Ta', 'To', 'rah', 'U_star', 'fc', 'LAI', 'Albedo']
# Exporting the results as .csv file
df = pd.DataFrame(OUTPUT)
df.to_csv('OSEB_To_ETa_Results.csv', header=headerList)
