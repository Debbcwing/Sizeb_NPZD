import os
import numpy as np
import time
import xarray as xr
from scipy.interpolate import make_interp_spline

from SizeModLogspace_clean import SizeModLogspace

Pnum = 150      # specify the nos of phytoplankton size groups
Ynum = 10       # no. of modeling year
mod_name = 'StdRun_clean'

# Import forcing data
wd = 'yourwd'
os.chdir(wd)

# Temperature and PAR data
LWST = np.genfromtxt('Forcing/LWST_13p_syn.csv', delimiter=',')
nSSI = np.genfromtxt('Forcing/nSSI_13p_syn.csv', delimiter=',')
for a in range(8):  # To synchonized the initial and end forcing value
    LWST[0, a + 1] = LWST[12, a + 1]
    nSSI[0, a + 1] = nSSI[12, a + 1]


# Converting nSSI into PAR
def nSSI2PAR(netSSI):
    """
    This function is to convert nSSI (W m-2) into PAR (Ein m-2 d-1/mol m-2 d-1)
    PAR is a term used to describle radiation in wavelengths usedful for photosynthesis [1]
    PAR is generally accepted to be wavelengths between 400 and 700 nm [1]
    Procedure:
    1. extr = extraction factor (only about 45% of the nSSI is in the 400-700nm range)
    2. conv = W m-2 converts into μmole m-2 s-1 = 4.57 [1,2]
    3. sec2day = converting from s-1 to d-1
    4. μmol2mol = converting from μmol to mol (1 mole = 1 Ein)

    [1] Sager and Farlane, 1997
    Sager, J.C., McFarlane, J.C. 1997. Radiation. In: Plant Growth Chamber Handbook. (Langhans, R.W. and
    Tibbitts, T.W., Eds.) Iowa State University: North Central Regional Research Publication No. 340, Iowa
    Agriculture and Home Economics Experiment Station Special Report no. 99, pp. 1-29.
    [2] Thimijan and Heins, 1983
    Thimijan, R.W., and R.D. Heins. 1983. Photometric, radiometric and quantum light units of measure:
    A review of procedures for interconversion. HortScience 18(6): 818-822
    """
    extr = 0.45
    conv = 4.57
    sec2day = 86400
    μmol2mol = 1e6
    coef = extr * conv * sec2day / μmol2mol  # 0.1776816
    return netSSI * coef


PAR = nSSI2PAR(nSSI)
PAR[:, 0] = nSSI[:, 0]

# Creating forcing curves
LWST_spl = make_interp_spline(LWST[:, 0], LWST[:, 5], k=2)
LWST_spl = LWST_spl(np.arange(365))
LWST_spl[364] = LWST_spl[0]
PAR_spl = make_interp_spline(PAR[:, 0], PAR[:, 5], k=2)
PAR_spl = PAR_spl(np.arange(365))
PAR_spl[364] = PAR_spl[0]

# Sinusoidal mixing forcing
"""
Mix depth of 82.125 is adopted here, the figure is obtained from the mean mixing depth from 8 lakes in Switzerland
     Lake name          Code    Zmix    ThermoZ
# 1. Greifensee         GR      30.0    4.30
# 2. Lake Zürich        LZ      135.0   5.80    
# 3. Hallwilersee       HA      48.0    2.15
# 4. Sempachersee       SE      87.0    0.5
# 5. Baldeggersee       BA      67.0    2.16
# 6. Vierwaldstätersee  VWS/LU  110.0   3.47
# 7. Upper Zürich       UZ      36.0    0.5
# 8. Walensee           WA      144.0   0.5
                                82.125  2.42
"""
# To simplify value for modeling: Mixing depth = 82m; thermocline depth = 2.5m
Zmix = 80
ThermoZ = 2.5
numMixRgm = 3
mld_sinu = np.zeros((365, numMixRgm))  # days; mixing regimes
mld_sinu[:, 0] = np.full((365,), Zmix)  # constant mixing depth
mld_sinu[:, 1] = (-((Zmix-ThermoZ)/2+ThermoZ) + (Zmix-ThermoZ)/2 * np.cos(np.arange(365) / 14.525)) * -1
    # medium mixing frequency
mld_sinu[:, 2] = (-((Zmix - ThermoZ) / 2 + ThermoZ) + (Zmix - ThermoZ) / 2 * np.cos(np.arange(365) / 4.825)) * -1
    # high mixing frequency
for i in range(numMixRgm):
    mld_sinu[364, i] = mld_sinu[0, i]       # synchonize the intial and end forcing value


# generate dmdt arrays
dmdt_sinu = np.zeros((365, numMixRgm))
for i in range(numMixRgm):
    dmdt_sinu[0:364, i] = np.diff(mld_sinu[:, i])  # dmdt (364,)
    dmdt_sinu[364, i] = dmdt_sinu[0, i]

# checking initial and end forcing values
LWST_spl[0], LWST_spl[364], PAR_spl[0], PAR_spl[364]
for i in range(numMixRgm):
    print(mld_sinu[0, i])
    print(mld_sinu[364, i])
    print(dmdt_sinu[0, i])
    print(dmdt_sinu[364, i])


mld_sinu = np.tile(mld_sinu, (Ynum, 1))  # multiplicating physical forcing arrays
dmdt_sinu = np.tile(dmdt_sinu, (Ynum, 1))
LWST_spl = np.tile(LWST_spl, Ynum)
PAR_spl = np.tile(PAR_spl, Ynum)

# Runs - loop
#   Nutrient schemes: Oligotrophic (1.0); Mesotrophic (15.0); Eutrophic (50.0)
#   Mixing regimes: Meromictic (just for observation); polymictic; dimictic; monomictic
N0_list = [1, 15, 50]
sol = np.zeros((Ynum*365, Pnum + 4, len(N0_list), numMixRgm))  # time; no. var; N0 schemes; mixing regimes

# Run solution
starttime = time.ctime()        # save it for documentation purposes
print(mod_name, str(Pnum), "start:", starttime)
start = time.time()             # start the timer

for i in range(len(N0_list)):   # nutrient level
    for j in range(numMixRgm):  # mixing regimes
        sol[:, :, i, j] = SizeModLogspace(mld=mld_sinu[:, j], par=PAR_spl, sst=LWST_spl, dmdt=dmdt_sinu[:, j],
                                          N0=N0_list[i], numP=Pnum, numYears=Ynum).solution

# Complete
end = time.time()  # stop the timer
endtime = time.ctime()  # save it for documentation purposes
print(mod_name, str(Pnum), "complete:", endtime)
runtime = (end - start) / 3600  # calculate run time (hour)
print(mod_name, str(Pnum), "run time(hour):", runtime)
print(mod_name, str(Pnum), "run time(day):", runtime/24)


# Storing output
out = xr.Dataset(
    {'Nut': (['time', 'n0s', 'MixRgm'], sol[:, 0, :, :]),
     'Zoo': (['time', 'ZooSize', 'n0s', 'MixRgm'], sol[:, 1:3, :, :]),
     'Phy': (['time', 'PhySize', 'n0s', 'MixRgm'], sol[:, 4:4+Pnum*1, :, :])},
    coords={'time': np.arange(Ynum*365),
            'n0s': N0_list,
            'MixRgm': np.arange(numMixRgm),
            'ZooSize': [5, 200],
            'PhySize': np.logspace(0, 2, Pnum)},
    attrs={'Title': 'This file contains sol output with different nutrient and mixing regimes',
           'No. phytoplankton size group': Pnum,
           'n0s': 'Nutrient regimes = 1, 15, 50',
           'LWST units': 'Degree celcius',
           'PAR units': 'Ein m$^{-2}$ d$^{-1}$',
           'MLD units': 'Meters',
           'Variable units': 'µmol NL$^{-1}$',
           'Size units': 'µm ESD',
           'Phyto size range': '1 - 100',
           'Zoo size range': '5 - 200',
           'Start run time': starttime,
           'Complete run time': endtime,
           'Total run time (hours)': runtime})
out.to_netcdf(wd + mod_name + '.nc')  # export
