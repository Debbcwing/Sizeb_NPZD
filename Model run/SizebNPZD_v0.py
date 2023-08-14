import numpy as np
from scipy.integrate import odeint


class SizeModLogspace:
    """
    Size-based model under the assumption of allometric relationships on nutrient
    uptake and grazing.
    -------
    Model has multiple groups of phytoplankton, 2 groups of zooplankton in different sizes,
    one nutrient source and a detritus state variable.
    -------
    Flexible parameters:
        - n0
        - prey tolerance
        - P size range
        - Z size range
        - no. of P size classes
        - no. of years for model running

    Fixed parameters:
        - Allometric parameters for Kn
        - Allometric parameters for Mu
        - Allometric parameters for Imax
        - Allometric parameters for grazing size preference
        - Model params
    """

    def __init__(self, mld, par, sst, dmdt, minP=0., maxP=2., numP=3, minZ=5, maxZ=200, numZ=2, numYears=5,
                 N0=15., Zsmall0=0.01, Zlarge0=0.01, D0=0.01,
                 Pmax=1.1, Kpar=0.1, alpha=0.15, Mz=0.1, Mp=0.2, Mz2=0.34,
                 Ieff=0.75, beta=0.69, Kp=3, Rem=0.6, k=0.1,
                 mu_alpha=-0.36, mu_beta=10 ** 0.69, Kn_alpha=0.52, Kn_beta=10 ** -0.71,
                 graz_alpha=-0.4, graz_beta=26, ops_alpha=0.56, ops_beta=0.65,
                 SizeTol=[0.2, 0.5]):
        self.mld = mld
        self.par = par  # photosynethically activated radiation (E*m^-2*d^-1)
        self.sst = sst  # sea surface temperature/lake surface water temperature (LSWT)
        self.dmdt = dmdt  # dmdt
        self.minP = minP  # min. size of P size range
        self.maxP = maxP  # max. size of P size range
        self.numP = numP  # number of phytoplankton size class
        self.minZ = minZ  # min. size of Z size range
        self.maxZ = maxZ  # max. size of Z size range
        self.numZ = numZ  # number of zooplankton size class
        self.numYears = numYears  # Number of years for model run
        self.N0 = N0  # Initial nutrient concentrations (µmol N L^-1)
        self.Zsmall0 = Zsmall0  # Initial Zsmall concentration (µmol N L^-1)
        self.Zlarge0 = Zlarge0  # Initial Zlarge concentration (µmol N L^-1)
        self.D0 = D0  # Initial detritus concentration (µmol N L^-1)
        self.Ps0 = np.repeat(0.01, self.numP, axis=0)  # Initial P concentration (µmol N L^-1)
        self.Pmax = Pmax  # Max photosynthesis rate
        self.Kpar = Kpar  # Water attenuation coefficient
        self.alpha = alpha  # Initial slope of P-I curve
        self.Mz = Mz  # Natural mortality - Zooplankton
        self.Mp = Mp  # Natural mortality - Phytoplankton
        self.Mz2 = Mz2  # Mortality from upper predators
        self.Ieff = Ieff  # Ingestion efficiency
        self.beta = beta  # Assimlation efficiency
        self.Kp = Kp  # Phytoplankton half-saturation constant (µmol N L^-1)
        self.Rem = Rem  # Remineralization rate (day^-1)
        self.k = k  # Cross thermocline mixing factor (m/day)
        self.mu_alpha = mu_alpha  # Slope of allometric max. uptk rate
        self.mu_beta = mu_beta  # Intercept of allometric max. uptk rate
        self.Kn_alpha = Kn_alpha  # Slope of allometric half-saturation constant
        self.Kn_beta = Kn_beta  # Intercept of allometric half-saturation constant
        self.graz_alpha = graz_alpha  # Slope of allometric max. grazing rate (day^-1)
        self.graz_beta = graz_beta  # Intercept of allometric max.grazing rate (day^-1)
        self.ops_alpha = ops_alpha  # Slope of allometric optimum prey size
        self.ops_beta = ops_beta  # Intercept of allometric optimum prey size
        self.SizeTol = SizeTol  # Prey size tolerance (should be an array)
        self.init_con = np.concatenate(([self.N0],  # N0
                                        [self.Zsmall0],  # Zsmall0
                                        [self.Zlarge0],  # Zlarge0
                                        [self.D0],  # D0
                                        [self.Ps0]), axis=None)  # Ps0
        self.t_array = np.arange(0, 365 * self.numYears)  # Time array
        self.Psize_array = np.logspace(self.minP, self.maxP, num=self.numP, base=10)  # Size array for phytoplankton
        self.Zsize_array = np.linspace(self.minZ, self.maxZ, self.numZ)  # Size array for zooplankton
        self.solution = odeint(self.AlloNPsZD, self.init_con, self.t_array, rtol=1.e-5, atol=1.e-5, hmax=1.e16,
                               hmin=1.e-8)

    def Vol2ESD(self, v):  # Vol2ESD(1) = 1.2407
        esd = ((v * 6) / np.pi) ** (1 / 3)
        return esd

    def Allo_mumax(self, PhySize):
        return self.mu_beta * ((PhySize / self.Vol2ESD(1)) ** self.mu_alpha)

    def Allo_NU(self, nut, PhySize):
        Kn_scal = self.Kn_beta * (
                (PhySize / self.Vol2ESD(1)) ** self.Kn_alpha)  # allometric function for half-saturation
        return nut / (nut + Kn_scal)

    def temp(self, time):
        """
        Eppley's formulation
        """
        return np.exp(0.063 * time)

    def LightLim(self, I, Z):
        """
        Light limitation function on photosynthesis based on Smith's function,
        illustrated by a P-I curve
        """
        I0 = I
        Iz = I0 * np.exp(-self.Kpar * Z)
        LL = (((self.Pmax / (self.Kpar * Z)) * np.log((self.alpha * I0 + np.sqrt((self.Pmax ** 2) +
                                                                                 (self.alpha * I0) ** 2)) / (
                                                              (self.alpha * Iz) + np.sqrt(self.Pmax ** 2 + (
                                                               self.alpha * Iz) ** 2)))))
        return LL

    def AlloImax(self, ZooSize):
        """
        Max. ingestion rate also depends on the size of the zooplankton
        ----------
        graz_alpha = -0.4  # Hansen et al. 1994, Banas, 2011
        graz_beta = 26  # Hansen et al. 1994, Banas, 2011 [day^-1]
        """
        return self.graz_beta * ((ZooSize / 1) ** self.graz_alpha)

    def AllogpP(self, PhySize, ZooSize, SizeTol):
        """
        Selected zooplankton size class of 10 & 200[µm]
        """
        OPS = self.ops_beta * ((ZooSize / 1) ** self.ops_alpha)  # µm ESD
        gpP = np.exp(-((np.log10(PhySize) - np.log10(OPS)) / SizeTol) ** 2)
        return gpP

    def AlloGraz(self, PhySize, ZooSize, Pi, SizeTol, denom):
        """
        PhySize = Size class/array of phytoplankton population;
        ZooSize = Size class/array of zooplankton population;
        Pi 	 	= Phytoplankton biomass
        """
        graz = self.AlloImax(ZooSize) * ((self.AllogpP(PhySize, ZooSize, SizeTol) * Pi) / (self.Kp + denom))
        return graz

    def mixing(self, time):
        """k refers to the cross thermocline mixing factor"""
        return (self.k + max(self.dmdt[time], 0.)) / self.mld[time]

    def AlloNPsZD(self, y, t):
        """
        A size-based trophic model with multiple phytoplankton and zooplankton size class
        ----------
        y: array of initial conditions of state variables (var)
        t: time array
        """
        # Initialization
        N = y[0]
        Zsmall = y[1]
        Zlarge = y[2]
        D = y[3]
        Ps = y[4:]
        var = np.zeros(len(y))
        t = int(t)

        # non-P dependent equations
        """
        Zsmall = specialist, narrower size tolerance factor;
        Zlarge = generalist, wider size tolerance factor.
        """
        Graz_denomZsmall = np.sum(self.AllogpP(self.Psize_array, self.Zsize_array[0], self.SizeTol[0]) * Ps)
        Graz_denomZlarge = np.sum(self.AllogpP(self.Psize_array, self.Zsize_array[1], self.SizeTol[1]) * Ps)

        # [dNdt]
        var[0] = (self.Rem * D  # Remineration
                  + self.mixing(t) * (self.N0 - N))  # MixingN

        # [dZ1dt]
        var[1] = (- self.Mz * Zsmall  # Natural mortality Z1
                  - self.Mz2 * (Zsmall ** 2)  # density-dependent top predator grazing Z1
                  - (self.dmdt[t] / self.mld[t]) * Zsmall)  # Mixing loss Z1

        # [dZ2dt]
        var[2] = (- self.Mz * Zlarge  # Natural mortality Z2
                  - self.Mz2 * (Zlarge ** 2)  # Top predator grazing Z2
                  - (self.dmdt[t] / self.mld[t]) * Zlarge)  # Mixing loss Z2

        # [dDdt]
        var[3] = (self.Mz * (Zsmall + Zlarge)  # MortalityZsmall+Zlarge
                  - self.Rem * D  # Remineration
                  - self.mixing(t) * D)  # Mixing loss D

        # P-dependent (loop)
        for i in range(len(Ps)):
            """
            Psize = Size class/array of phytoplankton population;
            """
            P = Ps[i]  # define P in the loop
            Psize = self.Psize_array[i]

            # [dNdt]
            var[0] = (var[0]
                      - self.Allo_mumax(Psize) * (self.temp(self.sst[t]) * self.Allo_NU(N, Psize) *
                                                  self.LightLim(self.par[t], self.mld[t])) * P
                      + (self.beta * (1. - self.Ieff) *
                         self.AlloGraz(Psize, self.Zsize_array[0], P, self.SizeTol[0], Graz_denomZsmall)) * Zsmall
                      # Excretion from grazing(Zsmall))
                      + (self.beta * (1. - self.Ieff) *
                         self.AlloGraz(Psize, self.Zsize_array[1], P, self.SizeTol[1], Graz_denomZlarge)) * Zlarge)
            # Excretion from grazing(Zlarge))

            # [dZ1dt]
            var[1] = (var[1]
                      + (self.beta * self.Ieff * self.AlloGraz(Psize, self.Zsize_array[0], P, self.SizeTol[0],
                                                               Graz_denomZsmall)) * Zsmall)  # GrazingP from Zsmall

            # [dZ2dt]
            var[2] = (var[2]
                      + (self.beta * self.Ieff * self.AlloGraz(Psize, self.Zsize_array[1], P, self.SizeTol[1],
                                                               Graz_denomZlarge)) * Zlarge)  # GrazingP from Zlarge

            # [dDdt]
            var[3] = (var[3]
                      + ((1. - self.beta) * self.AlloGraz(Psize, self.Zsize_array[0], P, self.SizeTol[0],
                                                          Graz_denomZsmall)) * Zsmall  # SloppyfeedP - Zsmall
                      + ((1. - self.beta) * self.AlloGraz(Psize, self.Zsize_array[1], P, self.SizeTol[1],
                                                          Graz_denomZlarge)) * Zlarge  # SloppyfeedP - Zlarge
                      + self.Mp * P)

            # [dPsdt]
            var[4 + i] = (var[4 + i]
                          + self.Allo_mumax(Psize) * (self.temp(self.sst[t]) * self.Allo_NU(N, Psize)
                                                      * self.LightLim(self.par[t], self.mld[t])) * P  # Uptk
                          - self.Mp * P  # MortalityP
                          - self.AlloGraz(Psize, self.Zsize_array[0], P, self.SizeTol[0],
                                          Graz_denomZsmall) * Zsmall  # GrazingP from Zsmall
                          - self.AlloGraz(Psize, self.Zsize_array[1], P, self.SizeTol[1],
                                          Graz_denomZlarge) * Zlarge  # GrazingP from Zlarge
                          - self.mixing(t) * P)  # Mixing

        return var
