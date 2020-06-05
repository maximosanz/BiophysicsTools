import numpy as np
from scipy.optimize import curve_fit

class CD_Unfolding():
    
    def __init__(self,T,Ellipticity,T_Units="K",Energy_Units="J"):
        Valid_T_Units = ["K","F","C"]
        Valid_E_Units = ["J","cal","kJ","kcal"]
        
        if T_Units not in Valid_T_Units:
            raise ValueError('Temperature units must be one of {}'.format(Valid_T_Units))
        if Energy_Units not in Valid_E_Units:
            raise ValueError('Energy units must be one of {}'.format(Valid_E_Units))
            
        self.Units = {"T" : T_Units,
                      "E" : Energy_Units}
        
        # Gas constant in SI units
        self.R = 8.31446261815324 # J * K^-1 * mol^-1
        
        self.T = np.array(T)
        self.Ellipticity = np.array(Ellipticity)
        if (self.T.shape != self.Ellipticity.shape
                or np.max([len(self.T.shape),len(self.Ellipticity.shape)]) > 1):
            raise ValueError('Temperature and Ellipticity must be 1-D arrays of the same length')
            
        # INITIAL ESTIMATES FOR EQUATION PARAMETERS:
        
        # Estimates for ThF and ThU
        ThF = np.min(self.Ellipticity)
        ThU = np.max(self.Ellipticity)
        
        # Estimate for Tm
        Fhalf = (ThU-ThF)/2. + ThF
        Tm = self.T[np.argmin(np.absolute(Fhalf-self.Ellipticity))]
        
        # Estimate for dH and dCp
        dH = -10.**5
        dCp = 0.
        
        self.Pars = [Tm,dH,dCp,ThU,ThF]
        
        self.dG = np.zeros(len(T))
        self.Th = np.zeros(len(T))
        
        return None
    
    def Calc_dG(self,T,Tm,dH,dCp,T_Unit="K",Tm_Unit="K"):
        if T_Unit != "K":
            T = self.T_in_K(T,T_Unit)
        if Tm_Unit != "K":
            Tm = self.T_in_K(Tm,Tm_Unit)
        return dH*(1 - T/Tm) - dCp * ((Tm - T) + T*np.log(T/Tm))
        
    def Gibbs_Helmholtz_Eq(self,T,Tm,dH,dCp,ThU,ThF):
        # Parameters:
        # Tm = melting temperature
        # dH = Enthalpy
        # dCp = Heat capacity
        # ThU = Ellipticity of the unfolded state
        # ThF = Ellipticity of the folded state

        T = self.T_in_K(T,self.Units["T"])
        Tm = self.T_in_K(Tm,self.Units["T"])
        
        # Compute deltaG (Gibbs' free energy)
        self.dG = self.Calc_dG(T,Tm,dH,dCp)
        
        K = np.exp(-self.dG/(self.R*T))
        alpha = K / (1. + K)
        
        # Calculate theoretical ellipticity
        Th = alpha*(ThF - ThU)+ThU
        self.Th = Th
        
        return Th
    
    def Fit(self):
        self.Pars = curve_fit(self.Gibbs_Helmholtz_Eq,
                                    self.T,
                                    self.Ellipticity,
                                    p0=self.Pars)[0]
        return self
    
    def Info(self,dG_Temp=None):
        if dG_Temp is None:
            # 25 deg celsius by default
            dG_TempD = {"C" : 25. ,
                        "K" : 25.+273. ,
                        "F" : 77.}
            
            dG_Temp = dG_TempD[self.Units["T"]]
            
        dG = self.Calc_dG(dG_Temp,*self.Pars[:3],T_Unit=self.Units["T"],Tm_Unit=self.Units["T"])
            
        E_Fc = 1.
        if self.Units["E"] in ["cal","kcal"]:
            E_Fc *= 0.239006
        if self.Units["E"][0] == "k":
            E_Fc /= 1000.
        print("Melting temperature = {:7.2f} {}".format(self.Pars[0],self.Units["T"]))
        print("Enthalpy            = {:10.4f} {} / mol".format(self.Pars[1]*E_Fc,self.Units["E"]))
        print("deltaG at {:7.2f} {} = {:10.4f} {} / mol".format(dG_Temp,self.Units["T"],dG*E_Fc,self.Units["E"]))
        return

    def T_in_K(self,T,input_unit):
        if input_unit == "C":
            T = np.copy(T)+273.
        elif input_unit == "F":
            T = (np.copy(T) + 459.67) * 5. / 9.
        return T

