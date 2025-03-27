# -*- coding: utf-8 -*-
"""
Created on Wed Feb 14 14:48:57 2024

@author: shang
"""

import pandas as pd
import matplotlib.pyplot as plt
import math
from scipy.interpolate import interp1d

def table_function_series(table_function_dict, s):
    # Erzeugen einer pandas.Series mit Zeit der Sim im Index
    # Stützpunkte durch Dictionary gegeben. Rest Interpoliert.
    tbf = pd.Series(index=s, dtype='float64')
    for key, value in table_function_dict.items():
        tbf[key] = value
    tbf.loc[sim_start:min(table_function_dict.keys())].backfill(inplace=True)
    tbf.loc[max(table_function_dict.keys()):sim_end].ffill(inplace=True)
    tbf.interpolate(method = 'index', inplace=True)
    return tbf

## Base Parameters
sim_start = 0
sim_end = 100
delta_t = 1
sim_length = sim_end - sim_start + delta_t # inklusive start and end
unit = "Years"
iterations = int(sim_length/delta_t)

## TIME-DEPENDENT PARAMETERS
## 7-YEAR CYCLE


## Species for crop rotation
class crop:
    def __init__(self, Y_max, N_uptake_unit, harvest_unit, straw, residue_max, CN_ratio, N_residue, N_fix, fert_org):
        # Max yield measured in OTS [kgOTS/(ha*a)]
        # Specific Nitrogen uptake [kgN/kgOTS]
        # Specific mass of harvested crop [kg/(kgOTS)]
        # Mass of straw and herb(?) [kgOTS/(ha*a)]
        # Max mass of crop residue [kgOTS/(ha*a)]
        # C/N ratio of straw and herb [-]
        # Max mass of nitrogen in crop residue [kgN/(ha*a)]
        # Specific Nitrogen fixation of crops [kgN/kgOTS]
        # Organic fertilizer pp.300-301 [kg/ha]
        self.Y_max = Y_max
        self.N_uptake_unit = N_uptake_unit
        self.harvest_unit = harvest_unit
        self.straw = straw
        self.residue_max = residue_max
        self.CN_ratio = CN_ratio
        self.N_residue = N_residue
        self.N_fix = N_fix
        self.fert_org = fert_org
grain =     crop(6500, 0.029, 1.15, 7700, 1700, 80, 17, 0,      13.3)
maize =     crop(7800, 0.032, 1.15, 8600, 2000, 55, 20, 0,      0) # no fert data
potato =    crop(9000, 0.031, 5.56, 6000, 1000, 30, 25, 0,      44)
beet =      crop(9800, 0.035, 7.14, 6000, 1000, 30, 25, 0,      0) # no fert data
peas =      crop(4400, 0.075, 1.14, 6000, 2300, 15, 53, 6.8E-02,6.7)
rapeseed =  crop(3600, 0.060, 1.11, 6900, 1300, 30, 30, 0,      0) # no fert data
clover =    crop(7000, 0.029, 1.14, 0,    4000, 15, 75, 0.025,  17)
grass =     crop(5200, 0.023, 1.15, 0,    4000, 50, 30, 0,      0) # no fert data


## Constant Simulation Parameters (Umweltdynamik, p302)
soildepth = 0.2         # SOIL [m]
OM_insoil = 0.04        # OGMTR [-]
HU_last_inOM = 0.75     # HFRAC [-]
C_OM_ini = 0.04
Navai_ini = 50          # NV [kg/ha]
CN_ratio = 20           # CN ratio of nutrious humus
NOx = 25                # NAT [kg/(ha*a)]
erosion = 5             # EROS, data from FARM model [t/(ha*a)]
planttime = 15          # GT, plant time of crops [week], line 106
harvesttime =32         # RP, harvest time [week], line 106, different from FEUCHT model
cut1 = 23               # time of the first cut, assumed
WILT = 1                # HUMUS model p.308 line 106
MANUR = 15              # assumed
COMPOST = 10            # assumed
Y_OTS_mk = 0            # mark value of yield, measured by OTS [kg/(ha*a)]
straw_remain = 0.1      # STRAW, share of straw remained in the field [-]

#fert = (clover.fert_org*2+potato.fert_org+grain.fert_org+peas.fert_org+6.7*2)/7
fert = 0

# YMAX NSP YSP SMAX RMAX CSTR MXRES NFIX
Y_max_total = (clover.Y_max * 2 + potato.Y_max + grain.Y_max * 3 + peas.Y_max) / 7
N_uptake_unit_total = (clover.N_uptake_unit * 2 + potato.N_uptake_unit + grain.N_uptake_unit * 3 + peas.N_uptake_unit) / 7
harvest_unit_total = (clover.harvest_unit * 2 + potato.harvest_unit + grain.harvest_unit * 3 + peas.harvest_unit) / 7
straw_total = (clover.straw * 2 + potato.straw + grain.straw * 3 + peas.straw) / 7
residue_max_total = (clover.residue_max * 2 + potato.residue_max + grain.residue_max * 3 + peas.residue_max) / 7
CN_ratio_total = (clover.CN_ratio * 2 + potato.CN_ratio + grain.CN_ratio * 3 + peas.CN_ratio) / 7
N_residue_total = (clover.N_residue * 2 + potato.N_residue + grain.N_residue * 3 + peas.N_residue) / 7
N_fix_total = (clover.N_fix * 2 + potato.N_fix + grain.N_fix * 3 + peas.N_fix) / 7
fert_org_total = (clover.fert_org * 2 + potato.fert_org + grain.fert_org * 3 + peas.fert_org) / 7


## Ordinary Table Functions
# Influence of C/N on decomposition, N/C on DECAYF
NC_x = [0, 0.01, 0.025, 0.033, 0.04, 0.05, 1] #  N/C
decay_y = [0, 0.05, 0.25, 0.5, 0.9, 1.0, 1.0] # DECAYF
decay_tbf = interp1d(NC_x, decay_y,'linear') 
# Influence of C/N on N transfer, N/C on TRANSF
NC_x = [0, 0.01, 0.025, 0.033, 0.04, 0.05, 1] # N/C
Ntrans_y = [1, 0.95, 0.75, 0.5, 0.1, 0, 0] # TRANSF
Ntrans_tbf = interp1d(NC_x, Ntrans_y,'linear') 
# Influence of N on crop growth
NIR_x = [0, 0.2, 0.35, 0.5, 1, 2, 3, 5, 10] # NIR
growth_y = [0, 0.2, 0.5, 0.8, 1, 0.9, 0.4, 0.1, 0] # HUNG
growth_tbf = interp1d(NIR_x, growth_y,'linear') 
# Influence of OM content on N leaching
OM_x = [0, 0.05, 0.1, 1] # OGMTR
Nleach_y = [1, 0.5, 0.2, 0.1] # RC
Nleach_tbf = interp1d(OM_x, Nleach_y,'linear') 


## Create simulation timeline and initialize integrators (stock variables)
sim_series = pd.Series(range(iterations)) # Range is exclusive of end
sim_series = pd.Series(range(iterations), index=sim_series)


## Initialize integrators
# CARB, C in OM [kg/ha]
# NOM, N in nutritious humus [kg/ha]
# NV, N available [kg/ha]
# RBIOM, relative biomass [-]
# CHUM, C in permanent humus [kg/ha]

C_OM = pd.Series(index=sim_series.index, dtype='O')
N_HU_nutri = pd.Series(index=sim_series.index, dtype='O')
Navai = pd.Series(index=sim_series.index, dtype='O')
biomass_rel = pd.Series(index=sim_series.index, dtype='O')
C_HU_perma = pd.Series(index=sim_series.index, dtype='O')

C_OM.iloc[0] = 0                        # assign value in the loop
N_HU_nutri.iloc[0] = 0                  # assign value in the loop
Navai.iloc[0] = Navai_ini
biomass_rel.iloc[0] = 0.01
C_HU_perma.iloc[0] = 0                  # assign value in the loop



## Step through timeline
for i in sim_series.index:
    # 1st stage variables (variable in stage n requires >=1 variable from stage n-1 to calculate)
    HU_nutri_insoil = OM_insoil*(1-HU_last_inOM)            # OMTR=OGMTR*(1-HFRAC) 
    C_OM[i] = HU_nutri_insoil*soildepth*10000*650*0.47      # CARB_1=OMTR*SOIL*10000*650*0.47
    C_HU_perma[i] = C_OM[i]*HU_last_inOM*(1-HU_last_inOM)   # CHUM=CARB_1*HFRAC/(1-HFRAC)
    N_HU_nutri[i] = C_OM[i]/CN_ratio                        # NOM_1=CARB_1/CN
    Navai[i] = Navai[i] + fert + NOx                        # NV=NV+FERT+NAT
    if i == cut1:
        Y_max_total = Y_max_total/2
        print("cut1 is done")
    Y_OTS = biomass_rel[i]*Y_max_total                      # YD_1=RBIOM*YMAX
    erosion_rel = erosion/(soildepth*10000*1600)            # ELOSS=EROS/(SOIL*10000*1600)
    CMANUR=MANUR*0.225*0.47
    NMANUR=CMANUR/20
    CCOMP=COMPOST*.75*.5*.47
    NCOMP=CCOMP/15
    C_OM[i] = C_OM[i]+CMANUR+CCOMP                          # CARB_2 = CARB_1+CMANUR+CCOMP
    N_HU_nutri[i] = N_HU_nutri[i]+NMANUR+NCOMP              # NOM_2=NOM_1+NMANUR+NCOMP
    
    
    # 2nd stage variables
    HU_last_loss = C_HU_perma[i]*(0.02+erosion_rel)                 # HLOSS=CHUM*(0.02+ELOSS)
    Navai_rel = Navai[i]/(Y_max_total*N_uptake_unit_total)          # NIR=NV/(YMAX*NSP)
    biomass = (Y_OTS/Y_max_total)*(Y_max_total+straw_total+residue_max_total) 
                                                                    # BIOM=(YD_1/YMAX)*(YMAX+SMAX+RMAX)
                                                                    
                                                                    
    # 3rd stage variables
    CN_ratio = C_OM[i]/N_HU_nutri[i]                                # CN=CARB_2/NOM_2
    OM_insoil = (C_OM[i]+C_HU_perma[i])/(soildepth*10000*650*0.47)  # OGMTR=(CARB_2+CHUM)/(SOIL*10000*650*0.47)
    HUNG = growth_tbf(Navai_rel)                                    # HUNG_tbf
    if i<planttime or i>harvesttime:
        cropgrowth_rel = 0                                          # RGPLANT
    else:
        cropgrowth_rel = (25*(20/52)/(harvesttime-planttime))*biomass_rel[i]*(1-biomass_rel[i])*WILT*HUNG
                                                                    # RGPLANT=(25*(20/52)/(RP-GT))*RBIOM*(1-RBIOM)*WILT*HUNG
    cropgrowth = cropgrowth_rel*Y_max_total                         # GPLANT=RGPLANT*YMAX
    N_uptake_crop = cropgrowth*(N_uptake_unit_total-N_fix_total)    # RNUP=GPLANT*(NSP-NFIX)
    
    
    # 4th stage variables
    DECAYF = decay_tbf(1/CN_ratio)                                  # DECAYF_tbf
    RDECAY=0.2*DECAYF*C_OM[i]*(1+0.5*math.sin(6.28*(i/100-0.25)))   # RDECAY=0.2*DECAYF*CARB_2*(1+0.5*SIN(6.28*(T-0.25)))
    TRANSF = Ntrans_tbf(1/CN_ratio)                                 # TRANSF_tbf
    N_transfer = TRANSF*Navai[i]                                    # RTRANSF=TRANSF*NV  
    RC = Nleach_tbf(OM_insoil)                                      # RC_tbf
    N_leach = 0.3*Navai[i]*RC                                       # RLEACH=0.3*NV*RC
    
    
    # Updates (set initial values for the next loop)
    HU_last_inOM = C_HU_perma[i]/(C_HU_perma[i]+C_OM[i])                # HFRAC_up=CHUM/(CHUM+CARB_2)
    C_OM[i+1] = C_OM[i]-delta_t*(RDECAY+erosion_rel*C_OM[i])            # CARB_up=CARB_2-DT*(RDECAY+ELOSS*CARB_2)
    C_HU_perma[i+1] = C_HU_perma[i]+delta_t*(RDECAY*0.25-HU_last_loss)  # CHUM_up=CHUM+DT*(RDECAY*0.25-HLOSS)
    biomass_rel[i+1] = biomass_rel[i]+delta_t*cropgrowth_rel            # RBIOM_up=RBIOM+DT*RGPLANT
    if cut1>0.1 and abs(i-cut1)<0.5*delta_t:
        biomass_rel[i] = 0.1                                            # RBIOM_up
        Y_OTS_mk1 = Y_OTS_mk
        Y_OTS_mk = 0.0
    if i>harvesttime:
        biomass_rel[i] = 0.01                                           # RBIOM_up
    Y_OTS = biomass_rel[i]*Y_max_total                                  # YD_2=RBIOM_up*YMAX
    if Y_OTS > Y_OTS_mk:
        Y_OTS_mk = Y_OTS                                                # YL=YD_2
    Navai[i+1] = Navai[i]+delta_t*(RDECAY/CN_ratio-N_transfer-N_uptake_crop-N_leach)
                                                                        # NV_up=NV+DT*(RDECAY/CN-RTRANSF-RNUP-RLEACH)
    N_HU_nutri[i+1] = N_HU_nutri[i]+delta_t*N_transfer-(RDECAY/CN_ratio)-erosion_rel*(C_OM[i]/CN_ratio)
                                                                        # NOM_up=NOM_2+DT*RTRANSF-(RDECAY/CN)-ELOSS*CARB_up/CN)            
    if i>harvesttime:
        C_OM[i+1] = C_OM[i+1]+0.47*(Y_OTS_mk/Y_max_total)*(straw_remain*straw_total+residue_max_total)
                                                                        #CARB_up=CARB_up+0.47*(YL/YMAX)*(STRAW*SMAX+RMAX)
        N_HU_nutri[i+1] = N_HU_nutri[i+1]+0.47*(Y_OTS_mk/Y_max_total)*(straw_remain*straw_total/N_residue_total)+(Y_OTS_mk/Y_max_total)*N_fix_total
                                                                        #NOM_up=NOM_3+0.47*(YL/YMAX)*(STRAW*SMAX/CSTR)+(YL/YMAX)*MXRES

    
## Static visualization
x_axis = sim_series
plt.plot(C_OM.index, C_OM)
plt.xlabel('Year')
plt.ylabel('C_OM') 
plt.title('C_OM')
plt.show()

plt.plot(N_HU_nutri.index, N_HU_nutri)
plt.xlabel('Year')
plt.ylabel('N_HU_nutri') 
plt.title('N_HU_nutri')
plt.show()

plt.plot(Navai.index, Navai)
plt.xlabel('Year')
plt.ylabel('Navai') 
plt.title('Navai')
plt.show()

plt.plot(biomass_rel.index, biomass_rel)
plt.xlabel('Year')
plt.ylabel('biomass_rel') 
plt.title('biomass_rel')
plt.show()

plt.plot(C_HU_perma.index, C_HU_perma)
plt.xlabel('Year')
plt.ylabel('C_HU_perma') 
plt.title('C_HU_perma')
plt.show()


