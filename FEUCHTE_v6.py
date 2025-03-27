# -*- coding: utf-8 -*-
"""
Created on Sat Feb 10 12:01:08 2024

@author: shang
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.interpolate import interp1d
import datetime
import random

# from matplotlib.animation import FuncAnimation
# from matplotlib.animation import PillowWriter
# import time


def step_function_series(step_function_dict, index_series):
    stpf = pd.Series(index=index_series, dtype='float64')
    for key, value in step_function_dict.items():
        stpf.loc[key] = value
    stpf.fillna(method='ffill', inplace=True)
    return stpf

## Base Parameters
# Case 1 in pp. 283
sim_start = 0
sim_end = 365
delta_t = 1

# Create simulation timeline
sim_index = pd.Index(np.arange(sim_start, sim_end, delta_t, dtype='float').round(2)) # float only ncecessary for delta_t < 1


## CODE BLOCK A
## Constant Simulation Parameters
watercap_unit = 0.15    # UCGD, capacity of groundwater storage [mWL/m]
soildepth = 0.3         # SOIL, depth of the plow layer [m]
OM = 0.08               # OGMTR, organic matter in soil in volume fraction [-]
clay = 0.3              # CLAY, clay in soil in volume fraction [-]

drought = 0             # DROUGHT, length of dry period [a]
harvest_loss = 0        # CX, harvest loss happens or not [0/1]
compression = 0         # CQ, soil compression happens or not [0/1]

planttime_week = 15     # GT, vegetation initial time (week)
planttime = int(planttime_week*7)
harvesttime_week = 30   # RP, harvest time (week)
harvesttime = int(harvesttime_week*7)

rootdepth = 1           # ROOT, [m]
biomass_max = 10000     # biomass in ODM [kg/ha]
fert_factor = 1          # FFAC, fertilization factor, 1-optimum fertilization, 0-no fert
tran_rate = 0.4         # Specific transpiration rate [m/(kg/ha)]
anevap = 0.444          # Annual evaporation [mWL/a]
rain_annual = 0
pump_rate = 150         #System Zoo 2 p64. pump_rate = 150 [1/year]
percolation_rate = 50   #System Zoo 2 p65. percolation_rate = 50 [1/year]

ground_covering = 1     # CVR, types of ground covering (0-no covering, 1-plants, 2-foil)
watering = 0            # IP, whether irrigation happens (1-yes, 2-no)


## CODE BLOCK B
# Ordinary Table Functions
# Influence of clay content[-] on groundwater level capacity[mWL/m], USFCAP=f(CLAY)
water_capacity = {0:0.01,   0.1:0.1,    0.2:0.16,   0.3:0.2, 
                  0.4:0.21, 0.5:0.15,   0.6:0.13,   1.0:0.05}
watercap_tbf = interp1d(list(water_capacity.keys()), 
                        list(water_capacity.values()), 
                        'linear')                           # USFCAP    
# Influence of clay content[-] on capillary rise of groundwater level[m], CPRISE=f(CLAY)
capillary_rise_from_clay_content ={0:0.1, 1:3}
capillary_rise_from_clay_content_tbf = interp1d(list(capillary_rise_from_clay_content.keys()), 
                                                list(capillary_rise_from_clay_content.values()), 
                                                'linear')   # CPRISE
# Influenc of moisture[-] on plant growth[-]
plant_growth_from_moisture = {0:0, 0.1:0.2, 0.5:1, 1:1, 2:1}
plant_growth_from_moisture_tbf = interp1d(list(plant_growth_from_moisture.keys()),
                                          list(plant_growth_from_moisture.values()), 
                                          'linear')         # WILT


## CODE BLOCK C
## Time-dependent Parameters (set-up as dict with key=time. Value=parametervalue)
# Daily precipitaion [mmWL/day]
rain = {0:2.21,     32:2.45,    60:1.52,    91:2.03, 
        121:2.03,   152:2.50,   182:2.77,   213:2.90, 
        244:2.20,   274:2.16,   305:2.28,   335:2.21}
rain = {key: value/1000*365 for key, value in rain.items()}
growth_time = {0:0, planttime:1, harvesttime+1:0}
## Set (time-dependent) parameters to timeline
rain = step_function_series(rain, sim_index)
growth_period = step_function_series(growth_time, sim_index)


## CODE BLOCK D
## Initialize integrators
moisture = pd.Series(index=sim_index, dtype='float64')
biomass_rel = pd.Series(index=sim_index, dtype='float64')
biomass = pd.Series(index=sim_index, dtype='float64')
water = pd.Series(index=sim_index, dtype='float64')
watertable = pd.Series(index=sim_index, dtype='float64')
irrigation_cumu = pd.Series(index=sim_index, dtype='float64')

moisture_j = 0.5          # INMOIST, initial moisture of soil [-]
biomass_rel_j = 0.01      # RBIOM, relative biomass [-]        
biomass_j = 100.0         # BIOM, biomass [kg/ha]
water_j = 0.0             # WATER, available groundwater [mWL], initial value assigned in the loop
watertable_j = -10.0      # WTBL, watertable [m], 10m underground
irrigation_cumu_j = 0.0   # IWTR, cumulative irrigation [mWL]


for i in sim_index: # Vorletzter Eintrag befüllt letzten Eintrag
    moisture[i] = moisture_j                    # INMOIST, initial moisture of soil [-]
    biomass_rel[i] = biomass_rel_j              # RBIOM, relative biomass [-]        
    biomass[i] = biomass_j                      # BIOM, biomass [kg/ha]
    water[i] = water_j                          # WATER, available groundwater [mWL], initial value assigned in the loop
    watertable[i] = watertable_j                # WTBL, watertable [m]
    irrigation_cumu[i] = irrigation_cumu_j

# 1. Waterholding capacity
## Lines 420-620
    # Waterholding capacity defined by organic matter [mWL]
    watercap_OM = soildepth * OM * 5                # soildepth * fraction of organic matter in soil * 5
    # Waterholding capacity corrected by soil compression
    watercap_soil = watercap_tbf(clay) * (0.94 if clay < 0.15 else 0.88) if compression != 0 else watercap_tbf(clay)
    # Waterholding capacity corrected by root depth and capillary rise
    watercap_root_CR = rootdepth + capillary_rise_from_clay_content_tbf(clay) # depth of attainable water level [mWL]
    # Waterholding capacity_final
    watercap_i = watercap_root_CR * watercap_soil if -watercap_root_CR > watertable[i] else watercap_OM + rootdepth * watercap_soil
    
    # Initialize soil water in the first loop
    water[i] = moisture[i] * watercap_i               # WATER=MOIST*WCAP

# 2. Precipitation OR Rain
    # later use rain[i] as a randomized version
    # e.g. rain = rain[i] * (0.5 + random.random()) /1000
    
# 3. Evaporation
    evaporation_i = moisture[i]*anevap*(1+0.95*math.sin(2*3.14*(i/365-0.25)))*(1-ground_covering/2)
                                                    # System Zoo 2 p64. evaporation

    
# 4. RGPLANT, Relative growth rate of plants [1/a]
    growthrate_rel_i = (25*(20/52)/(harvesttime_week-planttime_week)) * biomass_rel[i] * (1-biomass_rel[i])*plant_growth_from_moisture_tbf(moisture[i])
                                                                    # RGPLANT=(25*(20/52)/(RP-GT))*RBIOM*(1-RBIOM)*WILT
    growthrate_rel_i = growthrate_rel_i * growth_period[i]
                                                                    # growthrate_rel_factor, time-dependant factor
    growth_i = growthrate_rel_i*biomass_max                         # GPLANT=RGPLANT*BMAX
    
# 5. TSPIR, Transpiration [mWL/a]
    transpiration_i = (1/fert_factor)*tran_rate*growth_i/(1E+05)     # TSPIR= (1/FFAC)*400*GPLANT/1E+07
    
    
# 6. Updates
# 6.1 Biomass_relative, Biomass update
    # drought = 0                 if moisture[i]>0.3 else drought
    # drought = drought+delta_t   if moisture[i]<0.1 else drought
    # harvest_loss_i = 0
    # harvest_loss_i = growth_period[i] if drought>20/365 else harvest_loss_i
    # if harvest_loss_i == 1:
    #     biomass_rel_j = 0
    # else:
    #     if i > harvesttime:
    #         biomass_rel_j = 0                    
    #     else:
    #         biomass_rel_j = biomass_rel[i] + delta_t*growthrate_rel_i
    #                                                              # RBIOM=RBIOM+DT*RGPLANT
    biomass_rel_j = biomass_rel[i] + delta_t*growthrate_rel_i if i<harvesttime else 0   # assume no drought
    
    biomass_j = biomass_rel_j*biomass_max                        # BIOM=RBIOM*BMAX
         
# 6.2 Irrigation
    rain_annual = rain_annual + rain[i]*delta_t         # ANRAIN=ANRAIN+PRECIP*DT
    irrigation_i = 0.0                                           # IRRIGN=0 (irrigation amount at certain time [mWL])
    irrigation_i = 0.5 * watercap_i - water[i] if (watering == 1 and moisture[i] < 0.5 and i>planttime and i < harvesttime) else 0
    #
    irrigation_i = irrigation_i*pump_rate               # System Zoo 2 p64. Pump rate = 150 [1/year]; irrigation_i [mWL/year]
    #
    irrigation_cumu_j = irrigation_cumu[i] + irrigation_i        # IWTR=IWTR+IRRIGN
# 6.3 Watertable
    # down = irrigation_i/watercap_unit                            # DOWN=IRRIGN/UCGD, drop of water level [m/year]
    # watertable_j = watertable[i] + down                          # WTBL=WTBL+DOWN
# 6.4 Water, Moisture update
    water_excess = water_j - watercap_i                            # WXCESS=WATER-WCAP
    watertable_rise = water_excess/watercap_unit*percolation_rate if water_excess>0 else 0
    water_change = rain[i] - evaporation_i - transpiration_i + irrigation_i - watertable_rise
    water_j = water[i] + water_change                          # WATER=WATER+DT*(PRECIP-EVAP-TSPIR)
    water_j = watercap_i                        if water_j>watercap_i else water_j
    
    watertable_down = irrigation_i/watercap_unit
    watertable_change = watertable_rise - watertable_down
    waterttable_j = watertable[i] + watertable_change

    moisture_j = water_j/watercap_i                                # MOIST=WATER/WCAP



## Output results
x_axis = sim_index
# Biomass [kg/ha]
plt.plot(biomass.index, biomass)
plt.xlabel('Year')
plt.ylabel('Biomass [kg/ha]') 
plt.title('Model FEUCHTE - Biomass')
plt.show()

# Pricipitation [mm/a] (& Irrigation [mm/a])
plt.plot(rain.index, rain)
plt.xlabel('Year')
plt.ylabel('Precipitation [mm/a]') 
plt.title('Model FEUCHTE - Precipitation')
plt.show()

for i in irrigation_cumu.index:
    # Check if any value in the column is larger than zero, print diagram if so.
    if (irrigation_cumu[i] > 0).any():
        plt.plot(irrigation_cumu.index, irrigation_cumu)
        plt.xlabel('Year%')
        plt.ylabel('Irrigation [mm/a]') 
        plt.title('Model FEUCHTE - Irrigation')
        plt.show()
        break

# Available Groundwater [mWL]
plt.plot(water.index, water)
plt.xlabel('Year')
plt.ylabel('Available Groundwater [mWL]') 
plt.title('Model FEUCHTE - Available Groundwater')
plt.show()

# Groundwater Level [mWL]
plt.plot(watertable.index, watertable)
plt.xlabel('Year')
plt.ylabel('Groundwater Table [mWL]') 
plt.title('Model FEUCHTE - Groundwater Table')
plt.show()