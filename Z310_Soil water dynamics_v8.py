# -*- coding: utf-8 -*-
"""
Created on Wed Apr 17 11:07:59 2024

@author: shang
"""

import pandas as pd
import numpy as np
import math
from scipy.interpolate import interp1d
import random
import matplotlib.pyplot as plt


def step_function_series(step_function_dict, index_series):
    stpf = pd.Series(index=index_series, dtype='float64')
    for key, value in step_function_dict.items():
        stpf.loc[key] = value
    stpf.fillna(method='ffill', inplace=True)
    return stpf

def plot_diagram(data, title, ylim=None):
    plt.plot(data.index, data)
    plt.xlabel('1 Year')
    plt.ylabel(title.split(' - ')[-1])  # Extracting the y-label from the title
    plt.title(title)
    if ylim:
        plt.ylim(ylim)
    plt.show()

# Base Parameters
sim_start = 0
sim_end = 1
delta_t = 0.01      # TIME STEP
sim_index = pd.Index(np.arange(sim_start, sim_end, delta_t, dtype='float').round(2))

# Simulation Parameters
# Soil Parameters
DEPTH_OF_PLOW_LAYER = 0.3       # [m]
ORGANIC_MATTER_IN_SOIL = 0.08   # [1] wolumn fraction
CLAY_FRACTION_OF_SOIL = 0.3     # [m] column fraction
FIELD_CAPACITY_GROUNDWATER_AQUIFER = 0.15   # [mW/mSoil] m water per m soil
INITIAL_SOIL_MOISTURE = 0.5     # [1]
PERCOLATION_RATE = 50           # [1/Year] ********************Meaning?*********************
SOIL_COVER = 1                  # [1] 0: bare soil, 1: plant cover, 2: plastic foil cover
SOIL_COMPACTION_PARAMETER = 0   # [1] compaction occurs=1, doesn't occur=0
# Influence of clay fraction[1] on groundwater level capacity[mW/mSoil]
usable_field_capacity_of_soil = {0.0: 0.01, 0.1: 0.1,   0.2: 0.16,  0.3: 0.2, 
                                 0.4: 0.21, 0.5: 0.15,  0.6: 0.13,  1.0: 0.05}
usable_field_capacity_of_soil_tbf = interp1d(list(usable_field_capacity_of_soil.keys()), 
                                             list(usable_field_capacity_of_soil.values()), 'linear')
# Influence of clay fraction[1] on capillary rise of groundwater level[mW]
capillary_rise = {0: 0.1, 1: 3}
capillary_rise_tbf = interp1d(list(capillary_rise.keys()),
                              list(capillary_rise.values()), 'linear')
# Influenc of soil moisture[1] on plant growth[1]
moisture_effect = {0: 0, 0.1: 0.2, 0.5: 1, 1: 1, 2: 1}
moisture_effect_tbf = interp1d(list(moisture_effect.keys()),
                               list(moisture_effect.values()), 'linear')

# Crop Parameter
MAX_BIOMASS = 10000                 # [kg/ha] biomass in kgODM (organic dry matter) per hactare
SPECIFIC_TRANSPIRATION_RATE = 0.4   # [mW/(kg/ha)] m water transpired per biomass (biomass density in [kg/ha])
ROOT_DEPTH_OF_FIELD_CROP = 1        # [m]
BEGIN_GROWTH_PERIOD = 0.3           # 
HARVEST_DATE = 0.6                  #
FERTILIZATION_FACTOR = 1            # 1: optimal fertilization; 0: no fert, intermediates are valid for calculation
PUMP_RATE = 150                     # [1/Year] no irrigation=0 ********************Meaning?*********************

# Weather parameters
# Monthly average rainfall [mW/Month]
MONTHLY_RAINFALL  = {0: 0.068, 1: 0.068, 2: 0.068, 
                     3: 0.047, 4: 0.061, 5: 0.063, 
                     6: 0.075, 7: 0.086, 8: 0.09, 
                     9: 0.066, 10: 0.067, 11: 0.068, 12: 0.068}
MONTHLY_RAINFALL  = {key: value/1 for key, value in MONTHLY_RAINFALL.items()}
MONTHLY_RAINFALL  = step_function_series(MONTHLY_RAINFALL , sim_index)
MULTIPLIER_FOR_RAINFALL = 1         # 1: optimal rain, <1: less than optimum
ANNUAL_EVAPORATION = 0.444          # [mW/Year]
WEATHER_PERIOD = 3                  # [Day]

# Other parameters
growth_time = {0: 0, BEGIN_GROWTH_PERIOD: 1, HARVEST_DATE+delta_t: 0}
growth_period = step_function_series(growth_time, sim_index)
MONTHS_PER_YEAR = 12
PI = 3.14159
DELAY_TIME = 0.25/1                 # [Year]

# Initialize stock variables (integrators)
weather = pd.Series(index=sim_index, dtype='float64')
rain_amount = pd.Series(index=sim_index, dtype='float64')
precipitation = pd.Series(index=sim_index, dtype='float64')
soil_moisture = pd.Series(index=sim_index, dtype='float64')
groundwater_level = pd.Series(index=sim_index, dtype='float64')
relative_plant_biomass = pd.Series(index=sim_index, dtype='float64')
drought_period = pd.Series(index=sim_index, dtype='float64')
biomass =pd.Series(index=sim_index, dtype='float64')
soilwater =pd.Series(index=sim_index, dtype='float64')
irrigation_amount =pd.Series(index=sim_index, dtype='float64')

# Assign initial value to stock variables
weather_j = 0                   # [1]       Initial value obtained in the loop
rain_amount_j = 0               # [mW]      Accumulative precipitation amount, start from 0
precipitation_j = 0             # [mW/Year] Initial value obtained in the loop
soil_moisture_j = INITIAL_SOIL_MOISTURE     # [1]
groundwater_level_j = -10       # [mW]      negative for underground
relative_plant_biomass_j = 0.01 # [1]
drought_period_j = 0            # [1]       Initial value obtained in the loop
biomass_j = 0                   # [kg/ha]   Initial value obtained in the loop
soilwater_j = 0                 # [mW]      Initial value obtained in the loop, attainable_water_capacity required first
irrigation_amount_j = 0         # [mW]      Accumulative irrigation amount, start from 0


for i in sim_index:
# Update values for stock variables
    # value used in this iteration=[i], value updated for the next iteration=_j 
    # Uppercase=constants, Lowercase=intermediates
    weather[i] = weather_j
    rain_amount[i] = rain_amount_j
    soil_moisture[i] = soil_moisture_j
    groundwater_level[i] = groundwater_level_j
    relative_plant_biomass[i] = relative_plant_biomass_j
    drought_period[i] = drought_period_j
    biomass[i] = biomass_j
    soilwater[i] = soilwater_j
    irrigation_amount[i] = irrigation_amount_j
# Precipitation dynamics
    # Notice: set_to_zero[n+1] = new_value[n]
    # The difference between n_v and s_t_z is added to w_j, 
    #   which creates a random factor that, when multiplied, 
    #   induces fluctuations around the average value,
    #   while the area under the "rain curve" doesn't deviate far from the ave.
    new_value  = int(0.5 +random.random())/delta_t if abs(i*(365/WEATHER_PERIOD)-int(i*(365/WEATHER_PERIOD))) < delta_t*(365/WEATHER_PERIOD) else 0                # [1/Year]
    set_to_zero = weather[i]/delta_t if abs((i*delta_t/2)*(365/WEATHER_PERIOD)-int((i*delta_t/2)*(365/WEATHER_PERIOD))) < delta_t*(365/WEATHER_PERIOD) else 0   # [1/Year]
    weather_j = weather[i] + delta_t*(new_value - set_to_zero)          # [1]
    # M_F_R<1 if rain is less than optimal, otherwise M_F_R=1
    modified_rainfall = MONTHLY_RAINFALL[i] * MULTIPLIER_FOR_RAINFALL   # [mW/Month]
    # Why times 2?
    rain = MONTHS_PER_YEAR*modified_rainfall*weather_j*2                # [mW/Year]
    rain_amount_j = rain_amount[i] + delta_t*rain                       # [mW]
    precipitation[i] = rain                                             # [mW/Year]
# Waterholding capacity
    # Waterholding capacity defined by organic matter [mW]
    waterholding_capacity_of_organic_matter = DEPTH_OF_PLOW_LAYER * ORGANIC_MATTER_IN_SOIL * 5
    # Waterholding capacity corrected by soil compression [mW]
    reduced_field_capacity_from_compaction = usable_field_capacity_of_soil_tbf(CLAY_FRACTION_OF_SOIL) * (0.94 if CLAY_FRACTION_OF_SOIL<0.15 else 0.88) if SOIL_COMPACTION_PARAMETER>0 else usable_field_capacity_of_soil_tbf(CLAY_FRACTION_OF_SOIL)
    # Waterholding capacity corrected by root depth and capillary rise [mW]
    depth_of_accessible_water_horizon = ROOT_DEPTH_OF_FIELD_CROP + capillary_rise_tbf(CLAY_FRACTION_OF_SOIL)
    # Waterholding capacity_final [mW]
    groundwater_access = 1 if -depth_of_accessible_water_horizon < groundwater_level[i] else 0
    attainable_water_capacity = depth_of_accessible_water_horizon * reduced_field_capacity_from_compaction if groundwater_access>0 else waterholding_capacity_of_organic_matter + ROOT_DEPTH_OF_FIELD_CROP*reduced_field_capacity_from_compaction
# Evaporation [mW/Year]
    evaporation = soil_moisture[i]*ANNUAL_EVAPORATION*(1+0.95*math.sin(2*PI*(i/365-DELAY_TIME)))*(1-SOIL_COVER/2)
# Plant growth rate
    # Relative growth definitive function (Notice: tend to reach max when biomass_rel = 0.5) [1/Year]
    relative_growth = (25*(20/52)/(HARVEST_DATE -BEGIN_GROWTH_PERIOD)) * relative_plant_biomass[i] * (1-relative_plant_biomass[i]) * moisture_effect_tbf(soil_moisture[i])
    # Relative growth corrected by growth period, a time window for the process to happen [1/Year]
    relative_growth = relative_growth * growth_period[i]
    # Growth rate_final [kg/(ha*Year)]
    growth_rate = relative_growth * MAX_BIOMASS
# Transpiration [mW/Year]
    transpiration = (1/FERTILIZATION_FACTOR) * SPECIFIC_TRANSPIRATION_RATE*growth_rate/10000
# Biomass [kg/ha]
    # Drought happens if soil moisture is too low.
    drought = 0     if soil_moisture[i]>0.2 else 1                              # [1/Year]
    drought_end = drought_period[i]/delta_t     if soil_moisture[i]>0.2 else 0  # [1/Year]
    drought_period_j = drought_period[i] + delta_t*(drought-drought_end)        # [1]
    # If drought period is too long, or if harvesting occurs, then withering starts.
    #   Withering leads to sharp decrease of biomass. 
    harvest_loss = 1    if drought_period[i]>(20/365) and i>BEGIN_GROWTH_PERIOD and i<HARVEST_DATE else 0   # [1]
    withering = relative_plant_biomass[i]/delta_t      if i>HARVEST_DATE or harvest_loss == 1 else 0        # [1/Year]
    # Biomass calculated from biomass_rel
    relative_plant_biomass_j = relative_plant_biomass[i] + delta_t*(relative_growth-withering)              # [1]
    # biomass_final [kg/ha]
    biomass_j = relative_plant_biomass[i] * MAX_BIOMASS
# Soil moisture [1]
    # Assigning initial value to soilwater [mW (only as a boost in the first iteration)
    soilwater[i] = INITIAL_SOIL_MOISTURE*attainable_water_capacity  if i<0.01 else soilwater[i]    
    soil_moisture_j = soilwater[i]/attainable_water_capacity        # [1]
    soilwater_surplus = soilwater[i] - attainable_water_capacity    # [mW]
# Irrigation [mW/Year]
    irrigation = (0.5*attainable_water_capacity - soilwater[i])*PUMP_RATE if (soil_moisture[i]<0.5 and i>BEGIN_GROWTH_PERIOD and i<HARVEST_DATE) else 0
    irrigation_amount = irrigation_amount + delta_t*irrigation      # [mW]
# Groundwater level [mW]
    drop_of_groundwater_level = irrigation/FIELD_CAPACITY_GROUNDWATER_AQUIFER
    rise_of_groundwater_level = soilwater_surplus/FIELD_CAPACITY_GROUNDWATER_AQUIFER*PERCOLATION_RATE   if soilwater_surplus>0 else 0
    percolation = rise_of_groundwater_level
    groundwater_level_j = groundwater_level[i] + delta_t*(rise_of_groundwater_level-drop_of_groundwater_level)
    soilwater_j = soilwater[i] + delta_t*(irrigation +precipitation[i] -transpiration -evaporation -percolation)

## Output results
# x_axis = sim_index
plot_diagram(biomass, 'Model Soil Water Dynamics - Biomass [kg/ha]', (0, 20000))
plot_diagram(precipitation, 'Model Soil Water Dynamics - Precipitation [mW/Year]', (0, 2.5))
for i in irrigation_amount.index:
    if (irrigation_amount[i] > 0).any():
        plot_diagram(irrigation_amount, 'Model Soil Water Dynamics - Irrigation Amount [mW/Year]', (0, 10))
        break
plot_diagram(soilwater, 'Model Soil Water Dynamics - Soilwater [mW]', ylim=(0, 0.5))
plot_diagram(groundwater_level, 'Model Soil Water Dynamics - Groundwater Level [mW]', (-11, -9))