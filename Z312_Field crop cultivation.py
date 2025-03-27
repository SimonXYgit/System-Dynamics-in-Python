# -*- coding: utf-8 -*-
"""
Created on Thu May 30 13:56:39 2024

@author: shang
"""
# mWT/mmWT: The height of water table in meter/milimeter


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

# Convert a dictionary of monthly rainfall into a Pandas Series step function.
def convert_to_step_series(monthly_rainfall_dict):
    step_series = pd.Series(index=sim_index, dtype=float)
    for key, value in monthly_rainfall_dict.items():
        step_series[sim_index >= key] = value
    step_series.ffill(inplace=True)
    return step_series

# Plot individual lines
def plot_line(ax, x, y, color, ylabel, ylim):
    ax.plot(x, y, color=color)
    ax.set_ylabel(ylabel, color=color)
    ax.tick_params(axis='y', labelcolor=color)
    ax.set_ylim(ylim)

# Combine individual lines into a diagram - nutrient
def plot_combined_nutrient_diagram(total_biomass, fraction_of_organic_matter_in_soil, plant_available_nitrogen):
    fig, ax1 = plt.subplots()
    
    # Plot total biomass
    plot_line(ax1, total_biomass.index, total_biomass, 'tab:red', 'Total Biomass [kg/ha]', (0, 20000))

    # Create a second y-axis for fraction of organic matter in soil
    ax2 = ax1.twinx()
    plot_line(ax2, fraction_of_organic_matter_in_soil.index, fraction_of_organic_matter_in_soil, 'tab:blue', 'OM Fraction in Soil [1]', (0, 0.1))

    # Create a third y-axis for plant available nitrogen
    ax3 = ax1.twinx()
    ax3.spines['right'].set_position(('outward', 60))  # Move the third axis outward
    plot_line(ax3, plant_available_nitrogen.index, plant_available_nitrogen, 'tab:green', 'N Available to Plants [kg/ha]', (0, 250))
    
    fig.tight_layout()  # To ensure the right y-label is not clipped
    plt.title('Model Soil Nutrient Dynamics')
    plt.show()
    
# Combine individual lines into a diagram - soilwater
def plot_combined_soilwater_diagram(precipitation, irrigation_amount, soilwater, groundwater_level):
    fig, ax1 = plt.subplots()

    # Plot precipitation
    plot_line(ax1, precipitation.index, precipitation, 'tab:blue', 'Precipitation [mWT/Year]', (0, 2.5))

    # Create a second y-axis for irrigation amount
    ax2 = ax1.twinx()
    plot_line(ax2, irrigation_amount.index, irrigation_amount, 'tab:orange', 'Irrigation Amount [mWT/Year]', (0, 10))

    # Create a third y-axis for soil water
    ax3 = ax1.twinx()
    ax3.spines['right'].set_position(('outward', 60))  # Move the third axis outward
    plot_line(ax3, soilwater.index, soilwater, 'tab:green', 'Soilwater [mWT]', (0, 0.5))

    # Create a fourth y-axis for groundwater level
    ax4 = ax1.twinx()
    ax4.spines['right'].set_position(('outward', 120))  # Move the fourth axis further outward
    plot_line(ax4, groundwater_level.index, groundwater_level, 'tab:red', 'Groundwater Level [mWT]', (-11, -9))

    fig.tight_layout()  # To ensure the right y-label is not clipped
    plt.title('Model Soil Water Dynamics')
    plt.show()

# Base Parameters
sim_start = 0
sim_end = 1
delta_t = 0.01      # TIME STEP
sim_index = pd.Index(np.arange(sim_start, sim_end, delta_t, dtype='float').round(2))

# Simulation Parameters
## Soil Parameters - Soilwater model
DEPTH_OF_PLOW_LAYER = 0.3       # [m]
ORGANIC_MATTER_IN_SOIL = 0.08   # [1] wolumn fraction
CLAY_FRACTION_OF_SOIL = 0.3     # [m] column fraction
FIELD_CAPACITY_GROUNDWATER_AQUIFER = 0.15   # [mWT/mSoil] m water per m soil
INITIAL_SOIL_MOISTURE = 0.5     # [1]
PERCOLATION_RATE = 50           # [1/Year] ********************Meaning?*********************
SOIL_COVER = 1                  # [1] 0: bare soil, 1: plant cover, 2: plastic foil cover
SOIL_COMPACTION_PARAMETER = 0   # [1] compaction occurs=1, doesn't occur=0

## Soil parameters - Nutrient model
INITIAL_CN_IN_NUTRIENT_HUMUS = 20                           # [1]
ORGANIC_MATTER_IN_SOIL = 0.04                               # [1] volume fraction
INITIAL_FRACTION_PERMANENT_HUMUS_IN_ORGANIC_MATTER = 0.75   # [1] volume fraction
SOIL_LOSS_BY_EROSION = 5000                                 # [kg/(Year*ha)]
INITIAL_AVAILABLE_NITROGEN = 50                             # [kg/ha]
LEACHING_RATE = 0.3                                         # [1/Year]
SOIL_SPECIFIC_WEIGHT = 1600                                 # [kg/(m*m*m)]


## Fertilization parameters - Nutrient model
MANURE_APPLIED = 15000              # [kg/ha]
COMPOST_APPLIED = 10000             # [kg/ha]
ORGANIC_FERTILIZER_DATE = 0.2       # [Year] 0.20 = 10th week = beginning of March
MINERAL_FERTILIZER_DATE = 0.25      # [Year] 0.25 = 13th week = end of March
NITROGEN_FERTILIZER_APPLIED = 80    # [kg/ha]
CN_IN_COMPOST = 15                  # [1] ratio of C to N in compost
CN_IN_MANURE = 20                   # [1] ratio of C to N in manure
ODM_FRACTION_OF_MANURE = 0.25       # *1 [1] kg ODM per kg fresh matter
ODM_FRACTION_IN_COMPOST = 0.35      # *1 [1] kg ODM per kg fresh matter
NITROGEN_INPUT_FROM_ATMOSPHERE = 25 # kg/(ha*Year)
DECOMPOSITION_RATE = 0.2            # [1/Year]


## Crop Parameter - Soilwater model
MAX_CROP_YIELD = 6500               # [kg/ha] biomass in kgODM (organic dry matter) per hactare
SPECIFIC_TRANSPIRATION_RATE = 0.4   # [mWT/(kg/ha)] m water transpired per biomass (biomass density in [kg/ha])
ROOT_DEPTH_OF_FIELD_CROP = 1        # [m]
BEGIN_GROWTH_PERIOD = 0.3           # [1 Year]
HARVEST_TIME = 0.6                  # [1 Year]
FERTILIZATION_FACTOR = 1            # 1: optimal fertilization; 0: no fert, intermediates are valid for calculation
PUMP_RATE = 150                     # [1/Year] no irrigation=0 ********************Meaning?*********************

## Crop parameters - Nutrient model
SPECIFIC_NITROGEN_UPTAKE = 0.029        # [1] kg N per kg ODM harvest yield
HARVEST_SPECIFIC_WEIGHT = 1.15          # [1] kg fresh weight per kg ODM
MAX_CROP_RESIDUES = 1700                # [kg/ha]
AMOUNT_STRAW_AND_LEAVES = 7700          # [kg/ha]
CN_RATIO_IN_STRAW = 80                  # [1]
MAX_NITROGEN_IN_CROP_RESIDUES = 17      # [kg/ha]
NITROGEN_FIXATION = 0                   # [1] legumes only, otherwise = 0, in kg N per kg ODM yield
STRAW_AND_LEAVES_REMAIN_IN_FIELD = 1    # [1] straw or leaves remain in field = 1, are harvested = 0


## Weather parameters - Soilwater model
# Monthly average rainfall [mWT/Month]
MONTHLY_RAINFALL  = {0   : 0.068, 0.09: 0.068, 0.17: 0.068, 
                     0.25: 0.047, 0.34: 0.061, 0.42: 0.063, 
                     0.50: 0.075, 0.59: 0.086, 0.67: 0.09, 
                     0.75: 0.066, 0.84: 0.067, 0.92: 0.068}
MONTHLY_RAINFALL_series = convert_to_step_series(MONTHLY_RAINFALL)
MULTIPLIER_FOR_RAINFALL = 1             # 1: optimal rain, <1: less than optimum
ANNUAL_EVAPORATION = 0.444              # [mWT/Year]
WEATHER_PERIOD = 3                      # [Day]


## Other parameters - Soilwater model
growth_time = {0: 0, BEGIN_GROWTH_PERIOD: 1, HARVEST_TIME+delta_t: 0}
growth_period = step_function_series(growth_time, sim_index)
MONTHS_PER_YEAR = 12
PI = 3.14159
DELAY_TIME = 0.25                       # [Year]

## Other parameters - Nutrient model
ORGANIC_MATTER_SPECIFIC_WEIGHT = 650    # [kg/(m*m*m)] spec. weight, dry matter
sqm_per_ha = 10000                      # [m*m /ha]
C_per_OM = 0.47                         # [1] kg C per kg OM


## Table functions - Soilwater model
# Influence of clay fraction[1] on groundwater level capacity[mWT/mSoil]
usable_field_capacity_of_soil = {0.0: 0.01, 0.1: 0.1,   0.2: 0.16,  0.3: 0.2, 
                                 0.4: 0.21, 0.5: 0.15,  0.6: 0.13,  1.0: 0.05}
usable_field_capacity_of_soil_tbf = interp1d(list(usable_field_capacity_of_soil.keys()), 
                                             list(usable_field_capacity_of_soil.values()), 'linear')
# Influence of clay fraction[1] on capillary rise of groundwater level[mWT]
capillary_rise = {0: 0.1, 1: 3}
capillary_rise_tbf = interp1d(list(capillary_rise.keys()),
                              list(capillary_rise.values()), 'linear')
# Influenc of soil moisture[1] on plant growth[1]
moisture_effect = {0: 0, 0.1: 0.2, 0.5: 1.0, 1: 1, 2: 1}
moisture_effect_tbf = interp1d(list(moisture_effect.keys()),
                               list(moisture_effect.values()), 'linear')


## Table functions - Nutrient model
# Influence of 1/CN_ratio[1] on decomposition factor[1]
decomposition_factor = {0.0: 0.0, 0.01: 0.05, 0.025: 0.25, 0.033: 0.5, 0.04: 0.9, 0.05: 1.0, 1.0: 1.0}
decomposition_factor_tbf = interp1d(list(decomposition_factor.keys()), 
                                    list(decomposition_factor.values()), 'linear')
# Influence of 1/CN_ratio[1] on nitrogen transfer[1/Year]
N_transfer_function = {0.0: 0.0, 0.01: 0.95, 0.025: 0.75, 0.033: 0.5, 0.04: 0.2, 0.05: 0.0, 1.0: 0.0}
N_transfer_function_tbf = interp1d(list(N_transfer_function.keys()), 
                                   list(N_transfer_function.values()), 'linear')
# Influence of OM content[1] on N leaching factor[1]
nitrogen_leaching = {0.0: 1.0, 0.05: 0.5, 0.1: 0.2, 1.0: 0.1}
nitrogen_leaching_tbf = interp1d(list(nitrogen_leaching.keys()), 
                                 list(nitrogen_leaching.values()), 'linear')
# Influence of N availability[1] on crop growth factor[1]
nitrogen_effect = {0: 0, 0.2: 0.2, 0.35: 0.5, 0.5: 0.8, 1: 1, 2: 0.9, 3: 0.4, 5: 0.1, 10: 0}
nitrogen_effect_tbf = interp1d(list(nitrogen_effect.keys()), 
                               list(nitrogen_effect.values()), 'linear')


# Initialize stock variables (integrators) - Soilwater model
weather = pd.Series(index=sim_index, dtype='float64')
rain_amount = pd.Series(index=sim_index, dtype='float64')
precipitation = pd.Series(index=sim_index, dtype='float64')
soil_moisture = pd.Series(index=sim_index, dtype='float64')
groundwater_level = pd.Series(index=sim_index, dtype='float64')
relative_plant_biomass = pd.Series(index=sim_index, dtype='float64')
drought_period = pd.Series(index=sim_index, dtype='float64')
total_biomass =pd.Series(index=sim_index, dtype='float64')
soilwater =pd.Series(index=sim_index, dtype='float64')
irrigation_amount =pd.Series(index=sim_index, dtype='float64')

# Initialize stock variables (integrators) - Nutrient model
carbon_in_nutrient_humus = pd.Series(index=sim_index, dtype='float64')
nitrogen_in_nutrient_humus = pd.Series(index=sim_index, dtype='float64')
plant_available_nitrogen = pd.Series(index=sim_index, dtype='float64')
carbon_in_permanent_humus = pd.Series(index=sim_index, dtype='float64')
relative_crop_yield = pd.Series(index=sim_index, dtype='float64')
fraction_of_organic_matter_in_soil = pd.Series(index=sim_index, dtype='float64')


# Assign initial value to stock variables - Soilwater model
weather_j = 0                   # [1]       Initial value obtained in the loop
rain_amount_j = 0               # [mWT]      Accumulative precipitation amount, start from 0
precipitation_j = 0             # [mWT/Year] Initial value obtained in the loop
soil_moisture_j = INITIAL_SOIL_MOISTURE     # [1]
groundwater_level_j = -10       # [mWT]      negative for underground
relative_plant_biomass_j = 0.01 # [1]
drought_period_j = 0            # [1]       Initial value obtained in the loop
total_biomass_j = 0             # [kg/ha]   Initial value obtained in the loop
soilwater_j = 0                 # [mWT]      Initial value obtained in the loop, attainable_water_capacity required first
irrigation_amount_j = 0         # [mWT]      Accumulative irrigation amount, start from 0

# Assign initial value to stock variables - Nutrient model
carbon_in_nutrient_humus_j = 0.0                # Assign value in loop
nitrogen_in_nutrient_humus_j = 0.0              # Assign value in loop
plant_available_nitrogen_j = INITIAL_AVAILABLE_NITROGEN     # [kg/ha]
carbon_in_permanent_humus_j = ORGANIC_MATTER_IN_SOIL *(1-INITIAL_FRACTION_PERMANENT_HUMUS_IN_ORGANIC_MATTER) *DEPTH_OF_PLOW_LAYER *sqm_per_ha *ORGANIC_MATTER_SPECIFIC_WEIGHT *C_per_OM *INITIAL_FRACTION_PERMANENT_HUMUS_IN_ORGANIC_MATTER /(1-INITIAL_FRACTION_PERMANENT_HUMUS_IN_ORGANIC_MATTER) #[kg/ha]
relative_crop_yield_j = 0.0                     # [1]
fraction_of_organic_matter_in_soil_j = 0.0      # Assign value in loop



for i in sim_index:
# Update values for stock variables - Soilwater model
    # value used in this iteration=[i], value updated for the next iteration=_j 
    # Uppercase=constants, Lowercase=intermediates
    weather[i] = weather_j
    rain_amount[i] = rain_amount_j
    soil_moisture[i] = soil_moisture_j
    groundwater_level[i] = groundwater_level_j
    relative_plant_biomass[i] = relative_plant_biomass_j #
    drought_period[i] = drought_period_j
    total_biomass[i] = total_biomass_j #
    soilwater[i] = soilwater_j
    irrigation_amount[i] = irrigation_amount_j
# Update values for stock variables - Nutrient model
    carbon_in_nutrient_humus[i] = carbon_in_nutrient_humus_j
    nitrogen_in_nutrient_humus[i] = nitrogen_in_nutrient_humus_j
    plant_available_nitrogen[i] = plant_available_nitrogen_j
    carbon_in_permanent_humus[i] = carbon_in_permanent_humus_j
    relative_crop_yield[i] = relative_crop_yield_j
    fraction_of_organic_matter_in_soil[i] = fraction_of_organic_matter_in_soil_j

## Soilwater model
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
    modified_rainfall = MONTHLY_RAINFALL_series[i] * MULTIPLIER_FOR_RAINFALL   # [mWT/Month]
    # Why times 2?
    rain = MONTHS_PER_YEAR*modified_rainfall*weather_j*2                # [mWT/Year]
    rain_amount_j = rain_amount[i] + delta_t*rain                       # [mWT]
    precipitation[i] = rain                                             # [mWT/Year]
# Waterholding capacity
    # Waterholding capacity defined by organic matter [mWT]
    waterholding_capacity_of_organic_matter = DEPTH_OF_PLOW_LAYER * ORGANIC_MATTER_IN_SOIL * 5
    # Waterholding capacity corrected by soil compression [mWT]
    reduced_field_capacity_from_compaction = usable_field_capacity_of_soil_tbf(CLAY_FRACTION_OF_SOIL) * (0.94 if CLAY_FRACTION_OF_SOIL<0.15 else 0.88) if SOIL_COMPACTION_PARAMETER>0 else usable_field_capacity_of_soil_tbf(CLAY_FRACTION_OF_SOIL)
    # Waterholding capacity corrected by root depth and capillary rise [mWT]
    depth_of_accessible_water_horizon = ROOT_DEPTH_OF_FIELD_CROP + capillary_rise_tbf(CLAY_FRACTION_OF_SOIL)
    # Waterholding capacity_final [mWT]
    groundwater_access = 1 if -depth_of_accessible_water_horizon < groundwater_level[i] else 0
    attainable_water_capacity = depth_of_accessible_water_horizon * reduced_field_capacity_from_compaction if groundwater_access>0 else waterholding_capacity_of_organic_matter + ROOT_DEPTH_OF_FIELD_CROP*reduced_field_capacity_from_compaction
# Evaporation [mWT/Year]
    evaporation = soil_moisture[i]*ANNUAL_EVAPORATION*(1+0.95*math.sin(2*PI*(i/365-DELAY_TIME)))*(1-SOIL_COVER/2)

## Nutrient model
# Tier 1
# "C_input" for "C_input_in_nutrient_humus"
    # Carbon_in_nutrient_humus (CINH) inputs    
    C_in_crop_residue = C_per_OM*relative_plant_biomass[i] *(STRAW_AND_LEAVES_REMAIN_IN_FIELD *AMOUNT_STRAW_AND_LEAVES +MAX_CROP_RESIDUES)/delta_t  if abs(i-HARVEST_TIME) < delta_t/2  else 0
    C_input_from_crop_residues = C_in_crop_residue      #[kg/(ha*Year)]
    C_in_manure = MANURE_APPLIED * ODM_FRACTION_OF_MANURE * C_per_OM        #[kg/ha]
    C_in_compost = COMPOST_APPLIED * ODM_FRACTION_IN_COMPOST * C_per_OM     #[kg/ha]
    # CINH outputs
    relative_erosion_loss = SOIL_LOSS_BY_EROSION /(DEPTH_OF_PLOW_LAYER *sqm_per_ha *SOIL_SPECIFIC_WEIGHT)  #[1/Year]
    initial_C_in_nutrient_humus = ORGANIC_MATTER_IN_SOIL *(1 -INITIAL_FRACTION_PERMANENT_HUMUS_IN_ORGANIC_MATTER) *DEPTH_OF_PLOW_LAYER *sqm_per_ha *ORGANIC_MATTER_SPECIFIC_WEIGHT *C_per_OM   # [kg/ha]
    ## Assigning initial value in the first iteration
    carbon_in_nutrient_humus[i] = initial_C_in_nutrient_humus if i<0.01  else carbon_in_nutrient_humus[i]
    nitrogen_in_nutrient_humus[i] = initial_C_in_nutrient_humus/INITIAL_CN_IN_NUTRIENT_HUMUS if i<0.01  else nitrogen_in_nutrient_humus[i]
    N_in_crop_residues = (C_per_OM * relative_plant_biomass[i] * (STRAW_AND_LEAVES_REMAIN_IN_FIELD * AMOUNT_STRAW_AND_LEAVES / CN_RATIO_IN_STRAW + MAX_NITROGEN_IN_CROP_RESIDUES)) / delta_t if abs(i - HARVEST_TIME) < delta_t / 2 else 0 #[kg/(ha*year)]
    N_input_from_crop_residues = N_in_crop_residues     # [kg/(ha*Year)]
# Tier 2
    # C_inputs - completed
    C_input_from_organic_fertilizer = (C_in_manure + C_in_compost) /delta_t  if abs(i-ORGANIC_FERTILIZER_DATE) < delta_t/2  else 0  #[kg/(Year*ha)]
    # CINH outputs 
    C_loss_by_erosion = relative_erosion_loss * carbon_in_nutrient_humus[i] #[kg/(ha*Year)]
    CN_ratio = carbon_in_nutrient_humus[i] / nitrogen_in_nutrient_humus[i]
    N_input_from_organic_fertilizer = (C_in_manure/CN_IN_MANURE +C_in_compost/CN_IN_COMPOST)/delta_t if abs(i-ORGANIC_FERTILIZER_DATE) < delta_t/2  else 0
# Tier 3
    N_input_by_transfer = N_transfer_function_tbf(1/CN_ratio) * plant_available_nitrogen[i]    # [kg/(ha*Year)]
    decomposition_of_organic_matter = DECOMPOSITION_RATE * decomposition_factor_tbf(1/CN_ratio) * carbon_in_nutrient_humus[i] * (1+0.5*math.sin(2*PI*(i-DELAY_TIME))) #[kg/(ha*Year)]
# Tier 4
    N_loss_from_org_matter_decomposition = decomposition_of_organic_matter / CN_ratio   # [kg/(ha*Year)]
    N_loss_from_erosion = C_loss_by_erosion / CN_ratio  # [kg/(ha*Year)]  
# Tier 5.1
    N_input_from_org_matter_decomposition = N_loss_from_org_matter_decomposition #[kg/(ha*Year)]
    N_from_atmosphere = NITROGEN_INPUT_FROM_ATMOSPHERE #[kg/(Year*ha)]
    N_input_from_mineral_fertilizer = NITROGEN_FERTILIZER_APPLIED / delta_t if abs(i - MINERAL_FERTILIZER_DATE) < delta_t / 2 else 0 #[kg/(ha*Year)]
    relative_nitrogen_availability = plant_available_nitrogen[i] / (MAX_CROP_YIELD * SPECIFIC_NITROGEN_UPTAKE) #[1]

## Soilwater model
# Plant growth rate
    # Relative growth_rate definitive function (Notice: tend to reach max when biomass_rel = 0.5) [1/Year]
    relative_growth_rate = (25*(20/52)/(HARVEST_TIME -BEGIN_GROWTH_PERIOD)) *relative_plant_biomass[i] *(1-relative_plant_biomass[i]) *moisture_effect_tbf(soil_moisture[i]) *nitrogen_effect_tbf(relative_nitrogen_availability)
    # Relative growth corrected by growth period, a time window for the process to happen [1/Year]
    relative_growth_rate = relative_growth_rate * growth_period[i]
    # Growth rate_final [kg/(ha*Year)]
    growth_rate = relative_growth_rate * MAX_CROP_YIELD
    growth_rate_total_biomass = relative_growth_rate * total_biomass[i]
    
## Nutrient model
# Tier 5.2
    N_uptake_by_plants = growth_rate * (SPECIFIC_NITROGEN_UPTAKE - NITROGEN_FIXATION) #[kg/(ha*Year)]
    C_input_to_permanent_humus = 0.25 * decomposition_of_organic_matter #[kg/(ha*Year)]
    C_losses_of_permanent_humus = carbon_in_permanent_humus[i] * (0.2 + relative_erosion_loss) #[kg/(ha*Year)]
    
# Updates
# Carbon in permanent humus
    carbon_in_permanent_humus_j = carbon_in_permanent_humus[i] + delta_t*(C_input_to_permanent_humus - C_losses_of_permanent_humus)
# Fraction of organic matter in soil
    fraction_of_organic_matter_in_soil[i] = (carbon_in_nutrient_humus[i] + carbon_in_permanent_humus[i]) / (DEPTH_OF_PLOW_LAYER * sqm_per_ha * ORGANIC_MATTER_SPECIFIC_WEIGHT * C_per_OM) #[1]
# Plant available Nitrogen
    N_loss_by_leaching = LEACHING_RATE * plant_available_nitrogen[i] * nitrogen_leaching_tbf(fraction_of_organic_matter_in_soil[i]) #[kg/(ha*Year)]
    N_loss_by_transfer_to_org_matter = N_input_by_transfer #[kg/(ha*Year)] 
    plant_available_nitrogen_j = plant_available_nitrogen[i] + delta_t*(N_from_atmosphere +N_input_from_org_matter_decomposition +N_input_from_mineral_fertilizer -N_uptake_by_plants -N_loss_by_leaching -N_loss_by_transfer_to_org_matter)   # [kg/ha]
# Nitrogen in nutrient humus
    nitrogen_in_nutrient_humus_j = nitrogen_in_nutrient_humus[i] + delta_t*(N_input_from_organic_fertilizer +N_input_from_crop_residues +N_input_by_transfer +N_loss_from_org_matter_decomposition -N_loss_from_erosion)
# Carbon in nutrient humus 
    carbon_in_nutrient_humus_j = carbon_in_nutrient_humus[i] + delta_t*(C_input_from_organic_fertilizer + C_input_from_crop_residues - C_loss_by_erosion - decomposition_of_organic_matter)    # [kg/ha]
    fraction_permanent_humus_in_org_matter = carbon_in_permanent_humus[i] / (carbon_in_permanent_humus[i] + carbon_in_nutrient_humus[i]) #[1]
    harvest = relative_plant_biomass[i] /delta_t if abs(i-HARVEST_TIME) < delta_t/2 else 0 #[1/year]
# Relative crop yield
    relative_crop_yield_j = relative_crop_yield[i] + delta_t*harvest  # [1]
    harvest_yield_green_weight = HARVEST_SPECIFIC_WEIGHT * relative_crop_yield[i] * MAX_CROP_YIELD #[kg/ha]
# Carbon in crop residues
    C_in_crop_residues = (C_per_OM * relative_plant_biomass[i] * (STRAW_AND_LEAVES_REMAIN_IN_FIELD * AMOUNT_STRAW_AND_LEAVES + MAX_CROP_RESIDUES)) / delta_t if abs(i - HARVEST_TIME) < delta_t / 2 else 0 #[kg/(ha*year)]
    
## Soilwater model
# Transpiration [mWT/Year]
    transpiration = (1/FERTILIZATION_FACTOR) * SPECIFIC_TRANSPIRATION_RATE*growth_rate_total_biomass/10000
# Biomass [kg/ha]
    # Drought happens if soil moisture is too low.
    drought = 0     if soil_moisture[i]>0.2 else 1                              # [1/Year]
    drought_end = drought_period[i]/delta_t     if soil_moisture[i]>0.2 else 0  # [1/Year]
    drought_period_j = drought_period[i] + delta_t*(drought-drought_end)        # [1]
    # If drought period is too long, or if harvesting occurs, then withering starts.
    #   Withering leads to sharp decrease of biomass. 
    harvest_loss = 1    if drought_period[i]>(20/365) and i>BEGIN_GROWTH_PERIOD and i<HARVEST_TIME else 0   # [1]
    withering = relative_plant_biomass[i]/delta_t      if i>HARVEST_TIME or harvest_loss == 1 else 0        # [1/Year]
    # Biomass calculated from biomass_rel
    relative_plant_biomass_j = relative_plant_biomass[i] + delta_t*(relative_growth_rate - harvest - withering)              # [1]
    # biomass_final [kg/ha]
    total_biomass_j = relative_plant_biomass[i] * (MAX_CROP_YIELD + AMOUNT_STRAW_AND_LEAVES + MAX_CROP_RESIDUES)

# Soil moisture [1]
    # Assigning initial value to soilwater [mWT (only as a boost in the first iteration)
    soilwater[i] = INITIAL_SOIL_MOISTURE*attainable_water_capacity  if i<0.01 else soilwater[i]    # [mmWT]
    soil_moisture_j = soilwater[i]/attainable_water_capacity        # [1]
    soilwater_surplus = soilwater[i] - attainable_water_capacity    # [m]
# Irrigation [mWT/Year]
    irrigation = (0.5*attainable_water_capacity - soilwater[i])*PUMP_RATE if (soil_moisture[i]<0.5 and i>BEGIN_GROWTH_PERIOD and i<HARVEST_TIME) else 0
    irrigation_amount = irrigation_amount + delta_t*irrigation      # [mWT]
# Groundwater level [mWT]
    drop_of_groundwater_level = irrigation/FIELD_CAPACITY_GROUNDWATER_AQUIFER
    rise_of_groundwater_level = soilwater_surplus/FIELD_CAPACITY_GROUNDWATER_AQUIFER*PERCOLATION_RATE   if soilwater_surplus>0 else 0
    percolation = rise_of_groundwater_level
    groundwater_level_j = groundwater_level[i] + delta_t*(rise_of_groundwater_level-drop_of_groundwater_level)
    soilwater_j = soilwater[i] + delta_t*(irrigation +precipitation[i] -transpiration -evaporation -percolation)

## Output results
plot_combined_soilwater_diagram(precipitation, irrigation_amount, soilwater, groundwater_level)
plot_combined_nutrient_diagram(total_biomass, fraction_of_organic_matter_in_soil, plant_available_nitrogen)

