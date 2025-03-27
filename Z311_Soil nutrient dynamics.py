# -*- coding: utf-8 -*-
"""
Created on Fri Apr 19 12:07:51 2024

@author: shang
"""

import pandas as pd
import numpy as np
import math
from scipy.interpolate import interp1d
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
# Soil parameters
INITIAL_CN_IN_NUTRIENT_HUMUS = 20                           # [1]
ORGANIC_MATTER_IN_SOIL = 0.04                               # [1] volume fraction
INITIAL_FRACTION_PERMANENT_HUMUS_IN_ORGANIC_MATTER = 0.75   # [1] volume fraction
DEPTH_OF_PLOW_LAYER = 0.2                                   # [m]
SOIL_LOSS_BY_EROSION = 5000                                 # [kg/(Year*ha)]
INITIAL_AVAILABLE_NITROGEN = 50                             # [kg/ha]
LEACHING_RATE = 0.3                                         # [1/Year]
SOIL_SPECIFIC_WEIGHT = 1600                                 # [kg/(m*m*m)]

# Fertilization parameters
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

# Crop parameters
MAX_CROP_YIELD = 6500                   # [kg/ha]
SPECIFIC_NITROGEN_UPTAKE = 0.029        # [1] kg N per kg ODM harvest yield
HARVEST_SPECIFIC_WEIGHT = 1.15          # [1] kg fresh weight per kg ODM
MAX_CROP_RESIDUES = 1700                # [kg/ha]
AMOUNT_STRAW_AND_LEAVES = 7700          # [kg/ha]
CN_RATIO_IN_STRAW = 80                  # [1]
MAX_NITROGEN_IN_CROP_RESIDUES = 17      # [kg/ha]
NITROGEN_FIXATION = 0                   # [1] legumes only, otherwise = 0, in kg N per kg ODM yield
STRAW_AND_LEAVES_REMAIN_IN_FIELD = 1    # [1] straw or leaves remain in field = 1, are harvested = 0
BEGIN_GROWTH_PERIOD = 0.3               # [Year]
HARVEST_TIME = 0.6                      # [Year]
MOISTURE_EFFECT = 1                     # [1]

# Other parameters
ORGANIC_MATTER_SPECIFIC_WEIGHT = 650    # [kg/(m*m*m)] spec. weight, dry matter
PI = 3.14159                            # [1/Year]
sqm_per_ha = 10000                      # [m*m /ha]
DELAY_TIME = 0.25                       # [Year]
C_per_OM = 0.47                         # [1] kg C per kg OM


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



# Initialize stock variables (integrators)
relative_plant_biomass = pd.Series(index=sim_index, dtype='float64')
carbon_in_nutrient_humus = pd.Series(index=sim_index, dtype='float64')
nitrogen_in_nutrient_humus = pd.Series(index=sim_index, dtype='float64')
plant_available_nitrogen = pd.Series(index=sim_index, dtype='float64')
carbon_in_permanent_humus = pd.Series(index=sim_index, dtype='float64')
relative_crop_yield = pd.Series(index=sim_index, dtype='float64')
total_biomass = pd.Series(index=sim_index, dtype='float64')
fraction_of_organic_matter_in_soil = pd.Series(index=sim_index, dtype='float64')


# Assign initial value to stock variables
relative_plant_biomass_j = 0.01                 # [1]
carbon_in_nutrient_humus_j = 0.0                # Assign value in loop
nitrogen_in_nutrient_humus_j = 0.0              # Assign value in loop
plant_available_nitrogen_j = INITIAL_AVAILABLE_NITROGEN     # [kg/ha]
carbon_in_permanent_humus_j = ORGANIC_MATTER_IN_SOIL *(1-INITIAL_FRACTION_PERMANENT_HUMUS_IN_ORGANIC_MATTER) *DEPTH_OF_PLOW_LAYER *sqm_per_ha *ORGANIC_MATTER_SPECIFIC_WEIGHT *C_per_OM *INITIAL_FRACTION_PERMANENT_HUMUS_IN_ORGANIC_MATTER /(1-INITIAL_FRACTION_PERMANENT_HUMUS_IN_ORGANIC_MATTER) #[kg/ha]
relative_crop_yield_j = 0.0                     # [1]
total_biomass_j = 0.0                           # Assign value in loop
fraction_of_organic_matter_in_soil_j = 0.0      # Assign value in loop



for i in sim_index:
# Update values for stock variables
    relative_plant_biomass[i] = relative_plant_biomass_j 
    carbon_in_nutrient_humus[i] = carbon_in_nutrient_humus_j
    nitrogen_in_nutrient_humus[i] = nitrogen_in_nutrient_humus_j
    plant_available_nitrogen[i] = plant_available_nitrogen_j
    carbon_in_permanent_humus[i] = carbon_in_permanent_humus_j
    relative_crop_yield[i] = relative_crop_yield_j
    total_biomass[i] = total_biomass_j
    fraction_of_organic_matter_in_soil[i] = fraction_of_organic_matter_in_soil_j
    
# Tier 1
# C input for "C_input_in_nutrient_humus"
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
    # CINH inputs - completed
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
    
# Tier 5
    N_input_from_org_matter_decomposition = N_loss_from_org_matter_decomposition #[kg/(ha*Year)]
    
    

    N_from_atmosphere = NITROGEN_INPUT_FROM_ATMOSPHERE #[kg/(Year*ha)]
    N_input_from_mineral_fertilizer = NITROGEN_FERTILIZER_APPLIED / delta_t if abs(i - MINERAL_FERTILIZER_DATE) < delta_t / 2 else 0 #[kg/(ha*Year)]
    relative_nitrogen_availability = plant_available_nitrogen[i] / (MAX_CROP_YIELD * SPECIFIC_NITROGEN_UPTAKE) #[1]
    relative_growth_rate = 0 if (i<BEGIN_GROWTH_PERIOD) or (i>HARVEST_TIME) else (25*(20/52)/(HARVEST_TIME-BEGIN_GROWTH_PERIOD)) *relative_plant_biomass[i] *(1-relative_plant_biomass[i]) *MOISTURE_EFFECT *nitrogen_effect_tbf(relative_nitrogen_availability) #[1/year]
    growth_rate = relative_growth_rate * MAX_CROP_YIELD #[kg/(ha*year)]
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
# Total biomass
    relative_plant_biomass_j = relative_plant_biomass[i] + delta_t*(relative_growth_rate - harvest) # [1]
    C_in_crop_residues = (C_per_OM * relative_plant_biomass[i] * (STRAW_AND_LEAVES_REMAIN_IN_FIELD * AMOUNT_STRAW_AND_LEAVES + MAX_CROP_RESIDUES)) / delta_t if abs(i - HARVEST_TIME) < delta_t / 2 else 0 #[kg/(ha*year)]
    total_biomass[i] = relative_plant_biomass[i] * (MAX_CROP_YIELD + AMOUNT_STRAW_AND_LEAVES + MAX_CROP_RESIDUES) # [kg/ha]





## Output results
plot_diagram(total_biomass, 'Model Soil Nutrient Dynamics - Biomass [kg/ha]', (0, 20000))
plot_diagram(fraction_of_organic_matter_in_soil, 'Model Soil Nutrient Dynamics - OM Fraction in Soil [1]', (0, 0.1))
plot_diagram(plant_available_nitrogen, 'Model Soil Nutrient Dynamics - N Available to Plants [kg/ha]', (0, 250))