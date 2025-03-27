# -*- coding: utf-8 -*-
"""
Created on Tue Nov 26 12:56:26 2024

@author: Xinyuan
"""


import pandas as pd
import numpy as np
from scipy.interpolate import interp1d
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

# Summarize dimension data into DataFrame
def create_area_df(greenhouse_unit_number, area_of_outer_walls, floor_area, ceiling_area, area_of_south_windows, window_area):
    area_data = {
        "Calculation": [
            "Greenhouse Unit Number",
            "Area of Outer Walls",
            "Floor Area",
            "Ceiling Area",
            "Direct Solar-Exposed Area",
            "Area of Other Windows"
        ],
        "Result (m^2)": [
            greenhouse_unit_number,
            area_of_outer_walls,
            floor_area,
            ceiling_area,
            area_of_south_windows,
            window_area
        ]
    }
    # Create DataFrame and return it
    return pd.DataFrame(area_data)

# Disclimer: 
    # This model is adapted from System Zoo 1 with adjustments in parameter.
    # This model uses open-loop Proportional (P) control approach: control decision = f(heating_limit_temperature, outdoor_temperature)
    # For temperature control, better algorithm can be developed: control decision = f(T_target - T_actural)

# Decision Variables
DESIRED_ROOM_TEMPERATURE = 10  # [°C]
greenhouse_unit_number = 60

# Parameters changed:
# insolation_data adjusted to Champaign
    # https://www.solarenergylocal.com/states/illinois/champaign/
# K_INSULATION_XX adjusted to Ave Overall Heat Transfer Coeff = 4.0 [Wm−2C−1]
    # https://www.sciencedirect.com/science/article/pii/S019689040500275X
# COLLECTOR_EFFICIENCY_SOUTH_WINDOWS
# local T data adjusted to Champaign
# ENERGY_CONTENT_OF_FUEL = adjusted to natural gas
# EFFICIENCY_OF_FURNACE unchanged


# Base Parameters
sim_start = 0
sim_end = 1.0
#im_end = 1.02
delta_t = 0.02      # TIME STEP
sim_index = pd.Index(np.arange(sim_start, sim_end, delta_t, dtype='float').round(2))


# Time relationships
months_per_year = 12  # [Month/Year]
days_per_year = 365  # [Day/Year]
hours_per_year = 8760  # [Hour/Year]
month  = {0   : 1,  0.09: 2,  0.17: 3, 
          0.25: 4,  0.34: 5,  0.42: 6, 
          0.50: 7,  0.59: 8,  0.67: 9, 
          0.75: 10, 0.84: 11, 0.92: 12}
month_series = sim_index*12

# Simulation time parameters
INITIAL_TIME = 0  # [Year]
FINAL_TIME = 1  # [Year]
TIME_STEP = 0.02  # [Year]
SAVEPER = TIME_STEP  # [Year]



# Building parameters
HOUSE_LENGTH_EAST_WEST = 52  # [m]
HOUSE_WIDTH_NORTH_SOUTH = 3.2 * greenhouse_unit_number  # [m]
HEIGHT_TO_EAVES = 4.3  # [m]
WINDOW_FRACTION_SOUTH_WALL = 0.99  # [1]
WINDOW_FRACTION_OTHER_WALLS = 0.99  # [1]

# Calculations for areas
area_of_outer_walls = 2 * HOUSE_WIDTH_NORTH_SOUTH * HEIGHT_TO_EAVES + 2 * HOUSE_LENGTH_EAST_WEST * HEIGHT_TO_EAVES  # [m^2]
floor_area = HOUSE_WIDTH_NORTH_SOUTH * HOUSE_LENGTH_EAST_WEST  # [m^2]
ceiling_area = floor_area*(2*1.75/3.2)  # [m^2]   # Triangle roof
area_of_south_windows = WINDOW_FRACTION_SOUTH_WALL * HOUSE_LENGTH_EAST_WEST * HEIGHT_TO_EAVES + 0.5*ceiling_area  # [m^2]
window_area = WINDOW_FRACTION_OTHER_WALLS * (2 * HOUSE_WIDTH_NORTH_SOUTH * HEIGHT_TO_EAVES + HOUSE_LENGTH_EAST_WEST * HEIGHT_TO_EAVES)  # [m^2]
# Display the greenhouse dimensions
df = create_area_df(int(greenhouse_unit_number), area_of_outer_walls, floor_area, ceiling_area, area_of_south_windows, window_area)
print(df)



# Heat transfer parameters
# EVA copolymer (U-value = 5.62 W m−2 K−1)
# Glass tempered (U-value = 3.25 W m−2 K−1)
# ground floor has a U-value of 0.378 W m−2 K−1

K_INSULATION_WALLS = 5.62  # [W/(K*m^2)]   
K_INSULATION_ROOF = 5.62  # [W/(K*m^2)]
K_INSULATION_FLOOR = 0.378  # [W/(K*m^2)]
K_INSULATION_WINDOWS = 5.62  # [W/(K*m^2)]
COLLECTOR_EFFICIENCY_SOUTH_WINDOWS = 0.7  # [1]    # transmmisivity = 0.7
# Temperatures
GROUND_TEMPERATURE = 9  # [°C]
#DESIRED_ROOM_TEMPERATURE = 15  # [°C]

## Table functions
# Insolation on vertical wall [W/(m^2)]
# insolation_data = {0.0: 102, 0.5: 101, 1.5: 111, 2.5: 105, 3.5: 79, 
#                    4.5: 45, 5.5: 38, 6.5: 48, 7.5: 69, 8.5: 99, 
#                    9.5: 102, 10.5: 106, 11.5: 104, 12: 102}
insolation_data = {0.0: 129.2, 1: 166.7, 2: 208.3, 3: 212.5, 4: 212.5, 
                   5: 229.17, 6: 233.3, 7: 229.17, 8: 229.17, 9: 187.5, 
                   10: 145.8, 11: 129.2, 12: 129.2}
insolation_on_vertical_wall_tbf = interp1d(list(insolation_data.keys()), 
                                           list(insolation_data.values()), 'linear')
print(insolation_on_vertical_wall_tbf)

# Outdoor temperature [°C]
# temperature_data = {0.0: 17, 0.5: 17.8, 1.5: 17.3, 2.5: 14.2, 3.5: 9.1, 
#                     4.5: 4.9, 5.5: 1.5, 6.5: 0, 7.5: 0.8, 8.5: 4.6, 
#                     9.5: 8.8, 10.5: 13.2, 11.5: 16.4, 12: 17}
temperature_data = {0.0: -3.06, 1: -1.56, 2: 5.56, 3: 11.44, 4: 17.72, 
                      5: 22.83, 6: 23.94, 7: 23.22, 8: 20.17, 9: 12.83, 
                      10: 5.5, 11: -0.06, 12: -3.06}
outdoor_temperature_tbf = interp1d(list(temperature_data.keys()), 
                                   list(temperature_data.values()), 'linear')

# Furnace parameters
ENERGY_CONTENT_OF_FUEL = 10.75  # [kW*Hour/m3]   #10.75 [kWh/m3] Wiki
EFFICIENCY_OF_FURNACE = 0.85  # [1]


# Initialize stock variables (integrators) - Soilwater model
fuel_oil_consumption = pd.Series(index=sim_index, dtype='float64')
heating_degree_days = pd.Series(index=sim_index, dtype='float64')
heating_period = pd.Series(index=sim_index, dtype='float64')
solar_heating_gain = pd.Series(index=sim_index, dtype='float64')
furnace_power_demand = pd.Series(index=sim_index, dtype='float64')
outdoor_temperature = pd.Series(index=sim_index, dtype='float64')

# Assign initial value to stock variables - Soilwater model
fuel_oil_consumption_j = 0                   # [1]       Initial value obtained in the loop
heating_degree_days_j = 0
heating_period_j = 0
solar_heating_gain_j = 0
furnace_power_demand_j = 0
outdoor_temperature_j = 0 

for i in range(len(sim_index)):
# Update stock variables    
    fuel_oil_consumption[sim_index[i]] = fuel_oil_consumption_j
    heating_degree_days[sim_index[i]] = heating_degree_days_j
    heating_period[sim_index[i]] = heating_period_j
    furnace_power_demand[sim_index[i]] = furnace_power_demand_j

    outdoor_temperature[sim_index[i]] = outdoor_temperature_tbf(month_series[i])
    outdoor_temperature_j = outdoor_temperature[sim_index[i]]
# Heat gains and losses
    solar_heating_gain_j = area_of_south_windows * COLLECTOR_EFFICIENCY_SOUTH_WINDOWS * insolation_on_vertical_wall_tbf(month_series[i])  # [W]
    
    solar_heating_gain[sim_index[i]] = solar_heating_gain_j
    
    specif_heat_transfer_upper = (K_INSULATION_WALLS * (area_of_outer_walls - (window_area + area_of_south_windows)) +
                                  K_INSULATION_ROOF * ceiling_area + K_INSULATION_WINDOWS * (window_area + area_of_south_windows))  # [W/Â°C]
    specif_heat_transfer_to_ground = K_INSULATION_FLOOR * floor_area  # [W/Â°C]
    heat_loss_to_ground = specif_heat_transfer_to_ground * (DESIRED_ROOM_TEMPERATURE - GROUND_TEMPERATURE)  # [W]
# Power demand
    heating_limit_temperature = DESIRED_ROOM_TEMPERATURE + (heat_loss_to_ground - solar_heating_gain_j) / specif_heat_transfer_upper  # [Â°C]
    heating_power_demand = specif_heat_transfer_upper * (DESIRED_ROOM_TEMPERATURE - outdoor_temperature[sim_index[i]]) + heat_loss_to_ground - solar_heating_gain_j if heating_limit_temperature > outdoor_temperature_tbf(month_series[i]) else 0    #[W]
    furnace_power_demand_j = heating_power_demand / EFFICIENCY_OF_FURNACE/1000  # [kW]
    
    heating_parameter = 1 if furnace_power_demand_j > 0 else 0  # [1]
    heating = (heating_limit_temperature - outdoor_temperature[sim_index[i]]) * days_per_year * heating_parameter  # [C*Day/Year]
    heating_operational = heating_parameter * days_per_year  # [Day/Year]
    burn_rate = (hours_per_year / ENERGY_CONTENT_OF_FUEL) * furnace_power_demand_j  # [m3/Year]

    fuel_oil_consumption_j = fuel_oil_consumption[sim_index[i]] + delta_t*burn_rate
    heating_degree_days_j = heating_degree_days[sim_index[i]] + delta_t*heating
    heating_period_j = heating_period[sim_index[i]] + delta_t*heating_operational
print(furnace_power_demand[sim_index[-1]])


#%%
final_index = sim_index[-1]
final_furnace_power_demand = furnace_power_demand.sum()*hours_per_year*0.02
final_fuel_oil_consumption = fuel_oil_consumption[final_index]
print(f"Total Natural Gas Consumption: {final_furnace_power_demand:.2f} [kWh]")
print(f"Total Natural Gas Consumption: {final_fuel_oil_consumption:.2f} [E03 m3]")

#%%
# Function to plot individual lines with specific configurations
def plot_line(ax, x, y, color, ylabel, ylim):
    ax.plot(x, y, color=color)
    ax.set_ylabel(ylabel, color=color)
    ax.tick_params(axis='y', labelcolor=color)
    ax.set_ylim(ylim)

# Combined plot function for heating and energy dynamics
def plot_heating_and_energy_dynamics(sim_index, outdoor_temperature, solar_heating_gain, furnace_power_demand, fuel_oil_consumption):
    fig, ax1 = plt.subplots(figsize=(10, 6))

    # Plot Outdoor Temperature
    plot_line(ax1, sim_index, outdoor_temperature, 'red', 'Outdoor Temperature [°C]', (-5, 30))
    ax1.set_ylabel('Outdoor Temperature [°C]', color='red', fontsize=16)

    # Create second y-axis for Solar Heating Gain
    ax2 = ax1.twinx()
    plot_line(ax2, sim_index, solar_heating_gain / 1000, 'blue', 'Solar Heating Gain [kW]', (0, 500))
    ax2.set_ylabel('Solar Heating Gain [kW]', color='blue', fontsize=16)

    # Create third y-axis for Furnace Power Demand
    ax3 = ax1.twinx()
    ax3.spines['right'].set_position(('outward', 60))  # Move the axis outward
    plot_line(ax3, sim_index, furnace_power_demand, 'green', 'Furnace Power Demand [kW]', (0, 500))
    ax3.set_ylabel('Furnace Power Demand [kW]', color='green', fontsize=16)

    # Create fourth y-axis for Natural Gas Consumption
    ax4 = ax1.twinx()
    ax4.spines['right'].set_position(('outward', 120))  # Move the axis outward further
    plot_line(ax4, sim_index, fuel_oil_consumption / 1000, 'purple', 'Cumu Natural Gas Consumption [E03 m3]', (0, 100))
    ax4.set_ylabel('Cumu Natural Gas Consumption [E03 m3]', color='purple', fontsize=16)

    # Enable grid only for the rightmost y-axis
    ax4.grid(True, linestyle='--', linewidth=0.5)  # Enable grid lines for ax4
    ax1.grid(False)  # Disable grid lines for the main axis

    # Labels and title
    ax1.set_xlabel('Time [Years]', fontsize=18)
    fig.suptitle('Heating and Energy Dynamics Over Time', fontsize=20)

    # Adjust layout to avoid overlap
    fig.tight_layout()
    plt.show()

# Example usage with your data
plot_heating_and_energy_dynamics(sim_index, outdoor_temperature, solar_heating_gain, furnace_power_demand, fuel_oil_consumption)
