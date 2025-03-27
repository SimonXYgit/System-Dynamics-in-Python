# -*- coding: utf-8 -*-
"""
Created on Mon Dec  2 23:00:05 2024

@author: shang
"""

import pandas as pd
import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import os

# Define simulation and greenhouse parameters
sim_start = 0
sim_end = 1.0
delta_t = 0.02
sim_index = pd.Index(np.arange(sim_start, sim_end, delta_t, dtype='float').round(2))

# Fixed parameters
HOUSE_LENGTH_EAST_WEST = 52  # [m]
HEIGHT_TO_EAVES = 4.3  # [m]
WINDOW_FRACTION_SOUTH_WALL = 0.99  # [1]
WINDOW_FRACTION_OTHER_WALLS = 0.99  # [1]
COLLECTOR_EFFICIENCY_SOUTH_WINDOWS = 0.7  # [1]
ENERGY_CONTENT_OF_FUEL = 10.75  # [kW*Hour/m3]
EFFICIENCY_OF_FURNACE = 0.85  # [1]
# EVA copolymer (U-value = 5.62 W m−2 K−1)
# Glass tempered (U-value = 3.25 W m−2 K−1)
# ground floor has a U-value of 0.378 W m−2 K−1
K_INSULATION_WALLS = 5.62  # [W/(K*m^2)]
K_INSULATION_ROOF = 5.62  # [W/(K*m^2)]
K_INSULATION_FLOOR = 0.378  # [W/(K*m^2)]
K_INSULATION_WINDOWS = 5.62  # [W/(K*m^2)]
GROUND_TEMPERATURE = 9  # [°C]

# Interpolation data
insolation_data = {0.0: 129.2, 1: 166.7, 2: 208.3, 3: 212.5, 4: 212.5, 
                   5: 229.17, 6: 233.3, 7: 229.17, 8: 229.17, 9: 187.5, 
                   10: 145.8, 11: 129.2, 12: 129.2}
temperature_data = {0.0: -3.06, 1: -1.56, 2: 5.56, 3: 11.44, 4: 17.72, 
                    5: 22.83, 6: 23.94, 7: 23.22, 8: 20.17, 9: 12.83, 
                    10: 5.5, 11: -0.06, 12: -3.06}
insolation_on_vertical_wall_tbf = interp1d(list(insolation_data.keys()), 
                                           list(insolation_data.values()), 'linear')
outdoor_temperature_tbf = interp1d(list(temperature_data.keys()), 
                                   list(temperature_data.values()), 'linear')

# Simulation ranges for decision variables
greenhouse_unit_numbers = [1, 2, 3, 4, 5, 10, 15, 20, 40, 60]
desired_room_temperatures = range(10, 26)

# Initialize result storage
results = []

# Simulation loop
for greenhouse_unit_number in greenhouse_unit_numbers:
    for DESIRED_ROOM_TEMPERATURE in desired_room_temperatures:
        # Calculate areas
        HOUSE_WIDTH_NORTH_SOUTH = 3.2 * greenhouse_unit_number
        area_of_outer_walls = (2 * HOUSE_WIDTH_NORTH_SOUTH * HEIGHT_TO_EAVES + 
                               2 * HOUSE_LENGTH_EAST_WEST * HEIGHT_TO_EAVES)
        floor_area = HOUSE_WIDTH_NORTH_SOUTH * HOUSE_LENGTH_EAST_WEST
        ceiling_area = floor_area * (2 * 1.75 / 3.2)
        area_of_south_windows = (WINDOW_FRACTION_SOUTH_WALL * HOUSE_LENGTH_EAST_WEST * HEIGHT_TO_EAVES + 
                                 0.5 * ceiling_area)
        window_area = (WINDOW_FRACTION_OTHER_WALLS * 
                       (2 * HOUSE_WIDTH_NORTH_SOUTH * HEIGHT_TO_EAVES + 
                        HOUSE_LENGTH_EAST_WEST * HEIGHT_TO_EAVES))
        total_area = area_of_outer_walls + ceiling_area

        # Initialize stock variables
        fuel_oil_consumption_j = 0
        furnace_power_demand_j = 0
        
        # Time loop
        for t in sim_index:
            # Update outdoor temperature and solar heating gain
            month_series = t * 12
            outdoor_temperature = outdoor_temperature_tbf(month_series)
            solar_heating_gain = (area_of_south_windows * COLLECTOR_EFFICIENCY_SOUTH_WINDOWS * 
                                  insolation_on_vertical_wall_tbf(month_series))

            # Heat transfer calculations
            specif_heat_transfer_upper = (K_INSULATION_WALLS * (area_of_outer_walls - (window_area + area_of_south_windows)) +
                                          K_INSULATION_ROOF * ceiling_area + 
                                          K_INSULATION_WINDOWS * (window_area + area_of_south_windows))
            specif_heat_transfer_to_ground = K_INSULATION_FLOOR * floor_area
            heat_loss_to_ground = specif_heat_transfer_to_ground * (DESIRED_ROOM_TEMPERATURE - GROUND_TEMPERATURE)
            heating_limit_temperature = DESIRED_ROOM_TEMPERATURE + (heat_loss_to_ground - solar_heating_gain) / specif_heat_transfer_upper
            
            # Furnace power demand
            heating_power_demand = 0
            if heating_limit_temperature > outdoor_temperature:
                heating_power_demand = (specif_heat_transfer_upper * 
                                        (DESIRED_ROOM_TEMPERATURE - outdoor_temperature) + 
                                        heat_loss_to_ground - solar_heating_gain)
            furnace_power_demand_j = heating_power_demand / (EFFICIENCY_OF_FURNACE * 1000)
            burn_rate = (8760 / ENERGY_CONTENT_OF_FUEL) * furnace_power_demand_j
            fuel_oil_consumption_j += delta_t * burn_rate
        
        # Collect results for this scenario
        results.append([greenhouse_unit_number, DESIRED_ROOM_TEMPERATURE, 
                        area_of_outer_walls, ceiling_area, total_area, 
                        furnace_power_demand_j, fuel_oil_consumption_j])

# Convert results to DataFrame
columns = ["Greenhouse Unit Number", "Desired Room Temperature", 
           "Area of Outer Walls", "Ceiling Area", "Total Area", 
           "Final Furnace Power Demand (kWh)", "Final Fuel Oil Consumption (1000 m3)"]
result_df = pd.DataFrame(results, columns=columns)

# Export to CSV
output_path = r"C:\Users\shang\Desktop\2024 Fall\CEE 493 Sustainable Design\Simulation\Greenhouse python export"
os.makedirs(output_path, exist_ok=True)
result_df.to_csv(os.path.join(output_path, "simulation_results_EVA.csv"), index=False)

# Create a contour map
plt.figure(figsize=(10, 8))

# Create a meshgrid for the contour plot
x, y = np.meshgrid(greenhouse_unit_numbers, desired_room_temperatures)
z = result_df.pivot_table(
    index="Desired Room Temperature", 
    columns="Greenhouse Unit Number", 
    values="Final Furnace Power Demand (kWh)"
).values

# Plot the contour map
contour = plt.contourf(x, y, z, levels=20, cmap="viridis")
plt.colorbar(contour, label="Final Furnace Power Demand (kWh)")

# Label the axes
plt.xlabel("Greenhouse Unit Number")
plt.ylabel("Desired Room Temperature (°C)")
plt.title("Contour Map of Final Furnace Power Demand_EVA")

# Save the contour map
plt.savefig(os.path.join(output_path, "furnace_power_contour_map_EVA.png"))
plt.show()

# Create heatmap
heatmap_data = result_df.pivot_table(
    index="Desired Room Temperature", 
    columns="Greenhouse Unit Number", 
    values="Final Furnace Power Demand (kWh)"
)

plt.figure(figsize=(10, 8))
plt.imshow(heatmap_data, cmap="viridis", aspect="auto", origin="lower")
plt.colorbar(label="Final Furnace Power Demand (kWh)")
plt.xticks(ticks=np.arange(len(greenhouse_unit_numbers)), labels=greenhouse_unit_numbers, rotation=45)
plt.yticks(ticks=np.arange(len(desired_room_temperatures)), labels=desired_room_temperatures)
plt.xlabel("Greenhouse Unit Number")
plt.ylabel("Desired Room Temperature (°C)")
plt.title("Heatmap of Final Furnace Power Demand_EVA")
plt.savefig(os.path.join(output_path, "furnace_power_heatmap_EVA.png"))
plt.show()
