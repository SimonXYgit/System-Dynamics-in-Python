# -*- coding: utf-8 -*-
"""
Created on Sat Jan 27 11:36:26 2024

@author: shang
"""
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import time
# from matplotlib.animation import FuncAnimation
# from matplotlib.animation import PillowWriter
\

def print_with_delay(text, delay=0.01):
    for char in text:
        print(char, end='', flush=True)
        time.sleep(delay)
    print()

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

def visualization (display_until_year, deer_population, food, popchangerate, foodchangerate, FPDchangerate):
    # Filter for display
    filtered_deer_population = deer_population[deer_population.index <= display_until_year]
    filtered_food            = food           [food.index            <= display_until_year]

    # Plot for "deer-population + popchangerate + filtered_deer_population"
    fig, ax1 = plt.subplots()
    # Plot deer population on the left y-axis
    ax1.plot(deer_population.index, deer_population, color='b', label='Deer Population (Entire Simulation)')
    ax1.plot(filtered_deer_population.index, filtered_deer_population, color='r', label='Deer Population (Displayed Simulation)', linewidth=4)
    ax1.set_xlabel('Year')
    ax1.set_ylabel('Population Deer', color='b')
    ax1.tick_params('y', colors='b')
    ax1.set_title('KAIBAB: Year - Deer Population')
    # Create a second y-axis for the change rate
    ax2 = ax1.twinx()
    ax2.plot(popchangerate.index, popchangerate, color='r', linestyle='dashed', label='Change Rate of Population')
    ax2.set_ylabel('Change rate', color='r')
    ax2.tick_params('y', colors='r')
    # Display the legend
    fig.tight_layout()
    fig.legend(loc='upper right', bbox_to_anchor=(0.9, 0.9))
    # Set y-axis in scientific notation
    ax1.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
    ax2.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
    # Show the plot
    plt.show()

    # Plot for "food + foodchangerate"
    fig, ax1 = plt.subplots()
    # Plot food on the left y-axis
    ax1.plot(food.index, food, color='b', label='Food (Entire)')
    ax1.plot(filtered_food.index, filtered_food, color='r', label='Food (Displayed)', linewidth=4)
    ax1.set_xlabel('Year')
    ax1.set_ylabel('Food', color='b')
    ax1.tick_params('y', colors='b')
    ax1.set_title('KAIBAB: Year - Food')
    # Create a second y-axis for the food change rate
    ax2 = ax1.twinx()
    ax2.plot(foodchangerate.index, foodchangerate, color='r', linestyle='dashed', label='Change Rate of Food')
    ax2.set_ylabel('Change rate', color='r')
    ax2.tick_params('y', colors='r')
    # Display the legend
    fig.tight_layout()
    fig.legend(loc='upper right', bbox_to_anchor=(0.9, 0.9))
    # Set y-axis in scientific notation
    ax1.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    ax2.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    # Show the plot
    plt.show()

    # Plot for "food/deer_population + FPDchangerate"
    fig, ax3 = plt.subplots()
    # Plot food per deer on the left y-axis
    ax3.plot(deer_population.index, food/deer_population, color='b', label='Food per Deer (Entire)')
    ax3.plot(filtered_deer_population.index, filtered_food/filtered_deer_population, color='r', label='Food per Deer (Displayed)', linewidth=4)
    ax3.set_xlabel('Year')
    ax3.set_ylabel('Food per Deer', color='g')
    ax3.tick_params('y', colors='g')
    ax3.set_title('KAIBAB: Year - Food per Deer')
    # Create a second y-axis for the food per deer change rate
    ax4 = ax3.twinx()
    ax4.plot(FPDchangerate.index, FPDchangerate, color='m', linestyle='dashed', label='Change Rate of Food per Deer')
    ax4.set_ylabel('Change rate', color='m')
    ax4.tick_params('y', colors='m')
    # Display the legend
    fig.tight_layout()
    fig.legend(loc='upper right', bbox_to_anchor=(0.9, 0.9))
    # Set y-axis in scientific notation
    ax3.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    ax4.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    # Show the plot
    plt.show()
    return ()

### Initialization
## Base Parameters
sim_start = 1900 #(Z. 30)
sim_end = 1950
delta_t = 0.5
sim_length = sim_end - sim_start + delta_t # inklusive start and end
unit = "Years"
iterations = int(sim_length/delta_t)

## Constant Simulation Parameters
FMPD = 2000 # max. feed intake rate [kcal/(deer*d)]
MFOOD = 4.8E8 # max. feed capacity [kcal*a/d]
AREA = 800000 # area [acres]

## Ordinary Table Functions
# Feed rate per deer (Zeile 1150)
fpd_x = [0, 500, 1500, 2000, 1000000] # feed provision [kcal/(deer*d)]
fpd_y = [-0.5, -0.15, 0.15, 0.2, 0.2] # rel. net growth rate of deer [1/a]
fpd_tbf = interp1d(fpd_x, fpd_y,'linear') # As fpd goes up, growth rate approaches 0.2
# deer kills per predator (Zeile 1180, 1190)
deer_density = [0, 0.005, 0.01, 0.02, 0.025, 10] # [deer/acre]
dkpp_y = [0, 3, 13, 51, 56, 56] # prey rate per predator [deer/(predator*a)]
dkpp_tbf = interp1d(deer_density, dkpp_y, 'linear') # As deer_density goes up, dkpp aproaches 56
# time feed regeneration (Zeile 1210, 1220)
tfr_x = [0, 0.25, 0.5, 0.75, 1] # FOOD/MFOOD=avaiFeed/max.avaiFeed
tfr_y = [35, 15, 5, 1.5, 1] # regene time of feed [a]
tfr_tbf = interp1d(tfr_x, tfr_y,'linear') # As feed amount appoches its max.Capa, reg.Time goes down

## Create simulation timeline
sim_series = pd.Series(range(iterations)) # Range is exclusive of end
sim_series = pd.Series(range(iterations), index=sim_series*0.5+1900)
## Initialize integrators
deer_population = pd.Series(index=sim_series.index, dtype='O')
food = pd.Series(index=sim_series.index, dtype='O')
popchangerate = pd.Series(index=sim_series.index, dtype='O')
foodchangerate = pd.Series(index=sim_series.index, dtype='O')
FPDchangerate = pd.Series(index=sim_series.index, dtype='O')
deer_population.iloc[0] = 4000 
food.iloc[0] = 4.7E8 
popchangerate.iloc[0] = 0 
FPDchangerate.iloc[0] = 0
FPDchangerate.iloc[0] = 0

### Game Section
## Background information
#print_with_delay("Before 1900, the population of deer on the Kaibab Plateau was relatively stable. However, due to concerns about overgrazing and competition with livestock, predator control measures were implemented. The purpose was to increase the deer population for hunting purposes to increase the income of local companies")
print_with_delay("Deer population in 1900: 4000")
print_with_delay("Predator population in 1900: 266")
print_with_delay("Target: Reduce predator numbers to increase income from deers.")
## Input parameters and carry out simulation
# Until 1905
print_with_delay("------------------------------------(1905)------------------------------------")
ud_1 = int(input("Enter the number of predators you want to keep until 1905: "))
predator_population = {1900 : 266, 1905 : ud_1, 1920 : ud_1}
predator_population = table_function_series(predator_population, sim_series.index)
for i in sim_series.index: # Vorletzter Eintrag befüllt letzten Eintrag
   food_dp = food[i]/deer_population[i]
   NIR = fpd_tbf(food_dp)
   deer_density_i = deer_population[i] / AREA
   DNI = deer_population[i] * NIR
   DRP = predator_population[i]*dkpp_tbf(deer_density_i)
   FRR = (MFOOD-food[i])/tfr_tbf(food[i]/MFOOD)
   FCR = min(food[i], (deer_population[i]*FMPD))
   deer_population[i+0.5] = deer_population[i] + delta_t * (DNI-DRP)
   food[i+0.5] = food[i] + delta_t * (FRR-FCR)
   popchangerate[i] = (deer_population[i+0.5]-deer_population[i])/delta_t
   foodchangerate[i] = (food[i+0.5]-food[i])/delta_t
   FPDchangerate[i] = (food[i+0.5]/deer_population[i+0.5]-food[i]/deer_population[i])/delta_t
display_until_year = 1910
visualization (display_until_year, deer_population, food, popchangerate, foodchangerate, FPDchangerate)

# Until 1910
print("In 1905, you've ordered to keep ", ud_1, "predators alive.")
print_with_delay("------------------------------------(1910)------------------------------------")
ud_2 = int(input("Enter the number of predators you want to keep until 1910: "))
predator_population = {1900 : 266, 1905 : ud_1, 1910 : ud_2, 1920 : ud_2}
predator_population = table_function_series(predator_population, sim_series.index)
for i in sim_series.index: # Vorletzter Eintrag befüllt letzten Eintrag
   food_dp = food[i]/deer_population[i]
   NIR = fpd_tbf(food_dp)
   deer_density_i = deer_population[i] / AREA
   DNI = deer_population[i] * NIR
   DRP = predator_population[i]*dkpp_tbf(deer_density_i)
   FRR = (MFOOD-food[i])/tfr_tbf(food[i]/MFOOD)
   FCR = min(food[i], (deer_population[i]*FMPD))
   deer_population[i+0.5] = deer_population[i] + delta_t * (DNI-DRP)
   food[i+0.5] = food[i] + delta_t * (FRR-FCR)
   popchangerate[i] = (deer_population[i+0.5]-deer_population[i])/delta_t
   foodchangerate[i] = (food[i+0.5]-food[i])/delta_t
   FPDchangerate[i] = (food[i+0.5]/deer_population[i+0.5]-food[i]/deer_population[i])/delta_t
display_until_year = 1915
visualization (display_until_year, deer_population, food, popchangerate, foodchangerate, FPDchangerate)

# Until 1915
print("In 1905, you've ordered to keep ", ud_1, "predators alive.")
print("In 1910, you've ordered to keep ", ud_2, "predators alive.")
print_with_delay("------------------------------------(1915)------------------------------------")
ud_3 = int(input("Enter the number of predators you want to keep until 1915: "))
predator_population = {1900 : 266, 1905 : ud_1, 1910 : ud_2, 1915 : ud_3, 1920 : ud_3}
predator_population = table_function_series(predator_population, sim_series.index)
for i in sim_series.index: # Vorletzter Eintrag befüllt letzten Eintrag
   food_dp = food[i]/deer_population[i]
   NIR = fpd_tbf(food_dp)
   deer_density_i = deer_population[i] / AREA
   DNI = deer_population[i] * NIR
   DRP = predator_population[i]*dkpp_tbf(deer_density_i)
   FRR = (MFOOD-food[i])/tfr_tbf(food[i]/MFOOD)
   FCR = min(food[i], (deer_population[i]*FMPD))
   deer_population[i+0.5] = deer_population[i] + delta_t * (DNI-DRP)
   food[i+0.5] = food[i] + delta_t * (FRR-FCR)
   popchangerate[i] = (deer_population[i+0.5]-deer_population[i])/delta_t
   foodchangerate[i] = (food[i+0.5]-food[i])/delta_t
   FPDchangerate[i] = (food[i+0.5]/deer_population[i+0.5]-food[i]/deer_population[i])/delta_t
display_until_year = 1920
visualization (display_until_year, deer_population, food, popchangerate, foodchangerate, FPDchangerate)

# Until 1920
print("In 1905, you've ordered to keep ", ud_1, "predators alive.")
print("In 1910, you've ordered to keep ", ud_2, "predators alive.")
print("In 1910, you've ordered to keep ", ud_3, "predators alive.")
print_with_delay("------------------------------------(1920)------------------------------------")
ud_4 = int(input("Enter the number of predators you want to keep until 1920: "))
predator_population = {1900 : 266, 1905 : ud_1, 1910 : ud_2, 1915 : ud_3, 1920 : ud_4}
predator_population = table_function_series(predator_population, sim_series.index)
for i in sim_series.index: # Vorletzter Eintrag befüllt letzten Eintrag
   food_dp = food[i]/deer_population[i]
   NIR = fpd_tbf(food_dp)
   deer_density_i = deer_population[i] / AREA
   DNI = deer_population[i] * NIR
   DRP = predator_population[i]*dkpp_tbf(deer_density_i)
   FRR = (MFOOD-food[i])/tfr_tbf(food[i]/MFOOD)
   FCR = min(food[i], (deer_population[i]*FMPD))
   deer_population[i+0.5] = deer_population[i] + delta_t * (DNI-DRP)
   food[i+0.5] = food[i] + delta_t * (FRR-FCR)
   popchangerate[i] = (deer_population[i+0.5]-deer_population[i])/delta_t
   foodchangerate[i] = (food[i+0.5]-food[i])/delta_t
   FPDchangerate[i] = (food[i+0.5]/deer_population[i+0.5]-food[i]/deer_population[i])/delta_t
display_until_year = 1925
visualization (display_until_year, deer_population, food, popchangerate, foodchangerate, FPDchangerate)



