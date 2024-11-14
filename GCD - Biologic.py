# -*- coding: utf-8 -*-
"""
Created on Thu Apr 18 13:08:19 2024

@author: ok723

"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import glob
from natsort import natsorted
import os
from scipy.integrate import trapz
from matplotlib import cm
import cmocean

# My colours
# mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=['#c3121e', '#0348a1', '#ffb01c', '#027608', '#0193b0', '#9c5300', '#949c01', '#7104b5'])
#                                                      0sangre,   1neptune,  2pumpkin,  3clover,   4denim,    5cocoa,    6cumin,    7berry3
#Some colours
# mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=['#e41a1c','#377eb8', '#ff7f00', '#4daf4a','#f781bf', '#a65628', '#984ea3','#999999', '#dede00'])
# Nature colours
# mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=['0C5DA5', 'FF2C00', 'FF9500', '00B945', '845B97', '9c5300','474747', '9e9e9e','26d8f4','800080','#9c5300'])

# Font styles
plt.rcParams['axes.axisbelow'] = True

# plt.rcParams.update({
#     'font.size': 10,
#     'font.family': 'STIXGeneral',
#     'mathtext.fontset': 'stix' })

plt.rcParams.update({
    "font.family": "serif",   # specify font family here
    "font.serif": ["Garamond"],  # specify font here
    "font.size":11})          # specify font size here

plt.rcParams['axes.linewidth'] =0.5
plt.rcParams['grid.linewidth'] =0.5
plt.rcParams['lines.linewidth'] =1

# =============================================================================
# Creates a folder to store the graphs in the cwd
# cwd=os.getcwd()
# save_path = 'C:/Users/ok723/Desktop/GCD Test'
# if not os.path.exists(save_path):
#     os.makedirs(save_path)
# =============================================================================
# Mass of whatever you're normalising by, in mg, if you dont have mass in the datafiles.
# Make sure to comment this line out though: mass=extract_mass(file_path)
# mass=195.2
# =============================================================================

# Path files' folder.
path = 'C:/Users/ok723/Downloads/'

# Read files with *X* part in their name. Single * means read everything in the folder.
file_list = glob.glob(path + "*12_11_DEV4_MONO_02_GCPL.txt*")
file_list = natsorted(file_list)

# Colormap of your preference
colormap=cm.viridis

#Set the dictionaries to save data as global variables.
cycle_data_dict = {}
voltage_drops_dict = {}
specific_capacitance_dict = {}
specific_energy_dict = {}
specific_power_dict = {}
coulombic_efficiency_dict={}
energy_efficiency_dict = {}

av_cur_dict={}
# Not really needed unless you want to make sure which values are fitted and the slope values (E). Take a look later in the script at next 'y_stoixeia' occurence.
y_stoixeia=[]
E=[]

#Extracting current, Vmax, and mass
def extract_current(file_path):
    with open(file_path, 'r', encoding='latin1') as file:
        for line in file:
            if line.strip().startswith('Is'):
                parts = line.split()
                if len(parts) >= 3:
                    return float(parts[2])  #Extract the second number
    return None  #Return None if the line or number is not found

def extract_Vmax_value(file_path):
    with open(file_path, 'r', encoding='latin1') as file:
        for line in file:
            if line.strip().startswith("EM (V)"):
                parts = line.split()
                if len(parts) >= 4:
                    return float(parts[3])  #Extract the third character. (V) counts as one too for some reason.
    return None  #Return None if "EM (V)" is not found

def extract_mass(file_path):
    with open(file_path, "r", encoding="latin1") as file:
        for line in file:
            if "Characteristic mass" in line and "mg" in line:
                parts = line.split()
                for i, part in enumerate(parts):
                    if part == "mg" and i > 0:  #Ensure there is a value before 'mg'
                        try:
                            return float(parts[i-1])  #Convert the value before 'mg' to float
                        except ValueError:
                            return None  #Handle any conversion errors
    return None  #Return None if no 'Characteristic mass' is found

def voltage_drop(file_path):
    with open(file_path, 'r', encoding='latin1') as file:
        for line in file:
            if 'Nb header lines' in line:
                nb_header_lines = int(line.split(':')[-1].strip()) # Extract the number of header lines to skip
                break

    # Subtract 1 from the number of lines to skip. You need the headers of the columns.
    skip_rows = nb_header_lines - 1

    # Read the files
    df = pd.read_csv(file_path, delimiter='\t', skiprows=skip_rows, encoding='latin1')
    
    file_name = os.path.basename(file_path)
    num_cycles = len(df['cycle number'].unique())
    
    # Figure size and setting colors' range all across the colormap's spectrum
    fig = plt.figure(figsize=(3.5, 2.625), dpi=700)
    base_colors = [colormap(i / (num_cycles - 1)) for i in range(num_cycles)] if num_cycles > 1 else [colormap(0.5)]
    
    #Lists to pass to the global dictionaries
    voltage_drops_list = []
    specific_capacitance_list=[]
    specific_energy_list=[]
    specific_power_list=[]
    coulombic_efficiency_list=[]
    energy_efficiency_list=[]
    charge_list=[]
    discharge_list=[]
    av_cur_list=[]
    
    #Print datafiles' names
    print(f'Datafile name: {file_name}')
    
    #Analysis and plotting
    for cycle, color in zip(sorted(df['cycle number'].unique()), base_colors):
        cycle_data_original = df[df['cycle number'] == cycle]
        cycle_data = cycle_data_original.copy().reset_index(drop=True) #200 IQ.
        
        #Prints the original and reseted cycle data size, the Nth point of the cycle, and the entire cycle data respectively in case you need to validate anything.
        # print(f"Cycle {cycle}, Original Data Size: {cycle_data_original.shape}, Reset Data Size: {cycle_data.shape}")
        # print(cycle_data['Ewe/V'][0])
        # print(cycle_data)
        
        #Get current,V_max, mass and the transition points of charge and discharge for each cycle
        current=extract_current(file_path)
        v_max=extract_Vmax_value(file_path)
        mass=extract_mass(file_path)
        transition_points = cycle_data[cycle_data['ox/red'].diff() == -1]
        # print(transition_points)      
        # print(mass)
        # print(current)
        # print(v_max)
        
        #!!!!!!!!!!!!! Not accurate to use peak index as a method to get the transition points. It's better to use the current's 'ox/red' transition point approach above.
        # charge_data = cycle_data.iloc[:peak_index].copy()
        # discharge_data = cycle_data.iloc[peak_index:].copy()
        
        for index in transition_points.index:
            if index > 0 and (index + 1) in cycle_data.index:
                voltage_before = cycle_data.at[index - 1, 'Ewe/V'] #Voltage before uses -1
                voltage_after = cycle_data.at[index+1, 'Ewe/V'] #Voltage after uses +1 just in case there is a 'bigger' drop. Could potentially quantify what a bigger drop is.
                voltage_drop = voltage_before - voltage_after
                voltage_drops_list.append({'Cycle number': cycle, 'Voltage drop (V)': abs(voltage_drop),'ESR (Ω)':abs(voltage_drop*(10**3))/(2*current),'Voltage after (V)': abs(voltage_after)})
                voltage_drops_dict[file_name] = pd.DataFrame(voltage_drops_list)
                
                # Prints the voltage value after the drop and some random value that I can’t recall why I included here
                # print(cycle_data['Ewe/V'][index])
                # print(cycle_data.index)
                # print(voltage_drop)
                    
                # Plot the area under the curve to ensure everything runs properly.
                # discharge_data = cycle_data.iloc[index+1:].copy()
                # charge_data = cycle_data.iloc[:index-1].copy()   
                # fig = plt.figure(figsize=(3.5, 2.625),dpi=700)
                # plt.plot(discharge_data['time/s']-discharge_data['time/s'].iloc[0],discharge_data['Ewe/V'],color='#c3121e')
                # plt.fill_between(discharge_data['time/s']-discharge_data['time/s'].iloc[0],discharge_data['Ewe/V'])

        cycle_key = f'{os.path.splitext(file_name)[0]}_Cycle_{int(cycle)}'
        cycle_data_dict[cycle_key] = cycle_data
        
        # Need to make the dataframes into numpy arrays. [index+1] will fit starting at the point after the drop. Adjust accordingly.
        x=cycle_data['time/s']-cycle_data['time/s'].iloc[0]
        y=cycle_data['Ewe/V']
        x_data=x.iloc[index+1:].copy()
        y_data=cycle_data['Ewe/V'].iloc[index+1:].copy().to_numpy()
        
        #Not really needed unless you want to make sure which values are fitted. Take a look later in the script at next 'y_stoixeia' occurence.
        y_stoixeia.append(y_data)
        
        # fig = plt.figure(figsize=(3.5, 2.625),dpi=700) #Use if you want to plot each cycle in defferent figures (will need to move the figures style such as ticks etc a block to the right)
        plt.plot(cycle_data['time/s']-cycle_data['time/s'].iloc[0],cycle_data['Ewe/V'], label=f'Cycle {int(cycle)}', color=color,marker='o', markersize=2, mfc='w',linewidth=0.5)

# =============================================================================
#             # Capacitance linear fitt approach A
# =============================================================================
        # A = np.vstack([x_data, np.ones(len(x_data))]).T
        # m,c=np.linalg.lstsq(A,y_data,rcond=None)[0]
        # slope=round(m,5)
        # E.append(slope)# Not reallly needed, only here to view the slopes as a global variable.
        # plt.plot(x_data, m*x_data + c, 'r', label='Fitted line')
        # ax=plt.gca()
        # ax.text(0.8,0.1,f"Slope={slope}",va='center',ha='center',transform=ax.transAxes,fontsize=9, fontweight='bold')
# =============================================================================
#             # Capacitance linear fitt approach B
# =============================================================================
        # slope, intercept = np.polyfit(x_data,y_data, 1)
        # plt.plot(x,y, label=f'Cycle {int(cycle)}', color=color,marker='o', markersize=2, mfc='w',linewidth=0.5)
        # fitted_line = slope * x_data + intercept  # Calculate the fitted line
        # plt.plot(x_data,fitted_line, label='Fitted Line', linestyle='--')  # Plot the fitted line
        # plt.plot(x_data, y_data, label='Original Discharge Data')  # Plot original discharge data - u     seless
# ============================================================================= 
        # Chrge and discharge data using a point before and after the voltage drop    
        discharge_data = cycle_data.iloc[index+1:].copy()
        charge_data = cycle_data.iloc[:index-1].copy()      
        # print(charge_data)
        # print(discharge_data)
        
        # Area under the curve (energy) for charge and discharge.
        charge_auc = trapz(abs(charge_data['Ewe/V']), charge_data['time/s'])
        discharge_auc = trapz(abs(discharge_data['Ewe/V']), discharge_data['time/s'])
        # print(charge_auc)
        # print(discharge_auc)
        
        # Duration for charge and discharge
        charge_time_span = charge_data['time/s'].max() - charge_data['time/s'].min()
        discharge_time_span = discharge_data['time/s'].max() - discharge_data['time/s'].min()
        # print(f"Cycle {cycle}, Charge Time Span: {charge_time_span} seconds")
        # print(f"Cycle {cycle}, Discharge Time Span: {discharge_time_span} seconds"
# =============================================================================
#             # Discharge capacitance, energy, and power
# =============================================================================
        specific_energy=(discharge_auc*current)/(3.6*mass)
        specific_power=(discharge_auc*current*(10**3))/(mass*discharge_time_span)
        
        # specific_capacitance=current/(-slope*mass)
        specific_capacitance=(2*discharge_auc*current)/((voltage_after**2)*mass)

        #Coulombic and energy efficiency (need charge and discharge energy values for the latter one)
        charge_energy=(charge_auc*current*(10**3))
        discharge_energy=(discharge_auc*current*(10**3))
        energy_efficiency=discharge_energy/charge_energy
        coulombic_efficiency=discharge_time_span/charge_time_span
        
        # Testing stuff on how to calculated the average charge/discharge current.
        # av_ch_cur=charge_data['<I>/mA'].iloc[1:-1].mean()
        # av_dis_cur=discharge_data['<I>/mA'].iloc[1:-1].mean()

        # av_cur_list.append({'Cycle number': cycle, 'Average charge current': abs(av_ch_cur),'Cycle number': cycle, 'Average discharge current': abs(av_dis_cur)})
        # av_cur_dict[file_name] = pd.DataFrame(av_cur_list)
# =============================================================================
        #Passing specific capacitance, energy, and power values in the respective dictionaries.
        specific_capacitance_list.append({'Cycle number': cycle, 'Specific capacitance (F/g)': abs(specific_capacitance)})
        specific_capacitance_dict[file_name] = pd.DataFrame(specific_capacitance_list)
        
        specific_energy_list.append({'Cycle number': cycle, 'Specific energy (Wh/kg)': abs(specific_energy)})
        specific_energy_dict[file_name] = pd.DataFrame(specific_energy_list)
        
        specific_power_list.append({'Cycle number': cycle, 'Specific power (W/kg)': abs(specific_power)})
        specific_power_dict[file_name] = pd.DataFrame(specific_power_list)
        
        energy_efficiency_list.append({'Cycle number': cycle, 'Energy efficiency (%)':energy_efficiency})
        energy_efficiency_dict[file_name] = pd.DataFrame(energy_efficiency_list)
        
        coulombic_efficiency_list.append({'Cycle number': cycle, 'Coulombic efficiency (%)':coulombic_efficiency})
        coulombic_efficiency_dict[file_name] = pd.DataFrame(coulombic_efficiency_list)
# =============================================================================
        # Printing the capacitance per file per cycle. Irrelevant after addition of the next code cells.
        # print(f'Cycle {cycle}\nDischarge capacitance: {specific_capacitance:.3f} F/g \n\n')
# =============================================================================
        # Outdated.
        # with open(os.path.join(save_path, '1_data.txt'), 'a') as file:
        #                 file.write(f'{file_name} \n'f'Cycle {cycle} \n'
        #                         f'Discharge capacitance: {capacitance:.3f} F/g \n\n')
# =============================================================================
# Outdated as well.
    #     previous_filename = None         
    #     with open(os.path.join(save_path, '1_data.txt'), 'a') as file:
    #         if file_name != previous_filename:
    #             file.write(f'{file_name} \n')
    #             previous_filename = file_name

    #         file.write(f'Cycle {cycle} \n'
    #             f'Discharge capacitance: {specific_capacitance:.3f} F/g \n\n')
        
    #     discharge_data['Ewe/V'], discharge_data['time/s']
    # file_name = f'{os.path.splitext(file_name)[0]}' #Splits filane in name + format (eg .txt)
    # print(file_name)
# =============================================================================
    # Figure customisation
    plt.minorticks_on()
    plt.tick_params(direction='in', which='major', bottom=True, top=True, left=True, right=True,width=0.5,size=3)
    plt.tick_params(direction='in', which='minor', bottom=True, top=True, left=True, right=True,width=0.5,size=1.5)
    plt.grid(True, linestyle='dashed', linewidth='0.3', color='grey', alpha=0.8)
    plt.xlabel('Time (s)')
    plt.ylabel('Potential (V)')
    plt.legend(loc='best', frameon=False,fontsize=8)
    # plt.close()
    
    # Saving the figure
    # save_filename = os.path.join(save_path, f"{os.path.splitext(file_name)[0]}.png")
    # plt.savefig(save_filename,bbox_inches='tight', format='png', dpi=300)
    # plt.close(fig)  # Close the figure to free up memory
    
    # Saving the figure from DSC script
    # plt.savefig(mesa+'DSC plot above.png', dpi=900,bbox_inches="tight")
    # save_filename = os.path.join(save_path, f"{os.path.splitext(os.path.basename(file_path))[0]}_cycle_{cycle_number}.png")    
    return cycle_data_dict  # Return the dictionary containing DataFrames for each cycle

# Moment of truth
for file_path in file_list:
    cycle_data_dict = voltage_drop(file_path)
    
    
#Print run time
# end = time.time()
# total_time = end - start
# print("\n"+ 'Total run time: '+str(total_time))

#%% Plotting Nth cycle
def plot_specific_cycle(file_list, cycle_number):
    num_files = len(file_list)
    base_colors = [colormap(i / (num_files - 1)) for i in range(num_files)] if num_files > 1 else [colormap(0.5)]  # Use a mid-point for single file case
    # base_colors = [colormap(i / (num_files - 1)) for i in range(num_files)] if num_files > 1 else [colormap(0)]

    fig = plt.figure(figsize=(3.5, 2.625),dpi=700)
    for index, file_path in enumerate(file_list):
        with open(file_path, 'r', encoding='latin1') as file:
            for line in file:
                if 'Nb header lines' in line:
                    nb_header_lines = int(line.split(':')[-1].strip())
                    break
                
        skip_rows = nb_header_lines - 1
        df = pd.read_csv(file_path, delimiter='\t', skiprows=skip_rows, encoding='latin1')
        cycle_data = df[df['cycle number'] == cycle_number]
        if not cycle_data.empty:
            color = base_colors[index]
            # plt.plot(cycle_data['time/s'], cycle_data['Ewe/V'], label=os.path.basename(file_path), color=color) #Plots at time. Looks ugly. Could be useful though.
            plt.plot(cycle_data['time/s'] - cycle_data['time/s'].iloc[0], cycle_data['Ewe/V'], label = os.path.splitext(os.path.basename(file_path))[0], color=color,marker='o', markersize=2, mfc='w', linestyle='-')

    plt.minorticks_on()
    plt.tick_params(direction='in', which='major', bottom=True, top=True, left=True, right=True,width=0.5,size=3)
    plt.tick_params(direction='in', which='minor', bottom=True, top=True, left=True, right=True,width=0.5,size=1.5)
    plt.grid(True, linestyle='dashed', linewidth='0.3', color='grey', alpha=0.8)
    plt.xlabel('Time (s)')    
    plt.ylabel('Potential (V)')
    plt.legend(loc='best', frameon=False, fontsize=6)
    plt.ylim(-0.1,)

for i in range (5,6):
    plot_specific_cycle(file_list, i)
    plt.title(f'Cycle: {i}',fontsize=10)
#%% Specific capacitance
def plot_specific_capacitance(specific_capacitance_dict):
    fig = plt.figure(figsize=(3.5, 2.625),dpi=700)

    num_files = len(file_list)
    colors = [colormap(i / (num_files - 1)) for i in range(num_files)] if num_files > 1 else [colormap(0.5)]  # Use a mid-point for single file case
    for (file_name, df), color in zip(specific_capacitance_dict.items(), colors):
        plt.plot(df['Cycle number'], df['Specific capacitance (F/g)'], marker='o', markersize=2, mfc='w', linestyle='-', label=os.path.splitext(file_name)[0], color=color)
        
    plt.minorticks_on()
    plt.tick_params(direction='in', which='major', bottom=True, top=True, left=True, right=True,width=0.5,size=3)
    plt.tick_params(direction='in', which='minor', bottom=True, top=True, left=True, right=True,width=0.5,size=1.5)
    plt.grid(True, linestyle='dashed', linewidth='0.3', color='grey', alpha=0.8)
    plt.xlabel('Cycle number')
    plt.ylabel('Specific capacitance (F/g)')
    # plt.legend(loc='best', frameon=False, fontsize=6)

plot_specific_capacitance(specific_capacitance_dict)
#%% Ohmic drop/ESR
def plot_voltage_drops(voltage_drops_dict):
    fig = plt.figure(figsize=(3.5, 2.625),dpi=700)

    num_files = len(file_list)
    colors = [colormap(i / (num_files - 1)) for i in range(num_files)] if num_files > 1 else [colormap(0.5)]  # Use a mid-point for single file case
    for (file_name, df), color in zip(voltage_drops_dict.items(), colors):
        plt.plot(df['Cycle number'], df['Voltage drop (V)'], marker='o', markersize=2, mfc='w', linestyle='-', label=os.path.splitext(file_name)[0], color=color)
        
    plt.minorticks_on()
    plt.tick_params(direction='in', which='major', bottom=True, top=True, left=True, right=True,width=0.5,size=3)
    plt.tick_params(direction='in', which='minor', bottom=True, top=True, left=True, right=True,width=0.5,size=1.5)
    plt.grid(True, linestyle='dashed', linewidth='0.3', color='grey', alpha=0.8)
    plt.xlabel('Cycle number')
    plt.ylabel('Voltage drop (V)')
    # plt.legend(loc='best', frameon=False, fontsize=6)

plot_voltage_drops(voltage_drops_dict)
#%% ESR
def plot_ESR(voltage_drops_dict):
    fig = plt.figure(figsize=(3.5, 2.625),dpi=700)

    num_files = len(file_list)
    colors = [colormap(i / (num_files - 1)) for i in range(num_files)] if num_files > 1 else [colormap(0.5)]  # Use a mid-point for single file case
    for (file_name, df), color in zip(voltage_drops_dict.items(), colors):
        plt.plot(df['Cycle number'], df['ESR (Ω)'], marker='o', markersize=2, mfc='w', linestyle='-', label=os.path.splitext(file_name)[0], color=color)
        # ESR (Ω)
        #Voltage drop (V)
        
    plt.minorticks_on()
    plt.tick_params(direction='in', which='major', bottom=True, top=True, left=True, right=True,width=0.5,size=3)
    plt.tick_params(direction='in', which='minor', bottom=True, top=True, left=True, right=True,width=0.5,size=1.5)
    plt.grid(True, linestyle='dashed', linewidth='0.3', color='grey', alpha=0.8)
    plt.xlabel('Cycle number')
    plt.ylabel('ESR (Ω)')
    # plt.legend(loc='best', frameon=False, fontsize=6)

plot_ESR(voltage_drops_dict)

#%%
def ret(specific_capacitance_dict):
    fig = plt.figure(figsize=(3.5, 2.625),dpi=700)


    num_files = len(file_list)
    colors = [colormap(i / (num_files - 1)) for i in range(num_files)] if num_files > 1 else [colormap(0.5)]  # Use a mid-point for single file case
    for (file_name, df), color in zip(specific_capacitance_dict.items(), colors):
        plt.plot(df['Cycle number'],100* df['Specific capacitance (F/g)']/df['Specific capacitance (F/g)'].max(), marker='o', markersize=2, mfc='w', linestyle='-', label=os.path.splitext(file_name)[0], color=color)

    plt.minorticks_on()
    plt.tick_params(direction='in', which='major', bottom=True, top=True, left=True, right=True,width=0.5,size=3)
    plt.tick_params(direction='in', which='minor', bottom=True, top=True, left=True, right=True,width=0.5,size=1.5)
    plt.grid(True, linestyle='dashed', linewidth='0.3', color='grey', alpha=0.8)
    plt.xlabel('Cycle number')
    plt.ylabel('Capacitance retention (%)')
    plt.legend(loc='best', frameon=False, fontsize=6)

ret(specific_capacitance_dict)
#%%
def energy_efficiency(efficiency_dict):
    fig = plt.figure(figsize=(3.5, 2.625),dpi=700)

    num_files = len(file_list)
    colors = [colormap(i / (num_files - 1)) for i in range(num_files)] if num_files > 1 else [colormap(0.5)]  # Use a mid-point for single file case
    for (file_name, df), color in zip(efficiency_dict.items(), colors):
        plt.plot(df['Cycle number'],100*df['Energy efficiency (%)'], marker='o', markersize=2, mfc='w', linestyle='-', label=os.path.splitext(file_name)[0], color=color)

    plt.minorticks_on()
    plt.tick_params(direction='in', which='major', bottom=True, top=True, left=True, right=True,width=0.5,size=3)
    plt.tick_params(direction='in', which='minor', bottom=True, top=True, left=True, right=True,width=0.5,size=1.5)
    plt.grid(True, linestyle='dashed', linewidth='0.3', color='grey', alpha=0.8)
    plt.xlabel('Cycle number')
    plt.ylabel('Energy efficiency (%)')
    plt.legend(loc='best', frameon=False, fontsize=6)

energy_efficiency(energy_efficiency_dict)
#%%
def coulombic_efficiency(efficiency_dict):
    fig = plt.figure(figsize=(3.5, 2.625),dpi=700)


    num_files = len(file_list)
    colors = [colormap(i / (num_files - 1)) for i in range(num_files)] if num_files > 1 else [colormap(0.5)]  # Use a mid-point for single file case
    for (file_name, df), color in zip(efficiency_dict.items(), colors):
        plt.plot(df['Cycle number'],100*df['Coulombic efficiency (%)'], marker='o', markersize=2, mfc='w', linestyle='-', label=os.path.splitext(file_name)[0], color=color)

    plt.minorticks_on()
    plt.tick_params(direction='in', which='major', bottom=True, top=True, left=True, right=True,width=0.5,size=3)
    plt.tick_params(direction='in', which='minor', bottom=True, top=True, left=True, right=True,width=0.5,size=1.5)
    plt.grid(True, linestyle='dashed', linewidth='0.3', color='grey', alpha=0.8)
    plt.xlabel('Cycle number')
    plt.ylabel('Coulombic efficiency (%)')
    plt.legend(loc='best', frameon=False, fontsize=6)

coulombic_efficiency(coulombic_efficiency_dict)
#%%

mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=['#c3121e', '#0348a1', '#ffb01c', '#027608', '#0193b0', '#9c5300', '#949c01', '#7104b5'])
def plot_ragone(specific_energy_dict, specific_power_dict):
    fig = plt.figure(figsize=(3.5, 2.625),dpi=700)
    
    for file_name in specific_energy_dict:
        # Ensure the file exists in both dictionaries
        if file_name in specific_power_dict:
            energy_df = specific_energy_dict[file_name]
            power_df = specific_power_dict[file_name]

            # Merge dataframes based on 'Cycle number'
            merged_df = pd.merge(energy_df, power_df, on='Cycle number')

            plt.plot(merged_df['Specific energy (Wh/kg)'],merged_df['Specific power (W/kg)'],  marker='o', markersize=2, mfc='w', linestyle='-',label=os.path.splitext(file_name)[0])
    
    # # Add plot details
    plt.ylabel('Specific power (W/kg)')
    plt.xlabel('Specific energy (Wh/kg)')
    # ax.set_title('Specific Energy vs. Specific Power')
    plt.yscale('log')
    plt.xscale('log')
    plt.grid(True, linestyle='dashed', linewidth='0.3', color='grey', alpha=0.8)
    # ax.legend(loc='best', frameon=False, fontsize=6)
    plt.minorticks_on()
    plt.tick_params(direction='in', which='major', bottom=True, top=True, left=True, right=True, width=0.5, size=3)
    plt.tick_params(direction='in', which='minor', bottom=True, top=True, left=True, right=True, width=0.5, size=1.5)
    # plt.legend(loc='best', frameon=False)

plot_ragone(specific_energy_dict, specific_power_dict)