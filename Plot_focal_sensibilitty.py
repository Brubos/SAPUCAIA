#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 10 09:37:03 2023

@author: bruno.souza
"""

#%% Packages
import pandas as pd                    # Data manipulation and analysis
import numpy as np                     # Numerical operations and array handling
import matplotlib.pyplot as plt        # Visualization and plotting
from scipy.optimize import curve_fit   # Curve fitting in data analysis
from sklearn.metrics import r2_score   # Calculation of R^2
import re                              # Regular expressions for pattern matching

# %% Functions

# Function for reading the data file
def read_data_file():
    # Reading the file using pandas
    df = pd.read_csv(file_name, usecols=[0, 3, 4], delim_whitespace=True, skiprows=1,
                     names=["misalignment", "x_mean_position", "y_mean_position"])

    # Accessing data in specific columns
    misalignment = df["misalignment"].tolist()
    x_mean_position = df["x_mean_position"].tolist()
    y_mean_position = df["y_mean_position"].tolist()

    return misalignment, x_mean_position, y_mean_position

# Function for linear fit
def linear_fit(x, alpha, beta):
    x_array1 = np.array(x)       # Converting x to a numpy array
    return alpha * x_array1 + beta

# Function for parabolic fit
def parabolic_fit(x, a, b, c):
    x_array2 = np.array(x)       # Converting x to a numpy array
    return a * x_array2**2 + b * x_array2 + c

    
#%% Read the data file
file_name = "SPU Ver_f Tx=180 µm"
# file_name = "SPU Ver_f Ty=2.0 mm"
# file_name = "SPU Ver_f Tz=10.0 mm"
# file_name = "SPU Ver_f Rx=450 µrad"
# file_name = "SPU Ver_f Ry=90 µrad"
# file_name = "SPU Ver_f Rz=50 mrad"

# Regular expression to find 'Rx', 'Ry', 'Rz', 'Tx', 'Ty', or 'Tz'
match_deg = re.search(r'(Rx|Ry|Rz|Tx|Ty|Tz)', file_name)
if match_deg:
    deg_f = match_deg.group()

# Regular expression to find 'mm', 'um', 'urad', or 'mrad'
match_unit = re.search(r'(mm|µm|µrad|mrad)', file_name)
if match_unit:
    unit = match_unit.group()

misalignment, x_mean_position, y_mean_position = read_data_file()


#%% Fit type

fit_type_x='linear'
# fit_type_x='parabolic'

fit_type_y='parabolic'
# fit_type_y='linear'

if fit_type_x == 'linear':
    params_x, _ = curve_fit(linear_fit, misalignment, x_mean_position)
    fit_x = linear_fit(np.array(misalignment), *params_x)
    label_x = f'mx+n \nm={params_x[0]:.4f}, n={params_x[1]:.4f}'
else: 
    params_x, _ = curve_fit(parabolic_fit, misalignment, x_mean_position)
    fit_x = parabolic_fit(np.array(misalignment), *params_x)
    label_x = f'ax²+bx+c \na={params_x[0]:.4f}, b={params_x[1]:.4f}, c={params_x[2]:.4f}'
    
if fit_type_y == 'linear':
    params_y, _ = curve_fit(linear_fit, misalignment, y_mean_position)
    fit_y = linear_fit(np.array(misalignment), *params_y)
    label_y = f'mx+n \nm={params_y[0]:.4f}, n={params_y[1]:.4f}'
else: 
    params_y, _ = curve_fit(parabolic_fit, misalignment, y_mean_position)
    fit_y = parabolic_fit(np.array(misalignment), *params_y)
    label_y = f'ax²+bx+c \na={params_y[0]:.4f}, b={params_y[1]:.4f}, c={params_y[2]:.4f}'
    

#%% Plotting the graphs
plt.figure(figsize=(12, 5))



# Graph of x_mean_position versus misalignment
plt.subplot(1, 2, 1)
plt.scatter(misalignment, x_mean_position, label='Data')
plt.plot(misalignment, fit_x, color='red',label=label_x)
plt.title('Focal Sensitivity')
plt.xlabel(deg_f + ' [' + unit + ']')
plt.ylabel('X mean position [µm]')
plt.legend()
plt.grid(True)

# Graph of y_mean_position versus misalignment
plt.subplot(1, 2, 2)
plt.scatter(misalignment, y_mean_position, label='Data')
plt.plot(misalignment, fit_y, color='red',label=label_y)
plt.title('Focal Sensitivity')
plt.xlabel(deg_f + ' [' + unit + ']')
plt.ylabel('Y mean position [µm]')
plt.legend()
plt.grid(True)

# Adjust the space between subplots
plt.subplots_adjust(wspace=0.4)  # You can adjust the value as needed

# Save the generated figure
plt.savefig('Focal Sensitivity ' + deg_f,dpi=1000)  # Replace 'file_name.png' with the desired file name

# Display the graphs
plt.tight_layout()
plt.show()