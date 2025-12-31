#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 26 13:30:16 2023
Last Modified: 9 Feb, 2024 1:30 PM

@author: michaelrobert
"""

import numpy as np
from scipy.integrate import odeint #This is the ODE Solver
import matplotlib.pyplot as plt # This is for plotting


# Define your ODE Model 
# y are state variables (dependent on time and parameters)
# t is for time
# k is for parameters (rates, proportions, etc.)
def model(y,t,k):  
    # State Variables
    S = y[0]
    I = y[1]
    R = y[2]
    
    N = S+I+R
    
    # Parameters
    beta = k[0]
    gamma = k[1]
    
    # Differential Equations
    dsdt = -beta*S*I/N
    didt = beta*S*I/N - gamma*I
    drdt = gamma*I
    
    
    # Compile and return    
    dydt = [dsdt, didt, drdt]
    return dydt

# Initial Conditions
S0 = 499
I0 = 1
R0 = 0

N0 = S0+I0+R0

# Store initial conditions in a vector
y0 = [S0, I0, R0]

# Parameter Values
beta = .3
gamma = 1/10.0

# store parameter values in a vector
p = [beta, gamma]

# Setup Time
tB = 0  #Beginning Time
tF = 100 # Final Time
nTime = 100 # time points

# Create a vector of time points from beginning time to end time that includes nTime number of points
tspan = np.linspace(tB,tF,nTime)

# Initialize Storage for State Variables
# This creates zero vectors of length equal to that of tspan
#S = np.empty_like(tspan)
#I = np.empty_like(tspan)
#R = np.empty_like(tspan)



# Simulate the Model
ode_out = odeint(model, y0, tspan, args=(p,))

# Extract Solutions - Note you don't have to do this, but it makes things more clear sometimes
S = ode_out[:,0]
I = ode_out[:,1]
R = ode_out[:,2]


# Plot Results for 1 simulation 
plt.plot(tspan, S, 'b-', label="Susceptible") # These three lines plot each of the solutions
plt.plot(tspan, I, 'r:', label="Infectious")  # Note you need the labels for the plot legend
plt.plot(tspan, R, 'k--', label="Recovered")
plt.ylabel('Population') # Vertical axis label
plt.xlabel('time') # Horizontal Axis label
plt.legend(loc='best') # Show legend and put it in the "best" location
plt.show() # Show PLot


##################################################################

# Simulate the Model with a different value of gamma to compare
gamma1 = 1/5.0
p = [beta, gamma1]
ode_out = odeint(model, y0, tspan, args=(p,))
S2 = ode_out[:,0]
I2 = ode_out[:,1]
R2 = ode_out[:,2]



# Alternatively Make a Multi-plot
fig, axs = plt.subplots(2,2) # Plots a 2x2 plot
fig.tight_layout() # This creates necessary space between plots


axs[0,0].plot(tspan, S, 'b') # Plot 1
axs[0,0].plot(tspan, S2, 'b:') # Plot 1
axs[0,0].set(ylabel="Susceptible", xlabel="time", ylim=(0,1.1*N0))

axs[0,1].plot(tspan, I, 'r', label='$\gamma = %s$' %gamma)  # Plot 2
axs[0,1].plot(tspan, I2, 'r:', label='$\gamma = %s$' %gamma1)  # Plot 2
axs[0,1].set(ylabel="Infectious", xlabel="time", ylim=(0,1.1*N0))
axs[0,1].legend(loc='best') # Show legend and put it in the "best" location

axs[1,0].plot(tspan, R, 'k') # Plot 3
axs[1,0].plot(tspan, R2, 'k:') # Plot 3
axs[1,0].set(ylabel="Recovered", xlabel="time", ylim=(0,1.1*N0))

axs[1,1].axis('off') # Just clears the 4th plot since we don't need it here. 
plt.show()





