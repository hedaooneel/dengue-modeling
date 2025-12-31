import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

def wolbachia_enhanced_model(t, y, b, d, D, tau, q):
    I, U = y
    n = d + D  # Total death rate for infected individuals
    dI_dt = tau*b - n*I*(I+U)  # Adjusted for CI effect and fitness cost
    dU_dt = (1-tau)*b*I + (b*(1-q)*(I/(I+U))-d*(I+U))*U  # Uninfected with CI effect
    return [dI_dt, dU_dt]

# Parameters
b = 0.1  # Birth rate
d = 0.05  # Death rate of uninfected
D = 0.02  # Additional death rate due to infection
tau = 0.9  # Partial vertical transmission rate
q = 0.05  # Probability of cytoplasmic incompatibility

# Initial conditions
I0 = 1 # The inittial number of infected individuals
U0 = 1 # The initial number of uninfected individuals
initial_conditions = [I0, U0]

# Time span
t_span = [0, 100]
t_eval = np.linspace(*t_span, 200)

# Solve the ODE
sol = solve_ivp(wolbachia_enhanced_model, t_span, initial_conditions, args=(b, d, D, tau, q), t_eval=t_eval) # System to solve ODE, differential equation model
print (sol)

# Plot the results
plt.plot(sol.t, sol.y[0], label='Infected') # Plot infected individuals
plt.plot(sol.t, sol.y[1], label='Uninfected') # Plot uninfected individuals
plt.xlabel('Time') # Time as the independent variable
plt.ylabel('Population') # Population as the dependent variable
plt.legend()
plt.title('Enhanced Dynamics of Wolbachia Infection')
plt.show()
