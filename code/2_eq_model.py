import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

# Model function
def wolbachia_enhanced_model(t, y, b, d, D, tau, q):
    I, U = y
    n = d + D  # Total death rate for infected individuals
    dI_dt = tau * b - n * I * (I + U)  # Adjusted for CI effect and fitness cost
    dU_dt = (1 - tau) * b * I + (b * (1 - q) * (I / (I + U)) - d * (I + U)) * U  # Uninfected with CI effect
    return [dI_dt, dU_dt]

# Initial conditions and time span
initial_conditions = [0, 1]  # Initial numbers of infected and uninfected
t_span = [0, 100]
t_eval = np.linspace(*t_span, 200)

# Initial parameter values
initial_b = 0.1
initial_d = 0.05
initial_D = 0.02
initial_tau = 0.9
initial_q = 0.05


# Create the figure and axes
fig, ax = plt.subplots()
plt.subplots_adjust(left=0.1, bottom=0.35)

# Solve the ODE initially
sol = solve_ivp(wolbachia_enhanced_model, t_span, initial_conditions, args=(initial_b, initial_d, initial_D, initial_tau, initial_q), t_eval=t_eval)
lines = plt.plot(sol.t, sol.y[0], label='Infected'), plt.plot(sol.t, sol.y[1], label='Uninfected')

# Create sliders
axcolor = 'lightgoldenrodyellow'
sliders = {}
slider_positions = {'b': [0.25, 0.01, 0.65, 0.02], 'd': [0.25, 0.05, 0.65, 0.02],
                     'D': [0.25, 0.09, 0.65, 0.02], 'tau': [0.25, 0.13, 0.65, 0.02],
                     'q': [0.25, 0.17, 0.65, 0.02]}
slider_labels = {
    'b': 'Birth Rate (b)',
    'd': 'Death Rate of Uninfected (d)',
    'D': 'Additional Death Rate due to Infection (D)',
    'tau': 'Vertical Transmission Rate (tau)',
    'q': 'Probability of Cytoplasmic Incompatibility (q)'
}
for param, pos in slider_positions.items():
    ax_slider = plt.axes(pos, facecolor=axcolor)
    sliders[param] = Slider(ax_slider, slider_labels[param], 0.01, 1.0, valinit=eval(f'initial_{param}'))

# Update function for the sliders
def update(val):
    b = sliders['b'].val
    d = sliders['d'].val
    D = sliders['D'].val
    tau = sliders['tau'].val
    q = sliders['q'].val
    sol = solve_ivp(wolbachia_enhanced_model, t_span, initial_conditions, args=(b, d, D, tau, q), t_eval=t_eval)
    lines[0][0].set_ydata(sol.y[0])
    lines[1][0].set_ydata(sol.y[1])
    ax.relim()  # Recalculate limits
    ax.autoscale_view()  # Rescale the view based on the new data
    fig.canvas.draw_idle()

for slider in sliders.values():
    slider.on_changed(update)

# Set plot properties and show
plt.xlabel('Time')
plt.ylabel('Population')
plt.title('Dynamics of Wolbachia Infection with Interactive Parameters')
plt.legend()
plt.show()

## as uninfected gets fixed the infected rises 
## Add initial infection fix 
# infected clsoe to 1 uninfected close 0
