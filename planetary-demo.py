import rebound
import streamlit as st
import matplotlib.pyplot as plt

# Initialize REBOUND simulation
sim = rebound.Simulation()
sim.units = ('AU', 'yr', 'Msun') # Set units to astronomical units, years, and solar masses

# Add Sun
sim.add(m=1) # Sun mass in solar masses

# Predefined planets
PREDEFINED_PLANETS = [
    "Earth", "Jupiter", "Mars", "Venus", "Saturn", "Mercury", "Uranus", "Neptune"
]

# Streamlit Interface
st.title('Orbit of Planets - Gravity Simulation with REBOUND')

# Allow user to select planets and define initial states
selected_planets = st.multiselect('Select Planets:', options=PREDEFINED_PLANETS)

for planet in selected_planets:
    st.subheader(f"Initial Conditions for {planet}")
    mass = st.number_input(f'Mass of {planet} (solar masses)', key=f'mass_{planet}')
    semi_major_axis = st.number_input(f'Semi-Major Axis of {planet} (AU)', key=f'sma_{planet}')
    eccentricity = st.number_input(f'Eccentricity of {planet}', min_value=0.0, max_value=1.0, key=f'ecc_{planet}')
    sim.add(m=mass, a=semi_major_axis, e=eccentricity)

# Input for simulation time
simulation_time = st.number_input('Simulation Time (years):', min_value=0.01)

# Button to start the simulation
if st.button('Run Simulation'):
    # Lists to store paths
    paths_x = [[] for _ in range(sim.N)]
    paths_y = [[] for _ in range(sim.N)]

    # Run the simulation, saving the paths
    for t in range(100): # You can change this number to control the resolution of the path
        sim.integrate(t * simulation_time / 100.0)
        for i, particle in enumerate(sim.particles[1:]): # Skip the Sun
            paths_x[i].append(particle.x)
            paths_y[i].append(particle.y)

    # Plot the paths
    fig, ax = plt.subplots(figsize=(8,8))
    ax.axis('equal')
    ax.set_xlim(-2, 2)
    ax.set_ylim(-2, 2)
    for i, planet in enumerate(selected_planets):
        ax.plot(paths_x[i], paths_y[i], label=planet)
    ax.legend()
    st.pyplot(fig)
