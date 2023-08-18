import rebound
import streamlit as st
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Initialize REBOUND simulation
sim = rebound.Simulation()
sim.units = ('AU', 'yr', 'Msun') # Set units to astronomical units, years, and solar masses

# Add Sun
sim.add(m=1) # Sun mass in solar masses

# Predefined planets (name, mass in solar masses, semi-major axis, eccentricity)
PREDEFINED_PLANETS = {
    "Earth": [3e-6, 1, 0.017]
}

# Function to add planets to simulation
def add_planet(mass, semi_major_axis, eccentricity):
    sim.add(m=mass, a=semi_major_axis, e=eccentricity)

# Add predefined planets
for planet in PREDEFINED_PLANETS.values():
    add_planet(*planet)

# Streamlit Interface
st.title('Orbit of Planets - Gravity Simulation with REBOUND')

if st.button('Add Planet'):
    mass = st.number_input('Mass of the Planet (solar masses)', min_value=0.0)
    semi_major_axis = st.number_input('Semi-Major Axis (AU)', min_value=0.0)
    eccentricity = st.number_input('Eccentricity', min_value=0.0, max_value=1.0)
    add_planet(mass, semi_major_axis, eccentricity)
    st.write('Planet Added!')

# Create Animated Plot Function
def animate_planets():
    fig, ax = plt.subplots(figsize=(8,8))
    ax.axis('equal')
    ax.set_xlim(-2, 2)
    ax.set_ylim(-2, 2)

    lines = [ax.plot([], [], 'o')[0] for _ in range(sim.N)]

    def init():
        for line in lines:
            line.set_data([], [])
        return lines

    def update(frame):
        sim.integrate(sim.t + 0.01) # integrate by a small time step
        for i, line in enumerate(lines):
            particle = sim.particles[i]
            line.set_data(particle.x, particle.y)
        return lines

    ani = FuncAnimation(fig, update, frames=1000, init_func=init, blit=True)
    plt.show()

if st.button('Start Simulation'):
    animate_planets()
