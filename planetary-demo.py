import rebound
import numpy as np
import streamlit as st
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def planetaryorbit():
        # Initialize REBOUND simulation
        sim = rebound.Simulation()
        sim.units = ('AU', 'yr', 'Msun')  # Set units to astronomical units, years, and solar masses
        
        # Add Sun
        sim.add(m=1)  # Sun mass in solar masses
        
        # Predefined planets with default values (mass in solar masses, semi-major axis in AU, eccentricity)
        PREDEFINED_PLANETS = {
            "Earth": [3e-6, 1, 0.0167],
            "Jupiter": [9.54e-4, 5.20, 0.0489],
            "Mars": [3.21e-7, 1.52, 0.0934],
            "Venus": [2.45e-6, 0.72, 0.0068],
            "Saturn": [2.86e-4, 9.58, 0.0537],
            "Mercury": [1.65e-7, 0.39, 0.2056],
            "Uranus": [4.36e-5, 19.19, 0.0472],
            "Neptune": [5.15e-5, 30.07, 0.0086]
        }
        
        # Streamlit Interface
        st.title('Orbit of Planets - Gravity Simulation with REBOUND')
        
        # Allow user to select planets and define initial states
        selected_planets = st.multiselect('Select Planets:', options=list(PREDEFINED_PLANETS.keys()))
        
        for planet in selected_planets:
            st.subheader(f"Initial Conditions for {planet}")
            defaults = PREDEFINED_PLANETS[planet]
            mass = st.number_input(f'Mass of {planet} (solar masses)', value=defaults[0], key=f'mass_{planet}')
            semi_major_axis = st.number_input(f'Semi-Major Axis of {planet} (AU)', value=defaults[1], key=f'sma_{planet}')
            eccentricity = st.number_input(f'Eccentricity of {planet}', min_value=0.0, max_value=1.0, value=defaults[2], key=f'ecc_{planet}')
            sim.add(m=mass, a=semi_major_axis, e=eccentricity)
        
        # Input for simulation time
        simulation_time = st.number_input('Simulation Time (years):', min_value=0.01, value=1.0)
        
        # Inputs for controlling the view (zoom) of the plot
        x_range = st.slider('X Range (AU):', min_value=1, max_value=100, value=35)
        y_range = st.slider('Y Range (AU):', min_value=1, max_value=100, value=35)
        
        # Button to start the simulation
        if st.button('Run Simulation'):
            # Lists to store paths
            paths_x = [[] for _ in range(sim.N - 1)] # Excluding the Sun
            paths_y = [[] for _ in range(sim.N - 1)] # Excluding the Sun
        
            # Run the simulation, saving the paths
            for t in range(100):  # You can change this number to control the resolution of the path
                sim.integrate(t * simulation_time / 100.0)
                for i, particle in enumerate(sim.particles[1:]):  # Skip the Sun
                    paths_x[i].append(particle.x)
                    paths_y[i].append(particle.y)
        
            # Plot the paths
            fig, ax = plt.subplots(figsize=(10, 10))
            # Set plot with dark background
            # plt.style.use('dark_background')
            plt.style.use('seaborn-darkgrid')
            ax.axis('equal')
            ax.set_xlim(-x_range, x_range) # Controlled by the user
            ax.set_ylim(-y_range, y_range) # Controlled by the user
            colors = plt.cm.tab10.colors # Use a color map for the planets
        
            for i, planet in enumerate(selected_planets):
                ax.plot(paths_x[i], paths_y[i], label=planet, lw=2, color=colors[i % len(colors)]) # Line width and color
            ax.plot(0, 0, 'yo', markersize=20, label='Sun') # Sun's position
        
            ax.set_xlabel('X (AU)', fontsize=12)
            ax.set_ylabel('Y (AU)', fontsize=12)
            ax.set_title('Planetary Orbits in the Solar System', fontsize=16)
            ax.legend(fontsize=10, loc='upper right')
        
            st.pyplot(fig)

def gravitationalpotential():
        # Selection of celestial objects
        objects = {
            'Earth': 5.972e24,
            'Jupiter': 1.898e27,
            'Sun': 1.989e30,
            'Black Hole (10 solar masses)': 1.989e31
        }
        
        selected_object = st.selectbox('Select a celestial object:', list(objects.keys()))
        
        # Corresponding mass of selected object
        mass = objects[selected_object]
        
        # Gravitational constant
        G = 6.67430e-11
        
        # Speed of light
        c = 299792458
        
        # Sliders for x and y axis ranges
        x_range = st.slider('Select X-Axis Range:', -100, 100, (-10, 10))
        y_range = st.slider('Select Y-Axis Range:', -100, 100, (-10, 10))
        
        # Grid of x, y values based on selected ranges
        x = np.linspace(x_range[0], x_range[1], 100)
        y = np.linspace(y_range[0], y_range[1], 100)
        x, y = np.meshgrid(x, y)
        
        # Calculate r values (distances from the mass)
        r = np.sqrt(x**2 + y**2)
        
        # Avoid division by zero
        r[r == 0] = 1e-9
        
        # Compute GR potential V (Schwarzschild approximation)
        V = -G * mass / r + G * mass**2 / (c**2 * r**2)
        
        # Compute gradients of V with respect to x and y
        V_x = np.gradient(V, axis=1)
        V_y = np.gradient(V, axis=0)
        
        # Create 3D plot
        fig = plt.figure(figsize=(10, 10))
        ax = fig.add_subplot(111, projection='3d')
        plt.style.use('seaborn')
        
        # Plot the surface
        surface = ax.plot_surface(x, y, V, cmap='viridis', linewidth=0, antialiased=True, alpha=0.5)
        
        # Add arrows to represent the gradients (subsample for visualization)
        step = 5
        for i in range(0, V_x.shape[0], step):
            for j in range(0, V_x.shape[1], step):
                ax.quiver(x[i, j], y[i, j], V[i, j], V_x[i, j], V_y[i, j], 0, color='r', length=1e9, arrow_length_ratio=0.1)
        
        # Add color bar
        fig.colorbar(surface, ax=ax)
        
        # Labels and title
        ax.set_xlabel('X (m)')
        ax.set_ylabel('Y (m)')
        ax.set_zlabel('Potential (Joules)')
        ax.set_title(f'Gravitational Potential of {selected_object} (GR approximation)')
        
        # Display plot in Streamlit
        st.pyplot(fig)

######################### Navigation
st.sidebar.title('Maths-Demo')
page = st.sidebar.radio("Go to", ['Planetary Orbit', 'Gravitational Potential'])

if page == 'Planetary Orbit':
    planetaryorbit()
elif page == 'Gravitational Potential':
    gravitationalpotential()
