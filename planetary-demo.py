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
        AU_TO_M = 149597870700
        # Selection of celestial objects with corresponding mass and radius in AU
        objects = {
            'Earth': {'mass': 5.972e24, 'radius': 6371000 / AU_TO_M},
            'Jupiter': {'mass': 1.898e27, 'radius': 69911000 / AU_TO_M},
            'Sun': {'mass': 1.989e30, 'radius': 696340000 / AU_TO_M},
            'VY Canis Majoris': {'mass': 3.385e31, 'radius': 1.420e9 / AU_TO_M}
        }
        
        selected_object = st.selectbox('Select a celestial object:', list(objects.keys()))
        
        # Corresponding mass and radius of selected object
        mass = objects[selected_object]['mass']
        radius = objects[selected_object]['radius']
        
        # Gravitational constant
        G = 6.67430e-11
        
        # Slider for x and y axis ranges in AU
        x_range = st.slider('Select X-Axis Range (AU):', -1, 1, (-0.5, 0.5))
        y_range = st.slider('Select Y-Axis Range (AU):', -1, 1, (-0.5, 0.5))
        
        # Input for how many multiples of the radius to stop the plot
        radius_multiplier = st.number_input('Enter how many multiples of the radius to stop the plot:', min_value=1, value=2)
        
        # Convert selected ranges from AU to meters
        x_range_m = [r * AU_TO_M for r in x_range]
        y_range_m = [r * AU_TO_M for r in y_range]
        
        # Grid of x, y values based on selected ranges in meters
        x_m, y_m = np.meshgrid(np.linspace(x_range_m[0], x_range_m[1], 100),
                               np.linspace(y_range_m[0], y_range_m[1], 100))
        
        # Calculate r values (distances from the mass) in meters
        r_m = np.sqrt(x_m**2 + y_m**2)
        
        # Only consider distances greater than the specified multiples of the radius in meters
        mask = r_m < radius * AU_TO_M * radius_multiplier
        r_m[mask] = np.nan
        
        # Compute GR potential V (Schwarzschild approximation) in meters
        V_m = -G * mass / r_m + G * mass**2 / (c**2 * r_m**2)
        V_m[mask] = np.nan
        
        # Create 3D plot
        fig = plt.figure(figsize=(10, 10))
        ax = fig.add_subplot(111, projection='3d')
        plt.style.use('dark_background')
        
        # Plot the surface in AU
        surface = ax.plot_surface(x_m / AU_TO_M, y_m / AU_TO_M, V_m, cmap='viridis', linewidth=0, antialiased=True, alpha=0.5)
        
        # Add contours to the plot in AU
        contours = ax.contour(x_m / AU_TO_M, y_m / AU_TO_M, V_m, 10, colors='white', linestyles='solid', offset=np.nanmin(V_m))
        
        # Flip the z-axis
        ax.set_zlim(ax.get_zlim()[::-1])
        
        # Add color bar
        fig.colorbar(surface, ax=ax)
        
        # Labels and title in AU
        ax.set_xlabel('X (AU)')
        ax.set_ylabel('Y (AU)')
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
