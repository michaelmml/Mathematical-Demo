import rebound
import numpy as np
import streamlit as st
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.integrate import quad
from scipy.special import genlaguerre, factorial, sph_harm
from scipy.constants import epsilon_0, hbar, m_e, e

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
            'Star': {'mass': 1.989e30, 'radius': 19634000000 / AU_TO_M},
            'Object1': {'mass': 1.989e29, 'radius': 19634000000 / AU_TO_M},
            'Object2': {'mass': 1.989e29, 'radius': 19634000000 / AU_TO_M},
        }
        
        selected_object = st.selectbox('Select a celestial object:', list(objects.keys()))
        
        selected_object1 = st.selectbox('Select a second celestial object:', list(objects.keys()))
        mass_scaler1 = st.slider('Object 1 Mass Adjustment:', min_value=0.1, max_value=10.0, value=1.0)
        radius_scaler1 = st.slider('Object 1 Radius Adjustment:', min_value=0.1, max_value=10.0, value=1.0)
        dist_scaler1 = st.slider('Object 1 Distance Adjustment:', min_value=0.1, max_value=10.0, value=1.0)       
        
        selected_object2 = st.selectbox('Select a third celestial object:', list(objects.keys()))
        mass_scaler2 = st.slider('Object 2 Mass Adjustment:', min_value=0.1, max_value=10.0, value=1.0)
        radius_scaler2 = st.slider('Object 2 Radius Adjustment:', min_value=0.1, max_value=10.0, value=1.0)
        dist_scaler2 = st.slider('Object 2 Distance Adjustment:', min_value=0.1, max_value=10.0, value=1.0)        
        
        # Corresponding mass and radius of selected object
        mass = objects[selected_object]['mass']
        radius = objects[selected_object]['radius']
        mass1 = objects[selected_object1]['mass']*mass_scaler1
        radius1 = objects[selected_object1]['radius']*radius_scaler1
        mass2 = objects[selected_object2]['mass']*mass_scaler2
        radius2 = objects[selected_object2]['radius']*radius_scaler2       
        
        # Gravitational constant
        G = 6.67430e-11

        # Speed of light
        c = 299792458
        
        # Slider for x and y axis ranges in AU
        x_range = st.slider('Select X-Axis Range (AU):', -2.0, 2.0, (-0.5, 0.5), step=0.01)
        y_range = st.slider('Select Y-Axis Range (AU):', -2.0, 2.0, (-0.5, 0.5), step=0.01)
        
        # Input for how many multiples of the radius to stop the plot
        radius_multiplier = st.number_input('Enter how many multiples of the radius to stop the plot:', min_value=0.0, value=0.0, step=0.1)
        
        # Convert selected ranges from AU to meters
        x_range_m = [r * AU_TO_M for r in x_range]
        y_range_m = [r * AU_TO_M for r in y_range]
        
        # Grid of x, y values based on selected ranges in meters
        x_m, y_m = np.meshgrid(np.linspace(x_range_m[0], x_range_m[1], 200),
                               np.linspace(y_range_m[0], y_range_m[1], 200))
        
        # Calculate r values (distances from the mass) in meters
        r_m = np.sqrt(x_m**2 + y_m**2)
        
        # Only consider distances greater than the specified multiples of the radius in meters
        mask = r_m < radius * AU_TO_M * radius_multiplier
        r_m[mask] = np.nan
        
        # Compute GR potential V (Schwarzschild approximation) in meters
        # V_m = -G * mass / r_m + G * mass**2 / (c**2 * r_m**2)
        # V_m[mask] = np.nan

        ####### Compute GR potential V (Schwarzschild approximation) in meters and take into account of volume
        V_m0 = potential(x_m, y_m, mass, radius * AU_TO_M)
        V_m1 = potential(x_m - (AU_TO_M*dist_scaler1), y_m, mass1, radius1 * AU_TO_M) 
        V_m2 = potential(x_m, y_m - (AU_TO_M*dist_scaler2), mass2, radius2 * AU_TO_M)  
        V_m = np.log(np.abs(V_m0 + V_m1 + V_m2))
        V_m[mask] = np.nan
        
        # Create 3D plot
        # fig = plt.figure(figsize=(10, 10))
        # ax = fig.add_subplot(111, projection='3d')
        # plt.style.use('dark_background')

        # Create 3D plot
        fig = plt.figure(figsize=(10, 10))
        ax = fig.add_subplot(111, projection='3d')
        
        # Set black background for the plot
        ax.set_facecolor('black')
        fig.patch.set_facecolor('black')
        
        # Plot the surface in AU
        surface = ax.plot_surface(x_m / AU_TO_M, y_m / AU_TO_M, V_m, cmap='viridis', linewidth=0, antialiased=True, alpha=0.5)
        
        # Add contours to the plot in AU
        contours = ax.contour(x_m / AU_TO_M, y_m / AU_TO_M, V_m, 10, colors='white', linestyles='solid', offset=np.nanmin(V_m))
        
        # Flip the z-axis
        ax.set_zlim(ax.get_zlim()[::-1])
        
        # Remove axis labels
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.set_zticklabels([])
        
        # Add color bar
        fig.colorbar(surface, ax=ax)
        
        # Display plot in Streamlit
        st.pyplot(fig)


# Function to compute m(r)
def m(r, density):
        return (4/3) * np.pi * density * r**3

# Function to integrate
def integrand(r, density):
        G = 6.67430e-11  # m^3 kg^-1 s^-2, gravitational constant
        return G * m(r, density) / r

# Compute gravitational potential
def potential(x, y, mass, radius):
        # Constants
        G = 6.67430e-11  # m^3 kg^-1 s^-2, gravitational constant
        c = 2.998e8  # m/s, speed of light
        AU_TO_M = 1.496e11 # 1 Astronomical Unit in meters
        r = np.sqrt(x**2 + y**2)
        density = mass / ((4/3) * np.pi * radius**3)

        V = np.zeros_like(r)

        # Outside the sphere (Schwarzschild solution)
        r_s = 2 * G * mass / c**2  # Schwarzschild radius
        outside_mask = r >= radius
        V[outside_mask] = -G * mass / r[outside_mask] + (G * mass * r_s) / (2 * r[outside_mask]**2)

        # Inside the sphere (uniform density)
        inside_mask = r < radius
        V[inside_mask] = -(G * mass / (2 * radius**3)) * (3 * radius**2 - r[inside_mask]**2)

        # Inside the sphere (numerical integration)
        # for i in np.ndindex(r.shape):
        #         if r[i] < radius:
        #                 V[i], _ = quad(integrand, 0, r[i], args=(density,))

        return V

######################### Quantum Mechanics
# Define wavefunctions for different energy states of hydrogen atom

def wavefunction(n, l, m, r, theta, phi):
        a0 = 4 * np.pi * epsilon_0 * hbar**2 / (m_e * e**2)  # Bohr radius
        rho = 2 * r / (n * a0)
        L = genlaguerre(n-l-1, 2*l+1)
        Y = sph_harm(m, l, phi, theta)  # Note: sph_harm's arguments are m, l, phi, theta
        R = np.sqrt((2/n/a0)**3 * factorial(n-l-1)/(2*n*factorial(n+l))) * np.exp(-rho/2) * rho**l * L(rho)
        return R * Y


# Probability density
def probability_density(n, l, m, r, theta, phi):
        psi = wavefunction(n, l, m, r, theta, phi)
        return np.abs(psi)**2

def schrodinger():
        a0 = 4 * np.pi * epsilon_0 * hbar**2 / (m_e * e**2)  # Bohr radius
        
        # Streamlit Inputs for the first 3 states
        state_defaults = [(4, 0, 0), (4, 2, 0), (4, 2, 1)]
        states_inputs = []

        for i, default in enumerate(state_defaults, 1):
            n = st.selectbox(f'Select n for plot {i}', list(range(1, 5)), index=default[0])
            l = st.selectbox(f'Select l for plot {i}', list(range(0, n)), index=default[1])
            m = st.selectbox(f'Select m for plot {i}', list(range(0, l+1)), index=default[2])
            states_inputs.append((n, l, m))
        
        # Radial extent input
        radial_extent = st.slider('Select Radial Extent (Bohr Radii):', 1, 50, 20)
        
        # Create a grid of points in polar coordinates        
        r = np.linspace(0, radial_extent * a0, radial_extent*20)  # Increased extent and resolution
        theta = np.linspace(0, np.pi, radial_extent*20)  # Increased resolution for smooth reflection
        R, Theta = np.meshgrid(r, theta)
        
        # States to consider in the 3x3 grid
        states = states_inputs + [(2, 1, 1), (3, 0, 0), (3, 1, 0), (3, 1, 1), (3, 2, 0), (3, 2, 1)]
        
        # Create a 3x3 grid of contour plots
        fig, axs = plt.subplots(3, 3, figsize=(15, 15))
        
        for ax, (n, l, m) in zip(axs.flat, states):
            rho = probability_density(n, l, m, R, Theta, 0)
            X = R * np.sin(Theta)
            Y = R * np.cos(Theta)
            
            # Reflect the data
            ax.contourf(X, Y, rho, 100, cmap='viridis')
            ax.contourf(-X, Y, rho, 100, cmap='viridis')
            
            ax.set_title(f"n={n}, l={l}, m={m}")
            ax.axis('equal')
        
        plt.tight_layout()
        # Display in Streamlit
        st.pyplot(fig)

######################### Navigation
st.sidebar.title('Maths-Demo')
page = st.sidebar.radio("Go to", ['Planetary Orbit', 'Gravitational Potential', 'Quantum Probabilistic Space'])

if page == 'Planetary Orbit':
    planetaryorbit()
elif page == 'Gravitational Potential':
    gravitationalpotential()
elif page == 'Quantum Probabilistic Space':
    schrodinger()
