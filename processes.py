import numpy as np
import scipy as sc
import matplotlib.pyplot as plt
import cantera as ct
import ideal_thermodynamics as idt


# Returns moles of a species 'species' from a solution Instance or Mixture instance
def moles(mixture_or_solution, species):
    if isinstance(mixture_or_solution, ct.Mixture):
        return mixture_or_solution.species_moles[mixture_or_solution.species_index(0, species)]
    elif isinstance(mixture_or_solution, ct.Solution):
        return mixture_or_solution.mole_fraction_dict().get(species)


# Calculates change in enthalpy per mole for water, specified with temperature and pressure
def delta_h_water(T0, p0, T1, p1):
    water = ct.Water()
    water.TP = T0, p0
    h_0 = water.enthalpy_mole
    water.TP = T1, p1
    h_1 = water.enthalpy_mole
    return (h_1 - h_0) / 1000  # Enthalpy difference in J/mol


# Calculates change in enthalpy per mole for a Cantera Solution object,
# specified with temperature, pressure and species moles
def delta_h_solution(T0, p0, T1, p1, species_moles):
    solution = ct.Solution('gri30.yaml')
    solution.TPX = T0, p0, species_moles
    h0 = solution.enthalpy_mole
    solution.TP = T1, p1
    h1 = solution.enthalpy_mole
    return (h1 - h0) / 1000  # Enthalpy difference in J/mol


# Calculate reforming temperature for a given conversion extent of methane using the Cantera Mixture.equilibrate.
# A mixture object is used here, because the solution.equilibrate solver returns wrong values
def methane_reforming_temperature(T0, p0, conversion_extent_final, initial_species_moles):
    T_current = T0
    conversion_extent_current = 0

    temperatures = [T0]
    conversion_extents = [0]

    # Create an extensive state Mixture object, with a solution object and the total number of moles of all inputs
    solution = ct.Solution('gri30.yaml')
    solution.TPX = T0, p0, initial_species_moles
    mixture = ct.Mixture([(solution, sum(initial_species_moles.values()))])

    # Check for matching conversion extent with increasing temperature
    while conversion_extent_current <= conversion_extent_final:
        mixture.T = T_current
        mixture.equilibrate('TP')  # Find composition of the mixture @equilibrium in steady state
        conversion_extent_current = 1 - (moles(mixture, 'CH4') / initial_species_moles.get('CH4'))

        temperatures.append(T_current)
        conversion_extents.append(conversion_extent_current)

        T_current = T_current + 5  # Current reforming temperature

    # Format to numpy arrays and plot conversion extent over temperature
    npx = np.array(temperatures)
    npy = np.array(conversion_extents)

    fig, ax = plt.subplots()
    ax.plot(npx, npy)
    ax.set_ylabel('Conversion Extent, CH4')
    ax.set_xlabel('Process Temperature [$T$]')
    ax.set_title(f"Methane Reforming at {p0 / 100000} $Bar$, steady state")
    plt.show()
    return T_current  # K


# Calculates flow rates at the methane reactor outlet and returns them as dict object
def methane_reforming_flowrates(T0, p0, initial_species_moles):
    # Create an extensive state Mixture object, with a solution object and the total number of moles of all inputs
    solution = ct.Solution('gri30.yaml')
    solution.TPX = T0, p0, initial_species_moles
    mixture = ct.Mixture([(solution, sum(initial_species_moles.values()))])


def delta_h_mixture(T0, p0, initial_species_moles):
    # Create an extensive state Mixture object, with a solution object and the total number of moles of all inputs
    solution = ct.Solution('gri30.yaml')
    solution.TPX = T0, p0, initial_species_moles
    mixture = ct.Mixture([(solution, sum(initial_species_moles.values()))])

    h0 = solution.enthalpy_mole
    mixture.equilibrate('TP')
    h1 = solution.enthalpy_mole
    return (h1 - h0) / 1000


class Stream:
    # Create dictionaries with process data {'name of process': np.array[values during process]}
    # Convenient for plotting later
    temperatures = {}
    pressures = {}
    volumes = {}
    enthalpies = {}
    entropies = {}
