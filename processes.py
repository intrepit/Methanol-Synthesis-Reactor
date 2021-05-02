import numpy as np
import scipy as sc
import copy
import matplotlib.pyplot as plt
import cantera as ct
import ideal_thermodynamics as idt


# Returns moles of a species 'species' from a solution Instance or Mixture(one phase only) instance
def moles(mixture_or_solution, species):
    if isinstance(mixture_or_solution, ct.Mixture):
        return mixture_or_solution.species_moles[mixture_or_solution.species_index(0, species)]
    elif isinstance(mixture_or_solution, ct.Solution):
        return mixture_or_solution.mole_fraction_dict().get(species)


# Returns the moles of species of a mixture as dict object
def mole_dict(mixture):
    species_keys = mixture.species_names
    species_values = mixture.species_moles
    current_flows = {}
    i = 0
    for species in species_keys:
        current_flows.update({species: species_values[i]})
        i = i + 1
    return current_flows


# The Stream class is used to track process changes
# Inputs are a Solution object, (extensive) moles either as dict{'species_i': moles_i, 'species_k': moles_j} or as float
# The Stream class cannot use Mixture objects, because Mixtures cant modify its underlying thermophase object
# Some methods of Mixture would modify the calling solution object,
# which makes it dangerous to (re)use the solution object in main
class Stream:
    def __init__(self, Solution_object, initial_species_moles):
        # Track initial_moles and moles
        self._initial_moles = initial_species_moles
        self._flow_rates = self._initial_moles
        # Initialize extensive Mixture Object.
        # This will be used to keep track of changes for a stream over all its processes.
        # Must be a Solution object, otherwise throw exception
        if isinstance(Solution_object, ct.Solution):
            self._solution = Solution_object  # Keep it private!
            #self._solution.basis
            self._mixture = ct.Mixture([(self._solution, sum(self._initial_moles.values()))])  # Keep it private!
        else:
            raise TypeError("Stream(Solution_object,...) must be of type cantera.Solution")
        # Create dictionaries with process data {'name of process': np.array[values during process]}
        # Set resolution of storage
        # Convenient for plotting later
        self.temperatures = {}
        self.pressures = {}
        # self.volumes = {}  # per mole
        self.enthalpies = {}  # per mole
        self.entropies = {}  # per mole
        self.species_moles = {}
        self.total_moles = {}
        self.resolution = 20

    # Tracks changes in thermodynamic values between states (a process)
    # Writes changes in dictionaries with {'process name', array of thermodynamic values (len = resolution)},
    # when called from a process function
    def tracker(self, process_name):
        self.temperatures.setdefault(process_name, []).append(self._solution.T)
        self.pressures.setdefault(process_name, []).append(self._solution.P)
        # self.volumes.setdefault(process_name, []).append(self._solution.volume_mole) #macht das sinn?
        self.enthalpies.setdefault(process_name, []).append(self._solution.enthalpy_mole * self._mixture.phase_moles(0) / 1000)  # J/mol
        self.entropies.setdefault(process_name, []).append(self._solution.entropy_mole * self._mixture.phase_moles(0) / 1000)  # J/K
        self.species_moles.setdefault(process_name, []).append(mole_dict(self._mixture))
        self.total_moles.setdefault(process_name, []).append(self._mixture.phase_moles())

    # Returns moles of one calling species
    def get_species_moles(self, species):
        return moles(self._mixture, species)

    # Return species mole fraction of calling species
    def get_species_mole_fraction(self, species):
        return moles(self._solution, species)

    # Implementation of an isobaric process that steps to a different temperature. No chemical reactions
    def isobaric_change_of_temperature(self, T1, process_name):
        for temperature in np.linspace(self._solution.T, T1, self.resolution):
            self._solution.TP = temperature, self._solution.P
            self._mixture.T = temperature
            self.tracker(process_name)

    # Implementation of an isothermal process that steps to a different pressure. No chemical reactions
    def isothermal_change_of_pressure(self, p1, process_name):
        for pressure in np.linspace(self._solution.P, p1, self.resolution):
            self._solution.TP = self._solution.T, pressure
            self._mixture.P = pressure
            self.tracker(process_name)

    def isentropic_change_of_pressure(self, P1, process_name):
        for pressure in np.linspace(self._solution.P, P1, self.resolution):
            self._solution.SP = self._solution.s, pressure
            self._mixture = ct.Mixture([(self._solution, self._mixture.phase_moles(0))])
            self.tracker(process_name)

    def isentropic_change_of_volume(self, compression_ratio, process_name):
        pass

    def ignore_species_except(self):
        pass

    def ignore_species(self):
        pass

    # Equilibrates solution and mixture to new state by minimizing deltaGibbs
    # @isobaric and isothermal only
    def chemical_reaction(self, reaction_mechanism, process_name):
        # Chemical reaction until equilibrium
        self._mixture.equilibrate('TP')
        for i in np.linspace(0, self.resolution, self.resolution):
            self.tracker(process_name)

    def separate_gas(self, species, process_name, separated_stream_modeling_file='gri30.yaml'):
        current_flows = mole_dict(self._mixture)
        separated_flow = current_flows.pop(species)  # separated species moles
        remaining_flow = current_flows  # Flow was separated from current flow. Copy to new name for convenience

        self._solution.TPX = self._solution.T, self._solution.P, remaining_flow
        self._mixture = ct.Mixture([(self._solution, sum(remaining_flow.values()))])

        separated_solution = ct.Solution(separated_stream_modeling_file)
        separated_solution.TPX = self._solution.T, self._solution.P, current_flows[species]

        for i in np.linspace(0, self.resolution, self.resolution):
            self.tracker(process_name)
        return separated_solution, separated_flow

    def separate_water(self, process_name):
        current_flows = mole_dict(self._mixture)
        separated_water_flow = current_flows.pop('H2O')  # separated species moles
        remaining_flow = current_flows  # Flow was separated from current flow. Copy to new name for convenience

        self._solution.TPX = self._solution.T, self._solution.P, remaining_flow
        self._mixture = ct.Mixture([(self._solution, sum(remaining_flow.values()))])

        separated_water = ct.Water()
        separated_water.TP = self._solution.T, self._solution.P

        for i in np.linspace(0, self.resolution, self.resolution):
            self.tracker(process_name)
        return separated_water, separated_water_flow

    # Adds species with flow rates to the stream. Input argument as dict
    def add_ideal_gas(self, added_species_flows):
        current_flows = mole_dict(self._mixture)
        # Add new species flow rates to old flow rates
        for species in current_flows:
            if species in added_species_flows:
                current_flows[species] = (current_flows[species] + added_species_flows[species])
        self._solution.TPX = self._solution.T, self._solution.P, current_flows
        self._mixture = ct.Mixture([(self._solution, sum(current_flows.values()))])

    def set_reaction_mechanism(self, reaction_mechanism_file, initial_species_moles):
        # Initialize new Solution with new reaction mechanism
        reaction_solution = ct.Solution(reaction_mechanism_file)
        reaction_solution.TPX = self._solution.T, self._solution.P, initial_species_moles
        # Replace own solution and mixture with new reaction mechanism
        self._solution = reaction_solution
        self._mixture = ct.Mixture([(self._solution, sum(initial_species_moles.values()))])

    # Calculates temperature of a mixture for a given conversion extent of a species,
    # given as dict{'species': conversion_extent_final}
    def temperature_of_conversion_coefficient(self, conversion_extent):
        # Unpack argument to more readable variables
        species = [*conversion_extent][0]  # Species name as string
        conversion_extent_final = conversion_extent[species]  # Conversion extent as float

        # Create copies of solution and mixture object that only get used in this method
        # This is useful, because the properties self._solution and self._mixture should not be changed
        # Method returns temperature but doesnt set the class instance!
        copy_solution = ct.Solution(str(self._solution.ID) + '.yaml')  # Only works with yaml files. Legacy files .ct and .XMl not supported. Convert them to yaml
        copy_solution.TPX = self._solution.T, self._solution.P, self._solution.X

        # Mixture object must be created from copy_solution,
        # otherwise the _solution object will be set to new state!

        copy_mixture = ct.Mixture([(copy_solution, sum(self._initial_moles.values()))])  # Keep it private!

        # Initial conditions for while loop
        conversion_extent_current = 0
        T_current = copy_solution.T
        initial_mixture_moles = moles(copy_mixture, species)

        # Create arrays for plotting
        temperatures = [T_current]
        conversion_extents = [conversion_extent_current]

        # Check for matching conversion extent with increasing temperature
        while conversion_extent_current <= conversion_extent_final:
            copy_mixture.T = T_current
            copy_mixture.equilibrate('TP')  # Find composition of the mixture @equilibrium in steady state
            conversion_extent_current = 1 - (moles(copy_mixture, species) / initial_mixture_moles)

            temperatures.append(T_current)
            conversion_extents.append(conversion_extent_current)

            T_current = T_current + 2.5  # Current reforming temperature

        # Format to numpy arrays and plot conversion extent over temperature
        npx = np.array(temperatures)
        npy = np.array(conversion_extents)

        fig, ax = plt.subplots()
        ax.plot(npx, npy)
        ax.set_ylabel(f'Conversion Extent, {species}')
        ax.set_xlabel('Process Temperature [$T$]')
        ax.set_title(f'{species} Reforming at {round(copy_solution.P / 100000, 4)} $Bar$, steady state')
        plt.show()
        return T_current  # K


class WaterStream(Stream):
    pass
