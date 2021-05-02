import numpy as np
import copy
import matplotlib.pyplot as plt
import cantera as ct


# Im a Wrapper!
def moles(mixture_or_solution, species):
    """Returns moles of one species 'species' from a solution instance
       or Mixture (one phase only) instance"""
    if isinstance(mixture_or_solution, ct.Mixture):
        return mixture_or_solution.species_moles[mixture_or_solution.species_index(0, species)]
    elif isinstance(mixture_or_solution, ct.Solution):
        return mixture_or_solution.mole_fraction_dict().get(species)


# Im a Wrapper!
def mole_dict(mixture):
    """Returns moles of all species of a mixture as dict object"""
    species_keys = mixture.species_names
    species_values = mixture.species_moles
    current_flows = {}
    i = 0
    for species in species_keys:
        current_flows.update({species: species_values[i]})
        i = i + 1
    return current_flows


class Stream:
    """Tracks process changes between thermodynamic states

    Basis is per mole, units are SI.
    The Stream class is intended to work with solution objects
    but works with other Cantera.ThermoPhase objects such as pure fluids.
    Cannot use Mixtures, although the mixture solver is available.
    Use on your own risk

    Attributes
    ----------
        self._initial_moles : dict or float
        self._flow_rates : dict or float
        self._solution  : Cantera.Solution()
        self.temperatures : dict
        self.pressures : dict
        self.volumes : dict
        self.enthalpies : dict
        self.entropies : dict
        self.species_moles : dict
        self.total_moles : dict
        self.cp_mole : dict
        self._resolution: int

    Methods
    -------
        __init__()
        tracker()
        get_species_moles()
        get_species_mole_fraction()
        get_total_moles()
        get_temperature()
        get_pressure()
        get_volume()
        get_entropy()
        get_enthalpy()
        set_flowrates()
        isobaric_change_of_temperature()
        isothermal_change_of_pressure()
        isentropic_change_of_pressure
        isentropic_change_of_volume()
        ignore_species_except()
        ignore_species()
        chemical_reaction()
        separate_gases()
        separate_water()
        add_ideal_gas()
        set_reaction_mechanism()
        temperature_of_conversion_extent()

    """
    def __init__(self, thermo_stream: ct.Solution, initial_species_moles: dict):
        # Track initial_moles, current moles and solution
        self._initial_moles = initial_species_moles
        self._flow_rates = self._initial_moles
        self._solution = thermo_stream

        """
        # To DO modify type check to work with pure fluids
        if issubclass(Solution_object, ct._cantera.ThermoPhase):  
        else:
            raise TypeError("Stream(Solution_object,...) must be of type cantera.Solution")
        """

        # Create dictionaries with process data {'name of process': np.array[values during process]}
        # Set resolution of storage
        self.temperatures = {}
        self.pressures = {}
        self.volumes = {}  # per mole
        self.enthalpies = {}  # per mole
        self.entropies = {}  # per mole
        self.species_moles = {}
        self.total_moles = {}
        self.cp_mole = {}
        #  Only set before runtime.
        self._resolution = 200

    def tracker(self, process_name:str):
        """ Logs thermo properties between thermo states

        Method is called from thermo process methods and writes thermo properties
        into dictionaries with {'process name', list of thermodynamic values}

        :param process_name: str
        :return: 0
        """
        self.temperatures.setdefault(process_name, []).append(self._solution.T)
        self.pressures.setdefault(process_name, []).append(self._solution.P)
        # TO DO: Verify volume logging
        # self.volumes.setdefault(process_name, []).append(self._solution.volume_mole)
        if isinstance(self._flow_rates, dict):
            self.enthalpies.setdefault(process_name, []).append(
                self._solution.enthalpy_mole * sum(self._flow_rates.values()) / 1000)  # J
            self.entropies.setdefault(process_name, []).append(
                self._solution.entropy_mole * sum(self._flow_rates.values()) / 1000)  # J/K
            self.species_moles.setdefault(process_name, []).append(self._flow_rates)  # mol
            self.total_moles.setdefault(process_name, []).append(sum(self._flow_rates.values()))  # mol
            self.cp_mole.setdefault(process_name, []).append(self._solution.cp_mole /1000)  # J/ mol K
        elif isinstance(self._flow_rates, float):
            self.enthalpies.setdefault(process_name, []).append(
                self._solution.enthalpy_mole * self._flow_rates / 1000)  # J
            self.entropies.setdefault(process_name, []).append(
                self._solution.entropy_mole * self._flow_rates / 1000)  # J/K
            self.species_moles.setdefault(process_name, []).append(self._flow_rates)  # mol
            self.total_moles.setdefault(process_name, []).append(self._flow_rates)  # mol
            self.cp_mole.setdefault(process_name, []).append(self._solution.cp_mole / 1000)  # J/ mol K

    def get_species_moles(self, species: str) -> float:
        """Get moles of one calling species"""
        return self._flow_rates[species]

    def get_species_mole_fraction(self, species: str) -> float:
        """ Get species mole fraction of calling species"""
        return moles(self._solution, species)

    def get_total_moles(self) -> float:
        """Get total moles"""
        return sum(self._flow_rates.values())

    def get_temperature(self) -> float:
        """Get current temperature"""
        return self._solution.T

    def get_pressure(self) -> float:
        """Get current pressure"""
        return self._solution.P

    def get_volume(self) -> float:
        """Get current extensive volume"""
        # TO DO: Implement
        pass

    def get_entropy(self) -> float:
        """Get current entropy as extensive state"""
        if isinstance(self._flow_rates, float):
            return self._solution.entropy_mole * self._flow_rates / 1000  # J/K
        else:
            return self._solution.entropy_mole * sum(self._flow_rates.values()) / 1000  # J/K

    def get_enthalpy(self) -> float:
        """ Get current enthalpy as extensive state"""
        if isinstance(self._flow_rates, float):
            return self._solution.enthalpy_mole * self._flow_rates / 1000  # J
        else:
            return self._solution.enthalpy_mole * sum(self._flow_rates.values()) / 1000  # J

    def set_flowrates(self, flowrates: dict):
        """Set new flow rates"""
        self._flow_rates = flowrates

    def isobaric_change_of_temperature(self, T1: float, process_name: str):
        """Isobaric process that steps to a different temperature. No chemical reactions"""
        for temperature in np.linspace(self._solution.T, T1, self._resolution):
            self._solution.TP = temperature, self._solution.P
            self.tracker(process_name)

    def isothermal_change_of_pressure(self, p1: float, process_name: str):
        """Isothermal process that steps to a different pressure. No chemical reactions"""
        for pressure in np.linspace(self._solution.P, p1, self._resolution):
            self._solution.TP = self._solution.T, pressure
            self.tracker(process_name)

    def isentropic_change_of_pressure(self, P1: float, process_name: str):
        """Isentropic process that steps to a different pressure. No chemical reactions"""
        for pressure in np.linspace(self._solution.P, P1, self._resolution):
            self._solution.SP = self._solution.s, pressure
            self.tracker(process_name)

    def isentropic_change_of_volume(self, compression_ratio: float, process_name: str):
        """Isentropic process that steps to a different volume. No chemical reactions"""
        # TO DO: Implement
        pass

    def ignore_species_except(self, species: dict):
        # TO DO: Implement
        pass

    def ignore_species(self, species: dict):
        # TO DO: Implement
        pass

    def chemical_reaction(self, process_name):
        """Equilibrates solution and mixture to a new state by minimizing delta Gibbs energy
           Only for isobaric and isothermal chemical reactions"""

        # Set first state of of new process to old process for other calculations that depend on it
        self.tracker(process_name)

        # Chemical reaction until equilibrium
        mixture_chemical = ct.Mixture([(self._solution, sum(self._flow_rates.values()))])
        mixture_chemical.equilibrate('TP')
        self._flow_rates = mole_dict(mixture_chemical)
        self._solution.TPX = self._solution.T, self._solution.P, self._flow_rates

        # Call tracker
        for i in np.linspace(0, self._resolution, self._resolution - 1):
            self.tracker(process_name)

    def separate_gases(self, species_string, process_name, separated_stream_modeling_file='gri30.yaml'):
        # TO DO: Implement
        pass

    def separate_water(self, process_name: str):
        """Separate water from calling stream and return a water stream object

        :param process_name: str
        :return: seperated water: Stream
        """

        separated_water_flow = copy.copy(self._flow_rates['H2O'])

        # Create separated water return stream
        separated_water = ct.Water()
        separated_water.TP = self._solution.T, self._solution.P
        separated_water_stream = Stream(separated_water, 0)  # Initialized with 0 flow

        for i in np.linspace(1, 0, self._resolution):
            # Track changes with decreasing mole fraction of water in gas stream
            self._flow_rates['H2O'] = separated_water_flow * i
            self._solution.TPX = self._solution.T, self._solution.P, self._flow_rates
            self.tracker(process_name)
            # Track changes with increasing mole fraction of water in water stream
            separated_water_stream.set_flowrates(separated_water_flow * (1 - i))
            separated_water_stream.tracker(process_name)

        return separated_water_stream

    def add_ideal_gas(self, added_species_flows: dict):
        """Adds species with flow rates to the stream and models them as ideal gas"""
        for species in added_species_flows:
            # Add new species flow rates to old flow rates
            if species in self._flow_rates:
                self._flow_rates[species] = (self._flow_rates[species] + added_species_flows[species])
            else:
                self._flow_rates[species] = added_species_flows[species]
        self._solution.TPX = self._solution.T, self._solution.P, self._flow_rates

    def set_reaction_mechanism(self, reaction_mechanism_file, initial_species_moles):
        """Change modeling file of Solution"""
        # Initialize new Solution with new reaction mechanism
        reaction_solution = ct.Solution(reaction_mechanism_file)
        reaction_solution.TPX = self._solution.T, self._solution.P, initial_species_moles
        # Replace own solution and mixture with new reaction mechanism
        self._solution = reaction_solution

    def temperature_of_conversion_coefficient(self, conversion_extent: dict) -> float:
        """Calculates temperature of a mixture for a given conversion extent

        Method doesnt set any members of class,
        but creates a copy that is only used internally instead.
        :param conversion_extent: {'species', conversion_extent}
        :return: temperature
        """

        # Unpack argument to more readable variables
        species = [*conversion_extent][0]  # Species name as string
        conversion_extent_final = conversion_extent[species]  # Conversion extent as float

        # Create copies of solution that only gets used in this method
        # This is useful, because the property self._solution should not be changed by the mixture.
        copy_solution = ct.Solution(str(
            self._solution.ID) + '.yaml')  # Only works with yaml files. Legacy files .ct and .XMl not supported. Convert them to yaml
        copy_solution.TPX = self._solution.T, self._solution.P, self._solution.X
        reaction_mixture = ct.Mixture([(copy_solution, sum(self._initial_moles.values()))])  # Keep it private!

        # Initial conditions for while loop
        conversion_extent_current = 0
        T_current = copy_solution.T
        initial_species_moles = moles(reaction_mixture, species)

        # Create arrays for plotting
        temperatures = [T_current]
        conversion_extents = [conversion_extent_current]

        # Check for matching conversion extent with increasing temperature
        while conversion_extent_current <= conversion_extent_final:
            reaction_mixture.T = T_current
            reaction_mixture.equilibrate('TP')  # Find composition of the mixture @equilibrium in steady state
            conversion_extent_current = 1 - (moles(reaction_mixture, species) / initial_species_moles)

            temperatures.append(T_current)
            conversion_extents.append(conversion_extent_current)
            T_current = T_current + 2.5  # Current reforming temperature
        
        fig, ax = plt.subplots()
        ax.plot(np.array(temperatures), np.array(conversion_extents), color='#12719e')
        ax.set_facecolor('#f5f5f5')
        ax.set_ylabel(f'$Conversion$ $Extent$, {species}')
        ax.set_xlabel('$Process$ $Temperature$ $[T]$')
        ax.set_title(f'{species} $Reforming$ $at$ {round(copy_solution.P / 100000, 4)} $[Bar]$, $steady state$')
        ax.legend(['$F_{CH_4} = 1$    $mol$'  '\n' + '$F_{CO_2} = 0.3$ $mol$' + '\n' + '$F_{H_{2}} $  $= 2.5$ $mol$'])
        plt.savefig("Conversion_extent.png", dpi=1200)
        plt.show()

        return T_current  # K
