# Methanol Synthesis Reactor
import numpy as np
import matplotlib.pyplot as plt
import cantera as ct
from processes import Stream
import ideal_thermodynamics as idt
# Initial flow rates for all streams
F_CH4 = 1  # mol/s
F_CO2 = 0.3  # mol/s
F_H2O = 2.5  # mol/s

gas_flow_rates = {'CH4': F_CH4, 'CO2': F_CO2}
water_flow_rates = {'H2O': F_H2O}
flow_rates = {'CH4': F_CH4, 'CO2': F_CO2, 'H2O': F_H2O}  # Stored in dict and modified during the process chain. #stream1,2,3 flowrates

# Define state values @state 1
T1 = 300  # K
p1 = 15 * 100000  # Pa

# Initialize Cantera objects.
# One stream of water, water_in and another stream of gases CO2 + CH4 as gases_in
# Water
water_in = ct.Water()
water_in.TP = T1, p1
# Gases
gases_in = ct.Solution('gri30.yaml')
gases_in.TPX = T1, p1, gas_flow_rates
# Track gases in a process stream with the Stream class
reactor_stream = Stream(gases_in, gas_flow_rates)
resolution = 20
reactor_stream.resolution = resolution


# process 1-2, isobaric heating
# @state 2
T2 = 600  # K
p2 = p1  # Pa
# Heat gas stream
reactor_stream.isobaric_change_of_temperature(T2, '1-2')
# Boil water stream. I didnt write a class for this, thus ugly code
water_temperatures = {}
water_pressures = {}
water_entropies = {}
water_enthalpies = {}

for temperature in np.linspace(T1, T2, resolution):
    water_in.TP = temperature, p2
    water_enthalpies.setdefault('1-2', []).append(water_in.enthalpy_mole)
    water_entropies.setdefault('1-2', []).append(water_in.entropy_mole)
    water_temperatures.setdefault('1-2', []).append(water_in.T)
    water_pressures.setdefault('1-2', []).append(water_in.P)

print(water_temperatures)

# Mix streams @ T,P conditions of existing stream @ point 2
reactor_stream.add_ideal_gas(water_flow_rates)

# process 2-3, isobaric heating
# to specified conversion extent of methane
# @state 3
X_CH4 = 0.97  # Conversion extent methane
p3 = p2  # Pa
T3 = reactor_stream.temperature_of_conversion_coefficient({'CH4': X_CH4})  # K
reactor_stream.isobaric_change_of_temperature(T3, '2-3')

# process 3-4, isobaric isothermal catalytic methane reforming @equilibrium
# @state 4
T4 = T3
p4 = p3
# Uses the stream __init__(solution) mechanism as default, here 'gri30.yaml'
reactor_stream.chemical_reaction('gri30.yaml', '3-4')

# process 4-5, isobaric cooling
# @state 5
p5 = p4
T5 = idt.water_vapour_temperature(p4)  # solution.T_sat throws error. argument is partial pressure of water. Is this correct????
reactor_stream.isobaric_change_of_temperature(T5, '4-5')

# process 5-6, isobaric isothermal steam separation
# @state 6
T6 = T5
p6 = p5
# Remove H2O of stream instance, method returns solution object of separated stream, flowrate of separated stream
separated_water, separated_water_flow = reactor_stream.separate_water('5-6')  # Why do we track water enthalpy?????

# process 6-7, isobaric cooling
# @state 7
p7 = p6
T7 = 300  # K
reactor_stream.isobaric_change_of_temperature(T7, '6-7')

# process 7-8, isentropic compression
# @state 8
p8 = 60 * 100000  # Pa
reactor_stream.isentropic_change_of_pressure(p8, '7-8')

# process 8-9, isobaric change in heat
# @state 9
p9 = p8
T9 = 503  # K
reactor_stream.isobaric_change_of_temperature(T9, '8-9')

# process 9-10, isobaric isothermal catalytic methanol synthesis @equilibrium
p10 = p9
T10 = T9
# Create the catalytic reactor bed with the reaction mechanism 'methanol-synthesis.cti
# Ignore all species except H2, CO2, CO, thus setting new flow rates. This will cut on values like enthalpy etc (!)
# This is a wack hack.
# Better to define a custom .yaml, where species names match gri30.yaml but only catalytic reactions are allowed
reactor_stream.set_reaction_mechanism('methanol-synthesis.cti',
                                      {'H2': reactor_stream.get_species_moles('H2'),
                                       'CO': reactor_stream.get_species_moles('CO'),
                                       'CO2': reactor_stream.get_species_moles('CO2')})
reactor_stream.chemical_reaction('methanol-synthesis.cti', '9-10')


# Data Analysis
# Create numpy arrays for plotting
temperature_list = []
for key in reactor_stream.temperatures:
    temperature_list.extend(reactor_stream.temperatures[key])
temperatures = np.array(temperature_list)

pressure_list = []
for key in reactor_stream.pressures:
    pressure_list.extend(reactor_stream.pressures[key])
pressures = np.array(pressure_list)

entropy_list = []
add_water_entropy = False

for key in reactor_stream.entropies:  # Add water stream entropies starting from process 5-6 to the list
    if key == '5-6':
        add_water_entropy = True
    if not add_water_entropy:
        entropy_list.extend(reactor_stream.entropies[key])
    else:
        for entropies in reactor_stream.entropies[key]:
            entropy_list.append(entropies + (separated_water.entropy_mole * separated_water_flow / 1000))

entropies = np.array(entropy_list)

enthalpy_list = []
add_water_enthalpies = False

debug = []
debug_temperatures = []

for key in reactor_stream.enthalpies:
    if key == '5-6':
        add_water_enthalpies = True
    if not add_water_enthalpies:
        enthalpy_list.extend(reactor_stream.enthalpies[key])
    else:
        for enthalpies in reactor_stream.enthalpies[key]:
            enthalpy_list.append(enthalpies + (separated_water.enthalpy_mole * separated_water_flow / 1000))
            if key == '7-8':
                debug.append(enthalpies + (separated_water.enthalpy_mole * separated_water_flow / 1000))
                debug_temperatures = reactor_stream.temperatures['7-8']

enthalpies = np.array(enthalpy_list)

fig, ax = plt.subplots()
ax.plot(enthalpies, temperatures)
ax.plot(debug, debug_temperatures)
ax.set_ylabel('Temperature [$K$]')
ax.set_xlabel('Enthalpy [$J/(K)$]')
ax.set_title('T-H diagram')
plt.show()

fig, ax = plt.subplots()
ax.plot(entropies, temperatures)
ax.set_ylabel('Temperature [$K$]')
ax.set_xlabel('Entropy [$J/(K)$]')
ax.set_title('T-S diagram')
plt.show()



"""
# Calculate heat input per unit mole???
# Water passes through Boiler, CH4 and CO2 passes through heater, Both streams get mixed after @state 2
delta_h_water = process.delta_h_water(T1, p1, T2, p2)  # J/mol
delta_h_CH4_CO2 = process.delta_h_solution(T1, p1, T2, p2, {'CH4': F_CH4, 'CO2': F_CO2})  # J/mol
delta_h = (F_H2O * delta_h_water) + ((F_CH4 + F_CO2) * delta_h_CH4_CO2)  # J / s ?????????????
heat_1to2 = delta_h  # J / mol*s
"""


""""



# MethanolSynthesis from Methane, Water and CarbonDioxide
# Yeah!


# @1: Initiate Gas Solution Methane and CarbonDioxide
# CO2_CH4 = ct.Solution('gri30.yaml')
# CO2_CH4.X = {'CH4': FCH4, 'CO2': FCO2} #Important to set before .TP, otherwise pressure will change!
# CO2_CH4.TP = T1, p1

# Initiate Water as a Water Object
# water = ct.Water()
# water.TP = T1, p1

# create Values to keep track off


# @2: Thermodynamic state
T2 = 600  # K
p2 = p1

# CO2_CH4.TP = T2, p2 #1.1-2.1 Heat Gas Solution isobaric to 600K
# water.TP = T2, p2 #1.2-2.2 Heat Water isobaric to 600K

# New mixture to get extensive properties
temp = ct.Solution('gri30.yaml')
temp.X = {'CH4': F_CH4, 'CO2': F_CO2, 'H2O': F_H2O}  # Important to set before .TP, otherwise pressure will change!
temp.TP = T2, p2
CO2_CH4_H20 = ct.Mixture([(temp, (F_CH4 + F_CO2 + F_H2O))])

# @3: Thermodynamic state
T_reform = 600  # K, initial guess
p3 = p2
X_CH4 = 0.97  # =1 -(F_CH4@3/F_CH4@2) = 1-(n_CH4@3/n_CH4@2)

temperature_list = [T_reform]
conversionextent_list = [0]

i = 0  # start with 0 conversion of methane
while i <= X_CH4:
    CO2_CH4_H20.T = T_reform
    CO2_CH4_H20.equilibrate('TP')
    T_reform = T_reform + 5
    i = 1 - (CO2_CH4_H20.species_moles[13] / F_CH4)  # handwavy assumed to take F instead of n
    temperature_list.append(T_reform)
    conversionextent_list.append(i)

npx = np.array(temperature_list)
npy = np.array(conversionextent_list)

fig, ax = plt.subplots()
ax.plot(npx, npy)
ax.set_ylabel('Conversion Extent, CH4')
ax.set_xlabel('Process Temperature [$T$]')
ax.set_title(f"Methane Reforming at {p3 / 100000} $Bar$, steady state")
plt.show()

syngas_moles = CO2_CH4_H20.species_moles[0] + CO2_CH4_H20.species_moles[14] + CO2_CH4_H20.species_moles[
    15]  # Moles of H2, CO, CO2 at 4
x_H2 = CO2_CH4_H20.species_moles[0] / syngas_moles  # Molefraction H2 at 4
x_CO = CO2_CH4_H20.species_moles[14] / syngas_moles  # Molefraction CO at 4
x_CO2 = CO2_CH4_H20.species_moles[15] / syngas_moles  # Molefraction CO2 at 4


# @5 isobaric cooling to T = water gets liquid.

def clausius_clapeyron(p2):
    L = 40800  # J/mol
    R = 8.314462  # J/mol K
    T1 = 273.15 + 100  # K
    p1 = 101325  # Pa

    T2 = 1 / ((0 - np.log(p2 / p1) * (R / L)) + (1 / T1))
    return (T2)


T_H20 = clausius_clapeyron(p3)  # Water Vapor Temperature at p3
CO2_CH4_H20.T = T_H20

syngas_moles = CO2_CH4_H20.species_moles[0] + CO2_CH4_H20.species_moles[14] + CO2_CH4_H20.species_moles[
    15]  # Moles of H2, CO, CO2 at 4
print(syngas_moles)
x_H2 = CO2_CH4_H20.species_moles[0] / syngas_moles  # Molefraction H2 at 4
x_CO = CO2_CH4_H20.species_moles[14] / syngas_moles  # Molefraction CO at 4
x_CO2 = CO2_CH4_H20.species_moles[15] / syngas_moles  # Molefraction CO2 at 4
print(x_H2, x_CO, x_CO2)
S_module = (x_H2 - x_CO2) / (x_CO2 + x_CO)

T7 = 300  # K
p7 = 15 * 100000  # Pa
p8 = 60 * 100000  # Pa
methanolreactorphase = ct.Solution('methanol-synthesis.cti')
methanolreactorphase.X = {'H2': x_H2, 'CO': x_CO, 'CO2': x_CO2}
methanolreactorphase.TP = T7, p7


def adiabatic_compression(T0, p0, p1,
                          gamma):  # ideal isentropic compression. When is adiabaitc = reversible ok. Isnt isentropic = isotherm and adiabatic????
    T1 = T0 * (p0 / p1) ** ((1 - gamma) / gamma)
    return T1


def adiabatic_compression_real(p1, s0):
    while methanolreactorphase.s <= s0:
        methanolreactorphase.TP = (methanolreactorphase.T + 1), p1


T8 = adiabatic_compression(T7, p7, p8, methanolreactorphase.cp_mole / methanolreactorphase.cv_mole)

T8 = adiabatic_compression_real(p8, methanolreactorphase.s)
methanolreactorphase.TP = T8, p8

T9 = 503  # K
p9 = p8  # Pa
methanolreactorphase.TP = T9, p9
methanolreactorphase.equilibrate('TP')
a = 0

X_CO2_CO = methanolreactorphase.X[19] / (
            methanolreactorphase.X[13] + methanolreactorphase.X[14] + methanolreactorphase.X[19])


"""



