# Methanol Synthesis Reactor
import numpy as np
import matplotlib.pyplot as plt
import cantera as ct
import ideal_thermodynamics as idt
from operator import add
from Stream import Stream

# Modeling description

# Initial flow rates for the water and the gas stream
F_CH4 = 1  # mol/s
F_CO2 = 0.3  # mol/s
F_H2O = 2.5  # mol/s
gas_flows_in = {'CH4': F_CH4, 'CO2': F_CO2}  # mol/s
water_flows_in = {'H2O': F_H2O}  # mol/s

# Define state values @state 1
T1 = 300  # K
p1 = 15 * 100000  # Pa

# Initialize Cantera objects.
# Water to be modeled as pure fluid with a water EOS.
water_in = ct.Water()
water_in.TP = T1, p1

# Gases to be modeled as ideal gases with gri30
gases_in = ct.Solution('gri30.yaml')
gases_in.TPX = T1, p1, gas_flows_in

# Track both streams as a stream object.
water_in_stream = Stream(water_in, water_flows_in)
gas_stream = Stream(gases_in, gas_flows_in)

# process 1-2, isobaric heating
# @state 2
T2 = 600  # K
p2 = p1  # Pa
# Heat gas stream
gas_stream.isobaric_change_of_temperature(T2, '1-2')
# Boil water stream
water_in_stream.isobaric_change_of_temperature(T2, '1-2')
# Mix streams @ T,P conditions of existing ideal gas stream
gas_stream.add_ideal_gas(water_flows_in)

# process 2-3, isobaric heating
# to specified conversion extent of methane
# @state 3
X_CH4 = 0.97  # Conversion extent methane
p3 = p2  # Pa
T3 = gas_stream.temperature_of_conversion_coefficient({'CH4': X_CH4})  # K
gas_stream.isobaric_change_of_temperature(T3, '2-3')

# process 3-4, isobaric isothermal catalytic methane reforming @equilibrium
# @state 4
T4 = T3
p4 = p3
# Uses the stream __init__(solution) mechanism as default, here 'gri30.yaml'
gas_stream.chemical_reaction('3-4')

# process 4-5, isobaric cooling
# @state 5
p5 = p4
T5 = idt.water_vapour_temperature(p4)
gas_stream.isobaric_change_of_temperature(T5, '4-5')

# process 5-6, isobaric isothermal steam separation
# @state 6
T6 = T5
p6 = p5
# Remove H2O of stream instance, and track separated stream as a new Stream()
water_separated_stream = gas_stream.separate_water('5-6')

# process 6-7, isobaric cooling
# @state 7
p7 = p6
T7 = 300  # K
gas_stream.isobaric_change_of_temperature(T7, '6-7')

# process 7-8, isentropic compression
# @state 8
p8 = 60 * 100000  # Pa
gas_stream.isentropic_change_of_pressure(p8, '7-8')

# process 8-9, isobaric change in heat
# @state 9
p9 = p8
T9 = 503  # K
gas_stream.isobaric_change_of_temperature(T9, '8-9')

# process 9-10, isobaric isothermal catalytic methanol synthesis @equilibrium
p10 = p9
T10 = T9
# Create the catalytic reactor bed with the reaction mechanism 'methanol-synthesis.cti
# Ignore all species except H2, CO2, CO, thus setting new flow rates.
gas_stream.set_reaction_mechanism('methanol-synthesis.cti',
                                  {'H2': gas_stream.get_species_moles('H2'),
                                   'CO': gas_stream.get_species_moles('CO'),
                                   'CO2': gas_stream.get_species_moles('CO2')})
gas_stream.chemical_reaction('9-10')


# Data Analysis
def graphing():
    temperature_list = []
    for key in gas_stream.temperatures:
        temperature_list.extend(gas_stream.temperatures[key])

    pressure_list = []
    for key in gas_stream.pressures:
        pressure_list.extend(gas_stream.pressures[key])

    entropy_list = []
    for key in gas_stream.entropies:
        if key == '1-2':
            entropy_list.extend(map(add, gas_stream.entropies[key], water_in_stream.entropies[key]))
        elif key == '5-6':
            entropy_list.extend(map(add, gas_stream.entropies[key], water_separated_stream.entropies[key]))
        elif key in {'6-7', '7-8', '8-9', '9-10'}:
            # separated water isn't processed anymore, thus its current entropy is taken instead
            entropy_list.extend(map(lambda x: x + water_separated_stream.get_entropy(), gas_stream.entropies[key]))
        else:
            entropy_list.extend(gas_stream.entropies[key])

    debug_enthalpies = []
    debug_temps = []
    debug_key = '7-8'

    if debug_key == '1-2':
        debug_enthalpies.extend(map(add, gas_stream.enthalpies[debug_key], water_in_stream.enthalpies[debug_key]))
        debug_temps.extend(gas_stream.temperatures[debug_key])
    elif debug_key == '5-6':
        debug_enthalpies.extend(
            map(add, gas_stream.enthalpies[debug_key], water_separated_stream.enthalpies[debug_key]))
        debug_temps.extend(gas_stream.temperatures[debug_key])
    elif debug_key in {'6-7', '7-8', '8-9', '9-10'}:
        debug_enthalpies.extend(
            (map(lambda x: x + water_separated_stream.get_enthalpy(), gas_stream.enthalpies[debug_key])))
        debug_temps.extend(gas_stream.temperatures[debug_key])
    else:
        debug_enthalpies.extend(gas_stream.enthalpies[debug_key])
        debug_temps.extend(gas_stream.temperatures[debug_key])

    enthalpy_list = []
    for key in gas_stream.enthalpies:
        if key == '1-2':
            enthalpy_list.extend(map(add, gas_stream.enthalpies[key], water_in_stream.enthalpies[key]))
        elif key == '5-6':
            enthalpy_list.extend(map(add, gas_stream.enthalpies[key], water_separated_stream.enthalpies[key]))
        elif key in {'6-7', '7-8', '8-9', '9-10'}:
            # separated water isn't processed anymore, thus its current enthalpy is taken instead
            enthalpy_list.extend((map(lambda x: x + water_separated_stream.get_enthalpy(), gas_stream.enthalpies[key])))
        else:
            enthalpy_list.extend(gas_stream.enthalpies[key])

    # Plotting
    fig, ax = plt.subplots()
    ax.plot(np.array(enthalpy_list), np.array(temperature_list), color='#43a047')
    ax.plot(np.array(debug_enthalpies), np.array(debug_temps), color='#ff1744')
    ax.set_facecolor('#f5f5f5')
    ax.set_ylabel('Temperature [$K$]')
    ax.set_xlabel(f'Enthalpy [J]')
    ax.set_title('T-H diagram of  methanol synthesis')
    plt.text(-920000, 300, '1')
    plt.text(-770000, 600, '2')
    plt.text(-658000, 1200, '3')
    plt.text(-415000, 1200, '4')
    plt.text(-560000, 465, '5')
    plt.text(-618000, 434, '6')
    plt.text(-655000, 300, '7')
    plt.text(-640000, 436, '8')
    plt.text(-627000, 516, '9')
    plt.text(-693000, 516, '10')
    ax.legend(['$F_{CH_4} = 1$    $mol$'  '\n' + '$F_{CO_2} = 0.3$ $mol$' + '\n' + '$F_{H_{2}O} $  $= 2.5$ $mol$'])
    plt.savefig("T_H_diagram.png", dpi=1200)
    plt.show()

    fig, ax = plt.subplots()
    ax.plot(np.array(entropy_list), np.array(temperature_list), color='#e57373')
    ax.set_facecolor('#f5f5f5')
    ax.set_ylabel('Temperature [$K$]')
    ax.set_xlabel('Entropy [$J / K$]')
    ax.set_title('T-S diagram of  methanol synthesis')
    plt.text(383, 300, '1')
    plt.text(720, 600, '2')
    plt.text(868, 1200, '3')
    plt.text(1130, 1200, '4')
    plt.text(950, 465, '5')
    plt.text(797, 434, '6')
    plt.text(755, 300, '7')
    plt.text(720, 436, '8')
    plt.text(750, 516, '9')
    plt.text(610, 516, '10')
    ax.legend(['$F_{CH_4} = 1$    $mol$'  '\n' + '$F_{CO_2} = 0.3$ $mol$' + '\n' + '$F_{H_{2}O} $  $= 2.5$ $mol$'])
    plt.savefig("T_S_diagram.png", dpi=1200)
    plt.show()
    return


def thermal_efficiency():
    higher_heating_value_CH3OH = 0 - 726100  # J/mol

    efficiency = (gas_stream.get_species_moles('CH3OH') * higher_heating_value_CH3OH) / \
             ((gas_stream.enthalpies['1-2'][0] + water_in_stream.enthalpies['1-2'][0]) -
              (sum(filter(lambda x: x > 0, heats.values())) + sum(works.values())))
    return efficiency


def thermal_efficiency_limit():
    higher_heating_value_CH3OH = 0 - 726100  # J/mol

    def maximum_heat_exchange():
        if sum(heats.values()) - heats['3-4'] >= 0:
            return sum(heats.values())
        else:
            print(f'Anergy detected = {sum(heats.values()) - heats["3-4"]} J/s. Time to heat some buildings.')
            return heats['3-4']

    efficiency_recuperation = (gas_stream.get_species_moles('CH3OH') * higher_heating_value_CH3OH) / \
             ((gas_stream.enthalpies['1-2'][0] + water_in_stream.enthalpies['1-2'][0]) -
              (maximum_heat_exchange() + sum(works.values())))
    # print(gas_stream.get_species_moles('CH3OH') * higher_heating_value_CH3OH)
    # print((gas_stream.enthalpies['1-2'][0] + water_in_stream.enthalpies['1-2'][0]), (sum(heats.values()) - heats['3-4'] + sum(works.values())))
    return efficiency_recuperation


# Im a wrapper
def delta_enthalpy(stream, process_name):
    return stream.enthalpies[process_name][-1] - stream.enthalpies[process_name][0]


# Heat and Work calculation
heats = {'1-2': delta_enthalpy(gas_stream, '1-2') + delta_enthalpy(water_in_stream, '1-2'),
         # Add enthalpy change of both input streams
         '2-3': delta_enthalpy(gas_stream, '2-3'),
         '3-4': delta_enthalpy(gas_stream, '3-4'),
         '4-5': delta_enthalpy(gas_stream, '4-5'),
         '5-6': delta_enthalpy(gas_stream, '5-6') + (water_separated_stream.get_enthalpy() - 0),
         # Add enthalpies of separated water stream and gas stream
         '6-7': delta_enthalpy(gas_stream, '6-7'),
         '8-9': delta_enthalpy(gas_stream, '8-9'),
         '9-10': delta_enthalpy(gas_stream, '9-10')
         }
works = {'7-8': delta_enthalpy(gas_stream, '7-8')}


# Print results
print('Methane reforming temperature = ', T3, '\n')

# Ignore all species except syngas
x_CO = gas_stream.species_moles['5-6'][-1]['CO'] / sum(
    [gas_stream.species_moles['5-6'][-1]['CO'], gas_stream.species_moles['5-6'][-1]['CO2'], gas_stream.species_moles['5-6'][-1]['H2']])

x_CO2 = gas_stream.species_moles['5-6'][-1]['CO2'] / sum(
    [gas_stream.species_moles['5-6'][-1]['CO'], gas_stream.species_moles['5-6'][-1]['CO2'], gas_stream.species_moles['5-6'][-1]['H2']])

x_H2 = gas_stream.species_moles['5-6'][-1]['H2'] / sum(
    [gas_stream.species_moles['5-6'][-1]['CO'], gas_stream.species_moles['5-6'][-1]['CO2'], gas_stream.species_moles['5-6'][-1]['H2']])

S_module = (x_H2 - x_CO2) / (x_CO + x_CO2)
print(f"x_CO = {x_CO}, x_CO2 = {x_CO2}, x_H2 = {x_H2}", "\n"
      f'S_module ={S_module}\n')

print('@9 moles CO', gas_stream.species_moles['9-10'][0]['CO'], '\n'
      '@9 moles CO2', gas_stream.species_moles['9-10'][0]['CO2'], '\n'
      '@10 moles CO', gas_stream.species_moles['9-10'][-1]['CO'], '\n'
      '@10 moles CO2', gas_stream.species_moles['9-10'][-1]['CO2'], '\n'
      'X_CO+XO2 = ', 1 - ((gas_stream.species_moles['9-10'][-1]['CO'] + gas_stream.species_moles['9-10'][-1]['CO2']) / (gas_stream.species_moles['9-10'][0]['CO'] + gas_stream.species_moles['9-10'][0]['CO2'])), '\n')

print('heats = ', heats, '\nworks', works, '\n')

print(f'thermal efficiency limit = {thermal_efficiency_limit()}\n'
      f'thermal efficiency = {thermal_efficiency()}\n'
      f'methanol fraction of feed = {gas_stream.get_species_moles("CH3OH") / (F_CO2 + F_H2O + F_CH4)}', '\n'
      )

print(f"CO2   mfrac = {gas_stream.get_species_mole_fraction('CO2')}\n"
      f"CO    mfrac = {gas_stream.get_species_mole_fraction('CO')}\n"
      f"H2    mfrac = {gas_stream.get_species_mole_fraction('H2')}\n"
      f"CH3OH mfrac = {gas_stream.get_species_mole_fraction('CH3OH')}\n"
      )

cup_o_coffe_enthalpy = 4181 * 80 * 0.4  # Energy needed to produce one coffee
print(cup_o_coffe_enthalpy / -(sum(heats.values()) - heats["3-4"]))

graphing()