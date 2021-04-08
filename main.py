# MethanolSynthesis from Methane, Water and CarbonDioxide
# Yeah!

import math
import numpy as np
import scipy as sc
import matplotlib.pyplot as plt
import cantera as ct

T1 = 300  # K
p1 = 15 * 100000  # Pa
FCH4 = 1  # mol/s
FCO2 = 0.3  # mol/s
FH2O = 2.5  # mol/s

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
temp.X = {'CH4': FCH4, 'CO2': FCO2, 'H2O': FH2O}  # Important to set before .TP, otherwise pressure will change!
temp.TP = T2, p2
CO2_CH4_H20 = ct.Mixture([(temp, (FCH4 + FCO2 + FH2O))])

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
    i = 1 - (CO2_CH4_H20.species_moles[13] / FCH4)  # handwavy assumed to take F instead of n
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


def adiabatic_compression(T0, p0, p1, gamma):  # ideal isentropic compression. When is adiabaitc = reversible ok. Isnt isentropic = isotherm and adiabatic????
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

X_CO2_CO = methanolreactorphase.X[19]/(methanolreactorphase.X[13] + methanolreactorphase.X[14]+methanolreactorphase.X[19])
