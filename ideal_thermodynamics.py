#Ideal thermodynamic calculation tool base
#Calculations in SI units unless otherwise specified

import numpy as np

#Calculates the vapor temperature of an ideal solution with the Clausius-Clapeyron equation
def water_vapour_temperature(p1):
    L = 40800  # J/mol
    R = 8.314462  # J/mol K
    T0 = 273.15 + 100  # K
    p0 = 101325  # Pa

    T1 = 1 / ((0 - np.log(p1 / p0) * (R / L)) + (1 / T0))
    return (T1)


#Calculates temperature of an isentropic relationship between two pressure states
def isentropic_compression(T0, p0, p1, gamma):
    T1 = T0 * (p0 / p1) ** ((1 - gamma) / gamma)
    return T1

