# Ideal thermodynamic calculation tool base
# Calculations in SI units unless otherwise specified
import numpy as np
import math

# Calculates the vapor temperature of an ideal solution with the Clausius-Clapeyron equation
def water_vapour_temperature(p1):
    L = 40800  # J/mol
    R = 8.314462  # J/mol K
    T0 = 273.15 + 100  # K
    p0 = 101325  # Pa

    T1 = 1 / ((0 - np.log(p1 / p0) * (R / L)) + (1 / T0))
    return (T1)


# Calculates temperature of an isentropic relationship between two pressure states
def isentropic_compression(T0, p0, p1, gamma):
    T1 = T0 * (p0 / p1) ** ((1 - gamma) / gamma)
    return T1


def counter_flow_heat_exchange(massflow_a, cp_a, inlet_temperature_a, massflow_b, cp_b, inlet_temperature_b):
    # a is the initial hotter stream, be the initial colder stream
    # alpha =
    # efficiency = (1 - math.exp(-alpha)) / (1-(((massflow_a * cp_a) / (massflow_b * cp_b)) * math.exp(-alpha)))
    efficiency = 1
    outlet_temperature_stream_a = (efficiency * (inlet_temperature_b -inlet_temperature_a)) + inlet_temperature_a
    outlet_temperature_stream_b = (((massflow_a * cp_a) / (massflow_b * cp_b)) * efficiency * (inlet_temperature_a - inlet_temperature_b)) + inlet_temperature_b
    Q = massflow_a * cp_a * (outlet_temperature_stream_a - inlet_temperature_a)

    return outlet_temperature_stream_a, outlet_temperature_stream_b, Q
