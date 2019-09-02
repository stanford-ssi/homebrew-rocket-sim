import time

import numpy as np
# import quaternion
import matplotlib.pyplot as plt
from NRLMSISE00.nrlmsise_00_header import *
from NRLMSISE00.nrlmsise_00 import *
from DigitalDATCOM.datcom_lookup import lookup
import xml.etree.ElementTree as ET
from scipy.interpolate import griddata
import scipy.integrate as ode

def parse_thrustcurve(xmlroot):
    motor_masses = []
    motor_thrusts = []
    motor_times = []
    for child in xmlroot[0][0][1]:
        motor_masses.append(float(child.attrib['m']) / 1000)
        motor_thrusts.append(float(child.attrib['f']))
        motor_times.append(float(child.attrib['t']))
    burn_time = float(root[0][0].attrib['burn-time'])
    sample_times = time_step * np.arange(0, np.ceil(burn_time / time_step))
    motor_masses = np.interp(sample_times, motor_times, motor_masses)
    motor_thrusts = np.interp(sample_times, motor_times, motor_thrusts)
    return (motor_masses, motor_thrusts, sample_times)


def get_d_t_for_altitude(altitude):
    output = nrlmsise_output()
    model_input = nrlmsise_input()
    flags = nrlmsise_flags()
    aph = ap_array()
    flags.switches[0] = 0  # output in MKS
    model_input.alt = altitude / 1000  # convert to km
    model_input.g_lat = 32.990371  # Spaceport america
    model_input.g_long= -106.975116
    for i in range(7):
        aph.a[i] = 100
    for i in range(1, 24):
        flags.switches[i] = 1
    gtd7(model_input, flags, output)
    return output.d[5] * 1000, output.t[1]  # total air density in KG/M^3


def step_rocket(y, t):
    alt, vel = y

    thrust = 0  # Um... so are you just reseting this to 0 every time?
    mass = rocket_mass

    try:
        idx = np.where(np.isclose(sample_times, t))[0][0]
    except IndexError:
        idx = None

    if idx is not None:
        mass += motor_masses[idx]
        thrust += motor_thrusts[idx]

    weight = mass * g
    density, temperature = get_d_t_for_altitude(altitude)
    mach = vel / (20.05 * np.sqrt(temperature))

    # mach = 0.4
    # alpha = 1.0
    cg = 1.2

    # CD,CL,CM,CN,CA,XCP,CLA,CMA,CYB,CNB,CLB = coeff_for_conditions(mach, alpha, altitude, cg, mass)

    # We could look up the coefficients by mach, alpha, and altitude again,
    # but floating point error in DATCOM makes the values in the keys in coeffs
    # differ from mach, alpha, and altitude
    coeffs = list(lookup([mach], [alpha], [altitude], cg, mass).values())[0]
    drag_force = 0.5 * density * (vel ** 2) * area * coeffs['CD']

    net_force = thrust - weight - drag_force
    acc = net_force / mass

    #  doin' the integrating thing
    dydt = vel, acc

    # Bookkeeping
    times.append(t)
    positions.append(alt)
    velocities.append(vel)
    accelerations.append(acc)

    return dydt


# All SI units
rocket_mass = 16  # not including motor ... ? THOMAS: With motor rocket mass is 19kg
time_step = 3 # I have big questions about this timestep. The motor burns out in 3 seconds
c_d = 0.5
area = 0.00258064
altitude = 0 # So assuming a ground launch? Or something else?
g = 9.8 # NOTE: Upgrade to update based on altitude.
time, velocity, acc = 0, 0, 0
alpha = 0
positions, velocities, accelerations, times = [], [], [], []

# y0 = [altitude, velocity, acc, rocket_mass, time]
y0 = altitude, velocity

fname='Cesaroni_N5800.xml'
tree = ET.parse(fname)
root = tree.getroot()

motor_masses, motor_thrusts, sample_times = parse_thrustcurve(root)

t = np.linspace(0, 100)
ode.odeint(step_rocket, y0, t)
# loadDATCOM()

# while True:
#     times.append(time)
#     positions.append(altitude)
#     velocities.append(velocity)
#     accelerations.append(acc)
#     time += time_step
#
#     thrust = 0 #Um... so are you just reseting this to 0 every time?
#     mass = rocket_mass
#
#     try:
#         idx = np.where(np.isclose(sample_times, time))[0][0]
#     except IndexError:
#         idx = None
#
#     if idx is not None:
#         mass += motor_masses[idx]
#         thrust += motor_thrusts[idx]
#
#     weight = mass * g
#     density, temperature = get_d_t_for_altitude(altitude)
#     mach = velocity / (20.05 * np.sqrt(temperature))
#
#     mach = 0.4
#     alpha = 1.0
#     cg = 1.2
#
#     # CD,CL,CM,CN,CA,XCP,CLA,CMA,CYB,CNB,CLB = coeff_for_conditions(mach, alpha, altitude, cg, mass)
#
#     # We could look up the coefficients by mach, alpha, and altitude again,
#     # but floating point error in DATCOM makes the values in the keys in coeffs
#     # differ from mach, alpha, and altitude
#     coeffs = list(lookup([mach], [alpha], [altitude], cg, mass).values())[0]
#     drag_force = 0.5 * density * (velocity ** 2) * area * coeffs['CD']
#
#     net_force = thrust - weight - drag_force
#     acc = net_force / mass
#     delta_v = acc * time_step
#     velocity += delta_v
#     delta_x = velocity * time_step
#     altitude += delta_x
#     print(altitude)
#
#     if altitude < 0:
#         times.append(time)
#         positions.append(altitude)
#         velocities.append(velocity)
#         accelerations.append(acc)
#         break



plt.plot(times,positions)
plt.show()
