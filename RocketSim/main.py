import time

import numpy as np
import math
import quaternion
import matplotlib.pyplot as plt
from NRLMSISE00.nrlmsise_00_header import *
from NRLMSISE00.nrlmsise_00 import *
from DigitalDATCOM.datcom_lookup import lookup
import xml.etree.ElementTree as ET
from scipy.interpolate import griddata

# All SI units
rocket_mass = 16  # not including motor
time_step = 3
c_d = 0.5
area = 0.00258064
altitude = 0
g = 9.8
time, velocity, acc = 0, 0, 0
positions, velocities, accelerations, times, = [], [], [], []
mach_keys, alpha_keys, alt_keys =  [], [], []
coeff_data = None

def parse_thrustcurve(fname='Cesaroni_N5800.xml'):
    tree = ET.parse(fname)
    root = tree.getroot()
    motor_masses = []
    motor_thrusts = []
    motor_times = []
    for child in root[0][0][1]:
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

def loadDATCOM(fname='LookupTableTest.npz'):
    """ Load in the .npz file with the aerodynamic coefficients, and sort the keys for easy access during timestepping.
    """
    global mach_keys
    global alpha_keys
    global alt_keys
    with np.load(fname,allow_pickle=True) as data:
        global coeff_data
        coeff_data=data['arr_0'][()]
    mach_keys, alpha_keys, alt_keys = zip(*coeff_data.keys())

    mach_keys = np.unique(mach_keys)
    alpha_keys = np.unique(alpha_keys)
    alt_keys = np.unique(alt_keys)

def find_nearest(array,value):
    """ Helper method that finds the nearest value in a sorted array and returns the index.
    """
    idx = np.searchsorted(array, value, side="left")
    if idx > 0 and (idx == len(array) or math.fabs(value - array[idx-1]) < math.fabs(value - array[idx])):
        return idx-1,array[idx-1]
    else:
        return idx,array[idx]

def coeff_for_conditions(mach, alpha, alt):
    """ Helper method to retrieve coefficients. Given keys, it finds the nearest corresponding precalculated DATCOM value, and does a simple linear interpolation between the two.
    """
    nearest_mach_idx,nearest_mach = find_nearest(mach_keys,mach)
    nearest_alpha_idx,nearest_alpha = find_nearest(alpha_keys,alpha)
    nearest_alt_idx,nearest_alt = find_nearest(alt_keys,alt)

    if nearest_mach>mach:
        mach_1 = mach_keys[nearest_mach_idx-1]
        mach_2 = mach_keys[nearest_mach_idx]
    else:
        mach_1 = mach_keys[nearest_mach_idx]
        mach_2 = mach_keys[nearest_mach_idx+1]

    if nearest_alpha>alpha:
        alpha_1 = alpha_keys[nearest_alpha_idx-1]
        alpha_2 = alpha_keys[nearest_alpha_idx]
    else:
        alpha_1 = alpha_keys[nearest_alpha_idx]
        alpha_2 = alpha_keys[nearest_alpha_idx+1]

    if nearest_alt>alt:
        alt_1 = alt_keys[nearest_alt_idx-1]
        alt_2 = alt_keys[nearest_alt_idx]
    else:
        alt_1 = alt_keys[nearest_alt_idx]
        alt_2 = alt_keys[nearest_alt_idx+1]

    x = np.array([mach_1,mach_2])
    y = np.array([alpha_1,alpha_2])
    z = np.array([alt_1,alt_2])
    c1 = coeff_data[(mach_1,alpha_1,alt_1)]
    c2 = coeff_data[(mach_2,alpha_2,alt_2)]

    interpolated = {key: 0 for key in c1.keys()}

    for key in c1.keys():
        f = np.array([c1[key],c2[key]])
        dfdx = (f[1]-f[0])/(x[1]-x[0])
        dfdy = (f[1]-f[0])/(y[1]-y[0])
        dfdz = (f[1]-f[0])/(z[1]-z[0])
        print f[0], f[1], x[0],x[1], dfdx, mach
        interpolated[key]=f[0]+(dfdx*(mach-x[0]))+(dfdy*(alpha-y[0]))+(dfdz*(alt-z[0]))
    return interpolated

loadDATCOM()
print coeff_for_conditions(0.52,0.2,26500)
'''

motor_masses, motor_thrusts, sample_times = parse_thrustcurve()
loadDATCOM()


while True:
    times.append(time)
    positions.append(altitude)
    velocities.append(velocity)
    accelerations.append(acc)
    time += time_step

    thrust = 0
    mass = rocket_mass

    try:
        idx = np.where(np.isclose(sample_times, time))[0][0]
    except IndexError:
        idx = None

    if idx is not None:
        mass += motor_masses[idx]
        thrust += motor_thrusts[idx]

    weight = mass * g
    density, temperature = get_d_t_for_altitude(altitude)
    mach = velocity / (20.05 * np.sqrt(temperature))

    mach = 0.4
    alpha = 1.0
    cg = 1.2

    coeffs = coeff_for_conditions(mach, alpha, altitude, cg, mass)

    drag_force = 0.5 * density * (velocity ** 2) * area * coeffs['CD']

    net_force = thrust - weight - drag_force
    acc = net_force / mass
    delta_v = acc * time_step
    velocity += delta_v
    delta_x = velocity * time_step
    altitude += delta_x
    print(altitude)

    if altitude < 0:
        times.append(time)
        positions.append(altitude)
        velocities.append(velocity)
        accelerations.append(acc)
        break

plt.plot(times,positions)
plt.show()
'''
