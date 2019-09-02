import time

import numpy as np
import quaternion
import matplotlib.pyplot as plt
from NRLMSISE00.nrlmsise_00_header import *
from NRLMSISE00.nrlmsise_00 import *
from DigitalDATCOM.datcom_lookup import lookup
import xml.etree.ElementTree as ET
from scipy.interpolate import griddata


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


positions, velocities, accelerations, times = [], [], [], []
motor_masses, motor_thrusts, sample_times = parse_thrustcurve()

M = 16                              # rocket mass, kg
M_E = 5.974E24                      # Earth mass, kg
r_E = 6378100                       # Earth radius, m
m = np.diag([1.0, 1.0, 0.0])        # convenience matrix for computing tau_da 

Y_A0 = np.array([1.0, 0.0, 0.0])    # yaw axis
P_A0 = np.array([0.0, 1.0, 0.0])    # pitch axis
R_A0 = np.array([0.0, 0.0, 1.0])    # roll axis

t = 0.0                             # time, s
P = np.array([0.0, 0.0, 0.0])       # momentum, kg * m/s
L = np.array([0.0, 0.0, 0.0])       # angular momentum, kg * m^2/s
Q = np.quaternion(1, 0, 0, 0)       # rotation quaternion
                                    # ...(rotation relative to pointing
                                    # ...directly upward)
X = np.array([0.0, 0.0, 0.0])       # position relative to ground, m
I_0 = np.diag([1.0, 1.0, 1.0])      # moments of inertia, kg * m^2
W = 0.0                             # wind velocity, m/s

# TODO: define X_cp, X_cm, T

while True:
    X_dot = P / M                              # derivative of position
    R = Q.as_rotation_matrix()                 # rotation matrix
    R_A = np.dot(R, R_A0.T)                    # unit vector in roll axis 
    omega = np.linalg.multi_dot(
        R, np.linalg.inv(I_0), R.T, L.T)       # angular velocity
    s, v = Q.w, np.array([Q.x, Q.y, Q.z])      # components of quaternion
    s_dot = 0.5 * np.dot(omega, v)             # derivative of real part
    v_dot = 0.5 * (
        s * omega + np.cross(omega, v))        # derivative of the rest

    V_cm = X_dot + W                           # velocity of center of mass
    omega_hat = omega / np.linalg.norm(omega)  # normalized angular velocity
    X_bar = np.abs(X_cp - X_cm)
    V_omega = X_bar * np.sin(
        np.arccos(np.dot(R_A, omega_hat)) *
                  np.cross(R_A, omega))        # velocity of center of pressure
                                               # ...due to angular velocity
    V = V_cm + V_omega                         # total velocity
    V_hat = V / np.linalg.norm(V)              # normalized velocity
    alpha = np.arccos(np.dot(V_hat, R_A))      # angle of attack

    F_T = -T * R_A                             # force due to thrust
    g = M_E / (r_E + z) ** 2                   # gravitational acceleration
    F_g = np.array([0, 0, -M * g])             # force due to gravity
    F_A_mag = 0.5 * rho * V ** 2 * A_RB * C_A  # magnitude of axial
                                               # ...aerodynamic force
    F_A = -F_A_mag * R_A                       # axial aerodynamic force

    F_N_mag = 0.5 * rho * V ** 2 * A_RB * C_N  # magnitude of normal 
                                               # ...aerodynamic force
    F_N = F_N_mag * np.cross(
        R_A, np.cross(R_A, V_hat))             # normal aerodynamic force
    F = F_T + F_g + F_A + F_n                  # total force

    tau_N = (F_N_mag * X_bar *
             np.cross(R_A, V_hat))             # torque caused by normal force 
    tau_da = -C_da * np.linalg.multi_dot(
        R, m, np.linalg.inv(R), omega)         # torque due to thrust damping
                                               # ...(thrust slowing the
                                               # ...rocket's rotation)
    tau = tau_N + tau_da                       # total torque



'''
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

    # We could look up the coefficients by mach, alpha, and altitude again,
    # but floating point error in DATCOM makes the values in the keys in coeffs
    # differ from mach, alpha, and altitude
    coeffs = list(lookup([mach], [alpha], [altitude], cg, mass).values())[0]
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

'''

plt.plot(times,positions)
plt.show()
