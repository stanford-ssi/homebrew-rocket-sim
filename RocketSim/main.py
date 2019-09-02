import time

import numpy as np
import quaternion
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ET

from NRLMSISE00.nrlmsise_00_header import *
from NRLMSISE00.nrlmsise_00 import *
from DigitalDATCOM.datcom_lookup import lookup


def parse_thrust_curve(fname, time_step):
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
    return motor_masses, motor_thrusts, sample_times


def get_atmospheric_properties(altitude):
    output = nrlmsise_output()
    model_input = nrlmsise_input()
    flags = nrlmsise_flags()
    aph = ap_array()
    flags.switches[0] = 0  # output in MKS
    model_input.alt = altitude / 1000  # convert to km
    model_input.g_lat = 32.990371  # Spaceport America
    model_input.g_long = -106.975116
    for i in range(7):
        aph.a[i] = 100
    for i in range(1, 24):
        flags.switches[i] = 1
    gtd7(model_input, flags, output)
    return output.d[5] * 1000, output.t[1]  # total air density in KG/M^3


times, positions = [], []

dt = 0.1                            # timestep, s
curve_vals = parse_thrust_curve(
    'Cesaroni_N5800.xml', dt)
M_ms, T_ms, _ = curve_vals          # rocket masses, kg,
                                    # ... and rocket thrusts, N
curve_index = 0                     # index for previous two

M_r = 16                            # rocket mass, kg
X_cp = 0.25                         # rocket center of pressure
                                    # ...from nose tip, m
X_cm = 0.5                          # rocket center of gravity
                                    # ...from nose tip, m
M_E = 5.974E24                      # Earth mass, kg
r_E = 6378100                       # Earth radius, m

Y_A0 = np.array([1.0, 0.0, 0.0])    # yaw axis
P_A0 = np.array([0.0, 1.0, 0.0])    # pitch axis
R_A0 = np.array([0.0, 0.0, 1.0])    # roll axis
m = np.diag([1.0, 1.0, 0.0])        # convenience matrix for computing tau_da 

t = 0.0                             # time, s
P = np.array([0.0, 0.0, 0.0])       # momentum, kg * m/s
L = np.array([0.0, 0.0, 0.0])       # angular momentum, kg * m^2/s
Q = np.quaternion(1, 0, 0, 0)       # rotation quaternion
                                    # ...(rotation relative to pointing
                                    # ...directly upward)
X = np.array([0.0, 0.0, 0.0])       # position relative to ground, m
I_0 = np.diag([1.0, 1.0, 1.0])      # moments of inertia, kg * m^2
                                    # ...(in x, y, and z direction resp.)
W = 0.0                             # wind velocity, m/s

while True:
    times.append(t)
    positions.append(X)

    M = M_r                                    # total mass
    T = 0                                      # thrust
    if curve_index < len(M_ms):                # if motor isn't spent:
        M += M_ms[curve_index]                 # add motor mass
        T = T_ms[curve_index]                  # set thrust

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

    # TODO: calculate drag coefficients as in 1D case

    P += F * dt             # update momentum
    L += tau * dt           # update angular momentum
    Q.w += s_dot * dt       # update real part of quaternion
    Q.x += v_dot[0] * dt    # update rest of quaternion
    Q.y += v_dot[1] * dt
    Q.z += v_dot[2] * dt
    X += X_dot * dt         # update position
    t += dt                 # update time

    if X[2] < 0:  # if our z coordinate is underground
        times.append(t)
        positions.append(X)
        break


'''
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
'''

plt.plot(times, positions)
plt.show()
