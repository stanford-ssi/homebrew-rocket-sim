import time

import numpy as np
import quaternion
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ET

from NRLMSISE00.nrlmsise_00_header import *
from NRLMSISE00.nrlmsise_00 import *
from DigitalDATCOM.datcom_lookup import lookup


def parse_thrust_curve(fname, time_step):
    """Parses a thrust curve from an XML file.

    Args:
        fname: filename of XML file with motor data
        time_step: amount of time between time steps (s)
    Returns:
        motor_masses: list of motor masses at each time step (kg)
        motor_thrusts: list of thrusts at each time step (N)
        sample_times: times at which each sample is taken (s)
    """
    # Get the root of the XML file
    tree = ET.parse(fname)
    root = tree.getroot()

    # Create empty lists for data we're collecting
    motor_masses = []
    motor_thrusts = []
    motor_times = []

    # For each entry in our XML file
    for child in root[0][0][1]:
        # Append the data in that entry
        motor_masses.append(float(child.attrib['m']) / 1000)
        motor_thrusts.append(float(child.attrib['f']))
        motor_times.append(float(child.attrib['t']))

    # Get the burn time of the motor
    burn_time = float(root[0][0].attrib['burn-time'])

    # Calculate the times we want to sample at
    sample_times = time_step * np.arange(0, np.ceil(burn_time / time_step))

    # Interpolate motor masses and thursts to those time steps
    motor_masses = np.interp(sample_times, motor_times, motor_masses)
    motor_thrusts = np.interp(sample_times, motor_times, motor_thrusts)

    # Return it all
    return motor_masses, motor_thrusts, sample_times


def get_atmospheric_properties(altitude):
    """Parses a thrust curve from an XML file.

    Args:
        altitude: altitude to get data at (m)
    Returns:
        rho: air density at altitude (kg/m^3)
        temp: temperature (degrees C)
    """
    # Get output data, input data, and flags from NRLMSISE library
    output = nrlmsise_output()
    model_input = nrlmsise_input()
    flags = nrlmsise_flags()

    # Create an array to hold Ap-indexes for past (geomagnetic activity)
    aph = ap_array()

    # Ensure output is in MKS units
    flags.switches[0] = 0

    # Convert altitude to km
    model_input.alt = altitude / 1000

    # Set location to Spaceport America
    model_input.g_lat = 32.990371  # Spaceport America
    model_input.g_long = -106.975116

    # For simplicity, assume the Ap-index has been 100 for the last while
    for i in range(7):
        aph.a[i] = 100

    # Use all default settings
    for i in range(1, 24):
        flags.switches[i] = 1

    # Run the model
    gtd7(model_input, flags, output)

    # Get air density in kg/m^3
    rho = output.d[5] * 1000

    # Get temperature
    temp = output.t[1]
    return rho, temp


def safe_normalize(v):
    """Normalizes a vector, unless it's zero.
    In that case, it returns zero.

    Args:
        v: input vector
    Returns:
        v_hat: normalized vector
    """
    # Get the norm of the vector and divide it out unless zero
    norm = np.linalg.norm(v)
    return v / norm if norm > 0 else v



# Create lists for times and positions at each timestep
times, positions = [], []

# Set timestep and load thrust curve
dt = 0.1                            # timestep, s
curve_vals = parse_thrust_curve(
    'Cesaroni_N5800.xml', dt)
M_ms, T_ms, _ = curve_vals          # rocket masses, kg,
                                    # ... and rocket thrusts, N
curve_index = 0                     # index for previous two


# Set masses, distances, areas, etc.
M_r = 16                            # rocket mass, kg
X_cp = 0.25                         # rocket center of pressure
                                    # ...from nose tip, m
X_cm = 0.5                          # rocket center of gravity
                                    # ...from nose tip, m
A_RB = 0.0013                       # rocket cross-sectional area, m^2
M_E = 5.974E24                      # Earth mass, kg
r_E = 6378100                       # Earth radius, m
G = 6.673E-11                       # Gravitational constant, m^3/(kg * s)

# Set yaw axis as first coordinate, pitch axis as second, roll axis as third
Y_A0 = np.array([1.0, 0.0, 0.0])    # yaw axis
P_A0 = np.array([0.0, 1.0, 0.0])    # pitch axis
R_A0 = np.array([0.0, 0.0, 1.0])    # roll axis
m = np.diag([1.0, 1.0, 0.0])        # convenience matrix for computing tau_da

# Set current time, momentum, angular momentum, rotation to zero
t = 0.0                             # time, s
P = np.array([0.0, 0.0, 0.0])       # momentum, kg * m/s
L = np.array([0.0, 0.0, 0.0])       # angular momentum, kg * m^2/s
Q = np.quaternion(1, 0, 0, 0)       # rotation quaternion
                                    # ...(rotation relative to pointing
                                    # ...directly upward)

# Set initial position, moment of inertia, wind, aerodynamic damping
X = np.array([0.0, 0.0, 26000.0])   # position relative to ground, m
I_0 = np.diag([1.0, 1.0, 1.0])      # moments of inertia, kg * m^2
                                    # ...(in x, y, and z direction resp.)
W = 0.0                             # wind velocity, m/s
C_da = 0.0                          # damping coefficient (TODO: calculate)


# Loop until we hit the ground
while True:
    # Record current position and time
    times.append(t)
    positions.append(tuple(X))

    # Set current mass to rocket + remaining fuel, get thrust
    M = M_r                                    # total mass
    T = 0                                      # thrust
    if curve_index < len(M_ms):                # if motor isn't spent:
        M += M_ms[curve_index]                 # add motor mass
        T = T_ms[curve_index]                  # set thrust

    # Get velocity (just momentum over mass)
    X_dot = P / M                              # derivative of position

    # Get rotation and based on that, unit vector along which rocket is pointing
    R = quaternion.as_rotation_matrix(Q)       # rotation matrix
    R_A = np.dot(R, R_A0.T)                    # unit vector in roll axis

    # Get angular velocity and how fast we're rotating
    omega = np.linalg.multi_dot(
        (R, np.linalg.inv(I_0), R.T, L.T))     # angular velocity
    s, v = Q.w, np.array([Q.x, Q.y, Q.z])      # components of quaternion
    s_dot = 0.5 * np.dot(omega, v)             # derivative of real part
    v_dot = 0.5 * (
        s * omega + np.cross(omega, v))        # derivative of the rest

    # Compute how fast center of pressure is moving
    # (this is just how fast center of mass moves, plus a term
    # to account for rotation)
    V_cm = X_dot + W                           # velocity of center of mass
    omega_hat = safe_normalize(omega)          # normalized angular velocity
    X_bar = np.abs(X_cp - X_cm)
    V_omega = X_bar * np.sin(
        np.arccos(np.dot(R_A, omega_hat)) *
                  np.cross(R_A, omega))        # velocity of center of pressure
                                               # ...due to angular velocity
    V = V_cm + V_omega                         # total velocity

    # Find the angle of attack
    V_hat = safe_normalize(V)                  # normalized velocity
    alpha = np.arccos(np.dot(V_hat, R_A))      # angle of attack
                                               # TODO: fix reference angles to mesh with
                                               # DATCOM (default angle should be pi/2)

    # Compute thrust force, gravitational force
    F_T = T * R_A                              # force due to thrust
                                               # ...(note by differing
                                               # ...conventions we lose a
                                               # ...minus sign from the paper)

    z = X[2]                                   # z-coordinate of position
    g = G * M_E / (r_E + z) ** 2               # gravitational acceleration
                                               # ...(note typo in paper omits G)
    F_g = np.array([0.0, 0.0, 0.0])
    if z > 0:                                  # accounting for normal force
        F_g[2] = -M * g                        # force due to gravity

    # Get mach number
    rho, temp = get_atmospheric_properties(z)  # air density
    mach = np.linalg.norm(V) / (
        20.05 * np.sqrt(temp))

    # Round mach number up or down to avoid probelmatic region
    # which DATCOM can't model well (actually errors in the
    # regions we avoid!)
    if mach<=0.1:
        mach=0.1
    elif mach>=0.6 and mach <= 1:
        mach = 0.59
    elif mach >= 1 and mach <= 1.4:
        mach = 1.41                            # mach number (TODO: match paper)

    # If we're moving, get our new aerodynamic coefficients
    # as well as our new aerodynamic forces
    if mach > 0:                               # if we're moving
        angle_of_attack = 0.0
        # print(mach, angle_of_attack, z)
        lookup_results = lookup(
            [mach], [angle_of_attack], [z],
            X_cm, M)                               # DATCOM lookup results
        coeffs = list(lookup_results.values())[0]  # coefficients from DATCOM
        C_A, C_N = coeffs['CA'], coeffs['CN']      # axial and normal aerodynamic
                                                   # ...coefficients
                                                   # TODO: gracefully handle NaNs
                                                   # (altitude = 0)
        C_D = coeffs['CD']                         # Vehicle drag coefficient
        C_L = coeffs['CL']                         # Vehicle lift coefficient
        C_M = coeffs['CM']                         # Vehicle pitching-moment coeff. 
        C_L_A = coeffs['CLA']                      # Derivative of lift coeff. w.r.t. alpha
        C_M_A = coeffs['CMA']                      # Derivative of pitching moment coeff. w.r.t. alpha
        C_Y_B = coeffs['CYB']                      # Side force coeff. derivative w.r.t. sideslip angle
        C_N_B = coeffs['CNB']                      # Derivative of yawing-moment coeff. w.r.t. sideslip angle
        C_L_B = coeffs['CLB']                      # Derivative of rolling-moment coeff. w.r.t.sideslip angle




        #print(C_A, C_N)
        F_A_mag = 0.5 * rho * V ** 2 * A_RB * C_A  # magnitude of axial
                                                   # ...aerodynamic force
        F_A = -F_A_mag * R_A                       # axial aerodynamic force
        F_N_mag = 0.5 * rho * V ** 2 * A_RB * C_N  # magnitude of normal
                                                   # ...aerodynamic force
        F_N = F_N_mag * np.cross(
            R_A, np.cross(R_A, V_hat))             # normal aerodynamic force
    else:                                          # don't DATCOM at mach zero
        F_A = np.array([0.0, 0.0, 0.0])            # axial aerodynamic force
        F_N_mag = 0.0
        F_N = np.array([0.0, 0.0, 0.0])            # normal aerodynamic force

    # Compute total force (just add it up!)
    F = F_T + F_g + F_A + F_N                  # total force

    # Compute torques from aerodynamic force and thrust damping, and add them
    tau_N = (F_N_mag * X_bar *
             np.cross(R_A, V_hat))             # torque caused by normal force
    tau_da = -C_da * np.linalg.multi_dot(
        (R, m, np.linalg.inv(R), omega))       # torque due to thrust damping
                                               # ...(thrust slowing the
                                               # ...rocket's rotation)
    tau = tau_N + tau_da                       # total torque

    # Update all of our state variables by adding their derivative * dt
    P += F * dt             # update momentum
    L += tau * dt           # update angular momentum
    Q.w += s_dot * dt       # update real part of quaternion
    Q.x += v_dot[0] * dt    # update rest of quaternion
    Q.y += v_dot[1] * dt
    Q.z += v_dot[2] * dt
    X += X_dot * dt         # update position
    t += dt                 # update time
    curve_index += 1

    # If we're underground, end the simulation (this isn't "earthshot")
    z = X[2]                # get the z-coordinate
    if z < 0:               # if it's underground
        times.append(t)
        positions.append(tuple(X))
        break

# Plot everything!
plt.plot(times, positions)
plt.show()
