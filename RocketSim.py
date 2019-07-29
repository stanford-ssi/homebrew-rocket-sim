import numpy as np
import matplotlib.pyplot as plt
import time
from NRLMSISE00.nrlmsise_00_header import *
from NRLMSISE00.nrlmsise_00 import *
import xml.etree.ElementTree as ET
import quaternion
#All SI units

rocket_mass = 16 #not including motor
time_step = 0.05
c_d = 0.5
area = 0.00258064


altitude = 0
g= 9.8
time = 0
velocity = 0
positions = []
velocities = []
accelerations = []
times = []

def parse_thrustcurve():
    tree = ET.parse('Cesaroni_N5800.xml')
    root=tree.getroot()
    motor_masses = []
    motor_thrusts = []
    motor_times = []
    for child in root[0][0][1]:
        motor_masses.append(float(child.attrib['m'])/1000)
        motor_thrusts.append(float(child.attrib['f']))
        motor_times.append(float(child.attrib['t']))
    burn_time = float(root[0][0].attrib['burn-time'])
    sample_times = time_step*np.arange(0,np.ceil(burn_time/time_step))
    motor_masses = np.interp(sample_times,motor_times,motor_masses)
    motor_thrusts = np.interp(sample_times,motor_times,motor_thrusts)
    return (motor_masses, motor_thrusts, sample_times)

def get_density_for_altitude(altitude):
    output = nrlmsise_output()
    Input = nrlmsise_input()
    flags = nrlmsise_flags()
    aph = ap_array()
    flags.switches[0] = 0 #output in MKS
    Input.alt=altitude/1000 #convert to km
    Input.g_lat = 32.990371 #Spaceport america
    Input.g_long= -106.975116
    for i in range(7):
        aph.a[i]=100
    for i in range(1, 24):
        flags.switches[i]=1
    gtd7(Input, flags, output)
    return output.d[5]*1000 #total air density in KG/M^3

motor_masses, motor_thrusts, sample_times = parse_thrustcurve()
while True:
    time += time_step

    thrust=0
    mass = rocket_mass

    idx = -1
    try:
        idx = np.where(np.isclose(sample_times,time))[0][0]
    except IndexError:
        idx = -1

    if idx!=-1:
        mass += motor_masses[idx]
        thrust += motor_thrusts[idx]
    weight = mass*g

    drag_force = 0.5*get_density_for_altitude(altitude)*(velocity**2)*area*c_d

    net_force = thrust-weight-drag_force
    acc = net_force/mass
    delta_v = acc*time_step
    velocity += delta_v
    delta_x = velocity*time_step
    altitude += delta_x
    times.append(time)
    positions.append(altitude)
    velocities.append(velocity)
    accelerations.append(acc)
    if altitude<-80:
        break
plt.plot(times,positions)
plt.show()
