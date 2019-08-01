import subprocess
import numpy as np
import os


def lookup(mach, alpha, alt, cg, mass):
    with open('DigitalDATCOM/datcom_template.txt', 'r') as f:
        datcom_input = f.read()

    replacements = {
        'INSERT_MACH': str(mach),
        'INSERT_ALPHA': str(alpha),
        'INSERT_CG': str(cg),
        'INSERT_WEIGHT': str(mass)
    }
    for key, value in replacements.items():
        datcom_input.replace(key, value)

    with open('DigitalDATCOM/current.dcm', 'w') as f:
        f.write(datcom_input)

    command = 'echo DigitalDATCOM/current.dcm | DigitalDATCOM/datcom'
    with open('DigitalDATCOM/datcom_log.txt', 'w') as f:
        subprocess.call(command, shell=True, stdout=f)

    with open('DigitalDATCOM/datcom.out', 'r') as f:
        datcom_output = f.read()
    start, end = datcom_output.find('0\n'), datcom_output.find('0*** VEHICLE')
    table = np.array(datcom_output[start:end].split('\n'))
    str_coeffs = np.array(table[-2].split())

    coeffs = []
    for str_coeff in str_coeffs:
        try:
            coeffs.append(float(str_coeff))
        except ValueError:
            coeffs.append(0)
    return float_coeffs
