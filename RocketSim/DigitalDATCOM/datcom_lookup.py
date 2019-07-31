import subprocess
import numpy as np
import os
def lookup(mach, alpha, alt, cg, mass):
    file = open('DigitalDATCOM/datcom_template.txt','r')
    s = file.read()
    out = s.replace('INSERT_MACH',str(mach)).replace('INSERT_ALT',str(alt)).replace('INSERT_ALPHA',str(alpha)).replace('INSERT_CG',str(cg)).replace('INSERT_WEIGHT',str(mass))
    outfile = open('DigitalDATCOM/current.dcm','w')
    outfile.write(out)
    outfile.close()
    file.close()
    subprocess.call('cd DigitalDATCOM; echo current.dcm | ./datcom; cd ..',shell=True,stdout=open(os.devnull, 'wb'))
    datcom_out = open('DigitalDATCOM/datcom.out','r')
    s=datcom_out.read()
    start = s.find('0\n')
    end = s.find('0*** VEHICLE')
    table = np.array(s[start:end].split('\n'))
    coeffs = np.array(table[-2].split())

    try:
        CD = float(coeffs[1])
    except:
        CD = 0
    try:
        CL = float(coeffs[2])
    except:
        CL = 0
    try:
        CM = float(coeffs[3])
    except:
        CM = 0
    try:
        CN = float(coeffs[4])
    except:
        CN = 0
    try:
        CA = float(coeffs[5])
    except:
        CA = 0
    try:
        XCP  = float(coeffs[6])
    except:
        XCP = 0
    try:
        CLA  = float(coeffs[7])
    except:
        CLA = 0
    try:
        CMA  = float(coeffs[8])
    except:
        CMA = 0
    try:
        CYB   = float(coeffs[9])
    except:
        CYB = 0
    try:
        CNB = float(coeffs[10])
    except:
        CNB = 0
    try:
        CLB = float(coeffs[11])
    except:
        CLB = 0
    datcom_out.close()
    return (CD,CL,CM,CN,CA,XCP,CLA,CMA,CYB,CNB,CLB)
