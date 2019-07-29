import subprocess
import numpy as np

def lookup(mach, alpha, alt, cg, mass):
    file = open('DigitalDATCOM/datcom_template.txt','r')
    s = file.read()
    out = s.replace('INSERT_MACH',str(mach)).replace('INSERT_ALT',str(alt)).replace('INSERT_ALPHA',str(alpha)).replace('INSERT_CG',str(cg)).replace('INSERT_WEIGHT',str(mass))
    outfile = open('DigitalDATCOM/current.dcm','w')
    outfile.write(out)
    outfile.close()
    file.close()
    subprocess.call('cd DigitalDATCOM; echo current.dcm | ./datcom; cd ..',shell=True)
    datcom_out = open('DigitalDATCOM/datcom.out','r')
    s=datcom_out.read()
    start = s.find('0\n')
    end = s.find('0*** VEHICLE')
    table = np.array(s[start:end].split('\n'))
    coeffs = np.array(table[-2].split())
    try:
        CD = float(coeffs[1])
    except:
        CD = -1
    try:
        CL = float(coeffs[2])
    except:
        CL = -1
    try:
        CM = float(coeffs[3])
    except:
        CM = -1
    try:
        CN = float(coeffs[4])
    except:
        CN = -1
    try:
        CA = float(coeffs[5])
    except:
        CA = -1
    try:
        XCP  = float(coeffs[6])
    except:
        XCP = -1
    try:
        CLA  = float(coeffs[7])
    except:
        CLA = -1
    try:
        CMA  = float(coeffs[8])
    except:
        CMA = -1
    try:
        CYB   = float(coeffs[9])
    except:
        CYB = -1
    try:
        CNB = float(coeffs[10])
    except:
        CNB = -1
    try:
        CLB = float(coeffs[11])
    except:
        CLB = -1
    datcom_out.close()
    return (CD,CL,CM,CN,CA,XCP,CLA,CMA,CYB,CNB,CLB)
