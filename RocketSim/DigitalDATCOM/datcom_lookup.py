import subprocess
import numpy as np
import os

TEMPLATE_NAME = 'datcom_template.txt'
INPUT_NAME = 'current.dcm'
LOG_NAME = 'datcom_log.txt'
OUTPUT_NAME = 'datcom.out'
EXEC_NAME = 'datcom'
COLUMNS = ['ALPHA', 'CD', 'CL', 'CM', 'CN', 'CA', 'XCP',
           'CLA', 'CMA', 'CYB', 'CNB', 'CLB']

def lookup(machs, alphas, alts, cg, mass):
    current_dir = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(current_dir, TEMPLATE_NAME), 'r') as f:
        datcom_input = f.read()

    replacements = {
        'INSERT_MACHS': ','.join([str(mach) for mach in machs]),
        'INSERT_NMACH': str(len(machs)),
        'INSERT_ALPHAS': ','.join([str(alpha) for alpha in alphas]),
        'INSERT_NALPHA': str(len(alphas)),
        'INSERT_ALTS': ','.join([str(alt) for alt in alts]),
        'INSERT_NALT': str(len(alts)),
        'INSERT_CG': str(cg),
        'INSERT_WEIGHT': str(mass)
    }

    for key, value in replacements.items():
        datcom_input = datcom_input.replace(key, value)

    with open(os.path.join(current_dir, INPUT_NAME), 'w') as f:
        f.write(datcom_input)

    command = 'cd {}; echo {} | ./{}; cd ..'.format(
        current_dir, INPUT_NAME, EXEC_NAME)
    with open(os.path.join(current_dir, LOG_NAME), 'w') as f:
        # TODO this is the line that breaks things a bit: command's ./ ends up as "." and Windows doesn't recognise that as run
        subprocess.call(command, shell=True, stdout=f)

    with open(os.path.join(current_dir, OUTPUT_NAME), 'r') as f:
        datcom_output = f.read()

    coeffs = {}

    while True:
        card_start = datcom_output.find('FLIGHT CONDITIONS')
        if card_start == -1:
            break
        conds_start = card_start + datcom_output[card_start:].find('\n0') + 3
        conds_end = conds_start + datcom_output[(conds_start + 1):].find('\n0')
        diffs_start = conds_end + datcom_output[(conds_end + 1):].find('\n0\n') + 4
        diffs_end = diffs_start + datcom_output[diffs_start:].find('\n0*** VEHICLE')

        conds_text = datcom_output[conds_start:conds_end]
        diffs_text = datcom_output[diffs_start:diffs_end]
        conds = [float(cond) for cond in conds_text.split()]
        mach, alt = conds[:2]
        for diff_text in diffs_text.split('\n'):
            entries = [float(diff) for diff in diff_text.split()]
            values = dict(zip(COLUMNS[1:], entries[1:]))
            alpha = entries[0]
            coeffs[(mach, alpha, alt)] = values
        datcom_output = datcom_output[diffs_end:]

    return coeffs


# print(lookup([0.1, 0.2], [0.1, 0.2], [100, 200], 1, 1))
