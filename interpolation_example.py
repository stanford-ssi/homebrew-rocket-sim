from scipy.interpolate import griddata
import numpy as np

# this stands in for the call to DATCOM with some parameters a, b, c
def dummy_f(a, b, c):
    return a ** 2 + b * c

coeff_data = []

for _ in range(10000):
    a = np.random.normal()
    b = np.random.normal()
    c = np.random.normal()
    d = dummy_f(a, b, c)
    coeff_data.append( [a, b, c] + [d] )

coeff_data = np.array(coeff_data)


test_a = 0.2
test_b = 0.4
test_c = 0.6

real_d = dummy_f(0.2, 0.4, 0.6)
grid_d = griddata(coeff_data[:, :3], coeff_data[:, 3], np.array([test_a, test_b, test_c]))

print('Real:', real_d)
print('Interpolated:', grid_d)
