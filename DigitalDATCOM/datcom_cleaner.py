import numpy as np

#read in pre-cleaned file, where everything that can be cmd-f replaced is gone
file = open('datcom.txt')
s = file.read()

machs = [0.1,0.2,0.3,0.4,0.5,0.6,0.8,1.2,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.25,2.5]
alphas = [0.0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0]
alts =[26000.0,27000.0,28000.0,29000.0,30000.0,31000.0,32000.0,33000.0,34000.0,35000.0,36000.0,37000.0,38000.0,39000.0,40000.0,41000.0,42000.0,43000.0,44000.0]

c_ds = np.zeros((len(machs),len(alphas),len(alts)))

lastend = 0
start = s.find('ALPHA',lastend)
end = s.find('0*** VEHICLE')
table = np.array(s[start:end].split('\n'))
table = np.delete(table, len(table)-1)
table = np.delete(table, 1)
table = np.delete(table, 0)

for j in range(len(table)):
    l=table[j]
    coeffs = np.array(l.split())
    c_ds[0][j][0] = str(coeffs[1])
print c_ds

file.close()
