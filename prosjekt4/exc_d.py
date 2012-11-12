import numpy as np
import matplotlib.pyplot as mpl

infile = np.loadtxt("probability.txt")

S = sum(infile)
n = [i*4-(2*400) for i in range(len(infile))]
print S
mpl.plot(n,infile)
mpl.xlabel('energy')
mpl.ylabel('probability')
mpl.title('Probability distribution for different energies for T = 2.4')
mpl.show()