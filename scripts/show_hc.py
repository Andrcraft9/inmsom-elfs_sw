from hilbert_curve import d2xy
import numpy as np
import math as mth
import matplotlib.pyplot as plt

m = 3
print("Hilbert Curve index: %i" % m)
n = 2**m
print("Total points: %i" % n)

xplt = []
yplt = []
xmax, ymax = 0, 0
for d in range(n):
    x, y = d2xy(m, d)
    
    if x > xmax:
        xmax = x
    if y > ymax:
        ymax = y

    xplt.append(x)
    yplt.append(y)

print("xmax = %i and ymax = %i" % (xmax, ymax))
plt.plot(xplt, yplt, 'ks')
plt.plot(xplt, yplt, 'k')
    
plt.title("Hilbert curve, index: %i" % m)
plt.xticks(np.linspace(0, xmax, 4))
plt.yticks(np.linspace(0, ymax, 4))
plt.grid(True)
plt.show()