import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

data= np.loadtxt ('a0.csv',delimiter=',')
x= data[:,0]
y= data[:,1]
u= data[:,2]
v= data[:,3]


x_unique = np.unique(x)
y_unique = np.unique(y)

# Reshape into 2D grid
X, Y = np.meshgrid(x_unique, y_unique)

# Need to match shape of U, V
U = np.zeros_like(X)
V = np.zeros_like(Y)

for i, xv in enumerate(x_unique):
    for j, yv in enumerate(y_unique):
        idx = np.where((x == xv) & (y == yv))
        U[j, i] = u[idx][0]
        V[j, i] = v[idx][0]


plt.streamplot(X, Y, U, V)
plt.axis('equal')
plt.show()
plt.xlim([0,1])
plt.ylim([0,1])
