import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

data= np.loadtxt ('my-6um-box-zeta=0.csv',delimiter=',')
x= data[:,0]
y= data[:,1]
u= data[:,2]
v= data[:,3]

plt.quiver(x, y, u, v)
plt.axis('equal')

plt.show()

