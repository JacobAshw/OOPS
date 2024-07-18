import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os



X = ['5-cycle','10-cycle','Cube', 'Icosahedron', '2x3 Lattice'] 
Brute = [0.14,0.18,0.92,0.34,0.35] 
ES = [.97,.63,.45,.12,.83] 
  
X_axis = np.arange(len(X)) 
  
plt.bar(X_axis - 0.15, Brute, 0.3, label = 'Brute Force') 
plt.bar(X_axis + 0.15, ES, 0.3, label = 'Edge Swapping') 
  
plt.xticks(X_axis, X) 
plt.xlabel("Graph") 
plt.ylabel("Runtime(s)") 
plt.title("Runtime of verifying SRP for optimal known S3 pots") 
plt.yscale("log")
plt.legend() 

ax = plt.gca()
ax.set_ylim([0, 10000])

plt.show() 