import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os



X = ['Cube','Dodecahedron','2x4 Lattice', '2x5 Lattice', '3x5 Lattice'] 
Brute = [1.85,10000,8.31,2100,10000] 
ES = [.53,0.43,0.70,0.37,0.58] 
  
X_axis = np.arange(len(X)) 
  
plt.bar(X_axis - 0.15, Brute, 0.3, label = 'Brute Force') 
plt.bar(X_axis + 0.15, ES, 0.3, label = 'Edge Swapping') 
  
plt.xticks(X_axis, X) 
plt.xlabel("Graph") 
plt.ylabel("Runtime(s)") 
plt.title("Runtime of verifying SRP for optimal known S2 pots") 
plt.yscale("log")
plt.legend() 
plt.show() 