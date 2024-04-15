import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

folder_path = 'pot_output\graphing' #'C:\Users\ashworjs\OneDrive - Rose-Hulman Institute of Technology\Desktop\OOPS\pot_output\hq1\output'

all_files = os.listdir(folder_path)

# Filter out non-CSV files
csv_files = [f for f in all_files if f.endswith('.csv')]

print('a')
print(csv_files)

# Create a list to hold the dataframes
dfs = {}

for csv in csv_files:
    file_path = os.path.join(folder_path, csv)
    df = pd.read_csv(file_path)
    dfs.update({csv : df})

# print(df_list)

X = dfs.get('partition6c.csv').loc[:,"Graph"]

repls = [{'D~{': 'K_5'},
{'GhCGKC': 'C_8'},
{'F~CGG': 'L_4,3'},
{'E~`G': 'ConOpt'},
{'Gh`HGc': '4-Ladder'},
{'Gl_XIS': 'Cube'}
]
rep2 = {}
for s in repls:
    rep2.update(s)
repls = rep2
X2 = []
for x in X:
    if(x in repls.keys()):
        X2.append(repls.get(x))
    else:
        X2.append(x)
X = X2


hatchpatterns = ["/" , "-" , "O" , "x", "o", "+", ".", "*"]
  
X_axis = np.arange(len(X)) 

totalwidth=0.5
barwidth=0.5/len(dfs.keys())
num = 0

names = {}
names2 = [{'hq1': 'Hybrid'}, {'IPS': 'New'}, {'Qva': 'Original'}, {'par' : 'New'}, {'qva' : 'Original'}]
for name in names2:
    names.update(name)

totalnames = []
totaltimes = []

for csv in dfs.keys():
    name = names.get(csv[:3])
    totalnames.append(name)

    totaltime = 0
    for i in range(112):
        totaltime = totaltime + (dfs.get(csv).loc[:, "TileTime"])[i]
        totaltime = totaltime + (dfs.get(csv).loc[:, "BondTime"])[i]
    
    totaltimes.append(totaltime)

    # name = "asds"
    # plt.bar((X_axis - (totalwidth/2) + num*barwidth)[:112], (dfs.get(csv).loc[:, "TileTime"])[:112]+(dfs.get(csv).loc[:, "BondTime"])[:112], barwidth, label = name)
    # plt.bar(X_axis - (totalwidth/2) + num*barwidth, (dfs.get(csv).loc[:, "T_2 Value"])[:6], barwidth, label = csv, fill=False, hatch=hatchpatterns[num] )
    num = num + 1

for i in range(len(totalnames)):
    print(i)
    print(totalnames[i])
    print(totaltimes[i])
    plt.bar(totalnames[i], totaltimes[i])
  
# plt.bar(X_axis - 0.2, Ygirls, 0.4, label = 'hq1') 
# plt.bar(X_axis + 0.2, Zboys, 0.4, label = 'Boys') 
  
# plt.xticks(X_axis, X) 
plt.xlabel("Method") 
plt.ylabel("Total Time (seconds)") 
# plt.yscale('log')
plt.title("Total 6-vertex runtime") 
plt.legend() 
plt.show() 