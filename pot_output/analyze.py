from bz2 import compress
from re import L
import networkx as nx
import time
import csv
import tempfile
import sys
import math
import os
import matplotlib.pyplot as plt 
import pandas as pd

data = pd.read_csv("pot_output/7c.csv")

df = data.loc[data['SamePot']==False]
for g in df['Graph']:
    with tempfile.NamedTemporaryFile() as f:
        sstr = ">>graph6<<"+g+"\n"
        _ = f.write(sstr.encode())
        _ = f.seek(0)
        Graph = nx.read_graph6(f)
    print(g, ": ", Graph.edges())