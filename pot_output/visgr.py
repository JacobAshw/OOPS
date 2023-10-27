import chunk
import networkx as nx
import time
import csv
import tempfile
import sys
import math
import os
import matplotlib.pyplot as plt

graphname = "F?q`o"
with tempfile.NamedTemporaryFile() as f:
    sstr = ">>graph6<<"+graphname+"\n"
    _ = f.write(sstr.encode())
    _ = f.seek(0)
    Graph = nx.read_graph6(f)

nx.draw(Graph)
plt.show()



#Checking graph properties
examples = ["EQyw", "EQjo", "EEjo", "ECzw"]
