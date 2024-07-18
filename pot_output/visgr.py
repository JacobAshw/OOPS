import chunk
import networkx as nx
import time
import csv
import tempfile
import sys
import math
import os
import matplotlib.pyplot as plt

# graphname = "G?`fJ{"
# with tempfile.NamedTemporaryFile() as f:
#     sstr = ">>graph6<<"+graphname+"\n"
#     _ = f.write(sstr.encode())
#     _ = f.seek(0)
#     Graph = nx.read_graph6(f)

Graph = nx.hypercube_graph(4)
nx.draw(Graph)
# plt.show()

print("Num Nodes")
print(Graph.number_of_nodes())
print("Num Edges")
print(Graph.number_of_edges())



#Checking graph properties
examples = ["EQyw", "EQjo", "EEjo", "ECzw"]
