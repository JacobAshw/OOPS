import networkx as nx
import time
import sys
import tempfile
import math
import configparser
from matplotlib import pyplot as plt
import matplotlib.pyplot as plt

# * Below is an area sectioned off for building the target graph
# We use the networkx library (v3.1) for storing graphs (https://networkx.org/)
# If you are unfamiliar with graph construction using networkx, visit https://networkx.org/documentation/stable/tutorial.html
# Many graph generators already exist (https://networkx.org/documentation/stable/reference/generators.html)
# Consider using one of these generators if your are using a commonly studied type of graph
# ! ---------------------------------------------------------------------------------------------------

# Graph = nx
G = nx.cycle_graph(10)
# Graph = nx.lollipop_graph(4, 3)
# Graph = nx.cycle_graph(10)
# Graph = nx.dodecahedral_graph()
# Graph.add_edge(3, 1)
# Graph.add_node(4)
# Graph.add_node(5)
# Graph.add_edge(2,4)
# Graph.add_edge(4,1)
# Graph.add_edge(2,5)
# Graph = nx.ladder_graph(6)
#Graph = nx.cycle_graph(4)
# Graph = nx.octahedral_graph()
# Graph = nx.complete_graph(100)
# Graph = nx.cartesian_product(nx.petersen_graph(), nx.karate_club_graph())
# Graph = nx.truncated_tetrahedron_graph()

# ! ---------------------------------------------------------------------------------------------------

Graph = G

old_labels = {}
new_labels = {}
for index, node in enumerate(Graph.nodes()):
    new_labels.update({node: index})
    old_labels.update({index: node})
Graph = nx.relabel_nodes(Graph, new_labels)

Graphword = ""
for a in nx.generate_edgelist(Graph, data=['weight']):
    Graphword = Graphword + "\"" + a + "\"" + ","
print("Edgelist: " + Graphword[0:-1])

GraphwordG6 = nx.to_graph6_bytes(Graph)
print("G6: ", str(GraphwordG6)[12:-3])

pos = nx.kamada_kawai_layout(Graph)
#Draw the graph (with node labels)
nx.draw_networkx(Graph, pos, labels=old_labels, with_labels=True)
plt.show()