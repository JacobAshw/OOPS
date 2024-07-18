from numbers import Real
import networkx as nx
import time
import sys
import csv
import tempfile
import math
import configparser
from matplotlib import pyplot as plt
import matplotlib.pyplot as plt
from OOPS_files.algorithms import *
from OOPS_files.methods import *
from OOPS_files.gurobi import *

#C_4:
# pot = ['aa','AB','bb']
#Dodeca
# pot = ['aCC','ABC','bbc','ccc']
#Cube S3
# pot = ['aaa','eee','bbA','ccB','ddB','ACE','CDE','ADE']

# Fun loop (add a c to make it lame)
# pot = ['aa','bc','CCA','BA']
time_initial = time.perf_counter()
# ! ---------------------Put Target Pot Here-------------------------------
pot = ['aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa','A']
# ! -----------------------------------------------------------------------
print("Pot: ", pot)
display = False
ratios = getAllRatios(pot)
print("All ratios that can realiza a pot: ", ratios)


# Setup and Graph Input
# edgelist = ['0 2 {}', '0 4 {}', '0 5 {}', '1 3 {}', '1 5 {}', '1 4 {}', '2 3 {}']
# Graph = nx.parse_edgelist(edgelist)
Graph, pot, tile_assignments, orientations = getMinSizeRealizableNoGurobi(pot, ratios[0])

old_labels = {}
new_labels = {}
for index, node in enumerate(Graph.nodes()):
    new_labels.update({node: index})
    old_labels.update({index: node})
Graph = nx.relabel_nodes(Graph, new_labels)
pos = nx.kamada_kawai_layout(Graph)

display_labels = {}
for i in new_labels:
    display_labels.update({int(i) : str(i)})

# nx.draw_circular(Graph)
# plt.show() 

#Get Cantidate Pot
# pot, tile_assignments, orientations = S2_tiles_qvalue(Graph, False)

# display_orientation(Graph, pot, tile_assignments, orientations, pos, old_labels)


# Realization class
class Realization:
    def __init__(self, pot, orientations, labels):
        # Store this for later
        self.pot = pot
        self.labels = labels
        
        # Build the realization from the orientation
        self.Graph = nx.MultiDiGraph()

        # Add the right number of vertices
        self.num_nodes = len(orientations[0])
        for vertex in range(self.num_nodes):
            self.Graph.add_node(vertex)

        # Add in our edges
        self.bond_edge_types = len(orientations)
        for bond_edge_type in range(self.bond_edge_types):
            orientation = orientations[bond_edge_type]

            #Traverse each orientation matrix
            for start_node in range(len(orientation)):
                for end_node in range(len(orientation)):
                    #Add the correct number of edges corresponding to the matrix entry
                    for instances in range(orientation[start_node][end_node]):
                        self.Graph.add_edge(start_node, end_node, BET=bond_edge_type, key=bond_edge_type)
    
    def get_matrix(self):
        orientations = []

        for bond_edge_type in range(self.bond_edge_types):
            orientation = [[0 for i in range(self.num_nodes)] for j in range(self.num_nodes)]

            #Check our edges, add all of right BET into orientation
            for start_node in range(self.num_nodes):
                for end_node in range(self.num_nodes):
                    if(self.Graph.get_edge_data(start_node, end_node) != None):
                        for key, value in self.Graph.get_edge_data(start_node, end_node).items():
                            if(value['BET'] == bond_edge_type):
                                orientation[start_node][end_node] = orientation[start_node][end_node] + 1
            
            orientations.append(orientation)

        return orientations

    def gen_moves(self):
        moves = []

        edges = [e for e in self.Graph.edges.data()]
            
        #Get all possible moves
        for edge_1_index in range(len(edges)):
            #Avoid duplicate swaps
            for edge_2_index in range(edge_1_index + 1, len(edges)):
                edge1 = edges[edge_1_index]
                edge2 = edges[edge_2_index]
                if(edge1[2]['BET']==edge2[2]['BET']):
                    moves.append((edge1, edge2))
        
        return moves
    
    def copy(self):
        return Realization(self.pot, self.get_matrix(), self.labels)
    
    def realization_isomorphic(self, other):
        #NOT DONE
        return (self.get_matrix() == other.get_matrix())

    def graph_isomorphic(self, other):
        return (nx.is_isomorphic(self.Graph, other.Graph))

    def edge_swap(self, move):
        e1 = move[0]
        e2 = move[1]

        bond_edge_type = e1[2]['BET']
        if(e2[2]['BET'] != bond_edge_type):
            print('Edge swap not possible (mismatching BETs). Aborting.')
            return

        self.Graph.remove_edge(e1[0], e1[1], key=bond_edge_type)
        self.Graph.remove_edge(e2[0], e2[1], key=bond_edge_type)

        self.Graph.add_edge(e1[0], e2[1], BET=bond_edge_type, key=bond_edge_type)
        self.Graph.add_edge(e2[0], e1[1], BET=bond_edge_type, key=bond_edge_type)

    def display(self):
        display_orientation(self.Graph, pot, tile_assignments, self.get_matrix(), pos, self.labels)

R0 = Realization(pot, orientations, display_labels)
print("Displaying Initial")
print([e for e in R0.Graph.edges])
if(display): R0.display()

print("Beginning state spanning")

#List of realizations we have checked
checked_realizations = [R0]
#List of tuples (realization, move) we still need to check
moves = []

#Initial population of moves
for move in R0.gen_moves():
    moves.append((R0, move))

NRPStatus = True

#If there were several ratios, check them here
if(len(ratios) > 1):
    for i in range(2, len(ratios)):
        Graph, pot, tile_assignments, orientations = getMinSizeRealizableNoGurobi(pot, ratios[i])

        old_labels = {}
        new_labels = {}
        for index, node in enumerate(Graph.nodes()):
            new_labels.update({node: index})
            old_labels.update({index: node})
        Graph = nx.relabel_nodes(Graph, new_labels)
        pos = nx.kamada_kawai_layout(Graph)

        display_labels = {}
        for i in new_labels:
            display_labels.update({int(i) : str(i)})

        RX = Realization(pot, orientations, display_labels)

        if(not RX.graph_isomorphic(R0)):
            print("NRP Failed (on other ratio). Displaying nonisomorphic graph.")
            print([e for e in RX.Graph.edges])
            NRPStatus = False
            if(display): RX.display()
            break

         #If R is graph isomorphic to any checked realization, continue loop
        unique = True
        for other in checked_realizations:
            if(RX.realization_isomorphic(other)):
                unique=False
                break
        if(not unique):
            continue

        #If R is unique, add it to the checked realizations, and add its moves to the list
        checked_realizations.append(RX)
        for move in RX.gen_moves():
            moves.append((RX, move))

#The big loop
iters = 0
while(len(moves) != 0 and NRPStatus == True):
    iters = iters + 1
    if(iters % 100 == 1):
        print("Checked ", iters, " moves")
    #Grab a move
    move = moves.pop()

    #Execute it on a copy of R
    R = move[0].copy()
    R.edge_swap(move[1])

    #If R is not graph isomorphic to R0, we have failed the NRP. Stop.
    if(not R.graph_isomorphic(R0)):
        print("NRP Failed. Displaying nonisomorphic graph.")
        print([e for e in R.Graph.edges])
        NRPStatus = False
        if(display): R.display()
        break

    #If R is graph isomorphic to any checked realization, continue loop
    unique = True
    for other in checked_realizations:
        if(R.realization_isomorphic(other)):
            unique=False
            break
    if(not unique):
        continue

    #If R is unique, add it to the checked realizations, and add its moves to the list
    checked_realizations.append(R)
    for move in R.gen_moves():
        moves.append((R, move))

#if NRPStatus is still true, we have succeeded
if(NRPStatus == True):
    print("Succeeded the NRP")
print("Checked a total of ", iters, " moves")

time_final = time.perf_counter()
print("Took " + str(time_final-time_initial) + " seconds")