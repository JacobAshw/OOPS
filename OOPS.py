import networkx as nx
import time
import sys
import csv
import tempfile
import math
import configparser
from matplotlib import pyplot as plt
import matplotlib.pyplot as plt
from OOPS_files.methods import *
from OOPS_files.algorithms import *
from OOPS_files.gurobi import *

# ---------------------------------------------------------------------------------------------------------------------------------------------------
# The Optimal Oriented Pot Solver (OOPS)
version = "2.0"
# The goal of this code is to compute the optimal pot of tiles to construct a given graph
# Further detail of this code can be found here: TODO: Link to archive when paper is uploaded
# The current version only supports simple graphs, so self-loops and multi-edges may cause unsupported behavior
# ---------------------------------------------------------------------------------------------------------------------------------------------------
# Authors: Jacob Ashworth
# Other Contributors: Luca Grossmann, Fausto Navarro
# The work leading to this code was completed at ICERM, as part of Summer@ICERM 2023 (https://icerm.brown.edu/summerug/2023/)
# This codebase uses networkx to construct and build graphs (https://networkx.org/), v3.1
# Graphs are drawn using matplotlip's pyplot library (https://matplotlib.org/), v3.5.2
# The ILP solver used by this codebase is the Gurobi optimizer (https://www.gurobi.com/solutions/gurobi-optimizer/), v10.0
# ---------------------------------------------------------------------------------------------------------------------------------------------------

config_reader = configparser.ConfigParser()
config_path = 'config.ini'
if(len(sys.argv) > 1):
    config_path = sys.argv[1]

config_reader.read(config_path)
# print(config_reader.sections())

config = config_reader['Instance']

print(config['CLI_Enabled'])

CLI_Enabled = config.getboolean('CLI_Enabled')
display_graph = config.getboolean('display_graph') 
timer_enabled = config.getboolean('timer_enabled') 
gurobi_printiouts = config.getboolean('gurobi_printouts')
scenario_number = config.getint('scenario_number') 
bond_edges_verification = config.getboolean('bond_edges_verification')

# (1) qvalue
# (2) oneshot
# (3) partition
# (4) canonical
# (5) hybrid-qvalue
# (6) hybrid-canonical
method = ""
methodnum = config.getint('method')
if(methodnum==1):
    method = "qvalue"
if(methodnum==2):
    method = "oneshot"
if(methodnum==3):
    method = "partition"
if(methodnum==4):
    method = "canonical"
if(methodnum==5):
    method = "hybrid-quvalue"
if(methodnum==6):
    method = "hybrid-canonical"

graph = config['graph']
graphG6 = config.getboolean('graphG6')
outputfile = config['outputfile']

print(graph)

graph_name = graph

if(graphG6):
    # if(not 'graph6' in graph):
    #     graph = bytes(">>graphg6<<"+graph, 'ASCII')
    # else:
    graph = bytes(graph, 'ASCII')
    print(graph)
    Graph = nx.from_graph6_bytes(graph)
    # with tempfile.NamedTemporaryFile() as f:
    #     sstr = ">>graph6<<"+graph+"\n"
    #     _ = f.write(sstr.encode())
    #     _ = f.seek(0)
    #     Graph = nx.read_graph6(f)
else:
    Graph = nx.parse_edgelist(graph)

# TODO: Add Numerical Precision checks, output error if potential error is greater than 1
# TODO: Add ability to generate all optimal pots


print("OOPS solver version " + version)
# * Get our parameters from the command line interface, if enabled
if(CLI_Enabled):
    scenario_number, Graph = CLI_Setup(version)
    display_graph = True
    timer_enabled = True
    gurobi_printiouts = False
else:
    print("Using program parameters")

# * Relabel all of our nodes
# This is necessary due to non-integer node names. Certian generators (such as nx.grid_2D_graph)
# use different mathods of naming nodes, and our future code relies upon node names being integers
old_labels = {}
new_labels = {}
for index, node in enumerate(Graph.nodes()):
    new_labels.update({node: index})
    old_labels.update({index: node})
Graph = nx.relabel_nodes(Graph, new_labels)


# * Display the graph, if enabled
if(display_graph):
    #Get the layout to draw. If you want to display a differnet layout, check out the available layouts here: https://networkx.org/documentation/stable/reference/drawing.html
    pos = nx.kamada_kawai_layout(Graph)
    #Draw the graph (with node labels)
    nx.draw_networkx(Graph, pos, labels=old_labels, with_labels=True)
    plt.show()

# Do a final check to see if we can compute
valid_params = True
if(Graph.number_of_nodes() < 2):
    print("Graph has less than 2 nodes, aborting")
    valid_params = False
elif(not nx.is_connected(Graph)):
    print("Graph is disconnected, aborting")
    valid_params = False

if(valid_params):
    print("-----------------------------------------------------------------------")
    print("Beginning Computation...")
    if(timer_enabled): time_initial = time.perf_counter()
    print("")
    print("Gurobi:")

    data_headers = []
    data = []

    # If the chosen scenario is 1, we print out the information and close the program
    if(scenario_number == 1):
        # We now solve Scenario 1 for the graph, regardless of the scenario chosen by the user.
        # This is because the Scenario 1 ILP is quite fast, and it can be used to lower bound tile types in scenario 2.
        pot, tile_assignments, orientation = S1_optimal_pot(Graph, gurobi_printiouts)
        print("")

        print("Computation Complete")
        print("-----------------------------------------------------------------------")
        #Print out the final value for t1, and the pot
        print("T_1: " + str(len(pot)))
        print("pot: " + str(pot).replace("'",''))
        #Compute the tile usage, and print that out as well
        tile_usages = [len(tile_assignments.get(tile)) for tile in pot]
        print("tile usages: " + str(tile_usages))
        if(timer_enabled):
            time_final = time.perf_counter()
            print("Took " + str(time_final-time_initial) + " seconds")
            time_total = time_final - time_initial
        #If display_graph is enabled, draw the graph
        if(display_graph):
            # Build the digraph to display from the orientation matrix
            display_orientation(Graph, pot, tile_assignments, orientation, pos, old_labels)

        data_headers = ["Graph", "T_1", "T_1 pot", "T_1 ratios", "T_1 orientation", "T_1 time"]
        data = [graph_name, len(pot), pot, tile_assignments, orientation, time_total]


    if(scenario_number == 2):

        if(timer_enabled):
            time_initial = time.perf_counter()

        #Solve for tiles
        if(method == "oneshot"):
            t_pot, t_tile_assignments, t_orientations = S2_tiles_brute_force_partition(Graph, gurobi_printiouts)
        elif(method == "partition"):
            t_pot, t_tile_assignments, t_orientations = S2_tiles_partition_relaxation(Graph, gurobi_printiouts)
        elif(method == "qvalue"):
            t_pot, t_tile_assignments, t_orientations = S2_tiles_qvalue(Graph, gurobi_printiouts)
        elif(method == "hybrid-quvalue"):
            t_pot, t_tile_assignments, t_orientations = S2_tiles_hybrid_qvalue(Graph, gurobi_printiouts)
        
        #Timer code
        if(timer_enabled):
            time_final = time.perf_counter()
            print("Took " + str(time_final-time_initial) + " seconds")
            t_time = time_final - time_initial

        #Print results
        print("Tile type optimization:")
        optimal_tile_tiles, optimal_tile_bonds = print_pot(t_pot, t_tile_assignments)
        if(display_graph):
            display_orientation(Graph, t_pot, t_tile_assignments, t_orientations, pos, old_labels)

        #Bond code
        if(bond_edges_verification):
            print("Beginning bond edge optimization")

            #Timer code
            if(timer_enabled):
                time_initial = time.perf_counter()

            #Solve for bonds
            if(method == "oneshot"):
                b_pot, b_tile_assignments, b_orientations = S2_bonds_brute_force_partition(Graph, gurobi_printiouts, optimal_tile_tiles, optimal_tile_bonds)
            elif(method == "partition"):
                b_pot, b_tile_assignments, b_orientations = S2_bonds_partition_relaxation(Graph, gurobi_printiouts, optimal_tile_tiles, optimal_tile_bonds)
            elif(method == "qvalue"):
                b_pot, b_tile_assignments, b_orientations = S2_bonds_qvalue(Graph, gurobi_printiouts, optimal_tile_tiles, optimal_tile_bonds)
            elif(method == "hybrid-quvalue"):
                b_pot, b_tile_assignments, b_orientations = S2_bonds_hybrid_qvalue(Graph, gurobi_printiouts, optimal_tile_tiles, optimal_tile_bonds)

            #Timer code
            if(timer_enabled):
                time_final = time.perf_counter()
                print("Took " + str(time_final-time_initial) + " seconds")
                b_time = time_final - time_initial

            #Print results
            if(b_pot == False):
                print("Bond edge optimization is the same as tile optimization")
                t2_b2_match = True

                optimal_bond_bonds = optimal_tile_bonds
                optimal_bond_tiles = optimal_tile_tiles
            else:
                t2_b2_match = False
                #Print results
                print("Bond edge optimization")
                optimal_bond_tiles, optimal_bond_bonds = print_pot(b_pot, b_tile_assignments)

                #Counterexample detection
                if(optimal_bond_tiles != optimal_tile_tiles):
                    print("Optimal bond and tile mismatch")
                if(display_graph):
                    display_orientation(Graph, b_pot, b_tile_assignments, b_orientations, pos, old_labels)
        
        if(bond_edges_verification):
            data_headers = ["Graph","Methodnum","T_2 Value","T_2 Bonds","T_2 Pot","T_2 Pot Ratios","T_2 Orientation","B_2 Value","B_2 Tiles","B_2 Pot","B_2 Pot Ratios","B_2 Orientation","SamePot","TileTime","BondTime"]
            data = [graph_name, methodnum, optimal_tile_tiles, optimal_tile_bonds, t_pot, t_tile_assignments, t_orientations, optimal_bond_bonds, optimal_bond_tiles, b_pot, b_tile_assignments, b_orientations, t2_b2_match, t_time, b_time]
        else:
            data_headers = ["Graph","Methodnum","T_2 Value","T_2 Bonds","T_2 Pot","T_2 Pot Ratios","T_2 Orientation","TileTime"]
            data = [graph_name, methodnum, optimal_tile_tiles, optimal_tile_bonds, t_pot, t_tile_assignments, t_orientations, t_time]

    with open(outputfile, 'w', newline='\n') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(data_headers)
        writer.writerow(data)
