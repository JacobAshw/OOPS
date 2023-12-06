import networkx as nx
import time
import sys
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
config_reader.read('config.ini')
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

inputfile = config['inputfile']
outputfile = config['outputfile']

# * Below is an area sectioned off for building the target graph
# We use the networkx library (v3.1) for storing graphs (https://networkx.org/)
# If you are unfamiliar with graph construction using networkx, visit https://networkx.org/documentation/stable/tutorial.html
# Many graph generators already exist (https://networkx.org/documentation/stable/reference/generators.html)
# Consider using one of these generators if your are using a commonly studied type of graph

# TODO: Add Numerical Precision checks, output error if potential error is greater than 1
# TODO: Add ability to generate all optimal pots
# ! ---------------------------------------------------------------------------------------------------

# Graph = nx
# Graph = nx.cycle_graph(4)
# Graph.add_node(5)
# Graph.add_edge(1, 5)
# Graph.add_edge(2, 5)
Graph = nx.lollipop_graph(4, 3)
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
        #If display_graph is enabled, draw the graph
        if(display_graph):
            # Build the digraph to display from the orientation matrix
            display_orientation(Graph, pot, tile_assignments, orientation, pos, old_labels)


    if(scenario_number == 2):
        if(method == "oneshot"):
            #Timer code
            if(timer_enabled):
                time_initial = time.perf_counter()

            #Solve for tiles
            pot, tile_assignments, orientations = S2_tiles_brute_force_partition(Graph, gurobi_printiouts)
            
            #Timer code
            if(timer_enabled):
                time_final = time.perf_counter()
                print("Took " + str(time_final-time_initial) + " seconds")

            #Print results
            print("Tile type optimization:")
            optimal_tile_tiles, optimal_tile_bonds = print_pot(pot, tile_assignments)
            if(display_graph):
                display_orientation(Graph, pot, tile_assignments, orientations, pos, old_labels)

            #Bond code
            if(bond_edges_verification):
                print("Beginning bond edge optimization")

                #Timer code
                if(timer_enabled):
                    time_initial = time.perf_counter()

                #Solve for bonds
                pot, tile_assignments, orientations = S2_bonds_brute_force_partition(Graph, gurobi_printiouts, optimal_tile_tiles, optimal_tile_bonds)
                
                #Timer code
                if(timer_enabled):
                    time_final = time.perf_counter()
                    print("Took " + str(time_final-time_initial) + " seconds")

                #Print results
                if(pot == False):
                    print("Bond edge optimization is the same as tile optimization")
                else:
                    #Print results
                    print("Bond edge optimization")
                    optimal_bond_tiles, optimal_bond_bonds = print_pot(pot, tile_assignments)

                    #Counterexample detection
                    if(optimal_bond_tiles != optimal_tile_tiles):
                        print("Optimal bond and tile mismatch")
                    if(display_graph):
                        display_orientation(Graph, pot, tile_assignments, orientations, pos, old_labels)

        if(method == "partition"):
            #Timer code
            if(timer_enabled):
                time_initial = time.perf_counter()

            #Solve for tiles
            pot, tile_assignments, orientations = S2_tiles_partition_relaxation(Graph, gurobi_printiouts)
            
            #Timer code
            if(timer_enabled):
                time_final = time.perf_counter()
                print("Took " + str(time_final-time_initial) + " seconds")

            #Print results
            print("Tile type optimization:")
            optimal_tile_tiles, optimal_tile_bonds = print_pot(pot, tile_assignments)
            if(display_graph):
                display_orientation(Graph, pot, tile_assignments, orientations, pos, old_labels)

            #Bond code
            if(bond_edges_verification):
                print("Beginning bond edge optimization")

                #Timer code
                if(timer_enabled):
                    time_initial = time.perf_counter()

                #Solve for bonds
                pot, tile_assignments, orientations = S2_bonds_partition_relaxation(Graph, gurobi_printiouts, optimal_tile_tiles, optimal_tile_bonds)
                
                #Timer code
                if(timer_enabled):
                    time_final = time.perf_counter()
                    print("Took " + str(time_final-time_initial) + " seconds")
                
                #Counterexample detection
                if(pot == False):
                    print("Bond edge optimization is the same as tile optimization")
                else:
                    #Print results
                    print("Bond edge optimization")
                    optimal_bond_tiles, optimal_bond_bonds = print_pot(pot, tile_assignments)

                    #Counterexample detection
                    if(optimal_bond_tiles != optimal_tile_tiles):
                        print("Optimal bond and tile mismatch")
                    if(display_graph):
                        display_orientation(Graph, pot, tile_assignments, orientations, pos, old_labels)

        if(method == "qvalue"):
            #Timer code
            if(timer_enabled):
                time_initial = time.perf_counter()

            #Solve for tiles
            pot, tile_assignments, orientations = S2_tiles_qvalue(Graph, gurobi_printiouts)
            
            #Timer code
            if(timer_enabled):
                time_final = time.perf_counter()
                print("Took " + str(time_final-time_initial) + " seconds")

            #Print results
            print("Tile type optimization:")
            optimal_tile_tiles, optimal_tile_bonds = print_pot(pot, tile_assignments)
            if(display_graph):
                display_orientation(Graph, pot, tile_assignments, orientations, pos, old_labels)

            #Bond code
            if(bond_edges_verification):
                print("Beginning bond edge optimization")

                #Timer code
                if(timer_enabled):
                    time_initial = time.perf_counter()

                #Solve for bonds
                pot, tile_assignments, orientations = S2_bonds_qvalue(Graph, gurobi_printiouts, optimal_tile_tiles, optimal_tile_bonds)
                
                #Timer code
                if(timer_enabled):
                    time_final = time.perf_counter()
                    print("Took " + str(time_final-time_initial) + " seconds")
                    print("Bond edges optimization:")

                if(pot == False):
                    print("Bond edge optimization is the same as tile optimization")
                else:
                    #Print results
                    print("Bond edge optimization")
                    optimal_bond_tiles, optimal_bond_bonds = print_pot(pot, tile_assignments)

                    #Counterexample detection
                    if(optimal_bond_tiles != optimal_tile_tiles):
                        print("Optimal bond and tile mismatch")
                    if(display_graph):
                        display_orientation(Graph, pot, tile_assignments, orientations, pos, old_labels)

        if(method == "hybrid-qvalue"):
            #Timer code
            if(timer_enabled):
                time_initial = time.perf_counter()

            #Solve for tiles
            pot, tile_assignments, orientations = S2_tiles_hybrid_qvalue(Graph, gurobi_printiouts)
            
            #Timer code
            if(timer_enabled):
                time_final = time.perf_counter()
                print("Took " + str(time_final-time_initial) + " seconds")

            #Print results
            print("Tile type optimization:")
            optimal_tile_tiles, optimal_tile_bonds = print_pot(pot, tile_assignments)
            if(display_graph):
                display_orientation(Graph, pot, tile_assignments, orientations, pos, old_labels)

            #Bond code
            if(bond_edges_verification):
                print("Beginning bond edge optimization")

                #Timer code
                if(timer_enabled):
                    time_initial = time.perf_counter()

                #Solve for bonds
                pot, tile_assignments, orientations = S2_bonds_hybrid_qvalue(Graph, gurobi_printiouts, optimal_tile_tiles, optimal_tile_bonds)
                
                #Timer code
                if(timer_enabled):
                    time_final = time.perf_counter()
                    print("Took " + str(time_final-time_initial) + " seconds")
                    print("Bond edges optimization:")

                if(pot == False):
                    print("Bond edge optimization is the same as tile optimization")
                else:
                    #Print results
                    print("Bond edge optimization")
                    optimal_bond_tiles, optimal_bond_bonds = print_pot(pot, tile_assignments)

                    #Counterexample detection
                    if(optimal_bond_tiles != optimal_tile_tiles):
                        print("Optimal bond and tile mismatch")
                    if(display_graph):
                        display_orientation(Graph, pot, tile_assignments, orientations, pos, old_labels)

