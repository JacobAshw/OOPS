import networkx as nx
import time
import sys
import math
import matplotlib.pyplot as plt
from OOPS_methods import *
from OOPS_gurobi import *

# ---------------------------------------------------------------------------------------------------------------------------------------------------
# The Optimal Oriented Pot Solver (OOPS)
version = "1.2"
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

# * User Controlled Paramaters
# Parameters to control the overall program operation and data reporting
CLI_Enabled = False # * If true, all following parameters are ignored, and command line prompts are used
display_graph = True # If true, the graph will be displayed before and after computation begins
timer_enabled = True # If true, there will be intermittent progress reports and total time will be displayed
gurobi_printiouts = False # If true, gurobi (the ILP solver) will print out its progress. Many ILPs are run over the program, so this will be a lot
scenario_number = 2 # Which scenario to run (Scenario 1, 2, or 3). Currently, only scenario 1 and 2 are supported
one_shot = True

# * Parameters to control the search space of the program.
# You can use these to fine-tune the search space if you already have information about the graph
# For a full sweep of the search space, set both minimums to 1 and both maximums to INFINITY
# All bounds are inclusive of edge values
# Maximums will automatically be reduced to graph limits
bond_edge_type_minimum = 1
bond_edge_type_maximum = math.inf
tile_types_minimum = 1
tile_types_maximum = math.inf

# * Below is an area sectioned off for building the target graph
# We use the networkx library (v3.1) for storing graphs (https://networkx.org/)
# If you are unfamiliar with graph construction using networkx, visit https://networkx.org/documentation/stable/tutorial.html
# Many graph generators already exist (https://networkx.org/documentation/stable/reference/generators.html)
# Consider using one of these generators if your are using a commonly studied type of graph

# TODO: Add Numerical Precision checks, output error if potential error is greater than 1
# TODO: Add ability to generate all optimal pots
# ! ---------------------------------------------------------------------------------------------------

# Graph = nx
Graph = nx.grid_2d_graph(2, 4)
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
    # We now solve Scenario 1 for the graph, regardless of the scenario chosen by the user.
    # This is because the Scenario 1 ILP is quite fast, and it can be used to lower bound tile types in scenario 2.
    pot, tile_assignments, orientation = optimal_pot_s1(Graph, gurobi_printiouts)
    print("")

    # If the chosen scenario is 1, we print out the information and close the program
    if(scenario_number == 1):
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
        #Test code for one-shot
        if(one_shot):
            if(True):
                if(timer_enabled):
                    time_initial = time.perf_counter()
                partitions = []
                while(True):
                    print("part: ", partitions)
                    max_tiles = 5
                    max_bonds = 4
                    pot, tile_assignments, orientations = optimal_pot_s2_one_shot(Graph, max_bonds, max_tiles, False, partitions)
                    print(pot)
                    print(tile_assignments)
                    minsize, ratios = free_variable_solve(pot, max_bonds)
                    print("ms: ", minsize)
                    print("rt: ", ratios)
                    if(minsize < Graph.number_of_nodes()):
                        while(len(ratios) < max_tiles):
                            ratios.append(0)
                        partitions.append(ratios)
                        continue
                    break
                if(timer_enabled):
                    time_final = time.perf_counter()
                    print("Took " + str(time_final-time_initial) + " seconds")
                if(display_graph):
                    # Build the digraph to display from the orientation matrix
                    display_orientation(Graph, pot, tile_assignments, orientations, pos, old_labels)
            else:
                if(timer_enabled):
                    time_initial = time.perf_counter()
                bond_edge_types = 1
                tile_types = 1
                solved = False
                while not solved:
                    print("Checking B_2=" + str(bond_edge_types) + ", T_2=" + str(tile_types))
                    pot, tile_assignments, orientations = optimal_pot_s2_partition(Graph, bond_edge_types, tile_types, False)
                    if not len(pot) >= 1:
                        bond_edge_types = bond_edge_types + 1
                        if(bond_edge_types >= tile_types or bond_edge_types > bond_edge_type_maximum):
                            bond_edge_types = bond_edge_type_minimum
                            tile_types = tile_types + 1
                    else:
                        break
                print("Pot: ", pot)
                print("Tile_assign: ", tile_assignments)
                if(timer_enabled):
                    time_final = time.perf_counter()
                    print("Took " + str(time_final-time_initial) + " seconds")
                if(display_graph):
                    # Build the digraph to display from the orientation matrix
                    display_orientation(Graph, pot, tile_assignments, orientation, pos, old_labels)

        else:
            #We check if this is a full sweep, which lets us check for T_2, B_2 same-pot conjecture counterexamples later
            full_sweep = True
            # *Start by doing our best to improve the bounds given to us
            #We know T_2 > T_1, so improve that bound if possible
            if(tile_types_minimum <= len(pot)):
                print("Set tile types minimum through Scenario 1 check")
                tile_types_minimum = len(pot)
            else:
                full_sweep = False
            #We know T_2 < Size(G)
            if(tile_types_maximum >= Graph.number_of_nodes()):
                print("Set tile types maximum to number of nodes")
                tile_types_maximum = Graph.number_of_nodes()
            else:
                full_sweep = False
            #Since B_2 + 1 < T_2, and T_2<=n, this means B_2<=n-1
            if(bond_edge_type_maximum >= Graph.number_of_nodes()-1):
                print("Set bond edge types maximum to number of nodes - 1")
                bond_edge_type_maximum = Graph.number_of_nodes()-1
            else:
                full_sweep = False
            print("-----------------------------------------------------------------------")
            # The final check to see if we're doing a full sweep
            if(bond_edge_type_minimum != 1):
                full_sweep = False

            #Set our initial values
            tile_types = tile_types_minimum
            bond_edge_types = bond_edge_type_minimum

            print("Checking bond edges: [" + str(bond_edge_type_minimum) + "-" + str(bond_edge_type_maximum) + "], tile types: [" + str(tile_types_minimum) + "-" + str(tile_types_maximum) + "]")

            #Define the function to loop over all qvals in a specific tile and bond edge type
            def loop_through_qvals(Graph: nx.graph, tile_types: int, bond_edge_types: int, gurobi_print: bool):
                #Set our initial qvalue limits
                minqs = [0 for i in range(tile_types)]
                q_under_equal = -1

                initial = True
                
                #keep checking pots
                while(True):
                    #get our pot
                    pot, tile_assignments, orientations, qvals, savedata = optimal_pot_s2(Graph, bond_edge_types, tile_types, minqs, q_under_equal, gurobi_print)
                    
                    #if we got nothing and no equality is enforced, we cannot generate a pot under these constraints
                    if(q_under_equal<=0 and len(pot)==0):
                        return False, [], [], []

                    #if we got nothing and equailty is enforced, we move the threshold and try again
                    if(len(pot)==0):
                        minqs[q_under_equal] = 0
                        q_under_equal = q_under_equal - 1
                        minqs[q_under_equal] = minqs[q_under_equal] + 1
                        continue

                    #run the pot through the free variable checker
                    minsize, ratios = free_variable_solve(pot, bond_edge_types)

                    #if it passes, return the answer
                    if(minsize == Graph.number_of_nodes()):
                        return True, pot, tile_assignments, orientations
                    
                    #if it generates something bigger, an error has happened
                    elif(minsize > Graph.number_of_nodes()):
                        print("Pot cannot generate a graph. minsize: "+str(minsize))
                        #We want the program to crash now, so return nothing
                        return None

                    #if it generates something smaller we set the lower bound to these qvalues, plus one
                    else:
                        minqs = qvals
                        q_under_equal = len(minqs) - 1
                        minqs[q_under_equal] = minqs[q_under_equal] + 1
            
            #Now slowly increase our bounds until we find an answer
            solved = False
            while not solved:
                print("Checking B_2=" + str(bond_edge_types) + ", T_2=" + str(tile_types))
                solved, pot, tile_assignments, orientations = loop_through_qvals(Graph, tile_types, bond_edge_types, gurobi_printiouts)
                if not solved:
                    bond_edge_types = bond_edge_types + 1
                    if(bond_edge_types >= tile_types or bond_edge_types > bond_edge_type_maximum):
                        bond_edge_types = bond_edge_type_minimum
                        tile_types = tile_types + 1
                    if(tile_types > tile_types_maximum):
                        break
            print("-----------------------------------------------------------------------")
            if not solved:
                print("Solution not found within bounds. Consider increasing search space")
            else:
                print("Computation Complete")
                if(full_sweep):
                    print("Optimality Guaranteed up to T_2, B_2 same-pot conjecture")
                    print("T_2=" + str(tile_types))
                    print("B_2<=" + str(bond_edge_types))
                else:
                    print("Optimality Guaranteed within range, up to T_2, B_2 same-pot conjecture")
                    print("T_2<=" + str(tile_types))
                    print("B_2<=" + str(bond_edge_types))
                print("pot: " + str(pot).replace("'", ""))
                tile_usages = [len(tile_assignments.get(tile)) for tile in pot]
                print("tile usages: " + str(tile_usages))
                if(timer_enabled):
                    time_final = time.perf_counter()
                    print("Took " + str(time_final-time_initial) + " seconds")
                #If display_graph is enabled, draw the graph
                if(display_graph):
                    # Build the digraph to display from the orientation matrix
                    display_orientation(Graph, pot, tile_assignments, orientations, pos, old_labels)
                if(full_sweep):
                    print("-----------------------------------------------------------------------")
                    print("Would you like to verify optimality for B_2 (Y/N)? This may take a while.")
                    while True:
                        choice = input()
                        if(choice.lower()[0] == 'n'):
                            break
                        elif(choice.lower()[0] == 'y'):
                            print("Verifying optimailty...")
                            if(timer_enabled):
                                time_initial = time.perf_counter()
                            for bond_edges in range(bond_edge_type_minimum, bond_edge_types):
                                for tiles in range(tile_types+1, tile_types_maximum + 1):
                                    print("Checking B_2=" + str(bond_edges) + ", T_2=" + str(tiles))
                                    solved, pot, tile_assignments, orientations = loop_through_qvals(Graph, tiles, bond_edges, gurobi_printiouts)
                                    if(solved == True):
                                        print("T_2, B_2 same-pot conjecture counterexample detected!")
                                        print("T_2<=" + str(tiles))
                                        print("B_2=" + str(bond_edges))
                                        print("pot: " + str(pot).replace("'", ""))
                                        tile_usages = [len(tile_assignments.get(tile)) for tile in pot]
                                        print("tile usages: " + str(tile_usages))
                                        if(display_graph):
                                            # Build the digraph to display from the orientation matrix
                                            display_orientation(Graph, pot, tile_assignments, orientations, pos, old_labels)

                                        sys.exit()
                            if(timer_enabled):
                                time_final = time.perf_counter()
                                print("Took " + str(time_final-time_initial) + " seconds")
                            print("-----------------------------------------------------------------------")
                            print("Optimality verified. B_2=" +str(bond_edge_types))
                            break
                        else:
                            print("Please enter a single character, 'Y' or 'N'")