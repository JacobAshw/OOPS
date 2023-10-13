import networkx as nx
import time
import csv
import tempfile
import sys
import math
import os
import matplotlib.pyplot as plt
from OOPS_methods import *
from OOPS_gurobi import *
# from OOPS import loop_through_qvals


writedir = "pot_output"
graphdir = "noniso_graphs"
chunksize = 1
graphfilename = "graph7c.g6"
#---------------------------------------------------------------------------------
def loop_through_qvals(Graph: nx.graph, tile_types: int, bond_edge_types: int, gurobi_print: bool):
    #Set our initial qvalue limits
    minqs = [0 for i in range(tile_types)]
    q_under_equal = -1
    
    #keep checking pots
    while(True):
        #get our pot
        pot, tile_assignments, orientations, qvals, save_data = optimal_pot_s2(Graph, bond_edge_types, tile_types, minqs, q_under_equal, gurobi_print)
        
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
#--------------------------------------------------------------------------------------
graphfile = open(graphdir + "/" + graphfilename, 'r')
graphnames = [graphname[:len(graphname)-1] for graphname in graphfile.readlines()]
graphfile.close()
graphs = []
for graphname in graphnames:
    with tempfile.NamedTemporaryFile() as f:
        sstr = ">>graph6<<"+graphname+"\n"
        _ = f.write(sstr.encode())
        _ = f.seek(0)
        G = nx.read_graph6(f)
    graphs.append(G)

start_index = 0

if(not os.path.exists(writedir + "/" + graphfilename + "_c" + str(chunksize))):
    os.mkdir(writedir + "/" + graphfilename + "_c" + str(chunksize))

completedfiles = os.listdir(writedir + "/" + graphfilename + "_c" + str(chunksize))
for cf in completedfiles:
    cf2 = int(cf[cf.index("-")+1:cf.index(".")])
    if(start_index <= cf2):
        start_index = cf2 + 1

header = ["Graph","T_1 Value","T_1 Pot","T_1 Ratios","T_1 Time","T_1 Orientation","T_2 Value","T_2 Pot","T_2 Pot Ratios","T_2 Orientation","B_2 Value","B_2 Pot","B_2 Pot Ratios","B_2 Orientation","SamePot","TileTime","BondTime"]

while(start_index < len(graphs)):
    endindex = start_index + chunksize - 1
    if(endindex > len(graphs)-1):
        endindex = len(graphs)-1

    d_graphs = []#
    d_t1 = []#
    d_t1pot = []#
    d_t1ratios = []#
    d_t1time = []#
    d_t1orientation = []#
    d_t2vals = []#
    d_t2pots = []#
    d_t2potratios = []#
    d_t2orientation = []#
    d_b2vals = []#
    d_b2pot = []#
    d_b2potratios = []#
    d_b2orientation = []#
    d_samepots = []#
    d_tiletime = []#
    d_bondtime = []#

    for graphindex in range(start_index, endindex+1):
        #Save some starting data
        d_graphs.append(graphnames[graphindex])
        Graph = graphs[graphindex]
        # print(graphs)
        # print(list(Graph.edges()))

        old_labels = {}
        new_labels = {}
        for index, node in enumerate(Graph.nodes()):
            new_labels.update({node: index})
            old_labels.update({index: node})
        Graph = nx.relabel_nodes(Graph, new_labels)
        # print(Graph)
        print("Checking graph " + str(graphindex) + " in interval [" + str(start_index) + "-" + str(endindex) + "]" + ": " + str(Graph.edges()))

        #Run S1
        start_time = time.perf_counter()
        pot, tile_assignments, orientation = optimal_pot_s1(Graph, False)
        end_time = time.perf_counter()

        #Save S1 data
        tile_usages = [len(tile_assignments.get(tile)) for tile in pot]
        d_t1.append(len(pot)) 
        d_t1pot.append(str(pot).replace("'",""))
        d_t1time.append(end_time-start_time)
        d_t1ratios.append(tile_usages)
        d_t1orientation.append(orientation)

        #Set S2 bounds
        tile_types_minimum = len(pot)
        tile_types_maximum = Graph.number_of_nodes()
        bond_edge_type_minimum = 1
        bond_edge_type_maximum = Graph.number_of_nodes()-1

        tile_types = tile_types_minimum
        bond_edge_types = bond_edge_type_minimum

        #Start S2
        start_time = time.perf_counter()

        solved = False
        while not solved:
            solved, pot, tile_assignments, orientations = loop_through_qvals(Graph, tile_types, bond_edge_types, False)
            if not solved:
                bond_edge_types = bond_edge_types + 1
                if(bond_edge_types >= tile_types or bond_edge_types > bond_edge_type_maximum):
                    bond_edge_types = bond_edge_type_minimum
                    tile_types = tile_types + 1
                if(tile_types > tile_types_maximum):
                    break

        if not solved:
            print("Failed to solve: " + graphnames[graphindex])
            exit()
        
        end_time = time.perf_counter()

        #Save data
        d_t2vals.append(tile_types)
        d_t2pots.append(str(pot).replace("'",""))
        tile_usages = [len(tile_assignments.get(tile)) for tile in pot]
        d_t2potratios.append(tile_usages)
        d_t2orientation.append(orientations)
        d_tiletime.append(end_time-start_time)

        tilebonds = bond_edge_types

        #Run bond verification
        b2t2same = True
        start_time = time.perf_counter()
        for bond_edges in range(bond_edge_type_minimum, bond_edge_types):
            for tiles in range(tile_types+1, tile_types_maximum + 1):
                solved, pot, tile_assignments, orientations = loop_through_qvals(Graph, tiles, bond_edges, False)
                if(solved == True):
                    b2t2same = False
                    break
            if(solved == True):
                break
        end_time = time.perf_counter()
        #Save data
        d_bondtime.append(end_time-start_time)
        if(b2t2same):
            d_samepots.append("True")
            d_b2vals.append(tilebonds)
            d_b2pot.append("")
            d_b2orientation.append("")
            d_b2potratios.append("")
        else:
            d_samepots.append("False")
            d_b2vals.append(bond_edges)
            d_b2pot.append(str(pot).replace("'",""))
            d_b2orientation.append(orientations)
            tile_usages = [len(tile_assignments.get(tile)) for tile in pot]
            d_b2potratios.append(tile_usages)
    
    print("Recording interval [" + str(start_index) + "-" + str(endindex) + "]")
    csvname = writedir + "/" + graphfilename + "_c" + str(chunksize) + "/" + str(start_index) + "-" + str(endindex) + ".csv"
    with open(csvname, 'w', newline='\n') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(header)
        for i in range(len(d_graphs)):
            data = []
            data.append( str(d_graphs[i]) )
            data.append( str(d_t1[i]) )
            data.append( str(d_t1pot[i])) 
            data.append( str(d_t1ratios[i])) 
            data.append( str(d_t1time[i]))
            data.append( str(d_t1orientation[i])) 
            data.append( str(d_t2vals[i]) )
            data.append( str(d_t2pots[i]) )
            data.append( str(d_t2potratios[i])) 
            data.append( str(d_t2orientation[i])) 
            data.append( str(d_b2vals[i]) )
            data.append( str(d_b2pot[i]) )
            data.append( str(d_b2potratios[i])) 
            data.append( str(d_b2orientation[i])) 
            data.append( str(d_samepots[i]) )
            data.append( str(d_tiletime[i]) )
            data.append( str(d_bondtime[i]) )
            writer.writerow(data)

    start_index = start_index + chunksize