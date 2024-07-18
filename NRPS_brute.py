from re import L
import numpy as np
import gurobipy as gp
from gurobipy import GRB
from OOPS_files.methods import *
import time
import sys

time_initial = time.perf_counter()
# ! ---------------------Put Target Pot Here-------------------------------
pot = ['aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa','A']
# ! -----------------------------------------------------------------------
buildGraphs = False


print("Pot: " + str(pot))

possible_half_edges, possible_half_edges_hat = get_half_edge_labels()

matrix = []
for tile in pot:
    sums = [0 for i in range(len(possible_half_edges))]
    for item in tile:
        if item in possible_half_edges:
            sums[possible_half_edges.index(item)] = sums[possible_half_edges.index(item)] + 1
        elif item in possible_half_edges_hat:
            sums[possible_half_edges_hat.index(item)] = sums[possible_half_edges_hat.index(item)] - 1
    matrix.append(sums)
matrix = np.transpose(matrix)

M = matrix

bond_edge_types = 0
for row in M:
    nonzero = False
    for entry in row:
        if(entry != 0):
            bond_edge_types = bond_edge_types + 1
            nonzero = True
            break
    if(not nonzero):
        break
print("Bond edge types: " + str(bond_edge_types))

# Some checks to make sure the matrix is computable
for index in range(len(M)):
    # If the matrix is not rectangular, that means it was entered correctly
    if(len(M[index]) != len(M[0])):
        print("Matrix is not rectangular, row 0 is size " + str(len(M[0])) + " and row " + str(index) + " is length " + str(len(M[index])))
        quit()
    # If a row has no negative values there's no 'hatted' bond edges of that type, indicating an input error
    if(min(M[index]) >= 0):
        if(min(M[index] != 0)):
            print("Row " + str(index) + " has no negative values, so a graph can't be made")
            quit()
    # If a row has no positive values there's no 'non-hatted' bond edges of that type, indicating an input error
    if(max(M[index]) <= 0):
        if(max(M[index]) != 0):
            print("Row " + str(index) + " has no positive values, so a graph can't be made")
            quit()

# Make our model
m = gp.Model("free_variable_problem")

# Declare each variable with a lower bound of 0
variables = []
obj = gp.LinExpr()
for tile in range(len(M[0])):
    x = m.addVar(vtype=GRB.INTEGER, lb=0)
    obj = obj + x
    variables.append(x)

# Set our objective
m.setObjective(obj, GRB.MINIMIZE)

# Set the constraints (each sum of bond edge types has to be 0)
m.addConstr(obj >= 1)
for row in M:
    cstr = gp.LinExpr()
    for i, val in enumerate(row):
        cstr = cstr + val * variables[i]
    m.addConstr(cstr == 0)

# * To disable console output, uncomment this line
m.setParam(GRB.Param.LogToConsole, 0)

# * Find all solutions
m.setParam(GRB.Param.PoolSearchMode, 2)
m.setParam(GRB.Param.PoolSolutions, 1000000)
m.setParam(GRB.Param.PoolGap, 0.5)

# Optimize the model
m.optimize()

print("----------------------------------------------------------------------------")

# Check the status of the result and exit the program if an error happened
if(int(m.status) == 2):
    print("Optimal Solution Found")
elif(int(m.status) == 7):
    print("Operation Timed Out")
    quit()
elif(int(m.status) == 3):
    print("Problem is infeasible. Check pot input.")
    quit()
elif(int(m.status) == 4):
    print("Problem is either infeasible unbounded. Check pot input.")
    quit()
else:
    print("Other error.")
    print("Error code: " + str(m.status))
    quit()


# Print out the minimum size of a graph that can be formed, along with what tiles it uses
print("Minimum size of graph constructed from pot: " + str(m.ObjVal))

numSolutions = m.SolCount

print(str(numSolutions) + " solutions found")

tileusages = []
for sol in range(numSolutions):
    m.setParam(GRB.Param.SolutionNumber, sol)
    tileusage = {}
    for index, tile in enumerate(pot):
        tileusage.update({tile : int(round(m.Xn[index]))})
    tileusages.append(tileusage)

ratios = []
for tileusage in tileusages:
    print("Amounts of each tile type to build the graph: " + str(tileusage))
    ratio = []
    for tile in pot:
        ratio.append(tileusage.get(tile))
    ratios.append(ratio)


# Now that ratios are found, build graphs they can make ------------------------------------------------

Graphs = []
DisplayInfo = {}

num_verticies = sum(ratios[0])
edge_list = []
for i in range(num_verticies):
    for j in range(i+1, num_verticies):
        edge_list.append([i,j])


for ratio in ratios:
    def getVertexBetValue(vertex, bond_edge, hat):
        vertexrationum = 0
        vertexnum = vertex
        for value in ratio:
            vertexnum = vertexnum - ratio[vertexrationum]
            if(vertexnum < 0):
                vert = pot[vertexrationum]
                letter = ''
                if(hat == 0):
                    letter = possible_half_edges[bond_edge]
                else:
                    letter = possible_half_edges_hat[bond_edge]
                value = 0
                for bet in vert:
                    if letter == bet:
                        value = value + 1
                return value
            vertexrationum = vertexrationum + 1
        print("couldn't find vertex")
        return -999

    def other_node(node, edge):
        if(edge[0]==node):
            return edge[1]
        elif(edge[1]==node):
            return edge[0]

    tile_assignments = {}
    for tile in pot:
        tile_assignments.update({tile : []})
    currentvertexnum = 0
    for i in range(len(ratio)):
        cr = ratio[i]
        while(cr > 0):
            tile_assignments.get(pot[i]).append(currentvertexnum)
            currentvertexnum = currentvertexnum + 1
            cr = cr - 1
    print(tile_assignments)

    m = gp.Model("AllGraphsMade")

    vertex_bets_map = {}

    for vertex in range(num_verticies):
        for bond_edge in range(bond_edge_types):
            for hat in range(2):
                vertex_bets_map.update({(vertex, bond_edge, hat) : getVertexBetValue(vertex, bond_edge, hat)})

    # indexing: (vertex_from, vertex_to, bet, hat)
    edge_constraints_map = {}

    #Populate edge constraints map
    for edge in edge_list:
        for bond_edge in range(bond_edge_types):
            for hat in range(2):
                # normal orientation
                name = str(edge) + "b" + str(bond_edge) + ":" + str(hat)
                x = m.addVar(vtype=GRB.INTEGER, lb=0, name=name)
                edge_constraints_map.update({(edge[0], edge[1], bond_edge, hat) : x})
                # reverse orientation
                edge_swap = (edge[1], edge[0])
                name2 = str(edge_swap) + "b" + str(bond_edge) + ":" + str(hat)
                x2 = m.addVar(vtype=GRB.INTEGER, lb=0, name=name2)
                edge_constraints_map.update({(edge_swap[0], edge_swap[1], bond_edge, hat) : x2})

    # indexing: (vertex_from, vertex_to)
    edge_usage_map = {}

    for edge in edge_list:
        # normal orientation
        name = str(edge) + "bb" + str(bond_edge) + ":" + str(hat)
        x = m.addVar(vtype=GRB.BINARY, name=name)
        edge_usage_map.update({(edge[0], edge[1]) : x})
    
    # * Make our objective(s) (minimize number of edges)
    obj = gp.LinExpr()
    for k_key in list(edge_usage_map.keys()):
        k = edge_usage_map.get(k_key)
        obj = obj + k
    m.setObjective(obj, GRB.MINIMIZE) 

    # # Edges respect verticies
    for vertex in range(num_verticies):
        for bond_edge in range(bond_edge_types):
            for hat in range(2):
                cstr = gp.LinExpr()
                cstr = cstr + vertex_bets_map.get((vertex, bond_edge, hat))
                for edge in edge_list:
                    if vertex in edge:
                        other = other_node(vertex, edge)
                        cstr = cstr - edge_constraints_map.get((vertex, other, bond_edge, hat))
                m.addConstr(cstr == 0, name="edge_respect_verticies_"+str(vertex)+"b"+str(bond_edge)+":"+str(hat))

    # # Each (used) edge has a 2 bond edge exactly
    for edge in edge_list:
        cstr = gp.LinExpr()
        for bond_edge in range(bond_edge_types):
            for hat in range(2):
                cstr = cstr + edge_constraints_map.get((edge[0], edge[1], bond_edge, hat))
                cstr = cstr + edge_constraints_map.get((edge[1], edge[0], bond_edge, hat))
        m.addGenConstrIndicator(edge_usage_map.get((edge[0], edge[1])), True, cstr == 2, name="edge_"+str(edge)+"_2bet")
        m.addGenConstrIndicator(edge_usage_map.get((edge[0], edge[1])), False, cstr == 0, name="edge_"+str(edge)+"_2bet")

    # Edges can be built
    for edge in edge_list:
        for bond_edge in range(bond_edge_types):
            for hat in range(2):
                if(hat == 0):
                    hat2 = 1
                else:
                    hat2 = 0
                cstr = gp.LinExpr()
                cstr = cstr + edge_constraints_map.get((edge[0], edge[1], bond_edge, hat))
                cstr = cstr - edge_constraints_map.get((edge[1], edge[0], bond_edge, hat2))
                m.addConstr(cstr == 0, name="edge_"+str(edge)+"_bet_"+str(bond_edge)+"_built")

    m.setParam(GRB.Param.PoolSolutions, 1000000)
    m.setParam(GRB.Param.PoolGap, 0.5)
    m.setParam(GRB.Param.PoolSearchMode, 2)
    m.optimize()

    numSolutions = m.SolCount
    print(str(numSolutions) + " solutions found")
    tileusages = []
    for sol in range(numSolutions):
        m.setParam(GRB.Param.SolutionNumber, sol)

        orientations = []
        for bond_edge in range(bond_edge_types):
            orientation = []
            for vertex in range(num_verticies):
                this_row = []
                for othervertex in range(num_verticies):
                    #Since the edge list only has one copy of each edge, we must check both directions
                    if([vertex, othervertex] in edge_list or [othervertex, vertex] in edge_list):
                        this_row.append(int(0.1 + edge_constraints_map.get((vertex, othervertex, bond_edge, 0)).Xn))
                    else:
                        this_row.append(0)
                orientation.append(this_row)
            orientations.append(orientation)
        
        G = nx.MultiGraph()

        numnodes = sum(ratio)
        #Add each node
        for vertex in range(numnodes):
            G.add_node(vertex)

        #Add each (directed) edge
        edge_lables = {}
        half_edges, half_edges_hat = get_half_edge_labels()
        for index, orientation in enumerate(orientations):
            for vertex1 in range(numnodes):
                for vertex2 in range(numnodes):
                    if(orientation[vertex1][vertex2] == 1):
                        G.add_edge(vertex1, vertex2)

        Graphs.append(G)
        DisplayInfo.update({G : (pot, ratio, tile_assignments, orientations)})
        # for index, node 
        # display_orientation_noinfo(pot, ratio, tile_assignments, orientations)

print("Number of graphs realized: " + str(len(Graphs)))

nonisographs = []
for G in Graphs:
    iso = False
    for H in nonisographs:
        if nx.is_isomorphic(G, H):
            iso = True
            break
    if not iso:
        nonisographs.append(G)

G = nonisographs[0]
print('-')
print('Pot: ' + str(DisplayInfo.get(G)[0]))
print('-')
print('Tile_Assign: ' + str(DisplayInfo.get(G)[2]))
print('-')
print('Graph: ' + str([a for a in nx.generate_edgelist(G)]))

print("Number of distinct graphs realized: " + str(len(nonisographs)))
for G in nonisographs:
    info = DisplayInfo.get(G)
    display_orientation_noinfo(info[0], info[1], info[2], info[3])

time_final = time.perf_counter()
print("Took " + str(time_final-time_initial) + " seconds")