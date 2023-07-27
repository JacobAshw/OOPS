import numpy as np
import gurobipy as gp
from gurobipy import GRB
from OOPS_methods import *
import time
import sys

# The Subgraph Realization Problem Solver (SRPS)
# Given a pot, the SRPS aims to compute the smallest size graph the pot can produce.
# Further detail of this code can be found here: TODO: Link to archive when paper is uploaded
# ---------------------------------------------------------------------------------------------------------------------------------------------------
# Authors: Jacob Ashworth
# Other Contributors: Luca Grossmann, Fausto Navarro
# The work leading to this code was completed at ICERM, as part of Summer@ICERM 2023 (https://icerm.brown.edu/summerug/2023/)
# This codebase uses networkx to construct and build graphs (https://networkx.org/), v3.1
# Graphs are drawn using matplotlip's pyplot library (https://matplotlib.org/), v3.5.2
# The ILP solver used by this codebase is the Gurobi optimizer (https://www.gurobi.com/solutions/gurobi-optimizer/), v10.0
# ---------------------------------------------------------------------------------------------------------------------------------------------------

# ! ---------------------Put Target Pot Here-------------------------------
pot = {"aA", "aB", "bA"}
# ! -----------------------------------------------------------------------

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
# m.setParam(GRB.Param.LogToConsole, 0)

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

# Make tile usage dictionary
tileusage = {}
for index, tile in enumerate(pot):
    tileusage.update({tile : int(round(m.x[index]))})
# This output is in the same order as the construction matrix
# It will be an array of integers, specifying how many of the tile (represented by a column of M) are needed to build a smaller graph
print("Amounts of each tile type to build the graph: " + str(tileusage))