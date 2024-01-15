import networkx as nx
import gurobipy as gp
import numpy as np
from gurobipy import GRB
from OOPS_files.methods import *
import itertools

# This file holds the ILP solvers for OOPS.py
# All of our ILPs are solved using the Gurobi optimizer (https://www.gurobi.com/solutions/gurobi-optimizer/)

# This is the ILP solver for scenario 1
def optimal_pot_s1(G: nx.Graph, print_log: bool) -> list[str]:
    # Extract some graph information we need
    # The degree sequence of the graph
    degree_sequence = [d for n,d in G.degree()]
    degree_sequence.sort()
    # A list of each degree of the graph
    different_degrees = list(set(degree_sequence))
    # The number of verticies in the graph
    num_verticies = G.number_of_nodes()
    num_tiles = 2 * len(different_degrees)
    # The list of edges in the graph
    edge_list = list(G.edges())

    # This function returns the other node in an edge of the graph
    # Though it is a simple function, it is used quite often and defining it here saves space
    def other_node(node: int, edge: tuple[int]) -> int:
        if(edge[0]==node):
            return edge[1]
        elif(edge[1]==node):
            return edge[0]
    
    try:
        # Create a new model
        m = gp.Model("optimal_pot_s1")

        # Create dictionaries to store our decision variables in
        # Storing these variables in a dictionary is perferable, since in lists
        # certian sets of variables (primarily the edge_constraints) will have large
        # amounts of empty values, and it makes retriving answers more difficult.

        # Many of these dictionaries have multiple indexes, so they are standardized
        # as maps from tuples to decision variables
        # 0 always represents an unhatted edge, and 1 always represents a hatted edge

        # * Create decision variables and maps---------------------------------------------------------
        # vertex_bets_map holds the v_{i,h} decision variables, representing the
        # amount of each hat/non-hat in each vertex. The indexing tuple is (vertex, hat)
        vertex_bets_map = {}

        # populate vertex_bets_map
        for vertex in range(num_verticies):
            for hat in range(2):
                name = "node" + str(vertex) + ":" + str(hat) 
                x = m.addVar(vtype=GRB.INTEGER, lb=0, name=name)
                vertex_bets_map.update({(vertex, hat) : x})

        # tile_bets_map holds the t_{i,h} decision variables, representing the
        # amount of each hat/non-hat in each tile. The indexing tuple is (tile, hat)
        tile_bets_map = {}

        # populate tile_bets_map
        for tile in range(num_tiles):
            for hat in range(2):
                name = "tile" + str(tile) + ":" + str(hat) 
                x = m.addVar(vtype=GRB.INTEGER, lb=0, name=name)
                tile_bets_map.update({(tile, hat) : x})

        # vertex_tiles_decision_map holds the x_{i,j} decision variables, representing if
        # vertex i is tile type j. The indexing tuple is (vertex, tile)
        vertex_tiles_decision_map = {}

        # populate vertex_tiles_decision_map
        for vertex in range(num_verticies):
            for tile in range(num_tiles):
                name = "node" +  str(vertex) + "istile" + str(tile)
                x = m.addVar(vtype=GRB.BINARY, name=name)
                vertex_tiles_decision_map.update({(vertex, tile) : x})
        
        # k_map holds the k_{i} decision variables, representing if tile i is used at all.
        # The indexing tuple is: (tile) (for constiency, the single value is still a tuple)
        k_map = {}

        # populate k_map
        for tile in range(num_tiles):
            name = "k_tile" + str(tile)
            x = m.addVar(vtype=GRB.BINARY, name=name)
            k_map.update({(tile) : x})

        # edge_constraints_map holds the e_{i,j,b} decision variables, representing if the half
        # edge going from vertex i to vertex j is hat/non-hat type h indexing: (vertex_from, vertex_to, hat)
        edge_constraints_map = {}

        # populate edge_constraints_map
        for edge in edge_list:
            for hat in range(2):
                # normal orientation
                name = str(edge) + ":" + str(hat)
                x = m.addVar(vtype=GRB.INTEGER, lb=0, name=name)
                edge_constraints_map.update({(edge[0], edge[1], hat) : x})
                # reverse orientation
                # the reverse orientations are necessary because (0,2) could be a, and (2, 0) will then be A
                edge_swap = (edge[1], edge[0])
                name2 = str(edge_swap) + ":" + str(hat)
                x2 = m.addVar(vtype=GRB.INTEGER, lb=0, name=name2)
                edge_constraints_map.update({(edge_swap[0], edge_swap[1], hat) : x2})

        # * Make our objective---------------------------------------------------------------------------------------
        # Minimize over the sum of the k-values
        obj = gp.LinExpr()
        for k_key in list(k_map.keys()):
            k = k_map.get(k_key)
            obj = obj + k
        m.setObjective(obj, GRB.MINIMIZE) 

        # * Make our constraints----------------------------------------------------------------------------------
        # Add constraint: same number of as and As in the entire graph
        cstr = gp.LinExpr()
        for vertex in range(num_verticies):
            bet_no_hat = vertex_bets_map.get((vertex, 0))
            bet_hat = vertex_bets_map.get((vertex, 1))
            cstr = cstr + bet_no_hat - bet_hat
        m.addConstr(cstr == 0, "bond_edge_aA_match")

        # Add constraint: tiles are assigned to verticies of the correct degree
        for vertex in range(num_verticies):
            for tile in range(num_tiles):
                cstr = gp.LinExpr()
                for hat in range(2):
                    cstr = cstr + vertex_bets_map.get((vertex, hat))
                    cstr = cstr - tile_bets_map.get((tile, hat))
                m.addGenConstrIndicator(vertex_tiles_decision_map.get((vertex, tile)), True, cstr == 0, name="v"+str(vertex)+"t"+str(tile)+"chosen_degree_match")

        # Add constraint: verticies have bond edges equal to degree
        for vertex in range(num_verticies):
            cstr = gp.LinExpr()
            for hat in range(2):
                cstr = cstr + vertex_bets_map.get((vertex, hat))
            m.addConstr(cstr == G.degree[vertex], "vertex_"+str(vertex)+"_bets_match_degree")
        
        # Add constraint: each vertex must be exactly one tile type
        for vertex in range(num_verticies):
            cstr = gp.LinExpr()
            for tile in range(num_tiles):
                cstr = cstr + vertex_tiles_decision_map.get((vertex, tile))
            m.addConstr(cstr == 1, "v"+str(vertex)+"_one_tile_type")
        
        # Add constraint: Enforce tile selections
        for vertex in range(num_verticies):
            for tile in range(num_tiles):
                for hat in range(2):
                    cstr = gp.LinExpr()
                    cstr = cstr + vertex_bets_map.get((vertex, hat))
                    cstr = cstr - tile_bets_map.get((tile, hat))
                    for othertile in range(num_tiles):
                        if(othertile!=tile):
                            cstr = cstr - max(different_degrees) * vertex_tiles_decision_map.get((vertex, othertile))
                    m.addConstr(cstr <= 0, "enforce_tile_selection_v"+str(vertex)+"t"+str(tile)+":"+str(hat))

        # Add constraint: Ks count tiles used
        for tile in range(num_tiles):
            cstr = gp.LinExpr()
            for hat in range(2):
                cstr = cstr + tile_bets_map.get((tile, hat))
            cstr = cstr - max(different_degrees) * k_map.get((tile))
            m.addConstr(cstr <= 0, name="k_count_tile"+","+str(tile))

        # # Edges respect verticies
        for vertex in range(num_verticies):
            for hat in range(2):
                cstr = gp.LinExpr()
                cstr = cstr + vertex_bets_map.get((vertex, hat))
                for edge in edge_list:
                    if vertex in edge:
                        other = other_node(vertex, edge)
                        cstr = cstr - edge_constraints_map.get((vertex, other, hat))
                m.addConstr(cstr == 0, name="edge_respect_verticies_"+str(vertex)+":"+str(hat))

        # # Each edge has a 2 bond edge exactly
        for edge in edge_list:
            cstr = gp.LinExpr()
            for hat in range(2):
                cstr = cstr + edge_constraints_map.get((edge[0], edge[1], hat))
                cstr = cstr + edge_constraints_map.get((edge[1], edge[0], hat))
            m.addConstr(cstr == 2, name="edge_"+str(edge)+"_2bet")

        # Edges can be built
        for edge in edge_list:
            for hat in range(2):
                if(hat == 0):
                    hat2 = 1
                else:
                    hat2 = 0
                cstr = gp.LinExpr()
                cstr = cstr + edge_constraints_map.get((edge[0], edge[1], hat))
                cstr = cstr - edge_constraints_map.get((edge[1], edge[0], hat2))
                m.addConstr(cstr == 0, name="edge_"+str(edge)+"_built")

        # * Run out model ---------------------------------------------------------------------------
        if(not print_log):
            m.setParam(GRB.Param.LogToConsole, 0)

        # Optimize model
        m.optimize()

        # If our optimization is not optimal (status 2), print an error and return nothing
        if(m.status != 2):
            print("ERROR: optimal solution not found")
            print("Status of optimization: ", m.status)
            return None

        # * Construct our pot and orientation, and return it----------------------------------------------------------
        # Get our dictionary of tile assignments and build our pot
        pot = []
        tile_assignments = {}
        for tile_num in range(num_tiles):
            if(int(0.1 + k_map.get((tile_num)).X) == 1):
                tile = 'a' * int(0.1 + tile_bets_map.get((tile_num, 0)).X)
                tile = tile + 'A' * int(0.1 + tile_bets_map.get((tile_num, 1)).X)
                pot.append(tile)
                tile_assignments.update({tile : []})
                #add each vertex assigned to the tile
                for vertex in range(num_verticies):
                    if(int(0.1 + vertex_tiles_decision_map.get((vertex, tile_num)).X) == 1):
                        tile_assignments.get(tile).append(vertex)

        # Build the edges of the orientation of the graph
        orientations = []
        for bond_edge in range(1):
            orientation = []
            for vertex in range(num_verticies):
                this_row = []
                for othervertex in range(num_verticies):
                    #Since the edge list only has one copy of each edge, we must check both directions
                    if((vertex, othervertex) in edge_list or (othervertex, vertex) in edge_list):
                        this_row.append(int(0.1 + edge_constraints_map.get((vertex, othervertex, 0)).X))
                    else:
                        this_row.append(0)
                orientation.append(this_row)
            orientations.append(orientation)

        # Primarily for debugging: comment out to print all decision variable values
        # for v in m.getVars():
        #     print('%s %g' % (v.VarName, v.X))

        return pot, tile_assignments, orientations

    except gp.GurobiError as e:
        print('Error code ' + str(e.errno) + ': ' + str(e))

        u_map = {}

        for tile in range(num_tiles):
            name = "u_tile" + str(tile)
            x = m.addVar(vtype=GRB.BINARY, name=name)
            u_map.update({(tile):x})



def optimal_pot_s2_relaxation_partition(G, bond_edge_types, max_tiles, print_log, partitions):
    # Extract some graph info
    degree_sequence = [d for n,d in G.degree()]
    degree_sequence.sort()
    different_degrees = list(set(degree_sequence))
    num_verticies = G.number_of_nodes()
    edge_list = list(G.edges())

    def deg(node):
        deg = 0
        for edge in edge_list:
            if(node in edge):
                deg = deg + 1
        return deg

    def other_node(node, edge):
        if(edge[0]==node):
            return edge[1]
        elif(edge[1]==node):
            return edge[0]
    try:
        # Create a new model
        m = gp.Model("optimal_pot_s2")

        # Create the dictionaries of variables
        # indexing: (vertex, bet, hat)
        vertex_bets_map = {}
        # indexing: (tile, bet, hat)
        tile_bets_map = {}
        # indexing: (vertex, tile)
        vertex_tiles_decision_map = {}
        # indexing: (tile)
        k_map = {}
        # indexing: (tile)
        # u_map = {}
        # indexing: (vertex_from, vertex_to, bet, hat)
        edge_constraints_map = {}

        # * Create variables
        # make a variable for each bond edge of each vertex
        for vertex in range(num_verticies):
            for bond_edge in range(bond_edge_types):
                for hat in range(2):
                    name = "node" + str(vertex) + "bet" + str(bond_edge) + ":" + str(hat) 
                    x = m.addVar(vtype=GRB.INTEGER, lb=0, name=name)
                    vertex_bets_map.update({(vertex, bond_edge, hat) : x})

        # make a variable for each bond edge of each tile
        for tile in range(max_tiles):
            for bond_edge in range(bond_edge_types):
                for hat in range(2):
                    # name = "tile" + str(degree) + ":" + str(tile) + "bet" + str(bond_edge) + ":" + str(hat) 
                    name = "!" + "," + str(tile) + "," + str(bond_edge) + "," + str(hat)
                    x = m.addVar(vtype=GRB.INTEGER, lb=0, name=name)
                    tile_bets_map.update({(tile, bond_edge, hat) : x})

        # For each pair (vertex and tile), create a decision var
        for vertex in range(num_verticies):
            for tile in range(max_tiles):
                name = "node" +  str(vertex) + "istile" + str(tile)
                x = m.addVar(lb=0, ub=1, vtype=GRB.CONTINUOUS, name=name)
                vertex_tiles_decision_map.update({(vertex, tile) : x})

        # For each tile type, create a corresponding k var
        for tile in range(max_tiles):
            name = "k_tile" + str(tile)
            x = m.addVar(vtype=GRB.BINARY, name=name)
            k_map.update({(tile) : x})

        # For each side of each edge, create vars of each bet
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

        # * Make our objective(s) 
        obj = gp.LinExpr()
        for k_key in list(k_map.keys()):
            k = k_map.get(k_key)
            obj = obj + k
        m.setObjective(obj, GRB.MINIMIZE) 

        # * Make our constraints

        # Add constraint: same number of each bond edge type in all verticies (a graph can be made)
        for bond_edge in range(bond_edge_types):
            cstr = gp.LinExpr()
            for vertex in range(num_verticies):
                bet_no_hat = vertex_bets_map.get((vertex, bond_edge, 0))
                bet_hat = vertex_bets_map.get((vertex, bond_edge, 1))
                cstr = cstr + bet_no_hat - bet_hat
            m.addConstr(cstr == 0, "bond_edge_aA_match_"+str(bond_edge))

        # Add constraint: tiles are assigned to verticies of the correct degree
        for vertex in range(num_verticies):
            for tile in range(max_tiles):
                cstr = gp.LinExpr()
                for bond_edge in range(bond_edge_types):
                    for hat in range(2):
                        cstr = cstr + vertex_bets_map.get((vertex, bond_edge, hat))
                        cstr = cstr - tile_bets_map.get((tile, bond_edge, hat))
                m.addGenConstrIndicator(vertex_tiles_decision_map.get((vertex, tile)), True, cstr == 0, name="v"+str(vertex)+"t"+str(tile)+"chosen_degree_match")
        
        # # Add constraint: verticies have bond edges equal to degree
        for vertex in range(num_verticies):
            cstr = gp.LinExpr()
            for bond_edge in range(bond_edge_types):
                for hat in range(2):
                    cstr = cstr + vertex_bets_map.get((vertex, bond_edge, hat))
            m.addConstr(cstr == deg(vertex), "vertex_"+str(vertex)+"_bets_match_degree")
        
        # # Each vertex must be exactly one tile type
        for vertex in range(num_verticies):
            cstr = gp.LinExpr()
            for tile in range(max_tiles):
                cstr = cstr + vertex_tiles_decision_map.get((vertex, tile))
            m.addConstr(cstr == 1, "v"+str(vertex)+"_one_tile_type")
        
        # # Enforce tile selections
        for vertex in range(num_verticies):
            for tile in range(max_tiles):
                for bond_edge in range(bond_edge_types):
                    for hat in range(2):
                        cstr = gp.LinExpr()
                        cstr = cstr + vertex_bets_map.get((vertex, bond_edge, hat))
                        cstr = cstr - tile_bets_map.get((tile, bond_edge, hat))
                        for othertile in range(max_tiles):
                            if(othertile!=tile):
                                cstr = cstr - max(different_degrees) * vertex_tiles_decision_map.get((vertex, othertile))
                        m.addConstr(cstr <= 0, "enforce_tile_selection_v"+str(vertex)+"t"+str(tile)+"b"+str(bond_edge)+":"+str(hat))

        # # Ks count tiles used
        for tile in range(max_tiles):
            cstr = gp.LinExpr()
            for bond_edge in range(bond_edge_types):
                for hat in range(2):
                    cstr = cstr + tile_bets_map.get((tile, bond_edge, hat))
            cstr = cstr - max(different_degrees) * k_map.get((tile))
            m.addConstr(cstr <= 0, name="k_count_tile"+","+str(tile))

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

        # # Each edge has a 2 bond edge exactly
        for edge in edge_list:
            cstr = gp.LinExpr()
            for bond_edge in range(bond_edge_types):
                for hat in range(2):
                    cstr = cstr + edge_constraints_map.get((edge[0], edge[1], bond_edge, hat))
                    cstr = cstr + edge_constraints_map.get((edge[1], edge[0], bond_edge, hat))
            m.addConstr(cstr == 2, name="edge_"+str(edge)+"_2bet")

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

        # Use first k tiles
        for tile in range(max_tiles-1):
            m.addConstr(k_map.get((tile)) >= k_map.get((tile + 1)))

        # Use each tile
        for tile in range(max_tiles):
            cstr = gp.LinExpr()
            for vertex in range(num_verticies):
                cstr = cstr + vertex_tiles_decision_map.get((vertex, tile))
            m.addConstr(cstr >= k_map.get((tile)), name="t"+str(tile)+"_used")
        
        l = max(different_degrees) + 1
        for tile in range(max_tiles-1):
            cstr = gp.LinExpr()
            coeff = 0
            for bond_edge in range(bond_edge_types-1, -1, -1):
                for hat in range(1, -1, -1):
                    if(coeff == 0):
                        coeff = 1
                    else:
                        coeff = coeff * l
                    cstr = cstr + coeff * tile_bets_map.get((tile, bond_edge, hat))
                    cstr = cstr - coeff * tile_bets_map.get((tile+1, bond_edge, hat))
            m.addGenConstrIndicator(k_map.get((tile)), True, cstr >= 1, name="strictly_decreasing_a_t"+str(tile))

        for tile_perm in partitions:
            # print("ENFORCE: ", tile_perm)
            # print(tile_perm)
            nums = []
            for bond_edge in range(bond_edge_types):
                cstr = gp.LinExpr()
                for tile in range(max_tiles):
                    cstr = cstr + (int(tile_perm[tile]) * tile_bets_map.get((tile, bond_edge, 0)))
                    cstr = cstr - (int(tile_perm[tile]) * tile_bets_map.get((tile, bond_edge, 1)))
                valz = m.addVar(lb=-100, ub=100, vtype=GRB.INTEGER, name=str(tile_perm)+str(bond_edge)+str(tile))
                m.addConstr(cstr == valz)
                nums.append(valz)
            val2 = m.addVar(lb=0, ub=100)
            z = m.addGenConstrNorm(val2, nums, GRB.INFINITY)
            m.addConstr(val2 >= 1)

        # Optimize model
        if(not print_log):
            m.setParam(GRB.Param.LogToConsole, 0)

        m.optimize()

        # Turn our solutions into pots
        half_edges, half_edges_hat = get_half_edge_labels()
        if(m.status != 2):
            if(m.status == 3):
                return False, [], [], [] #
            else:
                print("Status of optimization: ", m.status)
                return None

        # for tile_num in range(max_tiles):
        #     print("Tile: " + str(tile_num))
        #     for bond_edge in range(bond_edge_types):
        #         print("Bond edge " + str(bond_edge) + " :" + str(tile_bets_map.get((tile_num, bond_edge, 0)).X))
        #         print("Bond edge " + str(bond_edge) + " hat :" + str(tile_bets_map.get((tile_num, bond_edge, 1)).X))

        # * Construct our pot and orientation, and return it----------------------------------------------------------
        # Get our dictionary of tile assignments and build our pot
        pot = []
        tile_assignments = {}
        for tile_num in range(max_tiles):
            if(int(0.1 + k_map.get((tile_num)).X) >= 1):
            # if(True):
                tile = ''
                for bond_edge in range(bond_edge_types):
                    tile = tile + half_edges[bond_edge] * int(0.1 + tile_bets_map.get((tile_num, bond_edge, 0)).X)
                    tile = tile + half_edges_hat[bond_edge] * int(0.1 + tile_bets_map.get((tile_num, bond_edge, 1)).X)
                pot.append(tile)
                tile_assignments.update({tile : []})
                #add each vertex assigned to the tile
                for vertex in range(num_verticies):
                    if(int( 0.1 + vertex_tiles_decision_map.get((vertex, tile_num)).X) >= 1):
                        tile_assignments.get(tile).append(vertex)

        # Build the edges of the orientation of the graph
        orientations = []
        for bond_edge in range(bond_edge_types):
            orientation = []
            for vertex in range(num_verticies):
                this_row = []
                for othervertex in range(num_verticies):
                    #Since the edge list only has one copy of each edge, we must check both directions
                    if((vertex, othervertex) in edge_list or (othervertex, vertex) in edge_list):
                        this_row.append(int(0.1 + edge_constraints_map.get((vertex, othervertex, bond_edge, 0)).X))
                    else:
                        this_row.append(0)
                orientation.append(this_row)
            orientations.append(orientation)

        return True, pot, tile_assignments, orientations


    except gp.GurobiError as e:
        print('Error code ' + str(e.errno) + ': ' + str(e))
    

def optimal_pot_s2_brute_force_partition(G, bond_edge_types, num_tiles, print_log):
    # Extract some graph info
    degree_sequence = [d for n,d in G.degree()]
    degree_sequence.sort()
    different_degrees = list(set(degree_sequence))
    num_verticies = G.number_of_nodes()
    edge_list = list(G.edges())

    def deg(node):
        deg = 0
        for edge in edge_list:
            if(node in edge):
                deg = deg + 1
        return deg

    def other_node(node, edge):
        if(edge[0]==node):
            return edge[1]
        elif(edge[1]==node):
            return edge[0]
    try:
        # Create a new model
        m = gp.Model("optimal_pot_s2")

        # Create the dictionaries of variables
        # indexing: (vertex, bet, hat)
        vertex_bets_map = {}
        # indexing: (tile, bet, hat)
        tile_bets_map = {}
        # indexing: (vertex, tile)
        vertex_tiles_decision_map = {}
        # indexing: (tile)
        k_map = {}
        # indexing: (tile)
        u_map = {}
        # indexing: (vertex_from, vertex_to, bet, hat)
        edge_constraints_map = {}

        # * Create variables
        # make a variable for each bond edge of each vertex
        for vertex in range(num_verticies):
            for bond_edge in range(bond_edge_types):
                for hat in range(2):
                    name = "node" + str(vertex) + "bet" + str(bond_edge) + ":" + str(hat) 
                    x = m.addVar(vtype=GRB.INTEGER, lb=0, name=name)
                    vertex_bets_map.update({(vertex, bond_edge, hat) : x})

        # make a variable for each bond edge of each tile
        for tile in range(num_tiles):
            for bond_edge in range(bond_edge_types):
                for hat in range(2):
                    # name = "tile" + str(degree) + ":" + str(tile) + "bet" + str(bond_edge) + ":" + str(hat) 
                    name = "!" + "," + str(tile) + "," + str(bond_edge) + "," + str(hat)
                    x = m.addVar(vtype=GRB.INTEGER, lb=0, name=name)
                    tile_bets_map.update({(tile, bond_edge, hat) : x})

        # For each pair (vertex and tile), create a decision var
        for vertex in range(num_verticies):
            for tile in range(num_tiles):
                name = "node" +  str(vertex) + "istile" + str(tile)
                x = m.addVar(vtype=GRB.BINARY, name=name)
                vertex_tiles_decision_map.update({(vertex, tile) : x})

        # For each tile type, create a corresponding k var
        for tile in range(num_tiles):
            name = "k_tile" + str(tile)
            x = m.addVar(vtype=GRB.BINARY, name=name)
            k_map.update({(tile) : x})

        # For each tile type, create a corresponding u var
        # for tile in range(num_tiles):
        #     name = "u_tile" + str(tile)
        #     x = m.addVar(vtype=GRB.BINARY, name=name)
        #     u_map.update({(tile) : x})

        # For each side of each edge, create vars of each bet
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

        # * Make our objective(s) 
        obj = gp.LinExpr()
        for k_key in list(k_map.keys()):
            k = k_map.get(k_key)
            obj = obj + k
        m.setObjective(obj, GRB.MINIMIZE) 

        # * Make our constraints

        # Add constraint: same number of each bond edge type in all verticies (a graph can be made)
        for bond_edge in range(bond_edge_types):
            cstr = gp.LinExpr()
            for vertex in range(num_verticies):
                bet_no_hat = vertex_bets_map.get((vertex, bond_edge, 0))
                bet_hat = vertex_bets_map.get((vertex, bond_edge, 1))
                cstr = cstr + bet_no_hat - bet_hat
            m.addConstr(cstr == 0, "bond_edge_aA_match_"+str(bond_edge))

        # Add constraint: tiles are assigned to verticies of the correct degree
        for vertex in range(num_verticies):
            for tile in range(num_tiles):
                cstr = gp.LinExpr()
                for bond_edge in range(bond_edge_types):
                    for hat in range(2):
                        cstr = cstr + vertex_bets_map.get((vertex, bond_edge, hat))
                        cstr = cstr - tile_bets_map.get((tile, bond_edge, hat))
                m.addGenConstrIndicator(vertex_tiles_decision_map.get((vertex, tile)), True, cstr == 0, name="v"+str(vertex)+"t"+str(tile)+"chosen_degree_match")
        
        # # Add constraint: verticies have bond edges equal to degree
        for vertex in range(num_verticies):
            cstr = gp.LinExpr()
            for bond_edge in range(bond_edge_types):
                for hat in range(2):
                    cstr = cstr + vertex_bets_map.get((vertex, bond_edge, hat))
            m.addConstr(cstr == deg(vertex), "vertex_"+str(vertex)+"_bets_match_degree")
        
        # # Each vertex must be exactly one tile type
        for vertex in range(num_verticies):
            cstr = gp.LinExpr()
            for tile in range(num_tiles):
                cstr = cstr + vertex_tiles_decision_map.get((vertex, tile))
            m.addConstr(cstr == 1, "v"+str(vertex)+"_one_tile_type")
        
        # # Enforce tile selections
        for vertex in range(num_verticies):
            for tile in range(num_tiles):
                for bond_edge in range(bond_edge_types):
                    for hat in range(2):
                        cstr = gp.LinExpr()
                        cstr = cstr + vertex_bets_map.get((vertex, bond_edge, hat))
                        cstr = cstr - tile_bets_map.get((tile, bond_edge, hat))
                        for othertile in range(num_tiles):
                            if(othertile!=tile):
                                cstr = cstr - max(different_degrees) * vertex_tiles_decision_map.get((vertex, othertile))
                        m.addConstr(cstr <= 0, "enforce_tile_selection_v"+str(vertex)+"t"+str(tile)+"b"+str(bond_edge)+":"+str(hat))

        # # Ks count tiles used
        for tile in range(num_tiles):
            cstr = gp.LinExpr()
            for bond_edge in range(bond_edge_types):
                for hat in range(2):
                    cstr = cstr + tile_bets_map.get((tile, bond_edge, hat))
            cstr = cstr - max(different_degrees) * k_map.get((tile))
            m.addConstr(cstr <= 0, name="k_count_tile"+","+str(tile))

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

        # # Each edge has a 2 bond edge exactly
        for edge in edge_list:
            cstr = gp.LinExpr()
            for bond_edge in range(bond_edge_types):
                for hat in range(2):
                    cstr = cstr + edge_constraints_map.get((edge[0], edge[1], bond_edge, hat))
                    cstr = cstr + edge_constraints_map.get((edge[1], edge[0], bond_edge, hat))
            m.addConstr(cstr == 2, name="edge_"+str(edge)+"_2bet")

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

        # Sum of us with bond edges is 0 (Linearize, someday)
        # for bond_edge in range(bond_edge_types):
        #     cstr = gp.LinExpr()
        #     for tile in range(num_tiles):
        #         num_bonds_used = tile_bets_map.get((tile, bond_edge, 0)) * u_map.get((tile))
        #         num_bonds_used_hat = tile_bets_map.get((tile, bond_edge, 0)) * u_map.get((tile))
        #         cstr = cstr + num_bonds_used - num_bonds_used_hat
        #     m.addConstr(cstr == 0)


        #!Integer Partition Madness
        #Get all the possible sequences of us
        def get_glist(number, listlen):
            def ruleAscLen(n, l):
                a = [0 for i in range(n + 1)]
                k = 1
                a[0] = 0
                a[1] = n
                while k != 0:
                    x = a[k - 1] + 1
                    y = a[k] - 1
                    k -= 1
                    while x <= y and k < l - 1:
                        a[k] = x
                        y -= x
                        k += 1
                    a[k] = x + y
                    yield a[:k + 1]

            n = number

            gs = [ruleAscLen(x, listlen) for x in range(1, n)]
            # g = ruleAscLen(n, 5)

            lsts = [[i for i in g] for g in gs]

            for lst in lsts:
                for elem in lst:
                    while(len(elem)!=listlen):
                        elem.append(0)

            blsts = []

            for lst in lsts:
                blsts.extend(lst)

            glst = []

            for blst in blsts:
                # print(blst)
                glst.extend(itertools.permutations(blst))

            return glst

        from math import factorial, gcd
        from functools import reduce
        def find_gcd(list):
            x = reduce(gcd, list)
            return x

        tile_perms = list(set(get_glist(G.number_of_nodes(), num_tiles)))
        p2 = [perm for perm in tile_perms if find_gcd(perm)==1]
        tile_perms = p2

        numin = 3
        numu = 0

        for tile_perm in tile_perms:
            # print(tile_perm)
            bigcstr = gp.LinExpr()
            for bond_edge in range(bond_edge_types):
                cstr = gp.LinExpr()
                for tile in range(num_tiles):
                    cstr = cstr + (int(tile_perm[tile]) * tile_bets_map.get((tile, bond_edge, 0)))
                    cstr = cstr - (int(tile_perm[tile]) * tile_bets_map.get((tile, bond_edge, 1)))
                gabs = m.addVar(lb=-1000, ub=1000, vtype=GRB.CONTINUOUS)
                valz = m.addVar(lb=-1000, ub=1000, vtype=GRB.CONTINUOUS, name=str(tile_perm)+str(bond_edge)+str(tile))
                m.addConstr(cstr == valz)
                m.addConstr(gabs == gp.abs_(valz))
                bigcstr = bigcstr + gabs
            m.addConstr(bigcstr >= 1)


        #!Madness over

        # Optimize model
        if(not print_log):
            m.setParam(GRB.Param.LogToConsole, 0)

        m.optimize()

        # for u in range(num_verticies):
        #     print("Uvalue" + str(u) + ": " + str(u_map.get((u)).X))

        # Turn our solutions into pots
        half_edges, half_edges_hat = get_half_edge_labels()
        if(m.status != 2):
            if(m.status == 3):
                return False, [], [] #
            else:
                print("Status of optimization: ", m.status)
                return None

        # * Construct our pot and orientation, and return it----------------------------------------------------------
        # Get our dictionary of tile assignments and build our pot
        pot = []
        tile_assignments = {}
        for tile_num in range(num_tiles):
            if(int(0.1 + k_map.get((tile_num)).X) == 1):
                tile = ''
                for bond_edge in range(bond_edge_types):
                    tile = tile + half_edges[bond_edge] * int(0.1 + tile_bets_map.get((tile_num, bond_edge, 0)).X)
                    tile = tile + half_edges_hat[bond_edge] * int(0.1 + tile_bets_map.get((tile_num, bond_edge, 1)).X)
                pot.append(tile)
                tile_assignments.update({tile : []})
                #add each vertex assigned to the tile
                for vertex in range(num_verticies):
                    if(int( 0.1 + vertex_tiles_decision_map.get((vertex, tile_num)).X) == 1):
                        tile_assignments.get(tile).append(vertex)

        # Build the edges of the orientation of the graph
        orientations = []
        for bond_edge in range(bond_edge_types):
            orientation = []
            for vertex in range(num_verticies):
                this_row = []
                for othervertex in range(num_verticies):
                    #Since the edge list only has one copy of each edge, we must check both directions
                    if((vertex, othervertex) in edge_list or (othervertex, vertex) in edge_list):
                        this_row.append(int(0.1 + edge_constraints_map.get((vertex, othervertex, bond_edge, 0)).X))
                    else:
                        this_row.append(0)
                orientation.append(this_row)
            orientations.append(orientation)

        return pot, tile_assignments, orientations


    except gp.GurobiError as e:
        print('Error code ' + str(e.errno) + ': ' + str(e))


#This is the ILP for scenario 2
def optimal_pot_s2(Graph, bond_edge_types, num_tiles, min_qs, q_and_under_equal, print_log, partitions=[]):
    G = Graph
    # Extract some graph info
    degree_sequence = [d for n,d in G.degree()]
    degree_sequence.sort()
    different_degrees = list(set(degree_sequence))
    num_verticies = G.number_of_nodes()
    edge_list = list(G.edges())

    def deg(node):
        deg = 0
        for edge in edge_list:
            if(node in edge):
                deg = deg + 1
        return deg

    def other_node(node, edge):
        if(edge[0]==node):
            return edge[1]
        elif(edge[1]==node):
            return edge[0]
    try:
        # Create a new model
        m = gp.Model("optimal_pot_s2")

        # Create the dictionaries of variables
        # indexing: (vertex, bet, hat)
        vertex_bets_map = {}
        # indexing: (tile, bet, hat)
        tile_bets_map = {}
        # indexing: (vertex, tile)
        vertex_tiles_decision_map = {}
        # indexing: (tile)
        k_map = {}
        # indexing: (vertex_from, vertex_to, bet, hat)
        edge_constraints_map = {}

        # * Create variables
        # make a variable for each bond edge of each vertex
        for vertex in range(num_verticies):
            for bond_edge in range(bond_edge_types):
                for hat in range(2):
                    name = "node" + str(vertex) + "bet" + str(bond_edge) + ":" + str(hat) 
                    x = m.addVar(vtype=GRB.INTEGER, lb=0, name=name)
                    vertex_bets_map.update({(vertex, bond_edge, hat) : x})

        # make a variable for each bond edge of each tile
        for tile in range(num_tiles):
            for bond_edge in range(bond_edge_types):
                for hat in range(2):
                    # name = "tile" + str(degree) + ":" + str(tile) + "bet" + str(bond_edge) + ":" + str(hat) 
                    name = "!" + "," + str(tile) + "," + str(bond_edge) + "," + str(hat)
                    x = m.addVar(vtype=GRB.INTEGER, lb=0, name=name)
                    tile_bets_map.update({(tile, bond_edge, hat) : x})

        # For each pair (vertex and tile), create a decision var
        for vertex in range(num_verticies):
            for tile in range(num_tiles):
                name = "node" +  str(vertex) + "istile" + str(tile)
                x = m.addVar(vtype=GRB.BINARY, name=name)
                vertex_tiles_decision_map.update({(vertex, tile) : x})

        # For each tile type, create a corresponding k var
        for tile in range(num_tiles):
            name = "k_tile" + str(tile)
            x = m.addVar(vtype=GRB.BINARY, name=name)
            k_map.update({(tile) : x})

        # For each side of each edge, create vars of each bet
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

        # * Make our objective(s) 
        # Minimize over each q, starting with tile 1
        l = max(different_degrees) + 1
        qs = []
        for index,tile in enumerate(range(num_tiles)):
            q = gp.LinExpr()
            coeff = pow(l, 2*bond_edge_types-1)
            for bet in range(bond_edge_types):
                for hat in range(2):
                    q = q + coeff * tile_bets_map.get((tile, bet, hat))
                    coeff = coeff / l
            qs.append(q)
            m.setObjectiveN(q, index, num_tiles-index)

        # * Make our constraints
        # Make sure we don't exceed min_pot_size
        obj = gp.LinExpr()
        for k_key in list(k_map.keys()):
            k = k_map.get(k_key)
            obj = obj + k
        m.addConstr(obj >= num_tiles, "min_pot_size_limit")

        qbounds = []
        # Make sure we don't exceed min_q
        for tile in range(num_tiles):
            # print(min_qs[index])
            q = gp.LinExpr()
            coeff = pow(l, 2*bond_edge_types-1)
            for bet in range(bond_edge_types):
                for hat in range(2):
                    q = q + coeff * tile_bets_map.get((tile, bet, hat))
                    coeff = coeff / l
            # print(q)
            # print(min_qs[tile])
            c = m.addConstr(q >= min_qs[tile], "minimum q-value " + str(tile))
            qbounds.append(c)

        # Force equality below q_and_under_equal
        for tile in range(q_and_under_equal):
            q = gp.LinExpr()
            coeff = pow(l, 2*bond_edge_types-1)
            for bet in range(bond_edge_types):
                for hat in range(2):
                    q = q + coeff * tile_bets_map.get((tile, bet, hat))
                    coeff = coeff / l
            c = m.addConstr(q == min_qs[tile], "force q-value " + str(tile))
            qbounds.append(c)

        # Add constraint: same number of each bond edge type in all verticies (a graph can be made)
        for bond_edge in range(bond_edge_types):
            cstr = gp.LinExpr()
            for vertex in range(num_verticies):
                bet_no_hat = vertex_bets_map.get((vertex, bond_edge, 0))
                bet_hat = vertex_bets_map.get((vertex, bond_edge, 1))
                cstr = cstr + bet_no_hat - bet_hat
            m.addConstr(cstr == 0, "bond_edge_aA_match_"+str(bond_edge))

        # Add constraint: tiles are assigned to verticies of the correct degree
        for vertex in range(num_verticies):
            for tile in range(num_tiles):
                cstr = gp.LinExpr()
                for bond_edge in range(bond_edge_types):
                    for hat in range(2):
                        cstr = cstr + vertex_bets_map.get((vertex, bond_edge, hat))
                        cstr = cstr - tile_bets_map.get((tile, bond_edge, hat))
                m.addGenConstrIndicator(vertex_tiles_decision_map.get((vertex, tile)), True, cstr == 0, name="v"+str(vertex)+"t"+str(tile)+"chosen_degree_match")
                            # m.addGenConstrIndicator(k_map.get((tile)), True, cstr >= 1, name="strictly_decreasing_a_t"+str(tile))

        
        # # Add constraint: verticies have bond edges equal to degree
        for vertex in range(num_verticies):
            cstr = gp.LinExpr()
            for bond_edge in range(bond_edge_types):
                for hat in range(2):
                    cstr = cstr + vertex_bets_map.get((vertex, bond_edge, hat))
            m.addConstr(cstr == deg(vertex), "vertex_"+str(vertex)+"_bets_match_degree")
        
        # # Each vertex must be exactly one tile type
        for vertex in range(num_verticies):
            cstr = gp.LinExpr()
            for tile in range(num_tiles):
                cstr = cstr + vertex_tiles_decision_map.get((vertex, tile))
            m.addConstr(cstr == 1, "v"+str(vertex)+"_one_tile_type")
        
        # # Enforce tile selections
        for vertex in range(num_verticies):
            for tile in range(num_tiles):
                for bond_edge in range(bond_edge_types):
                    for hat in range(2):
                        cstr = gp.LinExpr()
                        cstr = cstr + vertex_bets_map.get((vertex, bond_edge, hat))
                        cstr = cstr - tile_bets_map.get((tile, bond_edge, hat))
                        for othertile in range(num_tiles):
                            if(othertile!=tile):
                                cstr = cstr - max(different_degrees) * vertex_tiles_decision_map.get((vertex, othertile))
                        m.addConstr(cstr <= 0, "enforce_tile_selection_v"+str(vertex)+"t"+str(tile)+"b"+str(bond_edge)+":"+str(hat))

        # # Ks count tiles used
        for tile in range(num_tiles):
            cstr = gp.LinExpr()
            for bond_edge in range(bond_edge_types):
                for hat in range(2):
                    cstr = cstr + tile_bets_map.get((tile, bond_edge, hat))
            cstr = cstr - max(different_degrees) * k_map.get((tile))
            m.addConstr(cstr <= 0, name="k_count_tile"+","+str(tile))

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

        # # Each edge has a 2 bond edge exactly
        for edge in edge_list:
            cstr = gp.LinExpr()
            for bond_edge in range(bond_edge_types):
                for hat in range(2):
                    cstr = cstr + edge_constraints_map.get((edge[0], edge[1], bond_edge, hat))
                    cstr = cstr + edge_constraints_map.get((edge[1], edge[0], bond_edge, hat))
            m.addConstr(cstr == 2, name="edge_"+str(edge)+"_2bet")

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

        # * Limit our solution pool
        # Nonzero tiles are distinct (Limits {AAB, AAB, bba}-like solutions)
        for tile in range(num_tiles-1):
            cstr = gp.LinExpr()
            coeff = 0
            for bond_edge in range(bond_edge_types-1, -1, -1):
                for hat in range(1, -1, -1):
                    if(coeff == 0):
                        coeff = 1
                    else:
                        coeff = coeff * l
                    cstr = cstr + coeff * tile_bets_map.get((tile, bond_edge, hat))
                    cstr = cstr - coeff * tile_bets_map.get((tile+1, bond_edge, hat))
            m.addGenConstrIndicator(k_map.get((tile)), True, cstr >= 1, name="strictly_decreasing_a_t"+str(tile))

        # Use each bond edge type in a tile
        for bond_edge in range(bond_edge_types):
            for hat in range(2):
                cstr = gp.LinExpr()
                for tile in range(num_tiles):
                    cstr = cstr + tile_bets_map.get((tile, bond_edge, hat))
                m.addConstr(cstr >= 1, "bet_"+str(bond_edge)+":"+str(hat)+"_used")

        # Use each tile
        for tile in range(num_tiles):
            cstr = gp.LinExpr()
            for vertex in range(num_verticies):
                cstr = cstr + vertex_tiles_decision_map.get((vertex, tile))
            m.addConstr(cstr >= 1, name="t"+str(tile)+"_used")

        #Add partition stuff
        for tile_perm in partitions:
            # print("ENFORCE: ", tile_perm)
            # print(tile_perm)
            nums = []
            for bond_edge in range(bond_edge_types):
                cstr = gp.LinExpr()
                for tile in range(num_tiles):
                    cstr = cstr + (int(tile_perm[tile]) * tile_bets_map.get((tile, bond_edge, 0)))
                    cstr = cstr - (int(tile_perm[tile]) * tile_bets_map.get((tile, bond_edge, 1)))
                valz = m.addVar(lb=-100, ub=100, vtype=GRB.CONTINUOUS, name=str(tile_perm)+str(bond_edge)+str(tile))
                m.addConstr(cstr == valz)
                nums.append(valz)
            val2 = m.addVar(lb=0, ub=100)
            z = m.addGenConstrNorm(val2, nums, GRB.INFINITY)
            m.addConstr(val2 >= 1)
            

        # Optimize model
        if(not print_log):
            m.setParam(GRB.Param.LogToConsole, 0)

        m.optimize()



        # Turn our solutions into pots
        half_edges, half_edges_hat = get_half_edge_labels()
        if(m.status != 2):
            if(m.status == 3):
                return [], [], [], [], [] #
            else:
                print("Status of optimization: ", m.status)
                return None

        qvals = [q.getValue() for q in qs]

        # * Construct our pot and orientation, and return it----------------------------------------------------------
        # Get our dictionary of tile assignments and build our pot
        pot = []
        tile_assignments = {}
        for tile_num in range(num_tiles):
            if(int(0.1 + k_map.get((tile_num)).X) == 1):
                tile = ''
                for bond_edge in range(bond_edge_types):
                    tile = tile + half_edges[bond_edge] * int(0.1 + tile_bets_map.get((tile_num, bond_edge, 0)).X)
                    tile = tile + half_edges_hat[bond_edge] * int(0.1 + tile_bets_map.get((tile_num, bond_edge, 1)).X)
                pot.append(tile)
                tile_assignments.update({tile : []})
                #add each vertex assigned to the tile
                for vertex in range(num_verticies):
                    if(int(0.1 + vertex_tiles_decision_map.get((vertex, tile_num)).X) == 1):
                        tile_assignments.get(tile).append(vertex)

        # Build the edges of the orientation of the graph
        orientations = []
        for bond_edge in range(bond_edge_types):
            orientation = []
            for vertex in range(num_verticies):
                this_row = []
                for othervertex in range(num_verticies):
                    #Since the edge list only has one copy of each edge, we must check both directions
                    if((vertex, othervertex) in edge_list or (othervertex, vertex) in edge_list):
                        this_row.append(int(0.1 + edge_constraints_map.get((vertex, othervertex, bond_edge, 0)).X))
                    else:
                        this_row.append(0)
                orientation.append(this_row)
            orientations.append(orientation)

        # Primarily for debugging: comment out to print all decision variable values
        # for v in m.getVars():
        #     print('%s %g' % (v.VarName, v.X))
        savedata = [m, vertex_bets_map, tile_bets_map, vertex_tiles_decision_map, k_map, edge_constraints_map, qbounds, qs]

        return pot, tile_assignments, orientations, qvals, savedata


    except gp.GurobiError as e:
        print('Error code ' + str(e.errno) + ': ' + str(e))

def optimal_pot_s2_canonical(Graph, bond_edge_types, num_tiles, min_qs, q_and_under_equal, print_log, partitions=[]):
    G = Graph
    # Extract some graph info
    degree_sequence = [d for n,d in G.degree()]
    degree_sequence.sort()
    different_degrees = list(set(degree_sequence))
    num_verticies = G.number_of_nodes()
    edge_list = list(G.edges())

    def deg(node):
        deg = 0
        for edge in edge_list:
            if(node in edge):
                deg = deg + 1
        return deg

    def other_node(node, edge):
        if(edge[0]==node):
            return edge[1]
        elif(edge[1]==node):
            return edge[0]
    try:
        # Create a new model
        m = gp.Model("optimal_pot_s2")

        # Create the dictionaries of variables
        # indexing: (vertex, bet, hat)
        vertex_bets_map = {}
        # indexing: (tile, bet, hat)
        tile_bets_map = {}
        # indexing: (vertex, tile)
        vertex_tiles_decision_map = {}
        # indexing: (tile)
        k_map = {}
        # indexing: (vertex_from, vertex_to, bet, hat)
        edge_constraints_map = {}

        # * Create variables
        # make a variable for each bond edge of each vertex
        for vertex in range(num_verticies):
            for bond_edge in range(bond_edge_types):
                for hat in range(2):
                    name = "node" + str(vertex) + "bet" + str(bond_edge) + ":" + str(hat) 
                    x = m.addVar(vtype=GRB.INTEGER, lb=0, name=name)
                    vertex_bets_map.update({(vertex, bond_edge, hat) : x})

        # make a variable for each bond edge of each tile
        for tile in range(num_tiles):
            for bond_edge in range(bond_edge_types):
                for hat in range(2):
                    # name = "tile" + str(degree) + ":" + str(tile) + "bet" + str(bond_edge) + ":" + str(hat) 
                    name = "!" + "," + str(tile) + "," + str(bond_edge) + "," + str(hat)
                    x = m.addVar(vtype=GRB.INTEGER, lb=0, name=name)
                    tile_bets_map.update({(tile, bond_edge, hat) : x})

        # For each pair (vertex and tile), create a decision var
        for vertex in range(num_verticies):
            for tile in range(num_tiles):
                name = "node" +  str(vertex) + "istile" + str(tile)
                x = m.addVar(vtype=GRB.BINARY, name=name)
                vertex_tiles_decision_map.update({(vertex, tile) : x})

        # For each tile type, create a corresponding k var
        for tile in range(num_tiles):
            name = "k_tile" + str(tile)
            x = m.addVar(vtype=GRB.BINARY, name=name)
            k_map.update({(tile) : x})

        # For each side of each edge, create vars of each bet
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

        # * Make our objective(s) 
        # Minimize over each q, starting with tile 1
        l = max(different_degrees) + 1
        qs = []
        for index,tile in enumerate(range(num_tiles)):
            q = gp.LinExpr()
            coeff = pow(l, 2*bond_edge_types-1)
            for bet in range(bond_edge_types):
                for hat in range(2):
                    q = q + coeff * tile_bets_map.get((tile, bet, hat))
                    coeff = coeff / l
            qs.append(q)
            m.setObjectiveN(q, index, num_tiles-index)

        # * Make our constraints
        # Make sure we don't exceed min_pot_size
        obj = gp.LinExpr()
        for k_key in list(k_map.keys()):
            k = k_map.get(k_key)
            obj = obj + k
        m.addConstr(obj >= num_tiles, "min_pot_size_limit")

        qbounds = []
        # Make sure we don't exceed min_q
        for tile in range(num_tiles):
            # print(min_qs[index])
            q = gp.LinExpr()
            coeff = pow(l, 2*bond_edge_types-1)
            for bet in range(bond_edge_types):
                for hat in range(2):
                    q = q + coeff * tile_bets_map.get((tile, bet, hat))
                    coeff = coeff / l
            # print(q)
            # print(min_qs[tile])
            c = m.addConstr(q >= min_qs[tile], "minimum q-value " + str(tile))
            qbounds.append(c)

        # Force equality below q_and_under_equal
        for tile in range(q_and_under_equal):
            q = gp.LinExpr()
            coeff = pow(l, 2*bond_edge_types-1)
            for bet in range(bond_edge_types):
                for hat in range(2):
                    q = q + coeff * tile_bets_map.get((tile, bet, hat))
                    coeff = coeff / l
            c = m.addConstr(q == min_qs[tile], "force q-value " + str(tile))
            qbounds.append(c)

        # Add constraint: same number of each bond edge type in all verticies (a graph can be made)
        for bond_edge in range(bond_edge_types):
            cstr = gp.LinExpr()
            for vertex in range(num_verticies):
                bet_no_hat = vertex_bets_map.get((vertex, bond_edge, 0))
                bet_hat = vertex_bets_map.get((vertex, bond_edge, 1))
                cstr = cstr + bet_no_hat - bet_hat
            m.addConstr(cstr == 0, "bond_edge_aA_match_"+str(bond_edge))

        # Add constraint: tiles are assigned to verticies of the correct degree
        for vertex in range(num_verticies):
            for tile in range(num_tiles):
                cstr = gp.LinExpr()
                for bond_edge in range(bond_edge_types):
                    for hat in range(2):
                        cstr = cstr + vertex_bets_map.get((vertex, bond_edge, hat))
                        cstr = cstr - tile_bets_map.get((tile, bond_edge, hat))
                m.addGenConstrIndicator(vertex_tiles_decision_map.get((vertex, tile)), True, cstr == 0, name="v"+str(vertex)+"t"+str(tile)+"chosen_degree_match")
                            # m.addGenConstrIndicator(k_map.get((tile)), True, cstr >= 1, name="strictly_decreasing_a_t"+str(tile))

        
        # # Add constraint: verticies have bond edges equal to degree
        for vertex in range(num_verticies):
            cstr = gp.LinExpr()
            for bond_edge in range(bond_edge_types):
                for hat in range(2):
                    cstr = cstr + vertex_bets_map.get((vertex, bond_edge, hat))
            m.addConstr(cstr == deg(vertex), "vertex_"+str(vertex)+"_bets_match_degree")
        
        # # Each vertex must be exactly one tile type
        for vertex in range(num_verticies):
            cstr = gp.LinExpr()
            for tile in range(num_tiles):
                cstr = cstr + vertex_tiles_decision_map.get((vertex, tile))
            m.addConstr(cstr == 1, "v"+str(vertex)+"_one_tile_type")
        
        # # Enforce tile selections
        for vertex in range(num_verticies):
            for tile in range(num_tiles):
                for bond_edge in range(bond_edge_types):
                    for hat in range(2):
                        cstr = gp.LinExpr()
                        cstr = cstr + vertex_bets_map.get((vertex, bond_edge, hat))
                        cstr = cstr - tile_bets_map.get((tile, bond_edge, hat))
                        for othertile in range(num_tiles):
                            if(othertile!=tile):
                                cstr = cstr - max(different_degrees) * vertex_tiles_decision_map.get((vertex, othertile))
                        m.addConstr(cstr <= 0, "enforce_tile_selection_v"+str(vertex)+"t"+str(tile)+"b"+str(bond_edge)+":"+str(hat))

        # # Ks count tiles used
        for tile in range(num_tiles):
            cstr = gp.LinExpr()
            for bond_edge in range(bond_edge_types):
                for hat in range(2):
                    cstr = cstr + tile_bets_map.get((tile, bond_edge, hat))
            cstr = cstr - max(different_degrees) * k_map.get((tile))
            m.addConstr(cstr <= 0, name="k_count_tile"+","+str(tile))

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

        # # Each edge has a 2 bond edge exactly
        for edge in edge_list:
            cstr = gp.LinExpr()
            for bond_edge in range(bond_edge_types):
                for hat in range(2):
                    cstr = cstr + edge_constraints_map.get((edge[0], edge[1], bond_edge, hat))
                    cstr = cstr + edge_constraints_map.get((edge[1], edge[0], bond_edge, hat))
            m.addConstr(cstr == 2, name="edge_"+str(edge)+"_2bet")

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

        # * Limit our solution pool
        #------------------------ CANONICAL -----------------------------------------
        #! First, order tiles
        # Tiles are first ordered by size (decreasing, so largest first)
        for tile in range(num_tiles-1):
            cstr = gp.LinExpr()
            #Add this tile's size
            for bond_edge in range(bond_edge_types):
                for hat in range(2):
                    cstr = cstr + tile_bets_map.get((tile, bond_edge, hat))
            #Subtract the next tile's size
            for bond_edge in range(bond_edge_types):
                for hat in range(2):
                    cstr = cstr - tile_bets_map.get((tile+1, bond_edge, hat))
            #This must be at least 0
            m.addConstr(cstr >= 0)

        for vertex in range(num_verticies):
            for bond_edge in range(bond_edge_types):
                for hat in range(2):
                    name = "node" + str(vertex) + "bet" + str(bond_edge) + ":" + str(hat) 
                    x = m.addVar(vtype=GRB.INTEGER, lb=0, name=name)
                    vertex_bets_map.update({(vertex, bond_edge, hat) : x})

        # set up indicators for tiebreaker
        indicators = []
        for tile in range(num_tiles-1):
            x = m.addVar(vtype=GRB.BINARY)
            indicators.append(x)

            cstr = gp.LinExpr()
            #Add this tile's size
            for bond_edge in range(bond_edge_types):
                for hat in range(2):
                    cstr = cstr + tile_bets_map.get((tile, bond_edge, hat))
            #Subtract the next tile's size
            for bond_edge in range(bond_edge_types):
                for hat in range(2):
                    cstr = cstr - tile_bets_map.get((tile+1, bond_edge, hat))

            cstr = cstr - l * x

            #If they are equal, the indicator is off. If they are unequal, indicator is on.
            m.addConstr(cstr <= 0)
            m.addConstr(cstr >= (-1 * l) + 1)

        #Order the ties by number of As, then Bs, then Cs, ... (if indicator is off)
        for tile in range(num_tiles-1):
            thisscore = gp.LinExpr()
            nextscore = gp.LinExpr()

            mul = 1
            for bond_edge in range(bond_edge_types):
                thisscore = thisscore + mul * tile_bets_map.get((tile, bond_edge, 0))
                thisscore = thisscore + mul * tile_bets_map.get((tile, bond_edge, 1))

                nextscore = nextscore + mul * tile_bets_map.get((tile + 1, bond_edge, 0))
                nextscore = nextscore + mul * tile_bets_map.get((tile + 1, bond_edge, 1))
                mul = mul * (max(different_degrees) + 1)
            
            #if indicator off, this score has to be bigger
            m.addGenConstrIndicator(indicators[tile], False, thisscore >= nextscore)     

        #! Then, order bonds
        #Number of As > Number of Bs > Number of Cs > ...
        for bond_edge in range(bond_edge_types-1):
            cstr = gp.LinExpr()
            for tile in range(num_tiles):
                cstr = cstr + tile_bets_map.get((tile, bond_edge, 0))
                cstr = cstr + tile_bets_map.get((tile, bond_edge, 1))
            for tile in range(num_tiles):
                cstr = cstr - tile_bets_map.get((tile, bond_edge+1, 0))
                cstr = cstr - tile_bets_map.get((tile, bond_edge+1, 1))
            m.addConstr(cstr >= 0)

        #! ehh idk
        # for tile in range(num_tiles-1):
        #     cstr = gp.LinExpr()
        #     coeff = 0
        #     for bond_edge in range(bond_edge_types-1):
        #         for hat in range(1):
        #             if(coeff == 0):
        #                 coeff = 1
        #             else:
        #                 coeff = coeff * l
        #             cstr = cstr + coeff * tile_bets_map.get((tile, bond_edge, hat))
        #             cstr = cstr - coeff * tile_bets_map.get((tile+1, bond_edge, hat))
        #     m.addGenConstrIndicator(k_map.get((tile)), True, cstr >= 0, name="strictly_decreasing_a_t"+str(tile))

        #! Then, order hats
        # A >= hatA, B >= hatB, C >= hatC
        for bond_edge in range(bond_edge_types):
            cstr = gp.LinExpr()
            for tile in range(num_tiles):
                #Add number of unhatted
                cstr = cstr + tile_bets_map.get((tile, bond_edge, 0))
                #Subtract number hatted
                cstr = cstr - tile_bets_map.get((tile, bond_edge, 1))
            #Must be at least 0
            m.addConstr(cstr >= 0)



        #------------------------ CANONICAL -----------------------------------------

        # Use each bond edge type in a tile
        for bond_edge in range(bond_edge_types):
            for hat in range(2):
                cstr = gp.LinExpr()
                for tile in range(num_tiles):
                    cstr = cstr + tile_bets_map.get((tile, bond_edge, hat))
                m.addConstr(cstr >= 1, "bet_"+str(bond_edge)+":"+str(hat)+"_used")

        # Use each tile
        for tile in range(num_tiles):
            cstr = gp.LinExpr()
            for vertex in range(num_verticies):
                cstr = cstr + vertex_tiles_decision_map.get((vertex, tile))
            m.addConstr(cstr >= 1, name="t"+str(tile)+"_used")

        #Add partition stuff
        for tile_perm in partitions:
            # print("ENFORCE: ", tile_perm)
            # print(tile_perm)
            nums = []
            for bond_edge in range(bond_edge_types):
                cstr = gp.LinExpr()
                for tile in range(num_tiles):
                    cstr = cstr + (int(tile_perm[tile]) * tile_bets_map.get((tile, bond_edge, 0)))
                    cstr = cstr - (int(tile_perm[tile]) * tile_bets_map.get((tile, bond_edge, 1)))
                valz = m.addVar(lb=-100, ub=100, vtype=GRB.CONTINUOUS, name=str(tile_perm)+str(bond_edge)+str(tile))
                m.addConstr(cstr == valz)
                nums.append(valz)
            val2 = m.addVar(lb=0, ub=100)
            z = m.addGenConstrNorm(val2, nums, GRB.INFINITY)
            m.addConstr(val2 >= 1)
            

        # Optimize model
        if(not print_log):
            m.setParam(GRB.Param.LogToConsole, 0)

        m.setParam(GRB.Param.Symmetry, 2)

        m.optimize()



        # Turn our solutions into pots
        half_edges, half_edges_hat = get_half_edge_labels()
        if(m.status != 2):
            if(m.status == 3):
                return [], [], [], [], [] #
            else:
                print("Status of optimization: ", m.status)
                return None

        qvals = [q.getValue() for q in qs]

        # * Construct our pot and orientation, and return it----------------------------------------------------------
        # Get our dictionary of tile assignments and build our pot
        pot = []
        tile_assignments = {}
        for tile_num in range(num_tiles):
            if(int(0.1 + k_map.get((tile_num)).X) == 1):
                tile = ''
                for bond_edge in range(bond_edge_types):
                    tile = tile + half_edges[bond_edge] * int(0.1 + tile_bets_map.get((tile_num, bond_edge, 0)).X)
                    tile = tile + half_edges_hat[bond_edge] * int(0.1 + tile_bets_map.get((tile_num, bond_edge, 1)).X)
                pot.append(tile)
                tile_assignments.update({tile : []})
                #add each vertex assigned to the tile
                for vertex in range(num_verticies):
                    if(int(0.1 + vertex_tiles_decision_map.get((vertex, tile_num)).X) == 1):
                        tile_assignments.get(tile).append(vertex)

        # Build the edges of the orientation of the graph
        orientations = []
        for bond_edge in range(bond_edge_types):
            orientation = []
            for vertex in range(num_verticies):
                this_row = []
                for othervertex in range(num_verticies):
                    #Since the edge list only has one copy of each edge, we must check both directions
                    if((vertex, othervertex) in edge_list or (othervertex, vertex) in edge_list):
                        this_row.append(int(0.1 + edge_constraints_map.get((vertex, othervertex, bond_edge, 0)).X))
                    else:
                        this_row.append(0)
                orientation.append(this_row)
            orientations.append(orientation)

        # Primarily for debugging: comment out to print all decision variable values
        # for v in m.getVars():
        #     print('%s %g' % (v.VarName, v.X))
        savedata = [m, vertex_bets_map, tile_bets_map, vertex_tiles_decision_map, k_map, edge_constraints_map, qbounds, qs]

        return pot, tile_assignments, orientations, qvals, savedata


    except gp.GurobiError as e:
        print('Error code ' + str(e.errno) + ': ' + str(e))

def free_variable_solve(pot, bond_edge_types):
    possible_half_edges, possible_half_edges_hat = get_half_edge_labels()
    matrix = []
    for tile in pot:
        sums = [0 for i in range(bond_edge_types)]
        for item in tile:
            if item in possible_half_edges:
                sums[possible_half_edges.index(item)] = sums[possible_half_edges.index(item)] + 1
            elif item in possible_half_edges_hat:
                sums[possible_half_edges_hat.index(item)] = sums[possible_half_edges_hat.index(item)] - 1
        matrix.append(sums)
    matrix = np.transpose(matrix)

    m = gp.Model("free_variable_problem")

    variables = []
    obj = gp.LinExpr()
    for tile in range(len(pot)):
        x = m.addVar(vtype=GRB.INTEGER, lb=0)
        obj = obj + x
        variables.append(x)
    m.setObjective(obj, GRB.MINIMIZE)
    m.addConstr(obj >= 1)
    for row in matrix:
        cstr = gp.LinExpr()
        for i, val in enumerate(row):
            cstr = cstr + val * variables[i]
        m.addConstr(cstr == 0)
    m.setParam(GRB.Param.LogToConsole, 0)
    m.optimize()
    return m.ObjVal, m.x