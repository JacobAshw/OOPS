from OOPS_files.gurobi import *
from OOPS_files.methods import *

#Hybrid tuining params
num_to_partition = 1 + 0

def S1_optimal_pot(Graph, gurobi_flag):
    pot, tile_assignments, orientation = optimal_pot_s1(Graph, gurobi_flag)
    return pot, tile_assignments, orientation

def S2_tiles_brute_force_partition(Graph, gurobi_flag):
    bond_edges = nx.radius(Graph) + 1
    for num_tiles in range(1, Graph.number_of_nodes()):
        print("Checking number of tiles: ", num_tiles)
        pot, tile_assignments, orientations = optimal_pot_s2_brute_force_partition(Graph, bond_edges, num_tiles, gurobi_flag)
        if(pot != False):
            break
    return pot, tile_assignments, orientations

def S2_bonds_brute_force_partition(Graph, gurobi_flag, optimal_tiles, optimal_tile_bonds):
    tiles = Graph.number_of_nodes()
    for num_bonds in range(1, int(optimal_tile_bonds)):
        print("Checking number of bonds: ", num_bonds)
        pot, tile_assignments, orientations = optimal_pot_s2_brute_force_partition(Graph, num_bonds, tiles, gurobi_flag)
        if(pot != False):
            break
    return pot, tile_assignments, orientations

def S2_tiles_partition_relaxation(Graph, gurobi_flag):
    max_tiles = Graph.number_of_nodes()
    max_bonds = nx.radius(Graph) + 1
    partitions = []
    current_tiles = 3
    while(True):
        solved, pot, tile_assignments, orientations = optimal_pot_s2_relaxation_partition(Graph, max_bonds, current_tiles, gurobi_flag, partitions)
        if(not solved):
            for ratios in partitions:
                ratios.append(0)
            current_tiles = current_tiles + 1
            continue
        print(pot)
        print(orientations)
        print(tile_assignments)
        minsize, ratios = free_variable_solve(pot, max_bonds)
        if(minsize < Graph.number_of_nodes()):
            while(len(ratios) < max_tiles):
                ratios.append(0)
            partitions.append(ratios)
            continue
        break
    return pot, tile_assignments, orientations

def S2_bonds_partition_relaxation(Graph, gurobi_flag, optimal_tiles, optimal_tile_bonds):
    max_tiles = Graph.number_of_nodes()
    max_bonds = nx.radius(Graph) + 1
    partitions = []
    for num_bonds in range(1, int(optimal_tile_bonds)):
        while(True):
            solved, pot, tile_assignments, orientations = optimal_pot_s2_relaxation_partition(Graph, num_bonds, max_tiles, gurobi_flag, partitions)
            if(not solved):
                break
            minsize, ratios = free_variable_solve(pot, max_bonds)
            if(minsize < Graph.number_of_nodes()):
                while(len(ratios) < max_tiles):
                    ratios.append(0)
                partitions.append(ratios)
                continue
            else:
                return pot, tile_assignments, orientations
    return False, False, False

#Used by all the qvalue code
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

def S2_tiles_qvalue(Graph, gurobi_flag):
    #Solve S1
    pot, tile_assignments, orientation = S1_optimal_pot(Graph, gurobi_flag)

    #Set search space params
    tile_types_minimum = len(pot)
    tile_types_maximum = Graph.number_of_nodes()
    bond_edge_type_minimum = 1
    bond_edge_type_maximum = Graph.number_of_nodes()-1

    #Set our initial values
    tile_types = tile_types_minimum
    bond_edge_types = bond_edge_type_minimum

    #Define the function to loop over all qvals in a specific tile and bond edge type

    print("Checking bond edges: [" + str(bond_edge_type_minimum) + "-" + str(bond_edge_type_maximum) + "], tile types: [" + str(tile_types_minimum) + "-" + str(tile_types_maximum) + "]")

    #Now slowly increase our bounds until we find an answer
    solved = False
    while not solved:
        print("Checking B_2=" + str(bond_edge_types) + ", T_2=" + str(tile_types))
        solved, pot, tile_assignments, orientations = loop_through_qvals(Graph, tile_types, bond_edge_types, gurobi_flag)
        if not solved:
            bond_edge_types = bond_edge_types + 1
            if(bond_edge_types >= tile_types or bond_edge_types > bond_edge_type_maximum):
                bond_edge_types = bond_edge_type_minimum
                tile_types = tile_types + 1
            if(tile_types > tile_types_maximum):
                break
    return pot, tile_assignments, orientations

def S2_bonds_qvalue(Graph, gurobi_flag, tile_types, bond_edge_types):
    tile_types_maximum = Graph.number_of_nodes()
    bond_edge_type_minimum = 1
    for bond_edges in range(int(bond_edge_type_minimum), int(bond_edge_types)):
        for tiles in range(tile_types+1, tile_types_maximum + 1):
            print("Checking B_2=" + str(bond_edges) + ", T_2=" + str(tiles))
            solved, pot, tile_assignments, orientations = loop_through_qvals(Graph, tiles, bond_edges, gurobi_flag)
            if(solved):
                return pot, tile_assignments, orientations
    return False, False, False

def loop_through_qvals_hybrid(Graph: nx.graph, tile_types: int, bond_edge_types: int, gurobi_print: bool, partitions):
    #Set our initial qvalue limits
    minqs = [0 for i in range(tile_types)]
    q_under_equal = -1

    initial = True
    
    #keep checking pots
    while(True):
        #get our pot
        pot, tile_assignments, orientations, qvals, savedata = optimal_pot_s2(Graph, bond_edge_types, tile_types, minqs, q_under_equal, gurobi_print, partitions)
        
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


def S2_tiles_hybrid_qvalue(Graph, gurobi_flag):
    #Solve S1
    pot, tile_assignments, orientation = S1_optimal_pot(Graph, gurobi_flag)

    #Set search space params
    tile_types_minimum = len(pot)
    tile_types_maximum = Graph.number_of_nodes()
    bond_edge_type_minimum = 1
    bond_edge_type_maximum = Graph.number_of_nodes()-1

    #Set our initial values
    tile_types = tile_types_minimum
    bond_edge_types = bond_edge_type_minimum

    #Define the function to loop over all qvals in a specific tile and bond edge type

    print("Checking bond edges: [" + str(bond_edge_type_minimum) + "-" + str(bond_edge_type_maximum) + "], tile types: [" + str(tile_types_minimum) + "-" + str(tile_types_maximum) + "]")

    #Now slowly increase our bounds until we find an answer
    solved = False
    while not solved:
        print("Checking B_2=" + str(bond_edge_types) + ", T_2=" + str(tile_types))
        partitions = get_glist(num_to_partition, tile_types)
        solved, pot, tile_assignments, orientations = loop_through_qvals_hybrid(Graph, tile_types, bond_edge_types, gurobi_flag, partitions)
        if not solved:
            bond_edge_types = bond_edge_types + 1
            if(bond_edge_types >= tile_types or bond_edge_types > bond_edge_type_maximum):
                bond_edge_types = bond_edge_type_minimum
                tile_types = tile_types + 1
            if(tile_types > tile_types_maximum):
                break
    return pot, tile_assignments, orientations

def S2_bonds_hybrid_qvalue(Graph, gurobi_flag, tile_types, bond_edge_types):
    tile_types_maximum = Graph.number_of_nodes()
    bond_edge_type_minimum = 1
    for bond_edges in range(int(bond_edge_type_minimum), int(bond_edge_types)):
        for tiles in range(tile_types+1, tile_types_maximum + 1):
            print("Checking B_2=" + str(bond_edges) + ", T_2=" + str(tiles))
            partitions = get_glist(num_to_partition, tiles)
            solved, pot, tile_assignments, orientations = loop_through_qvals_hybrid(Graph, tiles, bond_edges, gurobi_flag, partitions)
            if(solved):
                return pot, tile_assignments, orientations
    return False, False, False