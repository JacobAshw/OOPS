import math
import random
from tkinter import E
import networkx as nx
import matplotlib.pyplot as plt
import itertools
from math import factorial, gcd
from functools import reduce

# This file contains methods used by OOPS.py
# Many functions are isolated here for decluttering of the main file
# See the documentation in OOPS.py for further information

# Make our list of half edges. Lowercase letters represent non-hatted half edges, and uppercase letters are the hatted versions.
# It is highly unlikely this code will ever solve a graph with more than 26 bond edge types, so alphabetical limitations are ignored
def get_half_edge_labels() -> tuple[list[chr],list[chr]]:
    return [chr(i) for i in range(97, 123)], [chr(i) for i in range(65, 91)]

def hat(bet):
    hel = get_half_edge_labels()
    if bet in hel[0]:
        return hel[1][hel[0].index(bet)]
    if bet in hel[1]:
        return hel[0][hel[1].index(bet)]

def find_gcd(list):
    x = reduce(gcd, list)
    return x

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

    print("Number: ", n)

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

    return [g for g in list(set(glst)) if find_gcd(g)==1]

# Get a list of n different colors, mainly used by OOPS to color the different tile types in the graph
def get_color_list(num_colors: int) -> list:
    #If num_colors is 20 or less, we do our best to get a bunch of visually distinct colors
    #Color pallet used from this thread: https://sashamaps.net/docs/resources/20-colors/
    colors = ['#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000']
    if(num_colors < len(colors)):
        return [colors[i] for i in range(num_colors)]
    #If there are too many colors, just generate random colors
    #Inspired by the approaches here: https://www.geeksforgeeks.org/create-random-hex-color-code-using-python/
    return [hex(random.randrange(0, 2**24)) for i in range(num_colors)]

def display_orientation(Graph, pot, tile_assignments, orientations, pos, old_labels):
    D = nx.DiGraph()
    #Add each node
    for vertex in range(Graph.number_of_nodes()):
        D.add_node(vertex)

    #Add each (directed) edge
    edge_lables = {}
    half_edges, half_edges_hat = get_half_edge_labels()
    for index, orientation in enumerate(orientations):
        for vertex1 in range(Graph.number_of_nodes()):
            for vertex2 in range(Graph.number_of_nodes()):
                if(orientation[vertex1][vertex2] == 1):
                    D.add_edge(vertex1, vertex2)
                    edge_lables.update({(vertex1, vertex2) : half_edges[index]})
    #Color our nodes
    color_list = get_color_list(len(tile_assignments))
    color_map = {}
    for index,assignment in enumerate(tile_assignments):
        for tile in tile_assignments.get(assignment):
            color_map.update({tile : color_list[index]})
    colors = [color_map.get(vertex) for vertex in range(Graph.number_of_nodes())]
    #Display the graph
    nx.draw_networkx(D, pos, labels=old_labels, node_color=colors, with_labels=True)
    nx.draw_networkx_edge_labels(D, pos, edge_lables)
    plt.show()


# This code runs the command line interface. It is largely uncommented, as the print statements provide sufficient information
def CLI_Setup(version: str):
    print("OOPS solver version " + version)
    print("Using command line interface")
    print("If you plan to use OOPS frequently, consider using a code editor and networkx for more efficient graph building")
    print("-----------------------------------------------------------------------")
    print("Please enter which scenario to find an optimal pot for ('1' or '2')")
    while True:
        try:
            scenario = int(input())
            if(scenario == 1 or scenario == 2):
                break
            if(scenario == 3):
                print("Unfortunarely, scenario 3 is not supported. Please enter '1' or '2'")
            else:
                print("Please input either '1' or '2'")
        except ValueError:
            print("Input not recognized as an integer. Please input either '1' or '2'")
    print("Scenario " + str(scenario) + " selected")
    print("Please enter your graph as a list of parenthesized edges")
    print("Only simple graphs are supported, so do not input self-loops or mutli-edges")
    print("Example: The three cycle would be '(0,1) (1,2) (2,0)'")
    while(True):
        # try:
        #Parse our input graph
        edges = input()
        edges = edges.replace(" ", "").replace(")", "")
        edges = edges.split("(")
        #Skip the first entry, since it should be empty
        edges = edges[1:]
        for i in range(len(edges)):
            edges[i] = edges[i].split(",")
        num_edges = len(edges)
        for edge_index in range(num_edges):
            if(len(edges[edge_index]) != 2):
                print("Edge " + str(edge_index) + " does not have 2 ends, please try again")
                break
            edge = (edges[edge_index][0], edges[edge_index][1])
            if(edge[0]==edge[1]):
                print("Self-loops are not allowed, please try again")
                continue
            edges[edge_index] = edge
        G = nx.Graph()
        G.add_edges_from(edges)
        break
        # except:
        #     print("Edge parsing failed. Please try again.")
    return scenario, G