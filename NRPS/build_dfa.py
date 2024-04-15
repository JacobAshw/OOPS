from email.utils import encode_rfc2231
from re import L
import networkx as nx
import time
import sys
import csv
import tempfile
import math
import configparser
from matplotlib import pyplot as plt
import matplotlib.pyplot as plt


pot = ['aaB', 'aab', 'AB', 'Ab']
tile_assign = {'aaB': [0, 1], 'aab': [], 'AB': [2], 'Ab': [3, 4, 5]}
edgelist = ['0 2 {}', '0 4 {}', '0 5 {}', '1 3 {}', '1 5 {}', '1 4 {}', '2 3 {}']
bond_edge_pairs = [('a', 'A'), ('b', 'B')]
orientation = {
    'a': [(0,2),(0,4),(1,5),(1,3)],
    'b': [(5,0),(4,1),(3,2)]
}
# ! ---------------------------------------------------------------------------------------------------
# G = nx.cycle_graph(4)
# G.add_node(4)
# G.add_edge(1, 4)
# G.add_edge(2, 4)

# G=nx.lollipop_graph(3,3)
# G.add_edge(5,3)

G = nx.parse_edgelist(edgelist)


# G = nx.dodecahedral_graph()
# ! ---------------------------------------------------------------------------------------------------
nx.draw_circular(G)
plt.show() 

# old_labels = {}
# new_labels = {}
# for index, node in enumerate(G.nodes()):
#     new_labels.update({node: index})
#     old_labels.update({index: node})
# G = nx.relabel_nodes(G, new_labels)

new_labels = {}
for index, node in enumerate(G.nodes()):
    new_labels.update({node: int(node)})
G = nx.relabel_nodes(G, new_labels)

class GraphPerm:
    def __init__(self, Graph : nx.MultiGraph):
        self.Graph = Graph

class DFA:
    def __init__(self, StartGraph : nx.MultiGraph):
        self.DFAGraph = nx.MultiDiGraph()
        self.nodenum = 0
        self.nodedict = {}

        self.addToDFA(StartGraph)
        self.level0 = [StartGraph]
        self.level1 = self.removedEdgeSet(self.level0)
        self.level2 = self.removedEdgeSet(self.level1)
        self.level3 = []
        self.level3R = nx.Graph()

        

        #Add nodes
        for graph in self.level1:
            self.addToDFA(graph)
        for graph in self.level2:
            self.addToDFA(graph)

        self.addToDFA(self.level3R)
        

        #Add out edges
        for edge in StartGraph.edges():
            newG = StartGraph.copy()
            newG.remove_edge(edge[0], edge[1])

            nextG = self.getIsoInLevel(1, newG)
            move = edge
            permutation = nx.vf2pp_isomorphism(newG, nextG)
            self.addEdgeDFA(StartGraph, nextG, move, permutation, 'r')

        for level1g in self.level1:
            for edge in level1g.edges():
                newG = level1g.copy()
                newG.remove_edge(edge[0], edge[1])

                nextG = self.getIsoInLevel(2, newG)

                move = edge
                permutation = nx.vf2pp_isomorphism(newG, nextG)
                self.addEdgeDFA(level1g, nextG, move, permutation, 'r')
        
        #Add in edges
        for level2g in self.level2:
            for i in range(level2g.number_of_nodes()):
                for j in range(level2g.number_of_nodes()):
                    newG = level2g.copy()
                    newG.add_edge(i, j)

                    nextG = self.getIsoInLevel(1, newG)
                    if(nextG == None):
                        nextG = self.level3R
                    
                    move = (i, j)
                    permutation = nx.vf2pp_isomorphism(newG, nextG)
                    self.addEdgeDFA(level2g, nextG, move, permutation, 'b')

        for level1g in self.level1:
            for i in range(level1g.number_of_nodes()):
                for j in range(level1g.number_of_nodes()):
                    newG = level1g.copy()
                    # print("-")
                    # print(newG.number_of_nodes())
                    newG.add_edge(i, j)
                    # print(newG.number_of_nodes())

                    nextG = self.getIsoInLevel(0, newG)
                    if(nextG == None):
                        nextG = self.level3R
                    
                    move = (i, j)
                    permutation = nx.vf2pp_isomorphism(newG, nextG)
                    self.addEdgeDFA(level1g, nextG, move, permutation, 'b')


    def get1Perms(self):
        # R1
        # print(l1outperms)
        # print(self.level0[0])
        for l1oute in self.DFAGraph.out_edges(0, data=True):
            if(l1oute[2]["color"] != "r"):
                continue
            perm1oute = l1oute[2]["permutation"]
            nextnode0 = l1oute[1]
            # print(perm1oute)
        # R2
            for l2oute in self.DFAGraph.out_edges(nextnode0, data=True):
                if(l2oute[2]["color"] != "r"):
                    continue
                perm2oute = l2oute[2]["permutation"]
                nextnode1 = l2oute[1]
                # print(perm2oute)
        # A1
                for l2ine in self.DFAGraph.out_edges(nextnode1, data=True):
                    if(l2ine[2]["color"] != "b"):
                        continue
                    perm2ine = l2ine[2]["permutation"]
                    nextnode2 = l2ine[1]
        # A2
                    for l1ine in self.DFAGraph.out_edges(nextnode2, data=True):
                        if(l1ine[2]["color"] != "b"):
                            continue
                        perm1ine = l1ine[2]["permutation"]
                        nextnode3 = l1ine[1]
                        # print(nextnode3)

    def movesFromOrient(self, orientation):
        moves = []
        for bond_edge_pair in bond_edge_pairs:
            edges = orientation.get(bond_edge_pair[0])
            for i in range(len(edges)):
                for j in range(i+1, len(edges)):
                    e1 = edges[i]
                    e2 = edges[j]
                    move = []
                    move.append(e1)
                    move.append(e2)
                    move.append((e1[0],e2[1]))
                    move.append((e2[0],e1[1]))
                    moves.append(move)
        return moves
        
    def assignOrientDigraph(self, assign, orient):
        G = nx.MultiDiGraph()

        currentnum = 1
        for key in assign.keys():
            for vertex in assign.get(key):
                G.add_vertex(vertex)
                for num in range(currentnum):
                    G.add_edge((vertex, vertex))
            currentnum = currentnum + 1
        
        currentnum = 1
        for bond_edge_pair in bond_edge_pairs:
            edges = orient.get(bond_edge_pair[0])
            for edge in edges:
                for num in range(currentnum):
                    G.add_edge(edge[0], edge[1])
            currentnum = currentnum + 1
        
        return G

    def assignOrientIsomorphic(self, a1, o1, a2, o2):
        G1 = self.assignOrientDigramp(a1, o1)
        G2 = self.assignOrientDigramp(a2, o2)
        return nx.is_isomorphic(G1, G2)

    def permOnMove(self, perm, move):
        if(perm == None):
            return None
        newMove = []
        for edge in move:
            newMove.append((perm.get(edge[0]), perm.get(edge[1])))
        # print(newMove)
        return newMove

    def finalMove(self, startmove):
        # print("MOVE", startmove)
        move = startmove
        perms = []
        for l1oute in self.DFAGraph.out_edges(0, data=True):
            # print(l1oute[2]["move"], " color ", l1oute[2]["color"])
            if(l1oute[2]["color"] != "r"):
                continue
            if(l1oute[2]["move"] != move[0] and l1oute[2]["move"] != (move[0][1],move[0][0])):
                continue
            perms.append(l1oute[2]["permutation"])
            move = self.permOnMove(l1oute[2]["permutation"], move)
            nextnode0 = l1oute[1]
            # print('cp1')
            break
        # print(nextnode0)
        # R2
        for l2oute in self.DFAGraph.out_edges(nextnode0, data=True):
            # print(move[1])
            # print(l2oute[2]["move"])
            # print(l2oute[2]["color"])
            if(l2oute[2]["color"] != "r"):
                continue
            if(l2oute[2]["move"] != move[1] and l2oute[2]["move"] != (move[1][1],move[1][0])):
                # print(l2oute[2]["move"], "!=", move[1])
                continue
            perms.append(l2oute[2]["permutation"])
            move = self.permOnMove(l2oute[2]["permutation"], move)
            nextnode1 = l2oute[1]
            break
        # A1
        for l2ine in self.DFAGraph.out_edges(nextnode1, data=True):
            if(l2ine[2]["color"] != "b"):
                continue
            if(l2ine[2]["move"] != move[2] and l2ine[2]["move"] != (move[2][1],move[2][0])):
                continue
            nextnode2 = l2ine[1]
            # print(l2ine[1])
            perms.append(l2ine[2]["permutation"])
            move = self.permOnMove(l2ine[2]["permutation"], move)
            # print('cp3')
            break
        # A2
        for l1ine in self.DFAGraph.out_edges(nextnode2, data=True):
            if(l1ine[2]["color"] != "b"):
                continue
            if(l1ine[2]["move"] != move[3] and l1ine[2]["move"] != (move[3][1],move[3][0])):
                continue
            nextnode3 = l1ine[1]
            perms.append(l1ine[2]["permutation"])
            move = self.permOnMove(l1ine[2]["permutation"], move)
            break

        finalperm = {}

        if None in perms:
            print("Nonisomorphic Graph Can Be Made")
            return None

        for i in perms[0].keys():
            ap1 = perms[0].get(i)
            ap2 = perms[1].get(ap1)
            ap3 = perms[2].get(ap2)
            ap4 = perms[3].get(ap3)
            finalperm.update({i : ap4})
        return finalperm

        

    def allNewOrientations(self, assign, orient):
        newOrientations = []

        print("#############################")
        print(orient)
        # print(self.movesFromOrient(orient))
        moves = self.movesFromOrient(orient)

        for move in moves:
            print(self.finalMove(move))


        print("#############################")

    def addToDFA(self, Graph):
        self.DFAGraph.add_node(self.nodenum, graph=Graph)
        self.nodedict.update({Graph : self.nodenum})
        self.nodenum = self.nodenum + 1

    def addEdgeDFA(self, Graph1, Graph2, Move, Permutation, Color):
        G1 = self.nodedict.get(Graph1)
        G2 = self.nodedict.get(Graph2)
        self.DFAGraph.add_edge(G1, G2, move=Move, permutation=Permutation, color=Color)

    def getIsoInLevel(self, levelint, Graph):
        level = []
        if levelint == 0:
            level = self.level0
        if levelint == 1:
            level = self.level1
        if levelint == 2:
            level = self.level2

        for other in level:
            if nx.is_isomorphic(other, Graph):
                return other
        return None

    def removedEdgeSet(self, Graphs):
        oneEdgeRemoves = []
        for Graph in Graphs:
            for edge in Graph.edges():
                newG = Graph.copy()
                newG.remove_edge(edge[0], edge[1])
                seen = False
                for other in oneEdgeRemoves:
                    if nx.is_isomorphic(newG, other):
                        seen = True
                        break
                if not seen:
                    oneEdgeRemoves.append(newG)
        return oneEdgeRemoves


# for edge in G.edges():
#     print(edge)
dfa = DFA(nx.MultiGraph(G))
# print("$")
# for edge in nx.MultiGraph(G).edges():
#     print(edge)

print('------------------------')
print('Start ')
print('TileAssign: ' + str(tile_assign))
print('Nodes: ' + str(list(G.nodes())))

# dfa.get1Perms(pot)
dfa.allNewOrientations(tile_assign, orientation)

# nx.draw_networkx(dfa.DFAGraph, with_labels=True)
colors = nx.get_edge_attributes(dfa.DFAGraph,'color').values()
nx.draw_circular(dfa.DFAGraph, with_labels=True,edge_color=colors)
plt.show()