import networkx as nx
import time
import csv
import tempfile
import configparser
import sys
import math
import os
import matplotlib.pyplot as plt

configwritedir = "pot_output/6cbenchmark/hc1"
graphdir = "noniso_graphs"
graphfilename = "graph6c.g6"

graphfile = open(graphdir + "/" + graphfilename, 'r')
graphnames = [graphname[:len(graphname)-1] for graphname in graphfile.readlines()]
graphfile.close()
# graphs = []
# for graphname in graphnames:
#     with tempfile.NamedTemporaryFile() as f:
#         sstr = ">>graph6<<"+graphname+"\n"
#         _ = f.write(sstr.encode())
#         _ = f.seek(0)
#         G = nx.read_graph6(f)
#     graphs.append(G)

default_config = "config.ini"
config = configparser.ConfigParser()
config.read(default_config)

print(config.sections())

# config = config_reader['Instance']

# print(config['CLI_Enabled'])
if(not os.path.exists(configwritedir + "/")):
    os.mkdir(configwritedir + "/")

if(not os.path.exists(configwritedir + "/" + "output" + "/")):
    os.mkdir(configwritedir + "/" + "output" + "/")

i = 0
for graphname in graphnames:
    config['Instance']['graph'] = graphname
    config['Instance']['display_graph'] = "False"
    config['Instance']['graphG6'] = "True"
    config['Instance']['outputfile'] = configwritedir + "/" + "output" + "/" + str(i) + ".csv"
    config['Instance']['method'] = "6"
    with open(configwritedir + "/" + str(i) + ".ini", 'w') as cfgfile:
        config.write(cfgfile)
    i = i + 1