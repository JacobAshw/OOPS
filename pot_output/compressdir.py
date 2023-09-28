from bz2 import compress
from re import L
import networkx as nx
import time
import csv
import tempfile
import sys
import math
import os
import matplotlib.pyplot as plt

compressdir = "pot_output/graph7c.g6_c1"
finish = "pot_output/7c.csv"

if not (os.path.exists(compressdir)):
    print("Not Found")
    exit()

files = os.listdir(compressdir)
header = ["Graph","T_1 Value","T_1 Pot","T_1 Ratios","T_1 Time","T_1 Orientation","T_2 Value","T_2 Pot","T_2 Pot Ratios","T_2 Orientation","B_2 Value","B_2 Pot","B_2 Pot Ratios","B_2 Orientation","SamePot","TileTime","BondTime"]


with open(finish, 'w', newline='\n') as wfile:
    writer = csv.writer(wfile)
    writer.writerow(header)
    for file in files:
        with open(compressdir + "/" + file) as rfile:
            reader = csv.reader(rfile)
            skip = 1
            for line in reader:
                if(skip > 0):
                    skip = skip - 1
                    continue
                writer.writerow(line)