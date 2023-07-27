# Overview
This repository contains computational tools for flexible tile models in DNA self-assembly.
The current tools provided are the OOPS and the SRPS, detailed more in their respective sections.

This software was developed as part of Summer@ICERM 2023 (https://icerm.brown.edu/summerug/2023/).
This codebase is the implementation of the algorithms outlined in (TODO: link paper once on archive).

## Setup
This code runs on Python 3. You can download the latest release here: https://www.python.org/downloads/
This codebase uses the python packages networkx, numpy, matplotlib, and gurobipy. You can install these packages by running the follow commands.
```
  pip install networkx
  pip install numpy
  pip install matplotlib
  pip install gurobipu
```

This code also uses the gurobi optimizer to solve mixed-integer programming problems. To use this optimizer, you will need to download and activate gurobi.
Gurobi is free for academic use, and you can follow the instructions here to download and activate gurobi: https://www.gurobi.com/features/academic-named-user-license/

# OOPS
Oriented Optimal Pot Solver for DNA Self-Assembly.
Given a target graph, the OOPS aims to compute the optimal (least tiles and bond-edge types) pot capable of realizing the graph, under various sets of restrictions.
Currently OOPS is capable of computing pots for graphs in Scenario 1 and Scenario 2.

## Usage
There is an included command-line interface, which is by default enabled. This will allow you to specify a scenario and a graph to run on by executing the code and following text prompts.

The reccomended way to use OOPS (if you have some coding experience) is to use a code editor to modify the parameters listed in the beginning of OOPS.py. 
This allows you to use all of the tools networkx provides to build target graphs: https://networkx.org/documentation/stable/tutorial.html

# SRPS
Subgraph Realization Problem Solver.
Given a pot, the SRPS aims to compute the smallest size graph the pot can produce. This is useful for helping verify pots in Scenario 2 and Scenario 3.

## Usage
Currently, the only way to use the SRPS is to use a code editor to modify the 'pot' parameter listed at the beginning.
