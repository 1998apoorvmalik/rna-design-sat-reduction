# RNA Design to SAT Reduction
This repository contains our project work for reducing the RNA design problem to SAT. It uses pysat (SAT solver) to find the boolean variables assignment which corresponds to the solution of the RNA design problem, if it exists. The solution to the RNA design problem is a sequence of nucleotides which uniqely folds into the input secondary structure. Any alternative structure that the sequence can fold into must contain strictly lesser number of base pairs than the input secondary structure.

The file `enumeration.py` contains the code for solving the RNA design problem using naive enumeration method. \
The file `sat.py` contains the code for reducing the RNA design problem to SAT and then solving it.

# Web Server
We have hosted our RNA Design SAT Solver at [RNA Design SAT Solver](http://rnadesignsat.herokuapp.com/). You can use this web server to design your own RNA sequences for a target secondary structure.