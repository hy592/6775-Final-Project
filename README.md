# MIMO_Ising

This repository contains code for using an Ising formulation to decode Multiple-Input-Multiple-Output decoding problems. The notebook "DI_MIMO.ipynb" contains code to generate instances, solve them and evaluate the bit error rates, using the Delta-Ising MIMO formulation described in [A Finite-Range Search Formulation of Maximum Likelihood MIMO Detection for Coherent Ising Machines (Abhishek et. al.)](https://arxiv.org/pdf/2205.05020.pdf). 

The instances generated are stored in the folder "di_mimo/instances". The ideal solutions are in "di_mimo/ideal_solutions", and the solutions from the Ising solver are stored in "di_mimo/solved_solutions"

## Trouble Shooting
```bash
# Set Up The Environment to Run Vitis
source /tools/Xilinx/Vitis/2022.2/settings64.sh
vitis_hls # open Vitis HLS GUI
```