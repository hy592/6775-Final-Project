# Sparse fpga make / MaxSAT FPGA
This directory contains the FPGA implementation that
1. initializes `num_trials` copies of the state vectors `s` and `a` with random initial conditions.
2. evolves all `num_trials` copies of the state vectors together, each doing sparse matrix-vector multiplication

This is a version of implementation of AnalogSat running on the FPGA that solves MaxSAT problems. It is an incomplete solver, and it only solves unweighted MaxSAT problems.

The only difference between this version and what's in `../sparse_parallel_threads` is the compilation flow. The other version in `../sparse_parallel_threads` uses the Vitis GUI to compile and run, which become not feasible after Vitis 2020.2 update because of the buggy new GUI.

## Build and Run

This version uses Makefile to build and run emulations.

To build:
`$ make all TARGET=sw_emu DEVICE=$AWS_PLATFORM`

For emulation builds, make sure it is a small build, otherwise it never finishes running.

`TARGET=sw_emu/hw_emu/hw`

To run:
`make run TARGET=sw_emu DEVICE=$AWS_PLATFORM`

The commandline arguments are specified in the Makefile. You need to modify the makefile if you are to change the commandline arguments.

There is no need to rebuild for each run.

## Adapting the Makefile for other version of code

### Case 1: still working on the ksat project
Please search for `FIXME` flags in the Makefile.

### Case 2: kernel with another name
In addition to case 1, you'll need to replace `ksat` in the makefile into whatever name you have for your kernel.

### Case 3: working in another repo
In addition to case 1+2, you'll need to copy over the `common/` directory in the top directory in this repo.

### Case 4: run with different command line argument interface
You'll need to make changes to what's under `run` or `test` in the Makefile

## File organization
#### test_data
This folder contains a couple of MaxSAT problems encoded in the WCNF format for sanity checks during development. 

#### .clang-format
CPP formatter configuration file to ensure the correct format across different development platforms.

#### host.cpp
Host code for MaxSAT FPGA. You need this file for compiling the host code. The host code is responsible for taking file I/O, orchestrate the FPGA(s), 
and solution checking. More information can be found in Owen's SP2021 report on slack.

The host code has the command line interface set up as follows. The commandline arguments are specified in the Makefile.
```
./ksat_host <problem_path> <number of variables> <number of clauses> <K> <number of steps> <A> <B> <run steps> <early stop> <timeout> <maxsat>

<problem_path> - path of a WCNF file
<number of variables> - the number of variables in the problem
<number of clauses> - the number of clauses in the problem
<K> - the maximum number of literals appearing in this problem
<number of steps> - number of steps before re-initialization
<A>, <B> - hyperparameters. Use A = 1.8 and B = 2.0 if don't know
<run steps> - number of steps between checking solutions. This number is usually << <number of steps>. <run steps> will be set equal to <number of steps> if pass in 0
<early stop> - 1 for True and 0 for False. Do we want to stop early if a solution is found?
<timeout> - in unit of microseconds. -1 if want to run indefinitely
<maxsat> - is this a maxsat problem. 1 for True and 0 for False
```

Note that maxsat must be set to 1 for now. You can also tune `MAX_VAR`, `MAX_CLAUSES`, `K_max`, `num_trials` before compilation to fit different needs. NOTE that you'll need small max_var, max_clauses, num_trials for running emulations...

#### ksat.cpp
FPGA kernel code for MaxSAT FPGA. Tunable parameters are all declared in the beginning of the file. Note that these tunable parameters need to agree between `host.cpp` and `ksat.cpp` otherwise the behavior is undefined. 

This cpp code will be transcribed into RTL language by Vitis HLS. Please read Owen's SP2021 report for more information.

#### Makefile
Build and run flow with Make

#### utils.mk
Not sure what this is doing but Makefile flow doesn't seem to work without this. It is specifying some build configs but I can't find any documentation around this file from Xilinx.

#### xrt.ini
Not sure what it does but it seems to help XRT.
