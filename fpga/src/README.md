# Final Project README

## Run C Simulation
- ecelinux

To run the software on ecelinux, first go to the `6775-Final-Project/fpga/src` directory, then run `make test`

The simulation result will be write into `test_output_low.txt` with paired `bestEnergy` and `spins`.

To check the functionality of dut and host program on software, run `make test_host`. `test_output_host.txt` should give exactly the same `bestEnergy` and `spins` as `test_output_low.txt` gives.

```bash
cd 6775-Final-Project/fpga/src # main path

source /classes/ece6775/setup-ece6775.sh # setup env

make test # do software simulation for functionality

make test_host # do software simulation with dut/host program
```
## Run C Synthesis
```bash
cd 6775-Final-Project/fpga/src # main path

source /classes/ece6775/setup-ece6775.sh # setup env

vivado_hls -f run_base.tcl # C Synthesis
```
<!-- You will implement and evaluate the performance for three designs:

- A baseline digitrec design that does not use any HLS optimization directives (`vivado_hls -f run_base.tcl`).
- An unrolled digitrec design which is similar to what you did in Lab 2 where unrolling and array partitioning are applied (`vivado_hls -f run_unroll.tcl`).
- A pipelined digitrec design which applies loop pipelining in addition to the previous optimizations (`vivado_hls -f run_pipeline.tcl`). -->

## Generating the Bitstream

1. The first step is to generate verilog from HLS
```bash
# In ecelinux
cd 6775-Final-Project/fpga/src # main path

source /classes/ece6775/setup-ece6775.sh # setup env

make test # builds and runs the csim

# generate Verilog
vivado_hls -f run_base.tcl
```

You can examine the directory cordic.prj/solution1/syn/verilog to check that the Verilog files have been generated.

2. The next step is to implement the hardware circuit described by the Verilog on an the FPGA.

```bash
source run_bitstream.sh # generate FPGA bitstream
```

## Programming the FPGA

Before you attempt to login, please first
check the occupancy of the available boards using the following link:
https://www.csl.cornell.edu/courses/ece6775/zedboard.html

```bash
zhang-zedboard-xx.ece.cornell.edu
```
where xx can be 01, 02, ..., 11.

All students have an account on each Zedboard with
their NetID as the username and their case-sensitive last name as their initial password.

```bash
# In ecelinux,
# Copy the bitstream from ecelinux to Zedboard
scp xillydemo.bit hy592@zhang-zedboard-09.ece.cornell.edu:~

# Login to the Zedboard
ssh hy592@zhang-zedboard-09.ece.cornell.edu

# In Zedboard
# Mount SD card. Try rebooting if there is a device busy error
mount /mnt/sd

# Copy the bitstream file to the SD card and reboot the Zedboard
cp xillydemo.bit /mnt/sd

sudo reboot # restart Zedboard
```
The Zedboard will restart. Wait about 30 seconds and login again — the FPGA will be fully programmed after the reboot

## Run FPGA Experiment
- Zedboard

To run the FPGA experiment on a Zedboard, first copy `6775-Final-Project` files to a Zedboard, 

```bash
# in ecelinux
zip -r final_project.zip 6775-Final-Project  # zip the file

scp final_project.zip <user>@zhang-zedboard-xx.ece.cornell.edu:~ # copy zip file to zedboard

ssh <user>@zhang-zedboard-xx.ece.cornell.edu # log in to a Zedboard
```

then go to the directory `6775-Final-Project/fpga/src`, finally run `make fpga`.

Notice, need to rename `Makefile_fpga` to `Makefile` so run on Zedboard environment.

```bash
# in zedboard
unzip final_project.zip # unzip the archive

cd 6775-Final-Project/fpga/src # main path

mv Makefile_fpga Makefile # for Zedboard Environment

make fpga # builds and runs FPGA experiment
```


## Testbench

To run the testbench, run ```make test``` in ```fpga/src/```. Currently, we have two implementations: ```AHC_low``` and ```AHC_new```. Modify 
- ```make test``` in ```fpga/src/Makefile```'s
- ```#include```s in ```fpga/src/AHC_test.cpp```'s

to change versions.

### Makefile

```cmake
test: AHC_test
	g++ ${CFLAGS} -o AHC_test AHC_test.cpp AHC_<low OR new>.cpp
	./AHC_test > test_output_<low OR new>.txt

AHC_test: AHC_test.cpp
	g++ ${CFLAGS} -o AHC_test AHC_test.cpp AHC_<low OR new>.cpp  
```