# README file for analysis

VMM3a/SRS Data Analysis Tool: 
vmm-sdat is the analysis software for VMM3a data, recorded with the SRS as PCAP or HDF5 files (GdGEM pipeline of the EFU). From the PCAP or HDF5 file, a root tree with the hits and clusters is created. 


## Further analysis with the root tree [`examples_read_data`]

The produced root tree can be further analyzed in different ways. Some general examples on how to read the data from the root tree are shown in `examples_read_data`.


### Analysis with Python
Python is becoming increasingly popular. The example tree-reader.py shows how with the help of uproot, a root file with the root tree inside is read. The data is then plotted with matplotlib. 
To select only specific entries of the tree, pandas can be used to apply a cut to the data.This is shown in the example tree-cutter.py.
For questions, please ask Lucian.Scharenberg@cern.ch


### Analysis with TSelector root script
The root framework provides the TSelector mechanism. Details can be found here:
https://root.cern.ch/developing-tselector
In short, the following steps are needed.
1. Open root from the command line
2. Open the root tree from the root prompt: `TFile f("example.root")`
3. Get the event tree from the file: `TTree *t = f.Get<TTree>("clusters_detector")`
4. Create the selector script: `t->MakeSelector("VMM3a_analysis")`
5. The last command creates two files, VMM3a_analysis.C, VMM3a_analysis.h
6. Add the plot that you want to create to the script. ATTENTION: In the example provided in the analysis folder, this has already been done!! 

So you need just to call the analysis script from your tree:
1. At the command line: `root example.root`
2. In root: `clusters_detector->Process("VMM3a_analysis.C++")`


### Analysis with compiled root script
In the CMakeList.txt in the main directory, the example executable accessTree is created.

`add_executable(accessTree ${CMAKE_CURRENT_SOURCE_DIR}/analysis/accessTree.cpp)`

This example can be edited to analyze and plot different data from the tree. In principle the same things can be done with the compiled root script and the TSelector root script. The compiled executable will just be faster.


### Analysis with CINT root script
This example creates a multitude of plots for the LET Multigrid detector that consists of grids and wires. It can be easily adapted to detector with an x/y readout by replacing wires with x and strips with y. To execute the script, the following strips are necessary:
1. start root by typing in a terminal: `root`
2. in root, type `.x plots_LET.C++("example_LET",".",0)`


### Normalization of histograms
The script normalizeHisto can be used to normalize histograms, or carry out a 
This example creates a multitude of plots for the LET Multigrid detector that consists of grids and wires. It can be easily adapted to detector with an x/y readout by replacing wires with x and strips with y. To execute the script, the following strips are necessary:
1. start root by typing in a terminal: `root`
2. in root, type `.x plots_LET.C++("example_LET",".",0)`


## Specific applications

In addition to the analysis methods that are introduced in general in `examples_read_data`, also specific analysis programs for one specific application are shown here.


### RD51 Beam Telescope [`application_rd51_beam_telescope`]

For the RD51 Beam Telescope (see [https://doi.org/10.1088/1748-0221/18/05/C05017](https://doi.org/10.1088/1748-0221/18/05/C05017)), which is read out with VMM3a/SRS, the data from the vmm-sdat root tree are taken and analysed by [anamicom](https://gitlab.physik.uni-muenchen.de/Jonathan.Bortfeldt/anamicom).
Anamicom takes the root tree from vmm-sdat and performs the event building and the track reconstruction for the charged particle tracks recorded with the beam telescope.
It saves the information in another root tree (not the one from vmm-sdat).
The information from the anamicom root tree are taken in the provided macro and used to calculate the spatial resolution and efficiency of a detector under test, following the method described [here](https://cds.cern.ch/record/2860765).
