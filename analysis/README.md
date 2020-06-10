# README file for analysis

VMM3a/SRS Data Analysis Tool: 
vmm-sdate is the analysis software for VMM3a data, recorded with the SRS as PCAP or HDF5 files (GdGEM pipeline of the EFU). From the PCAP or HDF5 file, a root tree with the hits and clusters is created. This root tree can be analyzed in different ways. For three ways, wich are explained here in more detail, example scripts can be found in the analysis folder 

## Analysis with Python
Python is becoming increasingly popular. The example tree-reader.py shows how with the help of uproot, a root file with the root tree inside is read. The data is then plotted with matplotlib. 
To select only specific entries of the tree, pandas can be used to apply a cut to the data.This is shown in the example tree-cutter.py.
For questions, please ask Lucian.Scharenberg@cern.ch

## Analysis with TSelector root script
The root framework provides the TSelector mechanism. Details can be found here:
https://root.cern.ch/developing-tselector
In short, the following steps are needed.
1. Open root from the command line
2. Open the root tree from the root prompt: TFile *f = TFile::Open("example.root")
3. Get the event tree from the file: TTree *t = (TTree *) f->Get("events")
4. Create the selector script: events->MakeSelector("VMM3a_analysis")
5. The last command creates two files, VMM3a_analysis.C, VMM3a_analysis.h
6. Add the plot that you want to create to the script. In the example provided in the analysis folder, this has already been done.
7. Execute the script from the root prompt: events->Process("VMM3a_analysis.C")

## Analysis with compiled root script
In the CMakeList.txt in the main directory, the example executable accessTree is created.
add_executable(accessTree ${CMAKE_CURRENT_SOURCE_DIR}/analysis/accessTree.cpp)
This example can be edited to plot different data from the tree.