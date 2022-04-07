# README file for run

VMM3a/SRS Data Analysis Tool: 
vmm-sdat is the analysis software for VMM3a data, recorded with the SRS or the ESS readout as PCAP files. From the PCAP file, a root tree with the hits and clusters is created. To create the pcapng file, the user has to directly get the UDP frames from the network card of the DAQ computer. This can be done with tcdump or Wireshark. The captured data has to be saved as pcapng file.

## Running the convertFile utility
The convertFile utility reads in the pcap file, and produces a root tree with the hits and clusters. The utility can be run directly from the command line, but since there are many parameters to be set, it is easier to use a script. In the analysis folder one finds the runConversion.py script. The script illustrates, how to run convertFile with the example.pcapng data file. The script shows how to store hits, plane clusters and detector clusters. The meaning of the command line parameters is explained in the main README file. 

## Geometry files
Instead of defining the geometry with the parameters -vmm, -axis, -map one can also use a JSON geometry file. The file example_geometry.json in the run directory illustrates the format. In the script runConversion_geometry.py shows how to use the geometry file. In case the detector readout is pads and not strips, the file example_geometry_pads.json illustrates how to define such a detector. The script runConversion_geometry_pad.py executes the pad example. The example geometry_LET.json shows a geometry implementation of a MultiGrid detector that consists of wires and grids. For the grids only some channels on one VMM on the hybrid are used, whereas for the wires two VMMs on one hybrid are used. The channels that are not connected on the VMM3a are marked with a -1. The data format of the example_LET.pcapng file is the ESS data format, whereas the other pcapng files are in the SRS data files.

## Calibration files
The VMM3a slow control tool https://gitlab.cern.ch/mguth/VMM-software-RD51 can generate calibration files that modify the ADC values and/or the time measurements of each channel. If a reliable calibration has been found, the JSON calibration file can be loaded into vmm-sdat to correct the ADC values. In the python script runConversion_calib.py exactly this is done, the example_calib.json file is applied to the data.

