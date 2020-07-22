# README file for run

VMM3a/SRS Data Analysis Tool: 
vmm-sdat is the analysis software for VMM3a data, recorded with the SRS as PCAP or HDF5 files (GdGEM pipeline of the EFU). From the PCAP or HDF5 file, a root tree with the hits and clusters is created. Normally the user would run the ESS DAQ (Event Formation Unit), which optionally can create an HDF5 file with the detector hits. As alternative, if very high rates are required, the user can also directly get the UDP frames that the SRS FEC sends. This can be done with tcdump or Wireshark. The captured data has to be saved as pcapng file.

## Running the convertFile utility
The convertFile utility reads in the hdf5 or pcap file, and produces a root tree with the hits and clusters. The utility can be run directly from the command line, but since there are many parameters to be set, it is easier to use a script. In the analysis folder one finds the runConvertFileH5.py and runConvertFilePcapng.py scripts. These scripts illustrate, how to run convertFile with the example.h5 or example.pcapng data file. The script for the hdf5 file stores only detector clusters, whereas the pcapng script shows how to store hits, plane clusters and detector clusters. The meaning of the command line parameters is explained in the main README file. 

## Calibration files
The VMM3a slow control tool https://gitlab.cern.ch/mguth/VMM-software-RD51 can generate calibration files that modify the ADC values and/or the time measurements of each channel. If a reliable calibration has been found, the JSON calibration file can be loaded into the ESS DAQ (event formation unit), and the data is automatically corrected. During the testing of the detector system, it is sometimes more convenient to let the DAQ store the raw data to file, and then use vmm-sdat to apply the calibration. In the python script runConvertFileH5_calib.py exactly this is done, the example_calib.json file is applied to the data.

