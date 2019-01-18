
# vmm-hdf5-to-root

Analysis software for VMM3a HDF5 files written by the EFU. From the HDF5 file, a root tree with the hits and clusters is created.
For more information about ROOT see [here](https://root.cern.ch/)

## Getting Started

### Prerequisites
- Boost system library [find_package( Boost REQUIRED COMPONENTS system)]
- HDF5 1.10 or newer [find_package(HDF5 1.10 REQUIRED)]
- h5cpp from https://github.com/ess-dmsc/h5cpp [find_package(h5cpp REQUIRED)]
- root6 from https://root.cern.ch/downloading-root [find_package(ROOT REQUIRED)]



## Installing

To build the program, start in the program directory vmm-hdf5-to-root
```
mkdir build
cd build
cmake .. or (to enable debug info) cmake -DENABLE_DTRACE=ON ..
```
## Built With
* [CMAKE](https://cmake.org/) - Cross platform makefile generation

## Deployment
How to run the program:
For help concerning the possible parameters, type:
```
./convert 
```

Complete command line:
```
./convert -f data.h5 -x 1,2,14,1,2,15,1,1,2,1,1,3,2,2,6,2,2,7,3,1,6,3,1,7 -y 1,1,0,1,1,1,1,2,0,1,2,1,2,2,4,2,2,5,3,1,4,3,1,5 -bc 20 -tac 100 -th 0 -cs 1 -cxys 2 -dt 200 -mst 2 -spc 500 -dp 200 -cha 0 -utpc 1 -hits 1 -n 100000
```
### Explanation of parameters
  
    -f: h5 data file with the extension .h5. The data file was created by ESS DAQ (EFU).
    -x: mapping of detectors, fecs and chips in x direction separated by comma (det,fec,chip, det,fec,chip etc)
    -y: mapping of detectors, fecs and chips in y direction separated by comma (det,fec,chip, det,fec,chip etc)
    -bc: bunch crossing clock. Optional argument (default 20 MHz)
    -tac: tac slope. Optional argument (default 100 ns)
    -th: threshold value in ADC counts. Optional argument
    -cs: minimum cluster size. Optional argument (default 3)
    -cxys: minimum cluster size in x and y together. Optional argument (default 6)
    -dt: maximum time difference between strips in time sorted vector. Optional argument (default 200)
    -mst: maximum missing strips in strip sorted vector. Optional argument (default 2);
    -spc: maximum time span of cluster in one dimension (determined by drift size and speed). Optional argument (default 500)
    -dp: maximum time between matched clusters in x and y. Optional argument (default 200);
    -cha: analyze TDC, BCID and ADC of all channels. Takes a long time. Optional argument (default 0)
    -utpc: use uTPC method to determine timestamp of cluster. Optional argument (default 1)
    -hits: store not only clusters but all hits (a hit is a VMM3 channel over threshold). Creates large files. Optional argument (default 1)
    -n: number of hits. Optional argument, if missing, all hits are analyzed.


## Accessing the produced ROOT tree
The script accessTree in the build directory shows how to access the produced trees. An example.root tree is provided in the source directory.


## Authors

* **Dorothea Pfeiffer** - *Initial work* - [dorotheapfeiffer](https://github.com/dorotheapfeiffer)

See also the list of [contributors](https://github.com/ess-dmsc/vmm-hdf5-to-root/contributors) who participated in this project.

## License

This project is licensed under the BSD 2-Clause "Simplified" License - see the [license](https://github.com/ess-dmsc/vmm-hdf5-to-root/contributors/LICENSE.md) file for details.
