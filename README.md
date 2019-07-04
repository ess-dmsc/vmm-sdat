
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
./convertFile 
```

Complete command line:
```
./convertFile -f data.h5 -vmm \"{1,0,2,0},{1,0,2,1},{1,0,2,2},{1,0,2,3},{1,1,2,6},{1,1,2,7},{1,1,2,8},{1,1,2,9}\" -axis \"{{1,0},0},{{1,1},0}\" -bc 40 -tac 60 -th 0 -cs 1 -ccs 2 -dt 200 -mst 1 -spc 500 -dp 200 -coin center-of-mass -ratio -hits 1 -json 0 -n 0

```
### Explanation of parameters
  
    -f: h5 data file with the extension .h5. The data file was created by ESS DAQ tool
    
    -vmm: mapping of detectors, plane, fecs and chips starting and ending with " and separated by brackets
        and comma {{det, plane, fec,chip}, {det, plane, fec, chip}, etc.}
        The tuples for the VMMs are defined as follows:
            detector (choose a number between 0 and 255)
            plane (0 or 1)
            fec (fecID set in firmware based on IP address, 10.0.0.1 is fecID 1, 10.0.0.2 is fecID 2 and so on)
            vmm (depends on connection of hybrid to FEC, FEC channel 1 equals VMMs 0 and 1, FEC channel 2 
                VMMs 2 and 3, FEC channel 8 VMMs 14 and 15)
        When looking at the detector, the following conventions are used:
            - top side of the hybrids is visible (if the hybrids are mounted in the readout plane)
            - side of the Hirose connector (bottom of the hybird) is visible (hybrids mounted on detector side)
            - plane 0 is at the bottom (HDMI cables go downwards)
            - plane 1 is at the right side (HDMI cables go to the right)
        If one looks at a VMM3a hybrid (connector to detector readout is on the bottom side), 
            the channel 0 of the VMM 0 is always where the HDMI cable is connected
        If the planes are correctly used as described above, the VMM IDs are always in icreasing order 
            PER HYBRID (e.g. 14, 15 or e.g. 0, 1)

    -axis: direction of axis. Detector, plane and direction flag (if direction flag = 1, axis direction is flipped). 
        Detector, plane and direction flag starting and ending with " and separated by bracket and comma 
            [[[det,plane],flag], [[det, plane],flag]]
        The tuples for the axes are defined as follows:
            - detector (choose a number between 0 and 255)
            - plane (0 or 1)
            - flip axis flag (0 or 1)
        Using the convention described above, if the plane axis is NOT FLIPPED:
            - plane 0 is at the bottom and goes from left (0) to right (255)
            - plane 1 is at the right and goes from bottom (0) to top (255)
        If the plane axis is FLIPPED:
            - plane 0 is at the bottom and goes from right (255) to left (0)
            - plane 1 is at the right and goes from top (0) to bottom (255)

    -bc: bunch crossing clock. Optional argument (default 40 MHz)

    -tac: tac slope. Optional argument (default 60 ns)

    -th: threshold value in ADC counts. Optional argument (default 0)

    -cs: minimum cluster size per plane. Optional argument (default 1)

    -ccs: minimum cluster size in plane 0 and plane 1 together. Optional argument (default 2)

    -dt: maximum time difference between strips in time sorted vector. Optional argument (default 200)

    -mst: maximum missing strips in strip sorted vector. Optional argument (default 2)

    -spc: maximum time span of cluster in one dimension (determined by drift size and speed). 
        Optional argument (default 500)

    -dp: maximum time between matched clusters in x and y. Optional argument (default 200)

    -coin: Valid clusters normally occur at the same time in plane 0 and plane 1 of a detctor. 
        The parameter -dp determines the permitted time difference between the planes
        The time can be calculated with the center-of-mass algorithm (center-of-mass), the uTPC method (utpc) 
        or the center-of-mass squared method (charge2). Optional argument (default center-of-mass)

    -ratio: Valid clusters normally have the same amount of charge in both detector planes 
        (ratio of (charge plane 0 / charge plane 1) is 100\% or 1
        The desired ratio for the matching can be set as optional argument, the default is 2 or 200\%, 
        i.e. the charge in plane 0 has to be between 50\% and 200\% of the charge in plane 1

    -hits: store not only clusters but all hits (a hit is a VMM3 channel over threshold). 
        Creates large files. Optional argument (default 1)

    -json: create a json file of the detector images. Optional argument (default 1)

    -n: number of hits to analyze. Optional argument (default 0, i.e. all hits)


  

## Accessing the produced ROOT tree
The script accessTree in the build directory shows how to access the produced trees. An example.root tree is provided in the source directory.


## Authors

* **Dorothea Pfeiffer** - *Initial work* - [dorotheapfeiffer](https://github.com/dorotheapfeiffer)

See also the list of [contributors](https://github.com/ess-dmsc/vmm-hdf5-to-root/contributors) who participated in this project.

## License

This project is licensed under the BSD 2-Clause "Simplified" License - see the [license](https://github.com/ess-dmsc/vmm-hdf5-to-root/contributors/LICENSE.md) file for details.
