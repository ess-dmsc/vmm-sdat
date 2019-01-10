
# vmm-hdf5-to-root

Analysis software for VMM3a HDF5 files written by the EFU. From the HDF5 file, a root tree with the hits and clusters is created.
For more information about ROOT see [ROOT](https://root.cern.ch/)

## Getting Started

### Prerequisites
- Boost system library [find_package( Boost REQUIRED COMPONENTS system)]
- HDF5 1.10 or newer [find_package(HDF5 1.10 REQUIRED)]
- h5cpp from https://github.com/ess-dmsc/h5cpp [find_package(h5cpp REQUIRED)]
- root6 from https://root.cern.ch/downloading-root [find_package(ROOT REQUIRED)]



### Installing

To build the program, start in the program directory vmm-hdf5-to-root
> mkdir build
> cd build
> cmake .. or (to enable debug info) cmake -DENABLE_DTRACE=ON ..



## Deployment
How to run the program:
For help concerning the possible parameters, type:
> ./convert 


Complete command line:
> ./convert -f ~/data/a00010.h5 -x 1,0,1,1 -y 2,14,2,15 -bc 20 -tac 100 -th 0 -cs 3 -cxys 6 -dt 200 -mst 2 -spc 500 -dp 200 -cha 0 -utpc 1 -hits 1 


## Built With
* [CMAKE](https://cmake.org/) - Cross platform makefile generation


## Contributing
TBD


## Authors

* **Dorothea Pfeiffer** - *Initial work* - [dorotheapfeiffer](https://github.com/dorotheapfeiffer)

See also the list of [contributors](https://github.com/ess-dmsc/vmm-hdf5-to-root/contributors) who participated in this project.

## License

TBD
