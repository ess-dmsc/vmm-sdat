
# vmm-sdat

VMM3a/SRS Data Analysis Tool: Analysis software for VMM3a data, recorded with the SRS as PCAP or HDF5 files (GdGEM pipeline of the EFU). From the PCAP or HDF5 file, a root tree with the hits and clusters is created.
For more information about ROOT see [here](https://root.cern.ch/)

## Getting Started

### Prerequisites
- Boost system library [find_package( Boost REQUIRED COMPONENTS system)]
- HDF5 1.10 or newer [find_package(HDF5 1.10 REQUIRED)]
- h5cpp from https://github.com/ess-dmsc/h5cpp [find_package(h5cpp REQUIRED)]
- root6 from https://root.cern.ch/downloading-root [find_package(ROOT REQUIRED)]
- Installation of ESS DAQ https://github.com/ess-dmsc/essdaq or event formation unit https://github.com/ess-dmsc/event-formation-unit


## Installing

To build the program, start in the program directory vmm-hdf5-to-root
```
mkdir build
cd build
cmake .. or (to enable debug info) cmake -DENABLE_DTRACE=ON ..

The program uses conan to install all the ESS dependencies.
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
./convertFile -f data.pcapng 
-vmm "[[1,0,2,2],[1,0,2,3],[1,0,2,0],[1,0,2,1],[1,1,2,8],[1,1,2,9],[1,1,2,6],[1,1,2,7]]" 
-axis "[[1,0],0],[[1,1],0]" -bc 40 -tac 60 -th 0 -cs 1 -ccs 2 -dt 200 -mst 1 -spc 500 
-dp 200 -coin center-of-mass -crl 0.5 -cru 2 -save 111 -json 0 -n 0 -algo 0 -swap 0 
-cal CalibrationFile.json

```

## Description of analysis program

The ESS DAQ https://github.com/ess-dmsc/essdaq provides in various detector pipelines the option to write hits 
to hdf5 files. The format of this hdf5 file is defined in the ESS DAQ in the Readout.h files. For GEM detectors 
and the SRS readout, the Readout.h file can be found in 
https://github.com/ess-dmsc/event-formation-unit/blob/master/src/gdgem/nmx/Readout.h

### Time calculation
The convertFile utility of the vmm-hdf5-to-root package analyses the hdf5 files of the gdgem/SRS pipeline. There 
are two time stamps in the hdf5 file, the srs_timestamp coming from the SRS front end card (FEC), and the 
chiptime coming from the VMM ASIC. The VMM measures time in BCID and TDC, with a 40 MHz BC clock the BCID has 
a 25 ns resolution. The BCID is a 12 bit number with values going from 0-4095, that means the bc_time is covering 
a range of 4096 * 25 ns = 102.4 us. The TDC gives information about the time between BCIDs. With a TAC slope of 
60 ns, the time resolution of the tdc_time is 60ns/256 bits = 0.23 ns. The tdc_time and the bc_time together are 
called chiptime. 

The formula is chiptime = bc_time - tdc_time = (BCID + 1) * 25 ns - TDC * 60ns / 256. 
The TDC time is subtracted, because the TDC starts counting when the peakfinder has found the hit, and is then 
stopped by the falling edge of the next BC clock pulse. A small TDC value means the hit occured shortly before 
the next BC clock, whereas a large TDC value means it appeared a long time before the next clock. 

Every 4096 * 25 ns the BCID overflows, these overflows are called offset in the FEC firmware. 
Example:
    Offset 24, BCID 100, TDC 80: 24 * 102.4 us + (100+1) * 25 ns - 80 * 60 ns/256 = 2460106.25 ns = 2.46 ms
But even with the offset the time range covered is only 32 * 102.4 us = 3278.8 us = 3.3 ms. Therefore every 
3.3 ms 42 bit time markers are send by the FEC card. All the offsets always refer to the previous time marker 
sent for the particular VMM. The EFU now takes the time markers and the offset and calculates the srs_timestamp 
from it in unit [ns]. 

### Calibration files
The ESS DAQ has the option to load JSON calibration files, that contain for each channel of each VMM ASIC 
an ADC (adc_offset and adc_slope) and a time (adc_offset and adc_slope) correction.
https://github.com/ess-dmsc/event-formation-unit/blob/master/src/gdgem/srs/CalibrationFile.cpp
The JSON calibration file can either be produced automatically with the VMM Slow Control program
https://gitlab.cern.ch/mguth/VMM-software-RD51
or by other means. If the calibration file is loaded during the data acquisition with the ESS DAQ, then the 
chip_time has been calculated using the following formula:
- new_chiptime = BCID * 25 ns + ( 25 ns - tdc_time - time_offset) * time_slope
For the ADC, the formula is: 
- new_adc = (adc - adc_offset) * adc_slope;
If the data has been saved to file without the use of calibration files in the ESS DAQ, then the calibration
can be loaded during the analysis by convertFile using the parameter -cal.

### Hits
As first step, the hits are stored in vectors (optionally also in a branch of the root tree), depending on the 
mapping of the VMM ASICs to the detectors and detector planes. Then, to create clusters for each detector plane, 
the hits are first sorted in time (which is the sum of srs_time and chiptime). 

### Time clusters
A time cluster contains all hits with subsequent timestamps, that have a time difference smaller or equal than 
the value defined in the -dt (delta t) parameter. Further, the time difference between the first and the last hit 
has to be smaller than the time span defined in the -spc (span cluster) parameter. 

### Plane clusters
The time clusters are then sorted by strip to create the plane clusters. For hits to belong to one cluster, the gap 
between neighbouring strips cannot be larger than the -mst (missing strips) parameter. For each of the plane clusters,
a position and a time is calculated using three different algorithms:
    - center-of-mass (charge as weight) 
    - center-of-mass2 (charge squared as weight)
    - utpc (latest time)
The user can define additional algorithms in the method Clusterer::AlgorithmUTPC() in
https://github.com/ess-dmsc/vmm-hdf5-to-root/blob/master/Clusterer.cpp

The additional algorithm is picked depending on the -algo parameter. At the moment, two additional algorithms are 
defined there 
    - 0 = utpc center-of-mass2 (center of mass squared of the strip with the latest time and its one or two neighbours)
    - 1 = utpc center-of-mass (center of mass of the strip with the latest time and its one or two neighbours)
If the number of strips in a plane cluster is equal or larger than the value specified in the -cs (cluster size) 
parameter, the cluster is valid. For X-ray data the cluster size parameter is usually 1 or 2, whereas for electron or
alpha tracks it can be 3 or larger. 
Valid plane clusters are stored in vectors (optionally also in a branch of the root tree). 

### Detector clusters
To determine now the detector clusters, the plane clusters in the two detector planes are matched depending on their
cluster times. The time difference between the cluster times of two plane clusters has to be smaller or equal than the 
value specified in the -dp (delta planes) value. The user has to choose in the -coin (coincidence) parameter, which of 
the four cluster times calculated with the different algorithms has to be used for the matching. The sum of all hists 
in the two planes of the detector cluster has to be larger or equal to the -ccs (coincident cluster size) parameter.
In a detector with equal charge sharing between the two planes, the charge of a detector cluster should be almost 
identical in the two detector planes. The -ratio (charge ratio) parameter defines the allowed charge ratio between
the two planes (charge plane 0 / charge plane 1) or (charge plane 1 / charge plane 0).


## Explanation of parameters
  
    -f: h5 data file with the extension .h5 or .pcapng. The data file was created by ESS DAQ tool
    
    -vmm: mapping of detectors, plane, fecs and chips starting and ending with " and separated by brackets
        and comma [[det, plane, fec,chip], [det, plane, fec, chip], etc.]
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
   	For a 1D detector,  only axis 0 is defined, for a 2D detector axis 0 and 1 have to be defined.
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
    
    -sc: Scale coordinates. Per detector a tuple with three values in mm, 
        e.g for two detectors [[s0,s1,s2], [s0,s1,s2]]
    
    -tl: Translate coordinates. Per detector a tuple with three values in mm, 
        e.g for two detectors [[t0,t1,t2], [t0,t1,t2]]
    
    -ro: Rotate around plane 0, plane 1, plane 2. Per detector a tuple with three angles in degrees, 
        e.g for two detectors [[r0,r1,r2], [r0,r1,r2]]
    
    -tr: Transform detector coordinates. 
        S=scale, T=translate, R0=rotation plane 0, R1=rotation plane1, R2=rotation plane2
        example (two detectors): -tr [[S,T,R2], [S,T,R0]]. 
        First detector scaling, then translation, then rotation around normal axis to plane0 and plane 1.
        Second detector scaling, then translation, then rotation around plane0 axis
        Normally beam in direction of plane 2 or z, hence: Plane 0 (x), plane 1 (y), plane 2 (z)
   
    -bc: bunch crossing clock. Optional argument (default 40 MHz)

    -tac: tac slope. Optional argument (default 60 ns)

    -th: threshold value in ADC counts. Optional argument (default 0, if -1, only hits with over threshold flag 1 are accepted)

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

    -crl: Valid clusters normally have the same amount of charge in both detector planes 
    	(ratio of charge plane 0 / charge plane 1 is 100% or 1). Depending on the readout, 
    	the charge sharing can be different, e.g. in a standard GEM strip readout the total 
    	charge is divided 60/40 between plane 0/plane 1.
		With -crl one sets the lower threshold for the plane0/plane1 charge ratio. 
		Optional argument (default 0.5).

	-cru: With -cru one sets the upper threshold for the plane0/plane1 charge ratio. 
		Optional argument (default 2).

	-swap: Same connectors on readout boards unintentionally swap odd and even channels. 
		With -swap 1 one can correct this. Optional parameter (default 0).

    -save: Select which data to store in root file. Input is a 3 bit binary number.
        bit 0 (LSB): hits (a hit is a VMM3a channel over threshold)
        bit 1      : clusters plane" << std::endl;
        bit 2 (MSB): clusters detector" << std::endl;
        Examples:
            0001: hits only 
            0100: clusters detector only 
            0111: everything, hits, clusters plane, clusters detector
   
    -stats: Show statistics of the run (default 1, show all stats)

    -json: create a json file of the detector images. Optional argument (default 0)

    -n: number of hits to analyze. Optional argument (default 0, i.e. all hits)

    -algo: There are three different algorithms implemented that calulate cluster times: 
        center-of-mass (charge as weight), 
        center-of-mass2 (charge squared as weight),
        utpc (latest time)
        The utpc algorithm works as follows:
            1. The strip with the largest time is chosen as position
            2. If there are several strips with identical largest times, the indices of the 
            outermost strips are retained
                Example: The cluster has 8 strips (strip 10-17, index 0-7), strips 12, 13 and 16
                have the same largest time, so strip 12 and 15 are retained. The indices of 
                strips 12 and 16 are 2 and 6 respectively
            3. The strip with the index that is the closest to the start or end of the track
            is retained. In the example, index 6 is one index away from the end of the track,
            and index 2 two indices, thus index 6 is retained.
            4. If the outermost indices have the same distance from the start/end of the track,
            the strip with the largest ADC value is retained.
            5. If the ADC values are the same, the center-of-mass2 is taken.

        The user can define an additional utpc algorithm in the method Clusterer::AlgorithmUTPC(),
        the algorithm is chosen depending on the -algo parameter. At the moment, two additional 
        algorithms are defined there, 
            0 = utpc center-of-mass2 (center of mass squared of the strip 
            with the latest time and its one or two neighbours)
            1 = utpc center-of-mass (center of mass of the strip 
            with the latest time and its one or two neighbours)

    -cal: Name of the calibration file. A calibration file is a JSON file containing an 
    	ADC and/or time correction in the form of a slope and an offset correction. 
    	The calibration file can be produced with the VMM slow control tool. Optional parameter. 

  

## Accessing the produced ROOT tree
The script accessTree in the build directory shows how to access the produced trees. An example.root tree is provided in the source directory.


## Authors

* **Dorothea Pfeiffer** - *Initial work* - [dorotheapfeiffer](https://github.com/dorotheapfeiffer)
* **Lucian Scharenberg** - *Debugging, Testing, Inspiration for new features, Python analysis*  - [lscharenberg](https://github.com/lscharenberg)

See also the list of [contributors](https://github.com/ess-dmsc/vmm-hdf5-to-root/contributors) who participated in this project.

## License

This project is licensed under the BSD 2-Clause "Simplified" License - see the [license](https://github.com/ess-dmsc/vmm-hdf5-to-root/contributors/LICENSE.md) file for details.
