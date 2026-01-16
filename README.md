
# vmm-sdat

VMM3a/SRS Data Analysis Tool: Analysis software for VMM3a data, recorded with the SRS or the ESS readout as PCAPNG file. From the PCAP file, a root tree with the hits and clusters is created.
For more information about ROOT see [here](https://root.cern.ch/)

## Getting Started

### Which branch do I need?
The main branch analyses PCAPNG files in the SRS or ESS data format. PCAPNG files can be generated with Wireshark or tcdump. The SRS branch analyses PCAPNG files, or HDF5 files that have been created with the GdGEM pipeline of the ESS DAQ/event formation unit. The old SRS data format (offsets 0-31) and the new SRS data format (offsets -1 to 15) are supported.

### Prerequisites
- Boost system library [find_package( Boost REQUIRED COMPONENTS system)]
- root6 from https://root.cern.ch/downloading-root [find_package(ROOT REQUIRED)]
- Installation of ESS DAQ https://github.com/ess-dmsc/essdaq or event formation unit https://github.com/ess-dmsc/event-formation-unit


## Installing

To build the program, start in the program directory vmm-sdat
```
mkdir build
cd build
cmake .. or (to enable debug info) cmake -DENABLE_DTRACE=ON ..
```
The program uses conan to install all the ESS dependencies.

## Built With
* [CMAKE](https://cmake.org/) - Cross platform makefile generation

## Deployment
How to run the program:
The executable is called convertFile, it has been created in the build directory. For help concerning the possible parameters, type:
```
./convertFile 
```

Example command line:
```
./convertFile -f data.pcapng 
-vmm "[[1,0,2,2],[1,0,2,3],[1,0,2,0],[1,0,2,1],[1,1,2,8],[1,1,2,9],[1,1,2,6],[1,1,2,7]]" 
-axis "[[1,0],0],[[1,1],0]" -bc 40 -tac 60 -th 0 -cs 1 -ccs 2 -dt 200 -mst 1 -spc 500 
-dp 200 -coin center-of-mass -crl 0.5 -cru 2 -save [[1],[1],[1]] -json 0 -n 0 -algo 0 -swap 0 
-cal CalibrationFile.json -df SRS

```

## Description of analysis program

### Time calculation
The convertFile utility of the vmm-sdat package analyses pcapng directly from the SRS FEC, or the hdf5 files of the 
gdgem/SRS pipeline. There are two time stamps in the hdf5 file, the srs_timestamp coming from the SRS front end card (FEC), and the 
chiptime coming from the VMM ASIC. To obtain the total time, these two timestamps are added, and the new field in the 
root file is just called time.

The VMM measures time in BCID and TDC, with a 44.444 MHz BC clock the BCID has a 22.5 ns resolution. The BCID is a 12 bit 
number with values going from 0-4095, that means the bc_time is covering a range of 4096 * 22.5 ns = 92.16 us. 
The TDC gives information about the time between BCIDs. With a TAC slope of 
60 ns, the time resolution of the tdc_time is 60ns/256 bits = 0.23 ns. The tdc_time and the bc_time together are 
called chiptime. 

The formula is chiptime = bc_time - tdc_time = (BCID + 1.5) * 22.5 ns - TDC * 60ns / 256. 
The TDC time is subtracted, because the TDC starts counting when the peakfinder has found the hit, and is then 
stopped by the falling edge of the next BC clock pulse. A small TDC value means the hit occured shortly before 
the next BC clock, whereas a large TDC value means it appeared a long time before the next clock. 

Every 4096 * 22.5 ns the BCID overflows, these overflows are called offset in the FEC firmware. 
Example:
    Offset 12, BCID 100, TDC 80: 12 * 92.16 us + (100+1.5) * 22.5 ns - 80 * 60 ns/256 = 1.18 ms
But even with the offset the time range covered is only 16 * 92.16 us = 1474.56 us = 1.47 ms. Therefore every 
1.47 ms 42 bit time markers are send by the FEC card. 

In the time marker we have th 42 bit srs clock counter, which is measured in 44 MHz clock cycles. So to convert it to a timestamp
in ns one has to multiply it with 22.5 ns. All the offsets always refer to the previous time marker 
sent for the particular VMM. vmm-sdat takes the time markers and the offset and calculates the srs_timestamp 
from it in unit [ns]. The chiptime is calculated from BCID and TDC, and the complete timestamp is the sum of 
srs_timestamp and chiptime as shown in the example above.

### Calibration files
The VMM slow control (https://gitlab.cern.ch/rd51-slow-control/vmmsc) can produce JSON calibration files, that 
contain for each channel of each VMM ASIC an ADC (adc_offset and adc_slope) and a time (adc_offset and adc_slope) correction.

If the calibration file is loaded during the analysis with vmm-sdat, then the 
chip_time has been calculated using the following formula:
- new_chiptime = BCID * 22.5 ns + ( 1.5*22.5 ns - tdc_time - time_offset) * time_slope
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
https://github.com/ess-dmsc/vmm-sdat/blob/master/src/Clusterer.cpp

The additional algorithm is picked depending on the -algo parameter. At the moment, four additional algorithms are 
defined there 
1: fit gaussian for position
2: center-of-mass over threshold only
3: center-of-mass2 over threshold only
4: position and time strip with highest ADC

If the number of strips in a plane cluster is equal or larger than the value specified in the -cs (cluster size) 
parameter, the cluster is valid. For X-ray data the cluster size parameter is usually 1 or 2, whereas for electron or
alpha tracks it can be 3 or larger. 
Valid plane clusters are optionally stored in a branch of the root tree. For pad detectors the clusters appear with the 
x-indices of the pads in plane 0 of the detector, and with the y-indices in plane 1 of the detector.  

### Detector clusters
To determine now the detector clusters, the plane clusters in the two detector planes are matched depending on their
cluster times. The time difference between the cluster times of two plane clusters has to be smaller or equal than the 
value specified in the -dp (delta planes) value. The user has to choose in the -coin (coincidence) parameter, which of 
the four cluster times calculated with the different algorithms has to be used for the matching. The sum of all hists 
in the two planes of the detector cluster has to be larger or equal to the -ccs (coincident cluster size) parameter.
In a detector with equal charge sharing between the two planes, the charge of a detector cluster should be almost 
identical in the two detector planes. The -ratio (charge ratio) parameter defines the allowed charge ratio between
the two planes (charge plane 0 / charge plane 1) or (charge plane 1 / charge plane 0). For pad detectors, the detector
clusters contain the x-index in the pos0 entry, and the y-index in the pos1 entry. Detector that are only 1D, contain the 
plane clusters also in the detector cluster tree. The entries of the missing plane (either plane 0 or plane 1) are set to zero.

## Explanation of parameters
  
    -f: data file with the extension .pcapng. The data file was created by Wireshark or tcdump
    
    -df: data format of pcapng file. Default is SRS format (valid offsets: -1 - 15), other option is data format TRIG 
   
    -geo:   the detector geometry has to be defined in a JSON file. Tow example geometry files (for strip and pad detectors) are in the run folder.
    
    -bc: bunch crossing clock. Optional argument (default 44.4444 MHz)

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

    -algo:  There are four different algorithms implemented that calulate cluster position/times: 
        1: fit gaussian for position
		2: center-of-mass over threshold only
		3: center-of-mass2 over threshold only
		4: position and time strip with highest ADC

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

    -save:  select which data to store in root file. Input is a list of lists of detectors, e.g. [[1,2],[1,2],[1,2,3]]
            first list : detectors for which to write the hits (hit is a VMM3a channel over threshold)
            second list : clusters plane
            third list : clusters detector
            Examples:
                [[1,2],[],[]]: hits for detectors 1 and 2 only 
                [[],[],[1,2]]: clusters detector for detector 1 and 2 only 
                [[2],[1],[1]]: hits for detector 2, clusters plane, clusters detector for detector 1 
   
    -stats: Show statistics of the run (default 1, show all stats)

    -json: create a json file of the detector images. Optional argument (default 0)

    -n: number of hits to analyze. Optional argument (default 0, i.e. all hits)

    -cal: Name of the calibration file. A calibration file is a JSON file containing an 
    	ADC and/or time correction in the form of a slope and an offset correction. 
    	The calibration file can be produced with the VMM slow control tool. Optional parameter. 

    -t0: SRS data format: Time correction in ns for each vmm in the format
        [[fec0, vmm0,correction0], [fec1, vmm1, correction1]]. The correction value is subtracted from all timestamps. 
        If instead of a number the word 'run' is put as correction, the first timestamp of the run is used as correction. For the ESS data format, the first timestamp is always substracted to get smaller timestamps.
        Optional argument for SRS data format (default 0)

## Running the program
In the run folder there is a script and example data that show how to use the convertFile utility. 
runConversion.py converts the Wireshark pcapng file example.pcapng into a root tree.
Building on this script you can construct your own scripts for file conversion. The meaning of the 
command line parameters is explained above. For sure you will have to change the mapping of FECs and VMMs
to the axes of your detector.
  

## Accessing the produced ROOT tree
The script accessTree in the build directory shows how to access the produced root trees. The data in the 
root trees can be further analyzed and plotted. An example.root tree is provided in the source directory.
The time of the first pcapng packet is stored in the file as date and as epochtime in seconds. You can retrieve the values like that:
TString t  = f->Get("unixtime")->GetTitle()
TString d  = f->Get("date")->GetTitle()

## Authors

* **Dorothea Pfeiffer** - *Initial work, implementation of new features* - [dorotheapfeiffer](https://github.com/dorotheapfeiffer)
* **Lucian Scharenberg** - *Debugging, testing, proposing new features, Python analysis*  - [lscharenberg](https://github.com/lscharenberg)

See also the list of [contributors](https://github.com/ess-dmsc/vmm-sdat/contributors) who participated in this project.

## License

This project is licensed under the BSD 2-Clause "Simplified" License - see the [license](https://github.com/ess-dmsc/vmm-sdat/contributors/LICENSE.md) file for details.
