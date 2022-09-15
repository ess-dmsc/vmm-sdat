# README file for monitoring and pcapng data storage 

## Monitoring a detector
The script monitorDetector.py takes data via tshark for a duration of a few seconds (duration and internet interface to be specified by user). Then the data in this pcapng file is analysed with vmm-sdat. The user has to edit the parameters for the analysis. The analysis creates a root tree, from which matplotlib plots are created. These plots are stored in png files. The html scripts monitorDetector.html continuously loadd the created images. In monitorDetector.html one can set the refresh time in seconds of the page in the meta tag: <meta http-equiv="refresh" content="10" />

To test the script with the provided test pcapng files, type
> python monitorDetector.py 
and open monitorDetector.html in a browser.

To use the monitorDetector script with real data from a network interface, edit the script according the explanation on top of the script.

## Monitoring the GdGEM
The GdGEM is divided in four quadrants, monitorGdGEM.html shows how to create a monitor with several tabs.

## Taking physics data and storing it
The python script storePcapngData.py illustrates how to take data using dumpcap either using a fixed duration per file or a fixed file size. The user has to edit the number of desired files, the file name and the internet interface.
