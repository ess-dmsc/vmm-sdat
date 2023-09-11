# README file for monitoring and pcapng data storage 

## How to monitor data
The script **_monitorData.py_** reads *pcapng* files from the disk and analyses them in *vmm-sdat*. The user has to edit the parameters for the analysis, in the same way as for their normal analysis. The analysis creates a *root* tree, from which *matplotlib* plots are created. These plots are stored in *png* files. The *html* script **_monitorData.html_** continuously refreshes and thus updates the created images. In **_monitorData.html_** one can set the refresh time in seconds at the top of the page in the meta tag: <meta http-equiv="refresh" content="10" />

To test the script with the provided example *pcapng* file **_example_monitoring_00000.pcapng_**, start it from the command line with 

> python3 monitorData.py 

and open **_monitorData.html_** in a browser.

To run the monitoring with real data, first edit the file **_storeMonitoringData.py_**. This script takes data via *dumpcap* until a certain file size is reached (file size and internet interface to be specified by user). It uses a ring buffer with the number of files specified by the user (10 is a reasonable choice for the buffer size).The number of the file that was last created is stored in the file **_monitoring_semaphore.txt_**. The script **_monitorData.py_** reads the file number from the file **_monitoring_semaphore.txt_**, and then subsequently analyses the file with the correct file number. 

**To read the correct semaphore file and data file, change the file name in monitorData.py  from "example" to "monitoring"!!** The open a second console window and start the program from the command line

> python3 storeMonitoringData.py 


## How to take physics data and store it
The python script **_storePcapngData.py_** illustrates how to take data using *dumpcap* either using a fixed duration per file or a fixed file size. The user has to edit the number of desired files, the file name and the internet interface.


##Interactive Monitoring
The folder interactive monitoring contains a similar type of monitoring, but with interactive plots.
