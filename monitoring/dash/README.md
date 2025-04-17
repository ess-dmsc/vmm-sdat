# README file for life monitoring

## Software to install for the monitoring tools

- pip install dash

- pip install plotly

- pip install numpy

- pip install pandas

- pip install matplotlib 

- pip install waitress 


## How to monitor data
The script **_monitor.py_** is a multi-process script. The first process acquires *pcapng* files from the network interface and stores them on disk. The second process analyses the oldest pcapng file in *vmm-essdat*, and produces a root tree. After analysing the pcapng the process deletes the file. The third process starts the waitress server and the dash application. The dash application consists of three pages in the page directory, hits, clusters and stats. The hit page loads the data from the oldest root file with uproot, and passes the data to the cluster and stats page.

For the monitoring to work, the hit page has to be always open. The user has to edit the parameters for the analysis in the config file, in the same way as for their normal analysis. The mapping/configuration file for the complete NMX detector is provided in the monitoring directory, as well as a calibration file.

To test the script with the provided example *pcapng* file **monitoring_00000.pcapng_**, *acquire_pcapng=0* has to be set in the script. Then start it from the command line with 

> python3 monitor.py 

and open **nmx_hits.html**, **nmx_clusters.html**  and **nmx_stats.html**  in a browser.

To run the monitoring with real data, *acquire_pcapng=1* has to be set, and the correct network interface has to be defined. If the user wants to keep the generated root trees as experimental data, *delete_root=0* has to be defined.

