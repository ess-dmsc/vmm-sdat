#!/usr/bin/python
import os
import subprocess
import re
import sys
import time as t

#Edit these fields to determine the acquisition
numberOfFiles = 100
fileName = "sipm_test"
interface="enp0s31f6"
#Choose either duration in seconds or filesize in kB
duration=60
filesize=1000000



cnt = 0
try:
    while cnt < numberOfFiles:
        theDate = t.strftime("%Y%m%d%H%M%S") 	
        name = fileName + "_" + f'{cnt:05}' + "_" + str(theDate) + ".pcapng"
        args = ['dumpcap', '-w', name, '-a', 'duration:' + str(duration), '-i', interface]
        #args = ['dumpcap', '-w', name, '-a', 'filesize:' + str(filesize), '-i', interface]
        subprocess.call(args)
        cnt = cnt + 1
except OSError:
    pass
