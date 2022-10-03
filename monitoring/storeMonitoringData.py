#!/usr/bin/python
import os,glob
import subprocess
import re
import sys
import time as t

#Edit these fields to determine the acquisition
#File sizes (specified in kB) between 10 MB and 20 MB seem to give the best results
filesize=20000
numberOfFiles = 10
fileName = "monitoring"
#interface="enp1s0f0"
interface="en0"





cnt = 0
try:
	fileList = glob.glob("./" + fileName + "*.root", recursive=True)
	for file in fileList:
    		try:
        		os.remove(file)
    		except OSError:
        		print("Error while root files..")

	fileList = glob.glob("./" + fileName + "*.pcapng", recursive=True)
	for file in fileList:
    		try:
        		os.remove(file)
    		except OSError:
        		print("Error while deleting pcapng files..")

	fileList = glob.glob("./" + fileName + "_semaphore.txt", recursive=True)
	for file in fileList:
    		try:
        		os.remove(file)
    		except OSError:
        		print("Error while deleting ./" + fileName + "_semaphore.txt file")

	while True:
		fileId = cnt%numberOfFiles
		name = fileName + "_" + f"{fileId:05}" + ".pcapng"
		args = ["dumpcap", "-w", name, "-a", "filesize:" + str(filesize), "-i", interface]
		subprocess.call(args)
		with open("./" + fileName + "_semaphore.txt", "w") as f:
			f.write(str(cnt%numberOfFiles))
		cnt = cnt + 1
except OSError:
    pass
