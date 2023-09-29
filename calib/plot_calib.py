#!/usr/bin/python
import os,glob
import subprocess
import re
import sys
import math
import uproot3 as uproot
import pandas as pd
import matplotlib.pyplot as plt
from timeit import default_timer as timer
import time as t
import numpy as np
from matplotlib.colors import LogNorm
import json
import argparse
import os


text_description =  "JSON calibration file plotter. The tool plot_calib_file.py plots one type of calibration at a time (adc, time, time walk)."
parser = argparse.ArgumentParser(description=text_description,
								 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("calib", help="calibration file")
parser.add_argument("-a", "--adc", action="store_true", help="plot ADC calib")
parser.add_argument("-t", "--time", action="store_true", help="plot time calib")
parser.add_argument("-w", "--timewalk", action="store_true", help="plot timewalk calib")


args = parser.parse_args()
config = vars(args)
calib_file = config["calib"]
is_adc = config["adc"]
is_time = config["time"]
is_timewalk = config["timewalk"]

h_channels = np.zeros(64)
color_1 = 'darkblue'
color_2 ='darkgreen'
color_3 = 'rebeccapurple'
color_4 = 'red'
line_width = 3

for i in range(64):
	h_channels[i] = i

hybrids = []
vmms = []
fecs = []
n=0
cols = 8
with open(calib_file, 'r') as file:
	the_calib  = json.load(file)
	for calib in the_calib["vmm_calibration"]:
		n=n+1



rows = int(n/cols)
if n%cols > 0:
	rows = rows + 1

	
fig, ax = plt.subplots(nrows=rows, ncols=cols, figsize=(30, 15))

n=0			
with open(calib_file, 'r') as file:
	the_calib  = json.load(file)
	for calib in the_calib["vmm_calibration"]:
		vmmID = calib["vmmID"]
		vmms.append(vmmID)
		fecID = calib["fecID"]
		fecs.append(fecID)
		hybridID = calib["hybridID"]
		hybrids.append(hybridID)
		if is_adc:
			slopes = calib["adc_slopes"]
			offsets = calib["adc_offsets"]
			for v in range(64):
				if offsets[v] > 15 or offsets[v] < -15:
					print(str(hybridID) + " - fec " + str(fecID) + " vmm" + str(vmmID) + " ch " + str(v) + ": offset " + str(offsets[v]))
			if rows > 1:
				ax[int(n/cols), int(n%cols)].plot(h_channels,offsets, color=color_1, linewidth=line_width, linestyle='solid')
				ax[int(n/cols), int(n%cols)].plot(h_channels,slopes, color=color_4, linewidth=line_width, linestyle='solid')
				ax[int(n/cols), int(n%cols)].title.set_text(hybridID[0:16]+"\n" + hybridID[16:32]+"\nfec"+str(fecID)+" vmm"+str(vmmID)+"\nadc offset/slope")
				ax[int(n/cols), int(n%cols)].set_xlabel('channel')
				#ax[int(n/cols), int(n%cols)].set_ylabel('slope/offset')
				ax[int(n/cols), int(n%cols)].set_ylim(-15, 15)
			else :
				ax[int(n%cols)].plot(h_channels,offsets, color=color_1, linewidth=line_width, linestyle='solid')
				ax[int(n%cols)].plot(h_channels,slopes, color=color_4, linewidth=line_width, linestyle='solid')
				ax[int(n%cols)].title.set_text(hybridID[0:16]+"\n" + hybridID[16:32]+"\nfec"+str(fecID)+" vmm"+str(vmmID)+"\nadc offset/slope")
				ax[int(n%cols)].set_xlabel('channel')
				#ax[int(n%cols)].set_ylabel('slope/offset')
				ax[int(n%cols)].set_ylim(-15, 15)			

	
		if is_time:
			slopes = calib["time_slopes"]
			offsets = calib["time_offsets"]
			for v in range(64):
				if offsets[v] > 20 or offsets[v] < -15:
					print(str(hybridID) + " - fec " + str(fecID) + " vmm" + str(vmmID) + " ch " + str(v) + ": offset " + str(offsets[v]))
			if rows > 1:
				ax[int(n/cols), int(n%cols)].plot(h_channels,offsets, color=color_1, linewidth=line_width, linestyle='solid')
				ax[int(n/cols), int(n%cols)].plot(h_channels,slopes, color=color_4, linewidth=line_width, linestyle='solid')
				ax[int(n/cols), int(n%cols)].title.set_text(hybridID[0:16]+"\n" + hybridID[16:32]+"\nfec"+str(fecID)+" vmm"+str(vmmID)+"\ntime offset/slope")
				ax[int(n/cols), int(n%cols)].set_xlabel('channel')
				#ax[int(n/cols), int(n%cols)].set_ylabel('slope/offset')
				ax[int(n/cols), int(n%cols)].set_ylim(-15, 20)
			else :
				ax[int(n%cols)].plot(h_channels,offsets, color=color_1, linewidth=line_width, linestyle='solid')
				ax[int(n%cols)].plot(h_channels,slopes, color=color_4, linewidth=line_width, linestyle='solid')
				ax[int(n%cols)].title.set_text(hybridID[0:16]+"\n" + hybridID[16:32]+"\nfec"+str(fecID)+" vmm"+str(vmmID)+"\ntime offset/slope")
				ax[int(n%cols)].set_xlabel('channel')
				#ax[int(n%cols)].set_ylabel('slope/offset')
				ax[int(n%cols)].set_ylim(-15, 20)			
		if is_timewalk:
			a = calib["timewalk_a"]
			b = calib["timewalk_b"]
			c = calib["timewalk_c"]
			d = calib["timewalk_d"]
			if rows > 1:
				ax[int(n/cols), int(n%cols)].plot(h_channels,a, color=color_1, linewidth=line_width, linestyle='solid')
				ax[int(n/cols), int(n%cols)].plot(h_channels,b, color=color_2, linewidth=line_width, linestyle='solid')
				ax[int(n/cols), int(n%cols)].plot(h_channels,c, color=color_3, linewidth=line_width, linestyle='solid')
				ax[int(n/cols), int(n%cols)].plot(h_channels,d, color=color_4, linewidth=line_width, linestyle='solid')
				ax[int(n/cols), int(n%cols)].title.set_text(hybridID[0:16]+"\n" + hybridID[16:32]+"\nfec"+str(fecID)+" vmm"+str(vmmID)+"\ntimewalk a/b/c/d")
				ax[int(n/cols), int(n%cols)].set_xlabel('channel')
				#ax[int(n/cols), int(n%cols)].set_ylabel('slope/offset')
				ax[int(n/cols), int(n%cols)].set_ylim(-15, 20)
			else :
				ax[int(n%cols)].plot(h_channels,a, color=color_1, linewidth=line_width, linestyle='solid')
				ax[int(n%cols)].plot(h_channels,b, color=color_2, linewidth=line_width, linestyle='solid')
				ax[int(n%cols)].plot(h_channels,c, color=color_3, linewidth=line_width, linestyle='solid')
				ax[int(n%cols)].plot(h_channels,d, color=color_4, linewidth=line_width, linestyle='solid')
				ax[int(n%cols)].title.set_text(hybridID[0:16]+"\n" + hybridID[16:32]+"\nfec"+str(fecID)+" vmm"+str(vmmID)+"\ntimewalk a/b/c/d")
				ax[int(n%cols)].set_xlabel('channel')
				#ax[int(n%cols)].set_ylabel('slope/offset')
				ax[int(n%cols)].set_ylim(-1000, 1000)			
		n=n+1

plt.show()