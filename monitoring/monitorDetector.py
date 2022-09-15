#!/usr/bin/python
import os
import subprocess
import re
import sys
import uproot3 as uproot
import pandas as pd
import matplotlib.pyplot as plt
import time as t
import numpy as np

#########################
#To take data via tshark, please uncomment line 113
#For the vmm-sdat analysis, edit parameters in line 115
#########################

#########################
#Edit the paramters below
duration=20
interface='en0'
detectors = 1 
channels_x = 256
channels_y = 256
file_name="example_monitoring"
##########################
    
def plot_data(fileName, ndet, binsx, binsy):
	# Get the tree in the ROOT file
	tree_detector = uproot.open(fileName)['clusters_detector']
	tree_hits = uproot.open(fileName)['hits']


	d_adc0 = tree_detector.array('adc0')
	d_pos0 = tree_detector.array('pos0')
	d_adc1 = tree_detector.array('adc1')
	d_pos1 = tree_detector.array('pos1')
	d_time0 = tree_detector.array('time0')
	d_time1 = tree_detector.array('time1')
	d_det = tree_detector.array('det')
	
	adc = tree_hits.array('adc')
	pos = tree_hits.array('pos')
	time = tree_hits.array('time')
	det = tree_hits.array('det')
	plane = tree_hits.array('plane')
	
	data_hits = {'adc': adc,'pos': pos,'det': det,'plane': plane}
	df_hits = pd.DataFrame(data_hits)

	data_clusters = {'adc0': d_adc0,'pos0': d_pos0,'adc1': d_adc1,'pos1': d_pos1,'det': d_det}
	df_clusters = pd.DataFrame(data_clusters)	


	theDate = t.strftime("%Y%m%d-%H%M%S") 	
	for i in range(0, ndet):
		fig, ax = plt.subplots(2, 4, gridspec_kw={'height_ratios': [1, 1]})
		hits0 = df_hits.query("plane == 0 and det == " + str(i))
		hits1 = df_hits.query("plane == 1 and det == " + str(i))
		ax[0, 0].hist(hits0['pos'], bins = binsx, range = [0.0, binsx], color='darkblue')
		ax[0, 1].hist(hits1['pos'], bins = binsy, range = [0.0, binsy], color='darkblue')		
		ax[0, 2].hist2d(hits0['pos'], hits0['adc'],bins =[binsx, 128],cmap=plt.cm.jet,range=np.array([(0, binsx), (0, 1024)]))
		ax[0, 3].hist2d(hits1['pos'], hits1['adc'],bins =[binsy, 128],cmap=plt.cm.jet, range=np.array([(0, binsy), (0,1024)]))		
		
		cl = df_clusters.query("det == " + str(i))
		ax[1, 0].hist(cl['pos0'], bins = binsx, range = [0.0, binsx], color='rebeccapurple')
		ax[1, 1].hist(cl['pos1'], bins = binsy, range = [0.0, binsy], color='rebeccapurple')
		ax[1, 2].hist(cl['adc0']+cl['adc1'], bins = 1000, range = [0.0, 5000],  color='rebeccapurple')
		ax[1, 3].hist2d(cl['pos1'], cl['pos0'],bins =[binsx, binsy],cmap=plt.cm.viridis,range=np.array([(0, binsx), (0, binsy)]))
		
		ax[0, 0].title.set_text("hits pos0")
		ax[0, 1].title.set_text("hits pos1")
		ax[0, 2].title.set_text("hits adc0")
		ax[0, 3].title.set_text("hits adc1")
		ax[1, 0].title.set_text("clusters pos0")
		ax[1, 1].title.set_text("clusters pos1")
		ax[1, 2].title.set_text("clusters adc0+adc1")
		ax[1, 3].title.set_text("clusters pos1:pos0")
		
		ax[0, 0].set_xlabel('x [pitch 0.4 mm]')
		ax[0, 0].set_ylabel('counts')
		ax[0, 1].set_xlabel('y [pitch 0.4 mm]')
		ax[0, 1].set_ylabel('counts')
		
		ax[0, 2].set_xlabel('x [pitch 0.4 mm]')
		ax[0, 2].set_ylabel('adc')
		ax[0, 3].set_xlabel('y [pitch 0.4 mm]')
		ax[0, 3].set_ylabel('adc')
		
		ax[1, 0].set_xlabel('x [pitch 0.4 mm]')
		ax[1, 0].set_ylabel('counts')
		ax[1, 1].set_xlabel('y [pitch 0.4 mm]')
		ax[1, 1].set_ylabel('counts')		
		
		ax[1, 2].set_xlabel('charge')
		ax[1, 2].set_ylabel('counts')
		ax[1, 3].set_xlabel('x [pitch 0.4 mm]')
		ax[1, 3].set_ylabel('y [pitch 0.4 mm]')		
		
		pdfName = "det_" + str(i) + "_" + theDate + ".png"
		fig.set_size_inches(17, 8)
		fig.tight_layout()
		#fig.savefig(pdfName, format="png")
		fig.savefig("det_" + str(i) + ".png", format="png")


try:
	cnt = 0
	limit = 20
	bufferSize=2
	while True:
		name = file_name + str(cnt) + ".pcapng"
		args = ['tshark', '-w', name, '-a', 'duration:'+str(duration), '-i', interface]
		#subprocess.call(args)
		
		args_vmmsdat = ['../build/convertFile', '-f', "./" + name, '-geo', 'example_monitoring.json', '-bc', '40', '-tac', '60', '-th','0', '-cs','1', '-ccs', '3', '-dt', '200', '-mst', '1', '-spc', '500', '-dp', '200', '-coin', 'center-of-masss', '-crl', '0.2', '-cru', '10', '-save', '[[0,1,2,3],[],[0,1,2,3]]', '-algo', '4', '-info', '', '-df','SRS']
		subprocess.call(args_vmmsdat)
		
		for file in os.listdir("."):
			if file.startswith(file_name + str(cnt)) and file.endswith(".root"):
				print(file)
				plot_data(file,detectors,channels_x,channels_y)
		
		if cnt == 0:
			cnt = 1
		else:
			cnt = 0
except OSError:
    pass
    
    
    
    
    