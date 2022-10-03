#!/usr/bin/python
import os,glob
import subprocess
import re
import sys
import uproot3 as uproot
import pandas as pd
import matplotlib.pyplot as plt
import time as t
import numpy as np
from matplotlib.colors import LogNorm

#########################
#To take data via tshark, please uncomment line 51
#For the vmm-sdat analysis, edit parameters in line 53
#########################

#########################
#Edit the paramters below
file_name="example"
#file_name="monitoring"
num_detectors = 1 
first_detector = 0 
channels_x = 256
channels_y = 256
max_charge = 10000
##########################
h1_total = [0] * channels_x    
h2_total = [0] * channels_y  
h3_total = np.zeros((channels_x, 128))
h4_total = np.zeros((channels_y, 128))
h5_total = [0] * channels_x    
h6_total = [0] * channels_y
h7_total = [0] * 1000
h8_total = np.zeros((channels_x, channels_y))    

h_cluster_rate = [0] * 1000  
h_hit_rate = [0] * 1000 
h_bit_rate = [0] * 1000 
h_time = [0] * 1000 
h_size0 = [0] * 64
h_size1 = [0] * 64 
h_size = [0] * 128
  
firstTime = True
cnt = 0
fileId=0
lastFileId = -1 

try:
	fileList = glob.glob('./det*.png', recursive=True)
	for file in fileList:
    		try:
        		os.remove(file)
    		except OSError:
        		print("Error while deleting monitoring .png files..")
	while True:
		if os.path.isfile('./' + file_name + '_semaphore.txt') :
			f = open('./' + file_name + '_semaphore.txt', "r") 
			fileId = int(f.read())
			print("Semaphore " + str(fileId))
		else:
			print("Cannot open file " + file_name + "_semaphore.txt! Waiting for data..")
			t.sleep(10)
			continue		
		name = file_name + "_" + f'{fileId:05}' + ".pcapng"
		if not os.path.isfile("./" + name) :
			print(name + " not found!")
			continue
		if lastFileId == fileId:
			print("File " + str(lastFileId) + " already analysed! Waiting for new data..")
			t.sleep(2)
			continue
		time_start = int( t.time() )	
		args_vmmsdat = ['../build/convertFile', '-f', "./" + name, '-geo', 'example_monitoring.json', '-bc', '40', '-tac', '60', '-th','0', '-cs','1', '-ccs', '3', '-dt', '200', '-mst', '1', '-spc', '500', '-dp', '200', '-coin', 'center-of-masss', '-crl', '0.2', '-cru', '10', '-save', '[[0],[],[0]]', '-algo', '4', '-info', '', '-df','SRS']
		#args_vmmsdat = ['../build/convertFile', '-f', "./" + name, '-cal', 'complete_calib_2022_09_23.json', '-vmm', '[[0,0,2,10],[0,0,2,11],[0,0,2,8],[0,0,2,9],[0,0,2,12],[0,0,2,13],[0,0,2,2],[0,0,2,3],[0,0,2,0],[0,0,2,1],[0,1,1,10],[0,1,1,11],[0,1,1,8],[0,1,1,9],[0,1,1,4],[0,1,1,5],[0,1,1,2],[0,1,1,3],[0,1,1,0],[0,1,1,1]]','-axis', '[[0,0],0],[[0,1],0]', '-bc', '44.44', '-tac', '60', '-th','0', '-cs','2', '-ccs', '5', '-dt', '200', '-mst', '1', '-spc', '500', '-dp', '500', '-coin', 'center-of-masss', '-crl', '0.2', '-cru', '5', '-save', '[[0],[],[0]]', '-json','0', '-algo', '0', '-info', 'monitor', '-df', 'SRS']			
		
		subprocess.call(args_vmmsdat)
		lastFileId = fileId
		cnt = cnt + 1
		time_middle = int( t.time() )

		print("Analysis " + str(cnt) + ": " + str(time_middle - time_start) + " s")
		
		for file in os.listdir("."):
			if file.startswith(file_name + "_" + f'{fileId:05}') and file.endswith(".root"):
				tree_hits = uproot.open(file)['hits']
				adc = tree_hits.array('adc')
				if adc.size == 0:
					continue
			
				pos = tree_hits.array('pos')
				time = tree_hits.array('time')
				det = tree_hits.array('det')
				plane = tree_hits.array('plane')
				
				tree_detector = uproot.open(file)['clusters_detector']
				d_adc0 = tree_detector.array('adc0')
				d_pos0 = tree_detector.array('pos0')
				d_adc1 = tree_detector.array('adc1')
				d_pos1 = tree_detector.array('pos1')
				d_time0 = tree_detector.array('time0')
				d_time1 = tree_detector.array('time1')
				d_size0 = tree_detector.array('size0')
				d_size1 = tree_detector.array('size1')
				d_det = tree_detector.array('det')
	
				data_hits = {'adc': adc,'pos': pos,'time': time, 'det': det,'plane': plane}
				df_hits = pd.DataFrame(data_hits)

				data_clusters = {'size0': d_size0,'size1': d_size1,'time0': d_time0,'time1': d_time1,'adc0': d_adc0,'pos0': d_pos0,'adc1': d_adc1,'pos1': d_pos1,'det': d_det}
				df_clusters = pd.DataFrame(data_clusters)	

				theDate = t.strftime("%Y%m%d-%H%M%S") 	
				for i in range(first_detector, num_detectors):
					hits0 = df_hits.query("plane == 0 and det == " + str(i))
					hits1 = df_hits.query("plane == 1 and det == " + str(i))
					cl = df_clusters.query("det == " + str(i))
					num_hits = hits0["time"].size +  hits1["time"].size
					if num_hits == 0:
						continue
					dt_hits = max(hits0["time"]) - min(hits0["time"])
					h_hit_rate.pop(0)
					h_bit_rate.pop(0)
					if dt_hits > 0:
						h_hit_rate.append(num_hits*1000000000.0/dt_hits)
						h_bit_rate.append(48*num_hits*1000000000.0/(dt_hits*1024*1024))
					else:
						h_hit_rate.append(0)
						h_bit_rate.append(0)						
										
					num_clusters = cl["time0"].size
					h_cluster_rate.pop(0)
					if num_clusters > 0:
						dt_cluster = max(cl["time0"]) - min(cl["time0"])
						h_cluster_rate.append(num_clusters*1000000000/dt_cluster)	
					else:
						h_cluster_rate.append(0)	
					h = plt.hist(cl['size0'], bins = 64, range = [0.0, 64], color='rebeccapurple')
					h_size0 = h_size0 + h[0]
					h = plt.hist(cl['size1'], bins = 64, range = [0.0, 64], color='rebeccapurple')
					h_size1 = h_size1 + h[0]
					h = plt.hist(cl['size0'] + cl['size1'], bins = 128, range = [0.0, 128], color='rebeccapurple')
					h_size = h_size + h[0]
					fig, ax = plt.subplots(2, 3, gridspec_kw={'height_ratios': [1, 2]})
					ax[0, 0].step(np.arange(0,64,1),h_size0, color='rebeccapurple')	
					ax[0, 1].step(np.arange(0,64,1),h_size1, color='rebeccapurple')			
					ax[0, 2].step(np.arange(0,128,1),h_size, color='rebeccapurple')	
					
					if firstTime == True:
						start_seconds = min(hits0["time"])/1000000000
						firstTime = False
					seconds = max(hits0["time"])/1000000000 - start_seconds
					h_time.pop(0)
					h_time.append(seconds)	
					ax[1, 0].plot(h_time,h_bit_rate, color='red')	
					ax[1, 1].plot(h_time,h_hit_rate, color='red')			
					ax[1, 2].plot(h_time,h_cluster_rate, color='red')	
					#print(h_time)
					ax[0, 0].title.set_text("clusters size0")
					ax[0, 1].title.set_text("clusters size1")
					ax[0, 2].title.set_text("clusters size0+size1")
					ax[0, 0].set_xlabel('size [strips]')
					ax[0, 0].set_ylabel('counts')
					ax[0, 1].set_xlabel('size [strips]')
					ax[0, 1].set_ylabel('counts')
					ax[0, 2].set_xlabel('size [strips]')
					ax[0, 2].set_ylabel('counts')
					ax[1, 0].title.set_text("bit rate")
					ax[1, 1].title.set_text("hit rate")
					ax[1, 2].title.set_text("clusters rate")		
					ax[1, 0].set_xlabel('time')
					ax[1, 0].set_ylabel('rate [Mbit/s]')
					ax[1, 1].set_xlabel('time')
					ax[1, 1].set_ylabel('rate [hits/s]')		
					ax[1, 2].set_xlabel('time')
					ax[1, 2].set_ylabel('rate [clusters/s]')			

					
					
					fig.set_size_inches(17, 8)
					fig.tight_layout()
					fig.savefig("det_stats_" + str(i) + ".png", format="png")
					
					
					fig, ax = plt.subplots(2, 4, gridspec_kw={'height_ratios': [1, 1]})
					h1 = ax[0, 0].hist(hits0['pos'], bins = channels_x, range = [0.0, channels_x], color='darkblue')
					h1_total = h1_total + h1[0]
					h2 = ax[0, 1].hist(hits1['pos'], bins = channels_y, range = [0.0, channels_y], color='darkblue')
					h2_total = h2_total + h2[0]		
					h3 = ax[0, 2].hist2d(hits0['pos'], hits0['adc'],bins =[channels_x, 128],cmap=plt.cm.jet,range=np.array([(0, channels_x), (0, 1024)]),norm=LogNorm())

					plt.colorbar(h3[3], ax=ax[0, 2], orientation='vertical')
					h3_total = h3_total + h3[0]
					h4 = ax[0, 3].hist2d(hits1['pos'], hits1['adc'],bins =[channels_y, 128],cmap=plt.cm.jet, range=np.array([(0, channels_y), (0,1024)]),norm=LogNorm())
					plt.colorbar(h4[3], ax=ax[0, 3], orientation='vertical')
					h4_total = h4_total + h4[0]
	
					
					h5 = ax[1, 0].hist(cl['pos0'], bins = channels_x, range = [0.0, channels_x], color='rebeccapurple')
					h5_total = h5_total + h5[0]
					h6 = ax[1, 1].hist(cl['pos1'], bins = channels_y, range = [0.0, channels_y], color='rebeccapurple')
					h6_total = h6_total + h6[0]
					h7 = ax[1, 2].hist(cl['adc0']+cl['adc1'], bins = 1000, range = [0.0, max_charge],  color='rebeccapurple')
					h7_total = h7_total + h7[0]
					h8 = ax[1, 3].hist2d(cl['pos0'], cl['pos1'],bins =[channels_x, channels_y],cmap=plt.cm.viridis,range=np.array([(0, channels_x), (0, channels_y)]))
					plt.colorbar(h8[3], ax=ax[1, 3], orientation='vertical')
					h8_total = h8_total + h8[0]
		
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

					fig.set_size_inches(17, 8)
					fig.tight_layout()
					fig.savefig("det_" + str(i) + ".png", format="png")

					
					
					fig, ax = plt.subplots(2, 4, gridspec_kw={'height_ratios': [1, 1]})
					ax[0, 0].step(np.arange(0,channels_x,1),h1_total, color='darkblue')	
					ax[0, 1].step(np.arange(0,channels_y,1),h2_total, color='darkblue')	
					xx = np.arange(0, channels_x+1, 1)  # len = 11
					yy = np.arange(0, 1032, 8)	
					im=ax[0, 2].pcolormesh(xx,yy,h3_total.transpose(), cmap=plt.cm.jet,norm=LogNorm())
					plt.colorbar(im, ax=ax[0, 2], orientation='vertical')

					xx = np.arange(0, channels_y+1, 1)  # len = 11
					yy = np.arange(0, 1032, 8)	
					im=ax[0, 3].pcolormesh(xx,yy,h4_total.transpose(), cmap=plt.cm.jet,norm=LogNorm())
					plt.colorbar(im, ax=ax[0, 3], orientation='vertical')
					
					ax[1, 0].step(np.arange(0,channels_x,1),h5_total, color='rebeccapurple')	
					ax[1, 1].step(np.arange(0,channels_y,1),h6_total, color='rebeccapurple')			
					ax[1, 2].step(np.arange(0,max_charge,10),h7_total, color='rebeccapurple')	
					im = ax[1, 3].imshow(h8_total.transpose(), cmap=plt.cm.viridis,origin='lower')
					plt.colorbar(im, ax=ax[1, 3], orientation='vertical')

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
		
					
					fig.set_size_inches(17, 8)
					fig.tight_layout()
					fig.savefig("det_total_" + str(i) + ".png", format="png")
					plt.close("all")
					
					time_end = int( t.time() )
					print("Plot " + str(cnt) + ": " + str(time_end-time_middle) + " s")
					cnt = cnt + 1


except OSError:
    pass
    
    
    
    
    
