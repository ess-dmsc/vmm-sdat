#!/usr/bin/python
import os,glob
import subprocess
import re
import sys
import uproot3 as uproot
import pandas as pd
import matplotlib.pyplot as plt
from timeit import default_timer as timer
import time as t
import numpy as np
from matplotlib.colors import LogNorm
from matplotlib import cm, ticker

#########################
#Edit the paramters below
file_name="tof"
#file_name="monitoring"
detector_id = 0 
channels_x = 256
channels_y = 256
max_charge = 10000
##########################
min_tof_s = 1e-5
max_tof_s = 10
tof_bins = 1000
min_energy_ev = 0.001
max_energy_ev = 1
energy_bins = 1000
min_wavelength = 0.01
max_wavelength = 10
wavelength_bins = 1000
min_thermal_tof = 0.0049*1.0e09
max_thermal_tof = 0.02*1.0e09


h1_total = [0] * channels_x    
h2_total = [0] * channels_y  
h3_total = np.zeros((channels_x, 128))
h4_total = np.zeros((channels_y, 128))
h5_total = [0] * channels_x    
h6_total = [0] * channels_y
h7_total = [0] * 1000
h8_total = np.zeros((channels_x, channels_y))

h_tof0 = [0] * tof_bins  
h_tof1 = [0] * tof_bins  
h_cl_tof = [0] * tof_bins 
h_energy0 = [0] * energy_bins  
h_energy1 = [0] * energy_bins  
h_cl_energy = [0] * energy_bins 
h_cl_wavelength = [0] * wavelength_bins 
h_thermal_image = np.zeros((channels_x, channels_y))

h_cluster_rate = [0] * 200  
h_hit_rate = [0] * 200 
h_time = [0] * 200 
h_percentage_x = [0] * 200 
h_percentage_y = [0] * 200 
h_size0 = [0] * 32
h_size1 = [0] * 32 
h_size01 = [0] * 64
h_deltaTime = [0] * 101 
h_max_missing_strip0 = [0] * 6 
h_max_missing_strip1 = [0] * 6 
h_span_cluster0 = [0] * 51 
h_span_cluster1 = [0] * 51 
h_max_delta_time0 = [0] * 51 
h_max_delta_time1 = [0] * 51 


cnt = 0
fileId=0
lastFileId = -1 

color_x = 'darkblue'
color_y ='darkgreen'
color_xy = 'rebeccapurple'
color_rate = 'red'
fig_w = 17
fig_h = 9
line_width = 3

start_time = timer()

def tofToEnergy(tof):
    distance = 19.4
    c = 299792458
    nmass = 939.56542052e6  # ev/c^2
    #correction = -5e-6
    correction = 0.0
    theTime = tof + correction

    gamma = np.sqrt(1.0 - pow(distance /theTime, 2) / (c * c))
    energy = (nmass / gamma - nmass)
    return energy

def tofToWavelength(tof):
    energy = tofToEnergy(tof)
	#energy in meV
    wavelength = 9.045 / np.sqrt(1000*energy)
    return wavelength

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
		args_vmmsdat = ['../build/convertFile', '-f', "./" + name, '-vmm', '[ [0,0,0,2],[0,0,0,3],[0,0,0,0],[0,0,0,1],[0,1,0,6],[0,1,0,7],[0,1,0,4],[0,1,0,5]]', '-axis', '[[0,0],1],[[0,1],1]', '-bc', '44.02', '-tac', '60', '-th','[0]', '-cs','[1]', '-ccs', '[2]', '-dt', '[200]', '-mst', '[1]', '-spc', '[500]', '-dp', '[500]', '-coin', 'center-of-masss', '-crl', '[0.1]', '-cru', '[10]','-save', '[[0],[0],[0]]', '-json','0', '-algo', '7', '-info', 'n_TOF', '-df', 'ESS', '-t0', '[[0,0,1710892800]]']
		#args_vmmsdat = ['../build/convertFile', '-f', "./" + name, '-vmm', '[[[0,0,0,8],[0,0,0,9],[0,0,0,6],[0,0,0,7],[0,0,0,4],[0,0,0,5],[0,0,0,2],[0,0,0,3],[0,0,0,0],[0,0,0,1],[0,1,1,0],[0,1,1,1],[0,1,1,2],[0,1,1,3],[0,1,1,4],[0,1,1,5],[0,1,1,6],[0,1,1,7],[0,1,1,8],[0,1,1,9]]','-axis', '[[0,0],1],[[0,1],0]', '-bc', '44.02', '-tac', '60', '-th','0', '-cs','1', '-ccs', '2', '-dt', '200', '-mst', '1', '-spc', '500', '-dp', '500', '-coin', 'center-of-masss', '-crl', '0.1', '-cru', '10', '-save', '[[0],[0],[0]]', '-json','0', '-algo', '0', '-info', 'monitor', '-df', 'ESS']
		
		subprocess.call(args_vmmsdat)
		lastFileId = fileId
		cnt = cnt + 1
		now_time = timer()

		print("Analysis " + str(cnt) + ": " + str(now_time - start_time) + " s")
		
		for file in os.listdir("."):
			if file.startswith(file_name + "_" + f'{fileId:05}') and file.endswith(".root"):
								
				####################################################################	
				# Hits
				####################################################################
				tree_hits = uproot.open(file)['hits']
				adc = tree_hits.array('adc')
				if adc.size == 0:
					continue
			
				pos = tree_hits.array('pos')
				time = tree_hits.array('time')
				tof = tree_hits.array('readout_time')
				det = tree_hits.array('det')
				plane = tree_hits.array('plane')
                
				data_hits = {'adc': adc,'pos': pos,'time': time, 'det': det,'plane': plane, 'tof': tof}
				df_hits = pd.DataFrame(data_hits)
				
				hits0 = df_hits.query("plane == 0 and det == " + str(detector_id))
				hits1 = df_hits.query("plane == 1 and det == " + str(detector_id))
						
				num_hits = hits0["time"].size +  hits1["time"].size
				if num_hits == 0:
					continue
		
				dt_hits = max(hits0["time"]) - min(hits0["time"])
				h_hit_rate.pop(0)
				if dt_hits > 0:
					h_hit_rate.append(num_hits*1000000.0/dt_hits)
				else:
					h_hit_rate.append(0)
				h_time.pop(0)
				h_time.append(now_time-start_time)
				
				####################################################################
				# Clusters
				####################################################################
				tree_detector = uproot.open(file)['clusters_detector']
				d_adc0 = tree_detector.array('adc0')
				d_pos0 = tree_detector.array('pos0')
				d_adc1 = tree_detector.array('adc1')
				d_pos1 = tree_detector.array('pos1')
				d_time0 = tree_detector.array('time0')
				d_time1 = tree_detector.array('time1')
				d_tof0 = tree_detector.array('time0_algo')
				d_tof1 = tree_detector.array('time1_algo')
				d_size0 = tree_detector.array('size0')
				d_size1 = tree_detector.array('size1')
				d_det = tree_detector.array('det')

				d_span_cluster0 = tree_detector.array('span_cluster0')
				d_span_cluster1 = tree_detector.array('span_cluster1')
				d_max_delta_time0 = tree_detector.array('max_delta_time0')
				d_max_delta_time1 = tree_detector.array('max_delta_time1')
				d_max_missing_strip0 = tree_detector.array('max_missing_strip0')
				d_max_missing_strip1 = tree_detector.array('max_missing_strip1')
				
				data_clusters = {'size0': d_size0,'size1': d_size1,'time0': d_time0,'time1': d_time1,'adc0': d_adc0,'pos0': d_pos0,'adc1': d_adc1,'pos1': d_pos1,'det': d_det, 'span_cluster0': d_span_cluster0,  'span_cluster1': d_span_cluster1,  'max_delta_time0': d_max_delta_time0,  'max_delta_time1': d_max_delta_time1,  'max_missing_strip0': d_max_missing_strip0,  'max_missing_strip1': d_max_missing_strip1, 'tof0': d_tof0, 'tof1': d_tof1}
				df_clusters = pd.DataFrame(data_clusters)
				
				thermal_clusters = {'pos0': d_pos0,'pos1': d_pos1,'det': d_det, 'tof0': d_tof0, 'tof1': d_tof1}
				df_thermal_clusters = pd.DataFrame(thermal_clusters)
				
				#Select clusters
				cl = df_clusters.query("det == " + str(detector_id))
				cl_thermal = df_thermal_clusters.query("det == " + str(detector_id) + " and tof0 >= " + str(min_thermal_tof) + " and tof0 <= " + str(max_thermal_tof))
				num_clusters = cl["time0"].size
				if num_clusters == 0:
					continue
				dt_cluster = max(cl["time0"]) - min(cl["time0"])
				h_cluster_rate.pop(0)
				if dt_cluster > 0:
					h_cluster_rate.append(num_clusters*1000000/dt_cluster)	
				else:
					h_cluster_rate.append(0)
						
				tree_plane= uproot.open(file)['clusters_plane']
				p_det = tree_plane.array('det')
				p_plane = tree_plane.array('plane')
				plane_clusters = {'det': p_det,'plane': p_plane}
				df_plane = pd.DataFrame(plane_clusters)
					 
				pl0 = df_plane.query("det == " + str(detector_id) + " and plane == 0")
				num_clusters_x = pl0["det"].size
				
				pl1 = df_plane.query("det == " + str(detector_id) + " and plane == 1")
				num_clusters_y = pl1["det"].size
				
				h_percentage_x.pop(0)
				h_percentage_y.pop(0)
				h_percentage_x.append(num_clusters*100/num_clusters_x)
				h_percentage_y.append(num_clusters*100/num_clusters_y)
	
				
                ####################################################################
				# ToF histograms
				####################################################################
				fig, ax = plt.subplots(2, 3, gridspec_kw={'height_ratios': [1,1]})
				h_temp = plt.hist(1e-09*hits0['tof'], bins = np.logspace(-5,1,1001))
				h_tof0 = h_tof0 + h_temp[0]
				ax[0, 0].stairs(h_tof0, h_temp[1],color=color_x)
				ax[0, 0].set_xscale("log", nonpositive='clip')
				ax[0, 0].set_yscale("log", nonpositive='clip')
				ax[0, 0].title.set_text("hits x ToF")
				ax[0, 0].set_xlabel('ToF [s]')
				ax[0, 0].set_ylabel('counts')

				h_temp = plt.hist(1e-09*hits1['tof'], bins = np.logspace(-5,1,1001))
				h_tof1 = h_tof1 + h_temp[0]
				ax[0, 1].stairs(h_tof1, h_temp[1],color=color_y)
				ax[0, 1].set_xscale("log", nonpositive='clip')
				ax[0, 1].set_yscale("log", nonpositive='clip')
				ax[0, 1].title.set_text("hits y ToF")
				ax[0, 1].set_xlabel('ToF [s]')
				ax[0, 1].set_ylabel('counts')

				h_temp = plt.hist(1e-09*cl["tof0"], bins = np.logspace(-5,1,1001))
				h_cl_tof = h_cl_tof + h_temp[0]
				ax[0, 2].stairs(h_cl_tof, h_temp[1],color=color_xy)
				ax[0, 2].set_xscale("log", nonpositive='clip')
				ax[0, 2].set_yscale("log", nonpositive='clip')
				ax[0, 2].title.set_text("clusters ToF")
				ax[0, 2].set_xlabel('ToF [s]')
				ax[0, 2].set_ylabel('counts')

				h_temp = plt.hist(tofToEnergy(1e-09*cl["tof0"]), bins = np.logspace(-4,1,1001))
				h_cl_energy = h_cl_energy + h_temp[0]
				ax[1, 0].stairs(h_cl_energy, h_temp[1],color=color_xy)
				ax[1, 0].set_xscale("log", nonpositive='clip')
				ax[1, 0].set_yscale("log", nonpositive='clip')
				ax[1, 0].title.set_text("clusters energy")
				ax[1, 0].set_xlabel('Energy [eV]')
				ax[1, 0].set_ylabel('counts')

				h_temp = plt.hist(tofToWavelength(1e-09*cl["tof0"]), bins = np.linspace(1.0,10,1001))
				h_cl_wavelength = h_cl_wavelength + h_temp[0]
				ax[1, 1].stairs(h_cl_wavelength, h_temp[1],color=color_xy)
				#ax[1, 1].set_xscale("log", nonpositive='clip')
				#ax[1, 1].set_yscale("log", nonpositive='clip')
				ax[1, 1].title.set_text("clusters wavelength")
				ax[1, 1].set_xlabel('Wavelength [A]')
				ax[1, 1].set_ylabel('counts')

				h_temp = plt.hist2d(cl_thermal['pos0'], cl_thermal['pos1'],bins =[channels_x, channels_y], cmap=plt.cm.magma_r,range=np.array([(0, channels_x), (0, channels_y)]))
				
				h_thermal_image = h_thermal_image + h_temp[0]
				im = ax[1, 2].imshow(h_thermal_image.transpose(), cmap=plt.cm.magma_r,origin='lower', norm=LogNorm())
				plt.colorbar(im, ax=ax[1, 2], orientation='vertical')
				ax[1, 2].title.set_text("thermal image")
				ax[1, 2].set_xlabel('x [0.4 mm]')
				ax[1, 2].set_ylabel('y [0.4 mm]')
				
				fig.set_size_inches(fig_w, fig_h)
				fig.tight_layout()
				fig.savefig("det_tof_" + str(detector_id) + ".png", format="png")
				
					
				####################################################################	
				# Stats histograms
				####################################################################
				h = plt.hist(cl['size0'], bins = 32, range = [-1.0, 31])
				h_size0 = h_size0 + h[0]
				h = plt.hist(cl['size1'], bins = 32, range = [-1.0, 31])
				h_size1 = h_size1 + h[0]
				h = plt.hist(cl['size0'] + cl['size1'], bins = 64, range = [-1.0, 63])
				h_size01 = h_size01 + h[0]
				h = plt.hist(cl['time1'] - cl['time0'], bins = 101, range = [-510, 500])
				h_deltaTime = h_deltaTime + h[0]
				
				h = plt.hist(cl['max_missing_strip0'], bins = 6, range = [-1.0, 5])
				h_max_missing_strip0 = h_max_missing_strip0 + h[0]
				h = plt.hist(cl['max_missing_strip1'], bins = 6, range = [-1.0, 5])
				h_max_missing_strip1 = h_max_missing_strip1 + h[0]
				h = plt.hist(cl['span_cluster0'], bins = 51, range = [-10, 500])
				h_span_cluster0 = h_span_cluster0 + h[0]
				h = plt.hist(cl['span_cluster1'], bins = 51, range = [-10, 500])
				h_span_cluster1 = h_span_cluster1 + h[0]
			
				h = plt.hist(cl['max_delta_time0'], bins = 51, range = [-10, 500])
				h_max_delta_time0 = h_max_delta_time0 + h[0]
				h = plt.hist(cl['max_delta_time1'], bins = 51, range = [-10, 500])
				h_max_delta_time1 = h_max_delta_time1 + h[0]
				
				####################################################################
				#Stats plots
				####################################################################
				fig, ax = plt.subplots(3, 4, gridspec_kw={'height_ratios': [1, 1,1]})
				ax[0, 0].step(np.arange(-0.5,31.5,1),h_size0, color=color_x, linewidth=line_width)
				ax[0, 1].step(np.arange(-0.5,31.5,1),h_size1, color=color_y, linewidth=line_width)
				ax[0, 2].step(np.arange(-0.5,63.5,1),h_size01, color=color_xy, linewidth=line_width)
				ax[0, 3].step(np.arange(-505,505,10),h_deltaTime, color=color_xy, linewidth=line_width)
					
				ax[1, 0].step(np.arange(-0.5,5.5,1),h_max_missing_strip0, color=color_x, linewidth=line_width)
				ax[1, 1].step(np.arange(-0.5,5.5,1),h_max_missing_strip1, color=color_y, linewidth=line_width)
				ax[1, 2].step(np.arange(-5,505,10),h_span_cluster0, color=color_x, linewidth=line_width)
				ax[1, 3].step(np.arange(-5,505,10),h_span_cluster1, color=color_y, linewidth=line_width)
				
				ax[2, 0].step(np.arange(-5,505,10),h_max_delta_time0, color=color_x, linewidth=line_width)
				ax[2, 1].step(np.arange(-5,505,10),h_max_delta_time1, color=color_y, linewidth=line_width)
				ax[2, 2].plot(h_time,h_percentage_x, color=color_x, linewidth=line_width)
				ax[2, 2].plot(h_time,h_percentage_y, color=color_y, linewidth=line_width)
				ax[2, 3].plot(h_time,h_hit_rate, color=color_rate, linewidth=line_width, linestyle='dashed')
				ax[2, 3].plot(h_time,h_cluster_rate, color=color_rate, linewidth=line_width)
			
				ax[0, 0].title.set_text("clusters size0")
				ax[0, 1].title.set_text("clusters size1")
				ax[0, 2].title.set_text("clusters size0+size1")
				ax[0, 3].title.set_text("delta time planes")
				ax[0, 0].set_xlabel('size [strips]')
				ax[0, 0].set_ylabel('counts')
				ax[0, 1].set_xlabel('size [strips]')
				ax[0, 1].set_ylabel('counts')
				ax[0, 2].set_xlabel('size [strips]')
				ax[0, 2].set_ylabel('counts')
				ax[0, 3].set_xlabel('time [ns]')
				ax[0, 3].set_ylabel('counts')
					
				ax[1, 0].title.set_text("max_missing_strip0")
				ax[1, 1].title.set_text("max_missing_strip1")
				ax[1, 2].title.set_text("span_cluster0")
				ax[1, 3].title.set_text("span_cluster1")
				ax[1, 0].set_xlabel('strips')
				ax[1, 0].set_ylabel('counts')
				ax[1, 1].set_xlabel('strips')
				ax[1, 1].set_ylabel('counts')
				ax[1, 2].set_xlabel('time [ns]')
				ax[1, 2].set_ylabel('counts')
				ax[1, 3].set_xlabel('time [ns]')
				ax[1, 3].set_ylabel('counts')
						
				ax[2, 0].title.set_text("max_delta_time0")
				ax[2, 1].title.set_text("max_delta_time1")
				ax[2, 2].title.set_text("common clusters: plane 0 (blue), plane 1 (green)")
				ax[2, 3].title.set_text("rate: hit (dashed), cluster (solid)")
				ax[2, 0].set_xlabel('ns')
				ax[2, 0].set_ylabel('counts')
				ax[2, 1].set_xlabel('time [ns]')
				ax[2, 1].set_ylabel('counts')
				ax[2, 2].set_xlabel('time since start of acq [s]')
				ax[2, 2].set_ylabel('percent [%]')	
				ax[2, 3].set_xlabel('time since start of acq [s]')
				ax[2, 3].set_ylabel('rate [kHz]')
												
				fig.set_size_inches(fig_w, fig_h)
				fig.tight_layout()
				fig.savefig("det_stats_" + str(detector_id) + ".png", format="png")
					

				####################################################################
				# Latest histograms
				####################################################################
				fig, ax = plt.subplots(2, 4, gridspec_kw={'height_ratios': [1, 1]})
				h1 = ax[0, 0].hist(hits0['pos'], bins = channels_x, range = [0.0, channels_x], color=color_x)
				h1_total = h1_total + h1[0]
				h2 = ax[0, 1].hist(hits1['pos'], bins = channels_y, range = [0.0, channels_y], color=color_y)
				h2_total = h2_total + h2[0]
				h3 = ax[0, 2].hist2d(hits0['pos'], hits0['adc'],bins =[channels_x, 128],cmap=plt.cm.jet,range=np.array([(0, channels_x), (0, 1024)]),norm=LogNorm())

				plt.colorbar(h3[3], ax=ax[0, 2], orientation='vertical')
				h3_total = h3_total + h3[0]
				h4 = ax[0, 3].hist2d(hits1['pos'], hits1['adc'],bins =[channels_y, 128],cmap=plt.cm.jet, range=np.array([(0, channels_y), (0,1024)]),norm=LogNorm())
				plt.colorbar(h4[3], ax=ax[0, 3], orientation='vertical')
				h4_total = h4_total + h4[0]

				
				h5 = ax[1, 0].hist(cl['pos0'], bins = channels_x, range = [0.0, channels_x], color=color_x)
				h5_total = h5_total + h5[0]
				h6 = ax[1, 1].hist(cl['pos1'], bins = channels_y, range = [0.0, channels_y], color=color_y)
				h6_total = h6_total + h6[0]
				h7 = ax[1, 2].hist(cl['adc0']+cl['adc1'], bins = 1000, range = [0.0, max_charge],  color=color_xy)
				h7_total = h7_total + h7[0]
				h8 = ax[1, 3].hist2d(cl['pos0'], cl['pos1'],bins =[channels_x, channels_y],cmap=plt.cm.magma_r,range=np.array([(0, channels_x), (0, channels_y)]))
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

				fig.set_size_inches(fig_w, fig_h)
				fig.tight_layout()
				fig.savefig("det_" + str(detector_id) + ".png", format="png")

					
				####################################################################
				# Total histograms
				####################################################################
				fig, ax = plt.subplots(2, 4, gridspec_kw={'height_ratios': [1, 1]})
				ax[0, 0].step(np.arange(0,channels_x,1),h1_total, color=color_x)
				ax[0, 1].step(np.arange(0,channels_y,1),h2_total, color=color_y)
				xx = np.arange(0, channels_x+1, 1)  # len = 11
				yy = np.arange(0, 1032, 8)
				im=ax[0, 2].pcolormesh(xx,yy,h3_total.transpose(), cmap=plt.cm.jet,norm=LogNorm())
				plt.colorbar(im, ax=ax[0, 2], orientation='vertical')

				xx = np.arange(0, channels_y+1, 1)  # len = 11
				yy = np.arange(0, 1032, 8)
				im=ax[0, 3].pcolormesh(xx,yy,h4_total.transpose(), cmap=plt.cm.jet,norm=LogNorm())
				plt.colorbar(im, ax=ax[0, 3], orientation='vertical')
				
				ax[1, 0].step(np.arange(0,channels_x,1),h5_total, color=color_x)
				ax[1, 1].step(np.arange(0,channels_y,1),h6_total, color=color_y)
				ax[1, 2].step(np.arange(0,max_charge,10),h7_total, color=color_xy)
				im = ax[1, 3].imshow(h8_total.transpose(), cmap=plt.cm.magma_r,origin='lower')
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
	
				
				fig.set_size_inches(fig_w, fig_h)
				fig.tight_layout()
				fig.savefig("det_total_" + str(detector_id) + ".png", format="png")
				
				####################################################################
				# clean up
				####################################################################
				plt.close("all")
				end_time = timer()
				print("Plot " + str(cnt) + ": " + str(end_time-start_time) + " s")
				cnt = cnt + 1


except OSError:
    pass
    
    
    
    
    
