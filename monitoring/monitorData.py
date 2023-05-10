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
import gc
from plotly.subplots import make_subplots
import plotly.graph_objects as go

#########################
#Edit the paramters below
file_name="example"
#file_name="monitoring"
detector_id = 0 
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
h7_total = [0] * int(max_charge*0.1)
h8_total = np.zeros((channels_x, channels_y))    

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

color_x = 'blue'
color_y ='green'
color_xy = 'purple'
color_rate = 'red'
fig_w = 17
fig_h = 10
line_width = 3

start_time = timer()

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
                args_vmmsdat = ['../build/convertFile', '-f', "./" + name, '-geo', 'example_monitoring.json', '-bc', '40', '-tac', '60', '-th','0', '-cs','1', '-ccs', '3', '-dt', '200', '-mst', '1', '-spc', '500', '-dp', '200', '-coin', 'center-of-masss', '-crl', '0.2', '-cru', '10', '-save', '[[0],[0],[0]]', '-algo', '4', '-info', '', '-df','SRS']
                #args_vmmsdat = ['../build/convertFile', '-f', "./" + name, '-vmm', '[[[0,0,0,8],[0,0,0,9],[0,0,0,6],[0,0,0,7],[0,0,0,4],[0,0,0,5],[0,0,0,2],[0,0,0,3],[0,0,0,0],[0,0,0,1],[0,1,1,0],[0,1,1,1],[0,1,1,2],[0,1,1,3],[0,1,1,4],[0,1,1,5],[0,1,1,6],[0,1,1,7],[0,1,1,8],[0,1,1,9]]','-axis', '[[0,0],1],[[0,1],0]', '-bc', '44.02', '-tac', '60', '-th','0', '-cs','1', '-ccs', '2', '-dt', '200', '-mst', '1', '-spc', '500', '-dp', '500', '-coin', 'center-of-masss', '-crl', '0.1', '-cru', '10', '-save', '[[0],[0],[0]]', '-json','0', '-algo', '0', '-info', 'monitor', '-df', 'ESS']                       
                
                subprocess.call(args_vmmsdat)
                lastFileId = fileId
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
                                det = tree_hits.array('det')
                                plane = tree_hits.array('plane')
                                data_hits = {'adc': adc,'pos': pos,'time': time, 'det': det,'plane': plane}
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
                                d_size0 = tree_detector.array('size0')
                                d_size1 = tree_detector.array('size1')
                                d_det = tree_detector.array('det')

                                d_span_cluster0 = tree_detector.array('span_cluster0')
                                d_span_cluster1 = tree_detector.array('span_cluster1')
                                d_max_delta_time0 = tree_detector.array('max_delta_time0')
                                d_max_delta_time1 = tree_detector.array('max_delta_time1')
                                d_max_missing_strip0 = tree_detector.array('max_missing_strip0')
                                d_max_missing_strip1 = tree_detector.array('max_missing_strip1')
                                
                                data_clusters = {'size0': d_size0,'size1': d_size1,'time0': d_time0,'time1': d_time1,'adc0': d_adc0,'pos0': d_pos0,'adc1': d_adc1,'pos1': d_pos1,'det': d_det, 'span_cluster0': d_span_cluster0,  'span_cluster1': d_span_cluster1,  'max_delta_time0': d_max_delta_time0,  'max_delta_time1': d_max_delta_time1,  'max_missing_strip0': d_max_missing_strip0,  'max_missing_strip1': d_max_missing_strip1}
                                df_clusters = pd.DataFrame(data_clusters)       

                                #Select clusters
                                cl = df_clusters.query("det == " + str(detector_id))
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
                                # Stats histograms
                                ####################################################################
                                h = np.histogram(cl['size0'], bins = 32, range = [-1.0, 31])
                                h_size0 = h_size0 + h[0]
                                h = np.histogram(cl['size1'], bins = 32, range = [-1.0, 31])
                                h_size1 = h_size1 + h[0]
                                h = np.histogram(cl['size0'] + cl['size1'], bins = 64, range = [-1.0, 63])
                                h_size01 = h_size01 + h[0]
                                h = np.histogram(cl['time1'] - cl['time0'], bins = 101, range = [-510, 500])
                                h_deltaTime = h_deltaTime + h[0]
                                
                                h = np.histogram(cl['max_missing_strip0'], bins = 6, range = [-1.0, 5])
                                h_max_missing_strip0 = h_max_missing_strip0 + h[0]
                                h = np.histogram(cl['max_missing_strip1'], bins = 6, range = [-1.0, 5])
                                h_max_missing_strip1 = h_max_missing_strip1 + h[0]
                                h = np.histogram(cl['span_cluster0'], bins = 51, range = [-10, 500])
                                h_span_cluster0 = h_span_cluster0 + h[0]
                                h = np.histogram(cl['span_cluster1'], bins = 51, range = [-10, 500])
                                h_span_cluster1 = h_span_cluster1 + h[0]                                        
                                
                                h = np.histogram(cl['max_delta_time0'], bins = 51, range = [-10, 500])
                                h_max_delta_time0 = h_max_delta_time0 + h[0]
                                h = np.histogram(cl['max_delta_time1'], bins = 51, range = [-10, 500])
                                h_max_delta_time1 = h_max_delta_time1 + h[0]
                                
                                ####################################################################    
                                #Stats plots
                                ####################################################################            
                                fig = make_subplots(rows = 3, cols = 4, subplot_titles = ("cluster size0", "clusters size1", "clusters size0+size1", "delta time planes",
                                                                                          "max_missing_strip0", "max_missing_strip1", "span_cluster0", "span_cluster1",
                                                                                          "max_delta_time0", "max_delta_time1", "common clusters: plane 0 (dashed), plane 1 (solid)", "rate: hit (dashed), cluster (solid)"))

                                fig.add_trace(go.Bar(x = np.arange(-0.5,31.5,1), y = h_size0, marker_color = color_x), row = 1, col = 1)
                                fig.add_trace(go.Bar(x = np.arange(-0.5,31.5,1), y = h_size1, marker_color = color_y), row = 1, col = 2)
                                fig.add_trace(go.Bar(x = np.arange(-0.5,63.5,1), y = h_size01, marker_color = color_xy), row = 1, col = 3)
                                fig.add_trace(go.Bar(x = np.arange(-505,505,10), y = h_deltaTime, marker_color = color_xy), row = 1, col = 4)

                                fig.add_trace(go.Bar(x = np.arange(-0.5,5.5,1), y = h_max_missing_strip0, marker_color = color_x), row = 2, col = 1)
                                fig.add_trace(go.Bar(x = np.arange(-0.5,5.5,1), y = h_max_missing_strip1, marker_color = color_y), row = 2, col = 2)
                                fig.add_trace(go.Bar(x = np.arange(-505,505,10), y = h_span_cluster0, marker_color = color_x), row = 2, col = 3)
                                fig.add_trace(go.Bar(x = np.arange(-505,505,10), y = h_span_cluster1, marker_color = color_y), row = 2, col = 4)

                                fig.add_trace(go.Bar(x = np.arange(-0.5,5.5,1), y = h_max_delta_time0, marker_color = color_x), row = 3, col = 1)
                                fig.add_trace(go.Bar(x = np.arange(-0.5,5.5,1), y = h_max_delta_time1, marker_color = color_y), row = 3, col = 2)
                                fig.add_trace(go.Scatter(x = h_time, y = h_percentage_x, marker_color = color_x, line = dict(dash = 'dot')), row = 3, col = 3)
                                fig.add_trace(go.Scatter(x = h_time, y = h_percentage_y, marker_color = color_y), row = 3, col = 3)
                                fig.add_trace(go.Scatter(x = h_time, y = h_hit_rate, marker_color = color_rate, line = dict(dash = 'dot')), row = 3, col = 4)
                                fig.add_trace(go.Scatter(x = h_time, y = h_cluster_rate, marker_color = color_rate), row = 3, col = 4)

                                fig.update_xaxes(title_text="size [strips]", row = 1, col = 1)
                                fig.update_xaxes(title_text="size [strips]", row = 1, col = 2)
                                fig.update_xaxes(title_text="size [strips]", row = 1, col = 3)
                                fig.update_xaxes(title_text="time [ns]", row = 1, col = 4)

                                fig.update_xaxes(title_text="strips", row = 2, col = 1)
                                fig.update_xaxes(title_text="strips", row = 2, col = 2)
                                fig.update_xaxes(title_text="time [ns]", row = 2, col = 3)
                                fig.update_xaxes(title_text="time [ns]", row = 2, col = 4)

                                fig.update_xaxes(title_text="ns", row = 3, col = 1)
                                fig.update_xaxes(title_text="time [ns]", row = 3, col = 2)
                                fig.update_xaxes(title_text="time since start of acq [s]", row = 3, col = 3)
                                fig.update_xaxes(title_text="time since start of acq [s]", row = 3, col = 4)

                                fig.update_yaxes(title_text="counts", row = 1, col = 1)
                                fig.update_yaxes(title_text="counts", row = 1, col = 2)
                                fig.update_yaxes(title_text="counts", row = 1, col = 3)
                                fig.update_yaxes(title_text="counts", row = 1, col = 4)

                                fig.update_yaxes(title_text="counts", row = 2, col = 1)
                                fig.update_yaxes(title_text="counts", row = 2, col = 2)
                                fig.update_yaxes(title_text="counts", row = 2, col = 3)
                                fig.update_yaxes(title_text="counts", row = 2, col = 4)

                                fig.update_yaxes(title_text="counts", row = 3, col = 1)
                                fig.update_yaxes(title_text="counts", row = 3, col = 2)
                                fig.update_yaxes(title_text="%", row = 3, col = 3)
                                fig.update_yaxes(title_text="rate [kHz]", row = 3, col = 4)

                                fig.update_layout(height = 1100, width = 1700, showlegend = False)
                                fig.write_html("det_stats_" + str(detector_id) + ".html")
                                fig.data = []
                                

                                ####################################################################    
                                # Latest histograms
                                ####################################################################
                                fig = make_subplots(rows = 2, cols = 4, subplot_titles = ("hits pos0", "hits pos1", "hits adc0", "hits adc1",
                                                                                          "clusters pos0", "clusters pos1", "clusters adc0+adc1", "clusters pos1:pos0"))

                                h1 = np.histogram(hits0['pos'], bins = channels_x, range = [0.0, channels_x])
                                h1_total = h1_total + h1[0]
                                h2 = np.histogram(hits1['pos'], bins = channels_y, range = [0.0, channels_y])
                                h2_total = h2_total + h2[0]             
                                h3 = np.histogram2d(hits0['pos'], hits0['adc'], bins = [channels_x, 128], range = np.array([(0, channels_x), (0, 1024)]))
                                h3_total = h3_total + h3[0]
                                h4 = np.histogram2d(hits1['pos'], hits1['adc'], bins = [channels_y, 128], range = np.array([(0, channels_y), (0,1024)]))
                                h4_total = h4_total + h4[0]

                                
                                h5 = np.histogram(cl['pos0'], bins = channels_x, range = [0.0, channels_x])
                                h5_total = h5_total + h5[0]
                                h6 = np.histogram(cl['pos1'], bins = channels_y, range = [0.0, channels_y])
                                h6_total = h6_total + h6[0]
                                h7 = np.histogram(cl['adc0']+cl['adc1'], bins = int(max_charge*0.1), range = [0.0, max_charge])
                                h7_total = h7_total + h7[0]
                                h8 = np.histogram2d(cl['pos0'], cl['pos1'], bins = [channels_x, channels_y], range = np.array([(0, channels_x), (0, channels_y)]))
                                h8_total = h8_total + h8[0]

                                fig.add_trace(go.Bar(x = h1[1][:-1], y = h1[0], marker_color = color_x), row = 1, col = 1)
                                fig.add_trace(go.Bar(x = h2[1][:-1], y = h2[0], marker_color = color_y), row = 1, col = 2)
                                fig.add_trace(go.Heatmap(z = np.transpose(h3[0]), colorbar_x = 0.63, colorbar_y = 0.5, colorbar_len = 0.2), row = 1, col = 3)
                                fig.add_trace(go.Heatmap(z = np.transpose(h4[0]), colorbar_x = 0.89, colorbar_y = 0.5, colorbar_len = 0.2), row = 1, col = 4)

                                fig.add_trace(go.Bar(x = h5[1][:-1], y = h5[0], marker_color = color_x), row = 2, col = 1)
                                fig.add_trace(go.Bar(x = h6[1][:-1], y = h6[0], marker_color = color_y), row = 2, col = 2)
                                fig.add_trace(go.Bar(x = h7[1][:-1], y = h7[0], marker_color = color_xy), row = 2, col = 3)
                                fig.add_trace(go.Heatmap(z = np.transpose(h8[0]), colorbar_x = 0.89, colorbar_y = -0.13, colorbar_len = 0.2), row = 2, col = 4)
                                
                                fig.update_xaxes(title_text=("x [pitch 0.4 mm]"), row = 1, col = 1)
                                fig.update_yaxes(title_text=("counts"), row = 1, col = 1)
                                fig.update_xaxes(title_text=("y [pitch 0.4 mm]"), row = 1, col = 2)
                                fig.update_yaxes(title_text=("counts"), row = 1, col = 2)
                                fig.update_xaxes(title_text=("x [pitch 0.4 mm]"), row = 1, col = 3)
                                fig.update_yaxes(title_text=("adc"), row = 1, col = 3)
                                fig.update_xaxes(title_text=("y [pitch 0.4 mm]"), row = 1, col = 4)
                                fig.update_yaxes(title_text=("adc"), row = 1, col = 4)
                                
                                fig.update_xaxes(title_text=("x [pitch 0.4 mm]"), row = 2, col = 1)
                                fig.update_yaxes(title_text=("counts"), row = 2, col = 2)
                                fig.update_xaxes(title_text=("y [pitch 0.4 mm]"), row = 2, col = 2)
                                fig.update_yaxes(title_text=("counts"), row = 2, col = 2)
                                fig.update_xaxes(title_text=("charge"), row = 2, col = 3)
                                fig.update_yaxes(title_text=("counts"), row = 2, col = 3)
                                fig.update_xaxes(title_text=("x [pitch 0.4 mm]"), row = 2, col = 4)
                                fig.update_yaxes(title_text=("y [pitch 0.4 mm]"), row = 2, col = 4)

                                fig.update_traces(colorbar_orientation='h', selector=dict(type='heatmap'))
                                fig.update_layout(height = 1100, width = 1700, showlegend = False)
                                fig.write_html("det_" + str(detector_id) + ".html")
                                fig.data = []
                                
                                ####################################################################    
                                # Total histograms
                                ####################################################################
                                fig = make_subplots(rows = 2, cols = 4, subplot_titles = ("hits pos0", "hits pos1", "hits adc0", "hits adc1",
                                                                                          "clusters pos0", "clusters pos1", "clusters adc0+adc1", "clusters pos1:pos0"))
                                
                                fig.add_trace(go.Bar(x = h1[1][:-1], y = h1_total, marker_color = color_x), row = 1, col = 1)
                                fig.add_trace(go.Bar(x = h2[1][:-1], y = h2_total, marker_color = color_y), row = 1, col = 2)
                                fig.add_trace(go.Heatmap(z = np.transpose(h3_total), colorbar_x = 0.63, colorbar_y = 0.5, colorbar_len = 0.2), row = 1, col = 3)
                                fig.add_trace(go.Heatmap(z = np.transpose(h4_total), colorbar_x = 0.89, colorbar_y = 0.5, colorbar_len = 0.2), row = 1, col = 4)

                                fig.add_trace(go.Bar(x = h5[1][:-1], y = h5_total, marker_color = color_x), row = 2, col = 1)
                                fig.add_trace(go.Bar(x = h6[1][:-1], y = h6_total, marker_color = color_y), row = 2, col = 2)
                                fig.add_trace(go.Bar(x = h7[1][:-1], y = h7_total, marker_color = color_xy), row = 2, col = 3)
                                fig.add_trace(go.Heatmap(z = np.transpose(h8_total), colorbar_x = 0.89, colorbar_y = -0.13, colorbar_len = 0.2), row = 2, col = 4)
                                
                                fig.update_xaxes(title_text=("x [pitch 0.4 mm]"), row = 1, col = 1)
                                fig.update_yaxes(title_text=("counts"), row = 1, col = 1)
                                fig.update_xaxes(title_text=("y [pitch 0.4 mm]"), row = 1, col = 2)
                                fig.update_yaxes(title_text=("counts"), row = 1, col = 2)
                                fig.update_xaxes(title_text=("x [pitch 0.4 mm]"), row = 1, col = 3)
                                fig.update_yaxes(title_text=("adc"), row = 1, col = 3)
                                fig.update_xaxes(title_text=("y [pitch 0.4 mm]"), row = 1, col = 4)
                                fig.update_yaxes(title_text=("adc"), row = 1, col = 4)
                                
                                fig.update_xaxes(title_text=("x [pitch 0.4 mm]"), row = 2, col = 1)
                                fig.update_yaxes(title_text=("counts"), row = 2, col = 2)
                                fig.update_xaxes(title_text=("y [pitch 0.4 mm]"), row = 2, col = 2)
                                fig.update_yaxes(title_text=("counts"), row = 2, col = 2)
                                fig.update_xaxes(title_text=("charge"), row = 2, col = 3)
                                fig.update_yaxes(title_text=("counts"), row = 2, col = 3)
                                fig.update_xaxes(title_text=("x [pitch 0.4 mm]"), row = 2, col = 4)
                                fig.update_yaxes(title_text=("y [pitch 0.4 mm]"), row = 2, col = 4)

                                fig.update_traces(colorbar_orientation='h', selector=dict(type='heatmap'))
                                fig.update_layout(height = 1100, width = 1700, showlegend = False)
                                fig.write_html("det_total_" + str(detector_id) + ".html")
                                fig.data = []
                                
                                ####################################################################    
                                # clean up
                                ####################################################################                    
                                gc.collect()

                                end_time = timer()
                                print("Plot " + str(cnt) + ": " + str(end_time-start_time) + " s")
                                cnt = cnt + 1


except OSError:
    pass
    
    
    
    
    
