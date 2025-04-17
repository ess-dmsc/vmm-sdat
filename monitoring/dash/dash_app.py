import os
import pandas as pd
import uproot
import time
from dash import Dash, html, dcc, page_container, Input, Output, State, no_update
from config import *

def create_dash_app():
	app = Dash(__name__, use_pages=True)
	app.title = "NMX monitoring"

	app.layout = html.Div([
		dcc.Location(id="url"),  # Required for use_pages
		dcc.Interval(id='refresh', interval=3000, n_intervals=0),
		dcc.Store(id='hit_data'),
		dcc.Store(id='plane_data'),
		dcc.Store(id='cluster_data'),
		page_container
	])
	return app


def register_callbacks(app):

	#print("DASH register_callbacks")		
	@app.callback(
		Output('hit_data', 'data'),
		Output('plane_data', 'data'),
		Output('cluster_data', 'data'),
		Input('refresh', 'n_intervals'),
		State("url", "pathname")
	)
	def send_data(n,pathname):
		try:
			print(pathname)
			if pathname == "/hits":
				if delete_root > 0:
					cleanup_old_root_files(".", max_files=10)
				root_files = count_files("*.root", file_name)
				if not root_files:
					raise ValueError("Waiting for root files")
				oldest_file = get_oldest_file(root_files)
				if not oldest_file:
					raise ValueError("No oldest root file")
				
				#print(f"Plotting {oldest_file}")
				tree_hits = uproot.open(oldest_file)['hits']
				tree_plane = uproot.open(oldest_file)['clusters_plane']
				tree_clusters = uproot.open(oldest_file)['clusters_detector']

				df_hits = tree_hits.arrays(['adc', 'ch', 'pos',  'vmm','plane', 'det'], library='pd')
				df_plane = tree_plane.arrays(['plane', 'det'], library='pd')
				if cluster_algorithm == "utpc":
					df_clusters = tree_clusters.arrays([
						'adc0', 'pos0_utpc', 'time0', 'size0', 'span_cluster0', 'max_delta_time0', 'max_missing_strip0',
						'adc1', 'pos1_utpc', 'time1','size1', 'span_cluster1', 'max_delta_time1', 'max_missing_strip1', 'det'
					], library='pd')
				elif cluster_algorithm == "charge2":
					df_clusters = tree_clusters.arrays([
						'adc0', 'pos0_charge2', 'time0', 'size0', 'span_cluster0', 'max_delta_time0', 'max_missing_strip0',
						'adc1', 'pos1_charge2', 'time1','size1', 'span_cluster1', 'max_delta_time1', 'max_missing_strip1', 'det'
					], library='pd')
				else:
					df_clusters = tree_clusters.arrays([
						'adc0', 'pos0', 'time0', 'size0', 'span_cluster0', 'max_delta_time0', 'max_missing_strip0',
						'adc1', 'pos1', 'time1','size1', 'span_cluster1', 'max_delta_time1', 'max_missing_strip1', 'det'
					], library='pd')			
				shared_data["hit_data"] = df_hits.to_dict('records')
				shared_data["plane_data"] = df_plane.to_dict('records')
				shared_data["cluster_data"] = df_clusters.to_dict('records')
				now = time.time()
				shared_data["hit_updated_at"] = now
				shared_data["plane_updated_at"] = now
				shared_data["cluster_updated_at"] = now
				if delete_root > 0:
					os.remove(oldest_file)
			return (
				shared_data.get("hit_data", []),
				shared_data.get("plane_data", []),
				shared_data.get("cluster_data", [])
			)
		except Exception as e:
			print(f"Message: {e}")
			return (
				shared_data.get("hit_data", []),
				shared_data.get("plane_data", []),
				shared_data.get("cluster_data", [])
			)