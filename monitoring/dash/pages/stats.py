import dash
from dash import html, dcc, Input, Output, State
import pandas as pd
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import plotly.io as pio
import plotly.express as px
import plotly.colors as pc
import numpy as np
from config import *
from timeit import default_timer as timer

last_seen_cluster_timestamp = 0.0
last_seen_plane_timestamp = 0.0

dash.register_page(__name__, path='/stats', name="Statistics")

layout = html.Div([
	html.Div([
		html.Div([
			# Log y-axis
			dcc.Checklist(
				id='logy_toggle',
				options=[{'label': '1D plots: y-axis log', 'value': 'logy'}],
				value=[],
				labelStyle={'display': 'inline-block', 'margin-right': '10px'}
			),
			html.Label("Color palette:", style={"margin-left": "20px", "margin-right": "2px"}),
			# Color palette dropdown
			dcc.Dropdown(
				id="color_palette",
				options=color_options,
				value="Config",
				clearable=False,
				style={"width": "200px", "margin-left": "5px", "verticalAlign": "middle"}
			),
			html.Div(id="color_preview", style={"display": "inline-block", "margin-left": "5px", "margin-right": "20px"}),
			], style={
			"display": "flex",
			"alignItems": "center",  # vertically align
			"justifyContent": "center",
			"margin-bottom": "20px"
		})
	]),
	dcc.Graph(id='statistics_graph'),
	dcc.Store(id='h_totals_stats')
])


	
@dash.callback(
	Output('statistics_graph', 'figure'),
	Output('h_totals_stats', 'data'),
	Input('plane_data', 'data'),
	Input('cluster_data', 'data'),
	Input('logy_toggle', 'value'),
	Input("color_palette", "value"),
	State('h_totals_stats', 'data')
)
def plot_data(plane_data,cluster_data, logy_toggle, color_palette, h_totals_stats):
	global last_seen_plane_timestamp
	global last_seen_cluster_timestamp
	h_totals_stats = h_totals_stats or {}
	cl_updated_at = shared_data.get("cluster_updated_at", 0.0)
	pl_updated_at = shared_data.get("plane_updated_at", 0.0)
	if not cluster_data or cl_updated_at == last_seen_cluster_timestamp or not plane_data or pl_updated_at == last_seen_plane_timestamp:
		return dash.no_update, h_totals_stats
	
	if color_palette == "Config":
		colors = color_config
	else:
		colors = palette_map.get(color_palette, px.colors.qualitative.Plotly)
	yaxis_type = 'log' if 'logy' in logy_toggle else 'linear'
		
	last_seen_plane_timestamp = pl_updated_at
	last_seen_cluster_timestamp = cl_updated_at
	df_plane = pd.DataFrame(plane_data)
	p00 = df_plane.query("plane == 0 and det == 0")
	p10 = df_plane.query("plane == 1 and det == 0")
	p01 = df_plane.query("plane == 0 and det == 1")
	p11 = df_plane.query("plane == 1 and det == 1")
	p02 = df_plane.query("plane == 0 and det == 2")
	p12 = df_plane.query("plane == 1 and det == 2")
	p03 = df_plane.query("plane == 0 and det == 3")
	p13 = df_plane.query("plane == 1 and det == 3")
		 	
	df_clusters = pd.DataFrame(cluster_data) 
	d0 = df_clusters.query("det == 0")
	d1 = df_clusters.query("det == 1")
	d2 = df_clusters.query("det == 2")
	d3 = df_clusters.query("det == 3")
	
	yaxis_type = 'log' if 'logy' in logy_toggle else 'linear'
	
	now_time = timer()
	h_start_time = np.array(h_totals_stats.get("h_start_time", [0]*1))
	if h_start_time[0] == 0:
		h_start_time[0] = now_time
	
	h_time_total = np.array(h_totals_stats.get("h_time_total", [0]*time_points))
	h_time_total = np.roll(h_time_total, -1)
	h_time_total[-1] =now_time-h_start_time[0]


	# size 0
	h1_0 = np.histogram(d0["size0"], bins = bins_size, range = [0.0, bins_size])
	h1_0_total = np.array(h_totals_stats.get("h1_0_total", [0]*bins_size))
	h1_0_total = h1_0_total + h1_0[0]
	
	
	h1_1 = np.histogram(d1["size0"], bins = bins_size, range = [0.0, bins_size])
	h1_1_total = np.array(h_totals_stats.get("h1_1_total", [0]*bins_size))
	h1_1_total = h1_1_total + h1_1[0]
	#y_step_1_1 = np.repeat(h1_1_total, 2)
	
	h1_2 = np.histogram(d2["size0"], bins = bins_size, range = [0.0, bins_size])
	h1_2_total = np.array(h_totals_stats.get("h1_2_total", [0]*bins_size))
	h1_2_total = h1_2_total + h1_2[0]
	#y_step_1_2 = np.repeat(h1_2_total, 2)
	
	h1_3 = np.histogram(d3["size0"], bins = bins_size, range = [0.0, bins_size])
	h1_3_total = np.array(h_totals_stats.get("h1_3_total", [0]*bins_size))
	h1_3_total = h1_3_total + h1_3[0]
	#y_step_1_3 = np.repeat(h1_3_total, 2)
	
	
	# size 1
	h2_0 = np.histogram(d0["size1"], bins = bins_size, range = [0.0, bins_size])
	h2_0_total = np.array(h_totals_stats.get("h2_0_total", [0]*bins_size))
	h2_0_total = h2_0_total + h2_0[0]
	
	h2_1 = np.histogram(d1["size1"], bins = bins_size, range = [0.0, bins_size])
	h2_1_total = np.array(h_totals_stats.get("h2_1_total", [0]*bins_size))
	h2_1_total = h2_1_total + h2_1[0]
	
	h2_2 = np.histogram(d2["size1"], bins = bins_size, range = [0.0, bins_size])
	h2_2_total = np.array(h_totals_stats.get("h2_2_total", [0]*bins_size))
	h2_2_total = h2_2_total + h2_2[0]
	
	h2_3 = np.histogram(d3["size1"], bins = bins_size, range = [0.0, bins_size])
	h2_3_total = np.array(h_totals_stats.get("h2_3_total", [0]*bins_size))
	h2_3_total = h2_3_total + h2_3[0]
	
	
	# size 0 + size 1
	h3_0 = np.histogram(d0["size0"]+d0["size1"], bins = 2*bins_size, range = [0.0, 2*bins_size])
	h3_0_total = np.array(h_totals_stats.get("h3_0_total", [0]*2*bins_size))
	h3_0_total = h3_0_total + h3_0[0]

	h3_1 = np.histogram(d1["size0"]+d1["size1"], bins = 2*bins_size, range = [0.0, 2*bins_size])
	h3_1_total = np.array(h_totals_stats.get("h3_1_total", [0]*2*bins_size))
	h3_1_total = h3_1_total + h3_1[0]
	
	h3_2 = np.histogram(d2["size0"]+d2["size1"], bins = 2*bins_size, range = [0.0, 2*bins_size])
	h3_2_total = np.array(h_totals_stats.get("h3_2_total", [0]*2*bins_size))
	h3_2_total = h3_2_total + h3_2[0]
	
	h3_3 = np.histogram(d3["size0"]+d3["size1"], bins = 2*bins_size, range = [0.0, 2*bins_size])
	h3_3_total = np.array(h_totals_stats.get("h3_3_total", [0]*2*bins_size))
	h3_3_total = h3_3_total + h3_3[0]
	
	
	# delta time plane
	h4_0 = np.histogram(d0["time1"] - d0["time0"], bins = bins_delta_plane, range = [-int(dp), int(dp)])
	h4_0_total = np.array(h_totals_stats.get("h4_0_total", [0]*bins_delta_plane))
	h4_0_total = h4_0_total + h4_0[0]
	
	h4_1 = np.histogram(d1["time1"] - d1["time0"], bins = bins_delta_plane, range = [-int(dp), int(dp)])
	h4_1_total = np.array(h_totals_stats.get("h4_1_total", [0]*bins_delta_plane))
	h4_1_total = h4_1_total + h4_1[0]
	
	h4_2 = np.histogram(d2["time1"] - d2["time0"], bins = bins_delta_plane, range = [-int(dp), int(dp)])
	h4_2_total = np.array(h_totals_stats.get("h4_2_total", [0]*bins_delta_plane))
	h4_2_total = h4_2_total + h4_2[0]
	
	h4_3 = np.histogram(d3["time1"] - d3["time0"], bins = bins_delta_plane, range = [-int(dp), int(dp)])
	h4_3_total = np.array(h_totals_stats.get("h4_3_total", [0]*bins_delta_plane))
	h4_3_total = h4_3_total + h4_3[0]
	
	
	# max missing strip 0
	h5_0 = np.histogram(d0["max_missing_strip0"], bins = bins_missing_strip, range = [0.0, bins_missing_strip])
	h5_0_total = np.array(h_totals_stats.get("h5_0_total", [0]*bins_missing_strip))
	h5_0_total = h5_0_total + h5_0[0]
	
	h5_1 = np.histogram(d1["max_missing_strip0"], bins = bins_missing_strip, range = [0.0, bins_missing_strip])
	h5_1_total = np.array(h_totals_stats.get("h5_1_total", [0]*bins_missing_strip))
	h5_1_total = h5_1_total + h5_1[0]
	
	h5_2 = np.histogram(d2["max_missing_strip0"], bins = bins_missing_strip, range = [0.0, bins_missing_strip])
	h5_2_total = np.array(h_totals_stats.get("h5_2_total", [0]*bins_missing_strip))
	h5_2_total = h5_2_total + h5_2[0]
	
	h5_3 = np.histogram(d3["max_missing_strip0"], bins = bins_missing_strip, range = [0.0, bins_missing_strip])
	h5_3_total = np.array(h_totals_stats.get("h5_3_total", [0]*bins_missing_strip))
	h5_3_total = h5_3_total + h5_3[0]			
	
	
	# max missing strio 1
	h6_0 = np.histogram(d0["max_missing_strip1"], bins = bins_missing_strip, range = [0.0, bins_missing_strip])
	h6_0_total = np.array(h_totals_stats.get("h6_0_total", [0]*bins_missing_strip))
	h6_0_total = h6_0_total + h6_0[0]
	
	h6_1 = np.histogram(d1["max_missing_strip1"], bins = bins_missing_strip, range = [0.0, bins_missing_strip])
	h6_1_total = np.array(h_totals_stats.get("h6_1_total", [0]*bins_missing_strip))
	h6_1_total = h6_1_total + h6_1[0]
	
	h6_2 = np.histogram(d2["max_missing_strip1"], bins = bins_missing_strip, range = [0.0, bins_missing_strip])
	h6_2_total = np.array(h_totals_stats.get("h6_2_total", [0]*bins_missing_strip))
	h6_2_total = h6_2_total + h6_2[0]
	
	h6_3 = np.histogram(d3["max_missing_strip1"], bins = bins_missing_strip, range = [0.0, bins_missing_strip])
	h6_3_total = np.array(h_totals_stats.get("h6_3_total", [0]*bins_missing_strip))
	h6_3_total = h6_3_total + h6_3[0]
	
	
	# max delta time 0
	h7_0 = np.histogram(d0["max_delta_time0"], bins = bins_max_delta_hits, range = [0.0, int(dt)])
	h7_0_total = np.array(h_totals_stats.get("h7_0_total", [0]*bins_max_delta_hits))
	h7_0_total = h7_0_total + h7_0[0]
	
	h7_1 = np.histogram(d1["max_delta_time0"], bins = bins_max_delta_hits, range = [0.0, int(dt)])
	h7_1_total = np.array(h_totals_stats.get("h7_1_total", [0]*bins_max_delta_hits))
	h7_1_total = h7_1_total + h7_1[0]	
	
	h7_2 = np.histogram(d2["max_delta_time0"], bins = bins_max_delta_hits, range = [0.0, int(dt)])
	h7_2_total = np.array(h_totals_stats.get("h7_2_total", [0]*bins_max_delta_hits))
	h7_2_total = h7_2_total + h7_2[0]
	
	h7_3 = np.histogram(d3["max_delta_time0"], bins = bins_max_delta_hits, range = [0.0, int(dt)])
	h7_3_total = np.array(h_totals_stats.get("h7_3_total", [0]*bins_max_delta_hits))
	h7_3_total = h7_3_total + h7_3[0]
	
	
	# max delta time 1	
	h8_0 = np.histogram(d0["max_delta_time1"], bins = bins_max_delta_hits, range = [0.0, int(dt)])
	h8_0_total = np.array(h_totals_stats.get("h8_0_total", [0]*bins_max_delta_hits))
	h8_0_total = h8_0_total + h8_0[0]
	
	h8_1 = np.histogram(d1["max_delta_time1"], bins = bins_max_delta_hits, range = [0.0, int(dt)])
	h8_1_total = np.array(h_totals_stats.get("h8_1_total", [0]*bins_max_delta_hits))
	h8_1_total = h8_1_total + h8_1[0]

	h8_2 = np.histogram(d2["max_delta_time1"], bins = bins_max_delta_hits, range = [0.0, int(dt)])
	h8_2_total = np.array(h_totals_stats.get("h8_2_total", [0]*bins_max_delta_hits))
	h8_2_total = h8_2_total + h8_2[0]
	
	h8_3 = np.histogram(d3["max_delta_time1"], bins = bins_max_delta_hits, range = [0.0, int(dt)])
	h8_3_total = np.array(h_totals_stats.get("h8_3_total", [0]*bins_max_delta_hits))
	h8_3_total = h8_3_total + h8_3[0]
	
	
	# span cluster 0
	h9_0 = np.histogram(d0["span_cluster0"], bins = bins_span, range = [0.0, int(spc)])
	h9_0_total = np.array(h_totals_stats.get("h9_0_total", [0]*bins_span))
	h9_0_total = h9_0_total + h9_0[0]

	h9_1 = np.histogram(d1["span_cluster0"], bins = bins_span, range = [0.0, int(spc)])
	h9_1_total = np.array(h_totals_stats.get("h9_1_total", [0]*bins_span))
	h9_1_total = h9_1_total + h9_1[0]
	
	h9_2 = np.histogram(d2["span_cluster0"], bins = bins_span, range = [0.0, int(spc)])
	h9_2_total = np.array(h_totals_stats.get("h9_2_total", [0]*bins_span))
	h9_2_total = h9_2_total + h9_2[0]
	
	h9_3 = np.histogram(d3["span_cluster0"], bins = bins_span, range = [0.0, int(spc)])
	h9_3_total = np.array(h_totals_stats.get("h9_3_total", [0]*bins_span))
	h9_3_total = h9_3_total + h9_3[0]


	# span cluster 1
	h10_0 = np.histogram(d0["span_cluster1"], bins = bins_span, range = [0.0, int(spc)])
	h10_0_total = np.array(h_totals_stats.get("h10_0_total", [0]*bins_span))
	h10_0_total = h10_0_total + h10_0[0]
	
	h10_1 = np.histogram(d1["span_cluster1"], bins = bins_span, range = [0.0, int(spc)])
	h10_1_total = np.array(h_totals_stats.get("h10_1_total", [0]*bins_span))
	h10_1_total = h10_1_total + h10_1[0]

	h10_2 = np.histogram(d2["span_cluster1"], bins = bins_span, range = [0.0, int(spc)])
	h10_2_total = np.array(h_totals_stats.get("h10_2_total", [0]*bins_span))
	h10_2_total = h10_2_total + h10_2[0]
	
	h10_3 = np.histogram(d3["span_cluster1"], bins = bins_span, range = [0.0, int(spc)])
	h10_3_total = np.array(h_totals_stats.get("h10_3_total", [0]*bins_span))
	h10_3_total = h10_3_total + h10_3[0]
	
	
	# cluster x/y percentage
	h11_00_total = np.array(h_totals_stats.get("h11_00_total", [0]*time_points))
	h11_00_total = np.roll(h11_00_total, -1)
	if p00["plane"].size > 0:
		h11_00_total[-1] = d0["size0"].size * 100 / p00["plane"].size
	else :
		h11_00_total[-1] = 0
	
	h11_01_total = np.array(h_totals_stats.get("h11_01_total", [0]*time_points))
	h11_01_total = np.roll(h11_01_total, -1)
	if p01["plane"].size > 0:
		h11_01_total[-1] = d1["size0"].size * 100 / p01["plane"].size
	else:
		h11_01_total[-1] = 0
	
	h11_02_total = np.array(h_totals_stats.get("h11_02_total", [0]*time_points))
	h11_02_total = np.roll(h11_02_total, -1)
	if p02["plane"].size > 0:
		h11_02_total[-1] = d2["size0"].size * 100 / p02["plane"].size
	else:
		h11_02_total[-1] = 0
	
	h11_03_total = np.array(h_totals_stats.get("h11_03_total", [0]*time_points))
	h11_03_total = np.roll(h11_03_total, -1)
	if p03["plane"].size > 0:
		h11_03_total[-1] = d3["size0"].size * 100 / p03["plane"].size
	else:
		h11_03_total[-1] = 0
				
	h11_10_total = np.array(h_totals_stats.get("h11_10_total", [0]*time_points))
	h11_10_total = np.roll(h11_10_total, -1)
	if p10["plane"].size > 0:
		h11_10_total[-1] = d0["size1"].size * 100 / p10["plane"].size
	else:
		h11_10_total[-1] = 0

	h11_11_total = np.array(h_totals_stats.get("h11_11_total", [0]*time_points))
	h11_11_total = np.roll(h11_11_total, -1)
	if p11["plane"].size > 0:
		h11_11_total[-1] = d1["size1"].size * 100 / p11["plane"].size
	else:
		h11_11_total[-1] = 0

	h11_12_total = np.array(h_totals_stats.get("h11_12_total", [0]*time_points))
	h11_12_total = np.roll(h11_12_total, -1)
	if p12["plane"].size > 0:
		h11_12_total[-1] = d2["size1"].size * 100 / p12["plane"].size
	else:
		h11_12_total[-1] = 0

	h11_13_total = np.array(h_totals_stats.get("h11_13_total", [0]*time_points))
	h11_13_total = np.roll(h11_13_total, -1)
	if p13["plane"].size > 0:
		h11_13_total[-1] = d3["size1"].size * 100 / p13["plane"].size
	else:
		h11_13_total[-1] = 0
	
	# cluster rate
	h12_0_total = np.array(h_totals_stats.get("h12_0_total", [0]*time_points))
	h12_0_total = np.roll(h12_0_total, -1)
	if d0["time0"].size > 0:
		if (max(d0["time0"]) -  min(d0["time0"])) > 0:
			h12_0_total[-1] = 1e+6*d0["time0"].size/(max(d0["time0"]) -  min(d0["time0"]))
		else:
			h12_0_total[-1] = 0
	else:
		h12_0_total[-1] = 0
	
	h12_1_total = np.array(h_totals_stats.get("h12_1_total", [0]*time_points))
	h12_1_total = np.roll(h12_1_total, -1)
	if d1["time0"].size > 0:
		if (max(d1["time0"]) -  min(d1["time0"])) > 0:
			h12_1_total[-1] = 1e+6*d1["time0"].size/(max(d1["time0"]) -  min(d1["time0"]))
		else:
			h12_1_total[-1] = 0
	else:
		h12_1_total[-1] = 0
	
	h12_2_total = np.array(h_totals_stats.get("h12_2_total", [0]*time_points))
	h12_2_total = np.roll(h12_2_total, -1)
	if d2["time0"].size > 0:
		if (max(d2["time0"]) -  min(d2["time0"])) > 0:
			h12_2_total[-1] = 1e+6*d2["time0"].size/(max(d2["time0"]) -  min(d2["time0"]))
		else:
			h12_2_total[-1] = 0
	else:
		h12_2_total[-1] = 0
	
	h12_3_total = np.array(h_totals_stats.get("h12_3_total", [0]*time_points))
	h12_3_total = np.roll(h12_3_total, -1)
	if d3["time0"].size > 1:
		if (max(d3["time0"]) -  min(d3["time0"])) > 0:
			h12_3_total[-1] = 1e+6*d3["time0"].size/(max(d3["time0"]) -  min(d3["time0"]))
		else:
			h12_3_total[-1] = 0
	else:
		h12_3_total[-1] = 0
	
	
	
	fig = make_subplots(rows = 3, cols = 4, horizontal_spacing=0.06, vertical_spacing=0.1, subplot_titles = ("clusters size0", "clusters size1","clusters size","delta time planes", "max missing strip0","max missing strip1", "max delta time0",  "max delta time1",  "span cluster0", "span cluster1", "common clusters x/y", "cluster rate"))
	fig.update_layout(
		plot_bgcolor='white',
		paper_bgcolor='white'
	)

	fig.add_trace(go.Bar(x = h1_0[1][:-1], y = h1_0_total, name='x q0', marker_color = colors[0]), row = 1, col = 1)
	fig.add_trace(go.Bar(x = h1_1[1][:-1], y = h1_1_total, name='x q1', marker_color = colors[1]), row = 1, col = 1)
	fig.add_trace(go.Bar(x = h1_2[1][:-1], y = h1_2_total, name='x q2', marker_color = colors[2]), row = 1, col = 1)
	fig.add_trace(go.Bar(x = h1_3[1][:-1], y = h1_3_total, name='x q3', marker_color = colors[3]), row = 1, col = 1)
	
	#fig.add_trace(go.Scatter(x=x_step_1,y=y_step_1_0,mode='lines',line=dict(shape='hv', color=color_x0, width=2)), row=1, col=1)
	#fig.add_trace(go.Scatter(x=x_step_1,y=y_step_1_1,mode='lines',line=dict(shape='hv', color=color_x1, width=2)), row=1, col=1)
	#fig.add_trace(go.Scatter(x=x_step_1,y=y_step_1_2,mode='lines',line=dict(shape='hv', color=color_x2, width=2)), row=1, col=1)
	#fig.add_trace(go.Scatter(x=x_step_1,y=y_step_1_3,mode='lines',line=dict(shape='hv', color=color_x3, width=2)), row=1, col=1)
	
	fig.add_trace(go.Bar(x = h2_0[1][:-1], y = h2_0_total, name='y q0', marker_color = colors[4]), row = 1, col = 2)
	fig.add_trace(go.Bar(x = h2_1[1][:-1], y = h2_1_total, name='y q1', marker_color = colors[5]), row = 1, col = 2)
	fig.add_trace(go.Bar(x = h2_2[1][:-1], y = h2_2_total, name='y q2', marker_color = colors[6]), row = 1, col = 2)
	fig.add_trace(go.Bar(x = h2_3[1][:-1], y = h2_3_total, name='y q3', marker_color = colors[7%len(colors)]), row = 1, col = 2)
	
	fig.add_trace(go.Bar(x = h3_0[1][:-1], y = h3_0_total, name='q0', marker_color = colors[8%len(colors)]), row = 1, col = 3)
	fig.add_trace(go.Bar(x = h3_1[1][:-1], y = h3_1_total, name='q1', marker_color = colors[9%len(colors)]), row = 1, col = 3)
	fig.add_trace(go.Bar(x = h3_2[1][:-1], y = h3_2_total, name='q2', marker_color = colors[10%len(colors)]), row = 1, col = 3)
	fig.add_trace(go.Bar(x = h3_3[1][:-1], y = h3_3_total, name='q3', marker_color = colors[11%len(colors)]), row = 1, col = 3)
	
	fig.add_trace(go.Bar(x = h4_0[1][:-1], y = h4_0_total, name='q0', marker_color = colors[8%len(colors)]), row = 1, col = 4)
	fig.add_trace(go.Bar(x = h4_1[1][:-1], y = h4_1_total, name='q1',marker_color = colors[9%len(colors)]), row = 1, col = 4)
	fig.add_trace(go.Bar(x = h4_2[1][:-1], y = h4_2_total, name='q2',marker_color = colors[10%len(colors)]), row = 1, col = 4)
	fig.add_trace(go.Bar(x = h4_3[1][:-1], y = h4_3_total, name='q3',marker_color = colors[11%len(colors)]), row = 1, col = 4)
	
	fig.add_trace(go.Bar(x = h5_0[1][:-1], y = h5_0_total,  name='x q0',marker_color = colors[0]), row = 2, col = 1)
	fig.add_trace(go.Bar(x = h5_1[1][:-1], y = h5_1_total,  name='x q1',marker_color = colors[1]), row = 2, col = 1)
	fig.add_trace(go.Bar(x = h2_2[1][:-1], y = h5_2_total,  name='x q2',marker_color = colors[2]), row = 2, col = 1)
	fig.add_trace(go.Bar(x = h2_3[1][:-1], y = h5_3_total,  name='x q3',marker_color = colors[3]), row = 2, col = 1)
	
	fig.add_trace(go.Bar(x = h6_0[1][:-1], y = h6_0_total,  name='y q0',marker_color = colors[4]), row = 2, col = 2)
	fig.add_trace(go.Bar(x = h6_1[1][:-1], y = h6_1_total,  name='y q1',marker_color = colors[5]), row = 2, col = 2)
	fig.add_trace(go.Bar(x = h6_2[1][:-1], y = h6_2_total,  name='y q2',marker_color = colors[6]), row = 2, col = 2)
	fig.add_trace(go.Bar(x = h6_3[1][:-1], y = h6_3_total,  name='y q3',marker_color = colors[7%len(colors)]), row = 2, col = 2)
	
	fig.add_trace(go.Bar(x = h7_0[1][:-1], y = h7_0_total,  name='x q0',marker_color = colors[0]), row = 2, col = 3)
	fig.add_trace(go.Bar(x = h7_1[1][:-1], y = h7_1_total,  name='x q1',marker_color = colors[1]), row = 2, col = 3)
	fig.add_trace(go.Bar(x = h7_2[1][:-1], y = h7_2_total,  name='x q2',marker_color = colors[2]), row = 2, col = 3)
	fig.add_trace(go.Bar(x = h7_3[1][:-1], y = h7_3_total,  name='x q3',marker_color = colors[3]), row = 2, col = 3)
	
	fig.add_trace(go.Bar(x = h8_0[1][:-1], y = h8_0_total,  name='y q3',marker_color = colors[4]), row = 2, col = 4)
	fig.add_trace(go.Bar(x = h8_1[1][:-1], y = h8_1_total,  name='y q3',marker_color = colors[5]), row = 2, col = 4)
	fig.add_trace(go.Bar(x = h8_2[1][:-1], y = h8_2_total,  name='y q3',marker_color = colors[6]), row = 2, col = 4)
	fig.add_trace(go.Bar(x = h8_3[1][:-1], y = h8_3_total,  name='y q3',marker_color = colors[7%len(colors)]), row = 2, col = 4)
	
	fig.add_trace(go.Bar(x = h9_0[1][:-1], y = h9_0_total,  name='x q0',marker_color = colors[0]), row = 3, col = 1)
	fig.add_trace(go.Bar(x = h9_1[1][:-1], y = h9_1_total,  name='x q1',marker_color = colors[1]), row = 3, col = 1)
	fig.add_trace(go.Bar(x = h9_2[1][:-1], y = h9_2_total,  name='x q2',marker_color = colors[2]), row = 3, col = 1)
	fig.add_trace(go.Bar(x = h9_3[1][:-1], y = h9_3_total,  name='x q3',marker_color = colors[3]), row = 3, col = 1)
	
	fig.add_trace(go.Bar(x = h10_0[1][:-1], y = h10_0_total,  name='y q0',marker_color = colors[4]), row = 3, col = 2)
	fig.add_trace(go.Bar(x = h10_1[1][:-1], y = h10_1_total,  name='y q1',marker_color = colors[5]), row = 3, col = 2)
	fig.add_trace(go.Bar(x = h10_2[1][:-1], y = h10_2_total,  name='y q2',marker_color = colors[6]), row = 3, col = 2)
	fig.add_trace(go.Bar(x = h10_3[1][:-1], y = h10_3_total,  name='y q3',marker_color = colors[7%len(colors)]), row = 3, col = 2)
	
	if d0["time0"].size > 0:					
		fig.add_trace(go.Scatter( x = h_time_total, y = h11_00_total, name='x q0', mode='markers', marker=dict(symbol='circle',size=16,color=colors[0])), row = 3, col = 3)
	if d1["time0"].size > 0:
		fig.add_trace(go.Scatter( x = h_time_total, y = h11_01_total, name='x q1', mode='markers', marker=dict(symbol='circle',size=16,color=colors[1])), row = 3, col = 3)
	if d2["time0"].size > 0:
		fig.add_trace(go.Scatter( x = h_time_total, y = h11_02_total, name='x q2', mode='markers', marker=dict(symbol='circle',size=16,color=colors[2])), row = 3, col = 3)
	if d3["time0"].size > 0:
		fig.add_trace(go.Scatter( x = h_time_total, y = h11_03_total, name='x q3', mode='markers', marker=dict(symbol='circle',size=16,color=colors[3])), row = 3, col = 3)
	
	if d0["time0"].size > 0:
		fig.add_trace(go.Scatter( x = h_time_total, y = h11_10_total, name='y q0', mode='markers', marker=dict(symbol='circle',size=16,color=colors[4])), row = 3, col = 3)
	if d1["time0"].size > 0:
		fig.add_trace(go.Scatter( x = h_time_total, y = h11_11_total, name='y q1', mode='markers', marker=dict(symbol='circle',size=16,color=colors[5])), row = 3, col = 3)
	if d2["time0"].size > 0:
		fig.add_trace(go.Scatter( x = h_time_total, y = h11_12_total, name='y q2', mode='markers', marker=dict(symbol='circle',size=16,color=colors[6])), row = 3, col = 3)
	if d3["time0"].size > 0:
		fig.add_trace(go.Scatter( x = h_time_total, y = h11_13_total, name='y q3', mode='markers', marker=dict(symbol='circle',size=16,color=colors[7%len(colors)])), row = 3, col = 3)
	
	if d0["time0"].size > 0:
		fig.add_trace(go.Scatter( x = h_time_total, y = h12_0_total, name='q0', mode='markers', marker=dict(symbol='circle',size=16,color=colors[8%len(colors)])), row = 3, col = 4)
	if d1["time0"].size > 0:
		fig.add_trace(go.Scatter( x = h_time_total, y = h12_1_total, name='q1', mode='markers', marker=dict(symbol='circle',size=16,color=colors[9%len(colors)])), row = 3, col = 4)
	if d2["time0"].size > 0:
		fig.add_trace(go.Scatter( x = h_time_total, y = h12_2_total, name='q2', mode='markers', marker=dict(symbol='circle',size=16,color=colors[10%len(colors)])), row = 3, col = 4)
	if d3["time0"].size > 0:
		fig.add_trace(go.Scatter( x = h_time_total, y = h12_3_total, name='q3', mode='markers', marker=dict(symbol='circle',size=16,color=colors[11%len(colors)])), row = 3, col = 4)


	fig.update_xaxes(title_text="size [strips]", row = 1, col = 1)
	fig.update_yaxes(type=yaxis_type,title_text="counts", row = 1, col = 1)
	fig.update_xaxes(title_text="size [strips]", row = 1, col = 2)
	fig.update_yaxes(type=yaxis_type,title_text="counts", row = 1, col = 2)
	fig.update_xaxes(title_text="size [strips]", row = 1, col = 3)
	fig.update_yaxes(type=yaxis_type,title_text="counts", row = 1, col = 3)
	

	fig.update_xaxes(title_text="time difference [ns]", row = 1, col = 4)
	fig.update_yaxes(type=yaxis_type,title_text="counts", row = 1, col = 4)
	fig.update_xaxes(title_text="missing [strips]", row = 2, col = 1)
	fig.update_yaxes(type=yaxis_type,title_text="counts", row = 2, col = 1)
	fig.update_xaxes(title_text="missing [strips]", row = 2, col = 2)
	fig.update_yaxes(type=yaxis_type,title_text="counts", row = 2, col = 2)

	fig.update_xaxes(title_text="time difference [ns]", row = 2, col = 3)
	fig.update_yaxes(type=yaxis_type,title_text="counts", row = 2, col = 3)
	fig.update_xaxes(title_text="time difference [ns]", row = 2, col = 4)
	fig.update_yaxes(type=yaxis_type,title_text="counts", row = 2, col = 4)
	
	fig.update_xaxes(title_text="time span [ns]", row = 3, col = 1)
	fig.update_yaxes(type=yaxis_type,title_text="counts", row = 3, col = 1)
	fig.update_xaxes(title_text="time span [ns]", row = 3, col = 2)
	fig.update_yaxes(type=yaxis_type,title_text="counts", row = 3, col = 2)
	fig.update_xaxes(title_text="time since start of run [s]", row = 3, col = 3)
	fig.update_yaxes(type=yaxis_type,title_text="percentage [%]", row = 3, col = 3)

	fig.update_xaxes(title_text="time since start of run [s]", row = 3, col = 4)
	fig.update_yaxes(type=yaxis_type,title_text="rate [kHz]", row = 3, col = 4)	
	
	fig.update_layout(yaxis=dict(type=yaxis_type), barmode='stack',uirevision='constant', height = the_height, width = the_width, showlegend = False)
	return fig, {
		"h1_0_total": h1_0_total.tolist(),
		"h1_1_total": h1_1_total.tolist(),
		"h1_2_total": h1_2_total.tolist(),
		"h1_3_total": h1_3_total.tolist(),		
		"h2_0_total": h2_0_total.tolist(),
		"h2_1_total": h2_1_total.tolist(),
		"h2_2_total": h2_2_total.tolist(),
		"h2_3_total": h2_3_total.tolist(),		
		"h3_0_total": h3_0_total.tolist(),
		"h3_1_total": h3_1_total.tolist(),
		"h3_2_total": h3_2_total.tolist(),
		"h3_3_total": h3_3_total.tolist(),		
		"h4_0_total": h4_0_total.tolist(),
		"h4_1_total": h4_1_total.tolist(),
		"h4_2_total": h4_2_total.tolist(),
		"h4_3_total": h4_3_total.tolist(),		
		"h5_0_total": h5_0_total.tolist(),
		"h5_1_total": h5_1_total.tolist(),
		"h5_2_total": h5_2_total.tolist(),
		"h5_3_total": h5_3_total.tolist(),		
		"h6_0_total": h6_0_total.tolist(),
		"h6_1_total": h6_1_total.tolist(),
		"h6_2_total": h6_2_total.tolist(),
		"h6_3_total": h6_3_total.tolist(),		
		"h7_0_total": h7_0_total.tolist(),
		"h7_1_total": h7_1_total.tolist(),
		"h7_2_total": h7_2_total.tolist(),
		"h7_3_total": h7_3_total.tolist(),
		"h8_0_total": h8_0_total.tolist(),
		"h8_1_total": h8_1_total.tolist(),
		"h8_2_total": h8_2_total.tolist(),
		"h8_3_total": h8_3_total.tolist(),
		"h9_0_total": h9_0_total.tolist(),
		"h9_1_total": h9_1_total.tolist(),
		"h9_2_total": h9_2_total.tolist(),
		"h9_3_total": h9_3_total.tolist(),
		"h10_0_total": h10_0_total.tolist(),
		"h10_1_total": h10_1_total.tolist(),
		"h10_2_total": h10_2_total.tolist(),
		"h10_3_total": h10_3_total.tolist(),		
		"h11_00_total": h11_00_total.tolist(),
		"h11_01_total": h11_01_total.tolist(),
		"h11_02_total": h11_02_total.tolist(),
		"h11_03_total": h11_03_total.tolist(),
		"h11_10_total": h11_10_total.tolist(),
		"h11_11_total": h11_11_total.tolist(),
		"h11_12_total": h11_12_total.tolist(),
		"h11_13_total": h11_13_total.tolist(),		
		"h12_0_total": h12_0_total.tolist(),
		"h12_1_total": h12_1_total.tolist(),
		"h12_2_total": h12_2_total.tolist(),
		"h12_3_total": h12_3_total.tolist(),
		"h_start_time": h_start_time.tolist(),
		"h_time_total": h_time_total.tolist()
	}
