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

last_seen_cluster_timestamp = 0.0

dash.register_page(__name__, path='/clusters', name="Clusters")

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
			
			# Log z-axis
			dcc.Checklist(
				id='logz_toggle',
				options=[{'label': '2D plots: z-axis log', 'value': 'logz'}],
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
				style={"width": "120px", "margin-left": "5px", "margin-right": "20px", "verticalAlign": "middle"}
			),
			html.Div(id="color_preview", style={"display": "inline-block", "margin-left": "5px", "margin-right": "20px"}),
			html.Label("Color map 2D:"),
			dcc.Dropdown(
				id="heatmap_color",
				options=get_colorscale_options(),
				value="Hot",
				style={"width": "200px", "margin-left": "5px", "verticalAlign": "middle"}
			)
			], style={
			"display": "flex",
			"alignItems": "center",  # vertically align
			"justifyContent": "center",
			"margin-bottom": "20px"
		})
	]),
	dcc.Graph(id='cluster_graph'),
	dcc.Store(id='h_totals_clusters')
])



@dash.callback(
	Output('cluster_graph', 'figure'),
	Output('h_totals_clusters', 'data'),
	Input('cluster_data', 'data'),
	Input('logy_toggle', 'value'),
	Input('logz_toggle', 'value'),
	Input("color_palette", "value"),
	Input("heatmap_color", "value"),
	State('h_totals_clusters', 'data')
)
def plot_data(cluster_data,logy_toggle,logz_toggle, color_palette, heatmap_color, h_totals_clusters):
	global last_seen_cluster_timestamp
	h_totals_clusters = h_totals_clusters or {}
	updated_at = shared_data.get("cluster_updated_at", 0.0)
	if not cluster_data or updated_at == last_seen_cluster_timestamp:
		return dash.no_update, h_totals_clusters
	last_seen_cluster_timestamp = updated_at
	df_clusters = pd.DataFrame(cluster_data) 
	d0 = df_clusters.query("det == 0")
	d1 = df_clusters.query("det == 1")
	d2 = df_clusters.query("det == 2")
	d3 = df_clusters.query("det == 3")

	if color_palette == "Config":
		colors = color_config
	else:
		colors = palette_map.get(color_palette, px.colors.qualitative.Plotly)
	yaxis_type = 'log' if 'logy' in logy_toggle else 'linear'
	zaxis_type = 'log' if 'logz' in logz_toggle else 'linear'

	if cluster_algorithm == "utpc":
		pos0x = d0["pos0_utpc"]
		pos1x = d1["pos0_utpc"]
		pos2x = 1280+d2["pos0_utpc"]
		pos3x = 1280+d3["pos0_utpc"]
		pos0y = d0["pos1_utpc"]
		pos1y = 1280+d1["pos1_utpc"]
		pos2y = d2["pos1_utpc"]
		pos3y = 1280+d3["pos1_utpc"]
		h4 = np.histogram2d(df_clusters["pos0_utpc"], df_clusters["pos1_utpc"], bins = [image_x, image_y], range = np.array([(0, image_x), (0, image_y)]))
		fig = make_subplots(rows = 2, cols = 4, horizontal_spacing=0.06, vertical_spacing=0.1, subplot_titles = ("clusters pos0_utpc", "clusters pos1_utpc", "clusters adc0+adc1", "clusters pos1_utpc:pos0_utpc","total clusters pos0_utpc", "total clusters pos1_utpc", "total clusters adc0+adc1", "total clusters pos1_utpc:pos0_utpc"))
		
	elif cluster_algorithm == "charge2":
		pos0x = d0["pos0_charge2"]
		pos1x = d1["pos0_charge2"]
		pos2x = 1280+d2["pos0_charge2"]
		pos3x = 1280+d3["pos0_charge2"]
		pos0y = d0["pos1_charge2"]
		pos1y = 1280+d1["pos1_charge2"]
		pos2y = d2["pos1_charge2"]
		pos3y = 1280+d3["pos1_charge2"]	
		h4 = np.histogram2d(df_clusters["pos0_charge2"], df_clusters["pos1_charge2"], bins = [image_x, image_y], range = np.array([(0, image_x), (0, image_y)]))
		fig = make_subplots(rows = 2, cols = 4, horizontal_spacing=0.06, vertical_spacing=0.1, subplot_titles = ("clusters pos0_charge2", "clusters pos1_charge2", "clusters adc0+adc1", "clusters pos1_charge2:pos0_charge2","total clusters pos0_charge2", "total clusters pos1_charge2", "total clusters adc0+adc1", "total clusters pos1_charge2:pos0_charge2"))
	
	else:
		pos0x = d0["pos0"]
		pos1x = d1["pos0"]
		pos2x = 1280+d2["pos0"]
		pos3x = 1280+d3["pos0"]
		pos0y = d0["pos1"]
		pos1y = 1280+d1["pos1"]
		pos2y = d2["pos1"]
		pos3y = 1280+d3["pos1"]	
		h4 = np.histogram2d(df_clusters["pos0"], df_clusters["pos1"], bins = [image_x, image_y], range = np.array([(0, image_x), (0, image_y)]))
		fig = make_subplots(rows = 2, cols = 4, horizontal_spacing=0.06, vertical_spacing=0.1, subplot_titles = ("clusters pos0", "clusters pos1", "clusters adc0+adc1", "clusters pos1:pos0","total clusters pos0", "total clusters pos1", "total clusters adc0+adc1", "total clusters pos1:pos0"))
	
	
	posx = np.concatenate([pos0x, pos1x, pos2x, pos3x])
	posy = np.concatenate([pos0y, pos1y, pos2y, pos3y])
	
	h1 = np.histogram(posx, bins = channels_x, range = [0.0, channels_x])
	h1_total = np.array(h_totals_clusters.get("h1", [0]*channels_x))
	h1_total += h1[0]
	
	h2 = np.histogram(posy, bins = channels_y, range = [0.0, channels_y])
	h2_total = np.array(h_totals_clusters.get("h2", [0]*channels_y))
	h2_total += h2[0]
	
	h3_0 = np.histogram(d0["adc0"]+d0["adc1"], bins = int(max_charge*charge_scale), range = [0.0, max_charge])
	h3_0_total = np.array(h_totals_clusters.get("h3_0", [0]*int(max_charge*charge_scale)))
	h3_0_total += h3_0[0]
	h3_1 = np.histogram(d1["adc0"]+d1["adc1"], bins = int(max_charge*charge_scale), range = [0.0, max_charge])
	h3_1_total = np.array(h_totals_clusters.get("h3_1", [0]*int(max_charge*charge_scale)))
	h3_1_total += h3_1[0]
	h3_2 = np.histogram(d2["adc0"]+d2["adc1"], bins = int(max_charge*charge_scale), range = [0.0, max_charge])
	h3_2_total = np.array(h_totals_clusters.get("h3_2", [0]*int(max_charge*charge_scale)))
	h3_2_total += h3_2[0]
	h3_3 = np.histogram(d3["adc0"]+d3["adc1"], bins = int(max_charge*charge_scale), range = [0.0, max_charge])
	h3_3_total = np.array(h_totals_clusters.get("h3_3", [0]*int(max_charge*charge_scale)))
	h3_3_total += h3_3[0]
	
	h4_total = np.array(
	h_totals_clusters.get("h4", np.zeros((image_x, image_y), dtype=np.float64).tolist()),
	dtype=np.float64
	)
	h4_total += h4[0]
	
		
	h1_patched = np.where(h1[0] == 0, 0.1, h1[0])
	h2_patched = np.where(h2[0] == 0, 0.1, h2[0])
	h4_patched = np.where(h4[0] == 0, 0.1, h4[0])
	h1_total_patched = np.where(h1_total == 0, 0.1, h1_total)
	h2_total_patched = np.where(h2_total == 0, 0.1, h2_total)
	h4_total_patched = np.where(h4_total == 0, 0.1, h4_total)
	
	if zaxis_type == 'log':
		h4_patched = np.log10(h4_patched)
		h4_total_patched = np.log10(h4_total_patched)	
	
		min_val_41 = np.floor(h4_patched.min())
		max_val_41 = np.ceil(h4_patched.max())
		min_val_42 = np.floor(h4_total_patched.min())
		max_val_42 = np.ceil(h4_total_patched.max())

		tickvals_41 = np.arange(min_val_41, max_val_41 + 1)
		ticktext_41 = [f"10^{int(v)}" for v in tickvals_41]
		tickvals_42 = np.arange(min_val_42, max_val_42 + 1)
		ticktext_42 = [f"10^{int(v)}" for v in tickvals_42]		
		colorbar_41 = dict(title='counts',tickvals=tickvals_41,ticktext=ticktext_41,x=1.0,y=0.79,len=0.49,thickness=20)
		colorbar_42 = dict(title='counts',tickvals=tickvals_42,ticktext=ticktext_42,x=1.0,y=0.24,len=0.49,thickness=20)

		fig.add_trace(go.Heatmap(z = np.transpose(h4_patched), colorbar=colorbar_41,colorscale=heatmap_color), row = 1, col = 4)
		fig.add_trace(go.Heatmap(z = np.transpose(h4_total_patched), colorbar=colorbar_42,colorscale=heatmap_color), row = 2, col = 4)	
	else:
		fig.add_trace(go.Heatmap(z = np.transpose(h4_patched), colorbar=dict(title="Counts",x=1.0,y=0.79,len=0.49,thickness=20),colorscale=heatmap_color), row = 1, col = 4)
		fig.add_trace(go.Heatmap(z = np.transpose(h4_total_patched), colorbar=dict(title="Counts",x=1.0,y=0.24,len=0.49,thickness=20),colorscale=heatmap_color), row = 2, col = 4)	


		
	fig.add_trace(go.Bar(x = h1[1][:-1], y = h1_patched,  name="x",marker_color = colors[0]), row = 1, col = 1)
	fig.add_trace(go.Bar(x = h2[1][:-1], y = h2_patched,  name="y",marker_color = colors[4]), row = 1, col = 2)
	fig.add_trace(go.Bar(x = h3_0[1][:-1], y = h3_0[0], name="q0", marker_color = colors[8%len(colors)]), row = 1, col = 3)
	fig.add_trace(go.Bar(x = h3_1[1][:-1], y = h3_1[0], name="q1",marker_color = colors[9%len(colors)]), row = 1, col = 3)
	fig.add_trace(go.Bar(x = h3_2[1][:-1], y = h3_2[0], name="q2",marker_color = colors[10%len(colors)]), row = 1, col = 3)
	fig.add_trace(go.Bar(x = h3_3[1][:-1], y = h3_3[0], name="q3",marker_color = colors[11%len(colors)]), row = 1, col = 3)

	fig.add_trace(go.Bar(x = h1[1][:-1], y = h1_total_patched,  name="x",marker_color = colors[0]), row = 2, col = 1)
	fig.add_trace(go.Bar(x = h2[1][:-1], y = h2_total_patched,  name="y",marker_color = colors[4]), row = 2, col = 2)
	fig.add_trace(go.Bar(x = h3_0[1][:-1], y = h3_0_total, name="q0", marker_color = colors[8%len(colors)]), row = 2, col = 3)
	fig.add_trace(go.Bar(x = h3_1[1][:-1], y = h3_1_total, name="q1",marker_color = colors[9%len(colors)]), row = 2, col = 3)
	fig.add_trace(go.Bar(x = h3_2[1][:-1], y = h3_2_total, name="q2",marker_color = colors[10%len(colors)]), row = 2, col = 3)
	fig.add_trace(go.Bar(x = h3_3[1][:-1], y = h3_3_total, name="q3",marker_color = colors[11%len(colors)]), row = 2, col = 3)
	
	fig.update_xaxes(title_text=("x [pitch 0.4 mm]"), row = 1, col = 1)
	fig.update_yaxes(type=yaxis_type,title_text=("counts"), row = 1, col = 1)
	
	fig.update_xaxes(title_text=("y [pitch 0.4 mm]"), row = 1, col = 2)
	fig.update_yaxes(type=yaxis_type,title_text=("counts"), row = 1, col = 2)
	
	fig.update_xaxes(title_text=("charge"), row = 1, col = 3)
	fig.update_yaxes(type=yaxis_type,title_text=("counts"), row = 1, col = 3)
	
	fig.update_xaxes(title_text=("x [pitch 0.4 mm]"), row = 1, col = 4)
	fig.update_yaxes(title_text=("y [pitch 0.4 mm]"), row = 1, col = 4)
	
	fig.update_xaxes(title_text=("x [pitch 0.4 mm]"), row = 2, col = 1)
	fig.update_yaxes(type=yaxis_type,title_text=("counts"), row = 2, col = 1)
	
	fig.update_xaxes(title_text=("y [pitch 0.4 mm]"), row = 2, col = 2)
	fig.update_yaxes(type=yaxis_type,title_text=("counts"), row = 2, col = 2)

	fig.update_xaxes(title_text=("charge"), row = 2, col = 3)
	fig.update_yaxes(type=yaxis_type,title_text=("counts"), row = 2, col = 3)
	
	fig.update_xaxes(title_text=("x [pitch 0.4 mm]"), row = 2, col = 4)
	fig.update_yaxes(title_text=("y [pitch 0.4 mm]"), row = 2, col = 4)
	
	fig.update_layout(height = the_height, width = the_width, showlegend = False, barmode='overlay',uirevision='constant',bargap=0,plot_bgcolor='white',paper_bgcolor='white')

	return fig, {
		"h1": h1_total.tolist(),
		"h2": h2_total.tolist(),
		"h3_0": h3_0_total.tolist(),
		"h3_1": h3_1_total.tolist(),
		"h3_2": h3_2_total.tolist(),
		"h3_3": h3_3_total.tolist(),
		"h4": h4_total.tolist()
	}
