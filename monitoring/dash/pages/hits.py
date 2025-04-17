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

dash.register_page(__name__, path='/hits', name="Hits")


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
				style={"width": "120px", "margin-left": "5px",  "margin-right": "5px", "verticalAlign": "middle"}
			),
			html.Div(id="color_preview", style={"display": "inline-block", "margin-left": "5px", "margin-right": "20px"}),
			html.Label("Color map 2D:"),
			dcc.Dropdown(
				id="heatmap_color",
				options=get_colorscale_options(),
				value="Jet",
				style={"width": "200px", "margin-left": "5px", "verticalAlign": "middle"}
			)
		], style={
			"display": "flex",
			"alignItems": "center",  # vertically align
			"justifyContent": "center",
			"margin-bottom": "20px"
		}),
		
	]),
	dcc.Graph(id='hit_graph'),
	dcc.Store(id='h_totals_hits')
])

@dash.callback(
	Output("color_preview", "children"),
	Input("color_palette", "value")
)
def update_preview(color_palette):
	if color_palette == "Config":
		colors = color_config
	else:
		colors = getattr(pc.qualitative, color_palette)
	return html.Div(
		[html.Div(style={
			"width": "15px", "height": "15px", "backgroundColor": c,
			"display": "inline-block", "marginRight": "3px", "border": "1px solid #ccc"
		}) for c in colors]
	)

@dash.callback(
	Output('hit_graph', 'figure'),
	Output('h_totals_hits', 'data'),
	Input('hit_data', 'data'),
	Input('logy_toggle', 'value'),
	Input('logz_toggle', 'value'),
	Input("color_palette", "value"),
	Input("heatmap_color", "value"),
	State('h_totals_hits', 'data'),
)
def plot_data(hit_data,logy_toggle, logz_toggle, color_palette, heatmap_color, h_totals_hits):
	h_totals_hits = h_totals_hits or {}
	if not hit_data:
		return dash.no_update, h_totals_hits
	df_hits = pd.DataFrame(hit_data)
	hits0 = df_hits.query("plane == 0 and (det >= 0 and det <=3)")
	hits1 = df_hits.query("plane == 1 and (det >= 0 and det <=3)")
	fig = make_subplots(rows = 2, cols = 4, horizontal_spacing=0.06, vertical_spacing=0.1, subplot_titles = ("hits pos0", "hits pos1", "hits adc0", "hits adc1", "total hits pos0", "total hits pos1", "total hits adc0", "total hits adc1"))

	if color_palette == "Config":
		colors = color_config
	else:
		colors = palette_map.get(color_palette, px.colors.qualitative.Plotly)
	yaxis_type = 'log' if 'logy' in logy_toggle else 'linear'
	zaxis_type = 'log' if 'logz' in logz_toggle else 'linear'

	
	h1 = np.histogram(hits0['pos'], bins = channels_x, range = [0.0, channels_x])
	h1_total = np.array(h_totals_hits.get("h1", [0]*channels_x))
	h1_total += h1[0]
	
	h2 = np.histogram(hits1['pos'], bins = channels_y, range = [0.0, channels_y])
	h2_total = np.array(h_totals_hits.get("h2", [0]*channels_y))
	h2_total += h2[0]

	h3 = np.histogram2d(hits0['pos'], hits0['adc'], bins = [channels_x, 128], range = np.array([(0, channels_x), (0,1024)]))
	h3_total = np.array(
	h_totals_hits.get("h3", np.zeros((channels_x, 128), dtype=np.float64).tolist()),
	dtype=np.float64
	)
	h3_total += h3[0]

	h4 = np.histogram2d(hits1['pos'], hits1['adc'], bins = [channels_y, 128], range = np.array([(0, channels_y), (0,1024)]))
	h4_total = np.array(
	h_totals_hits.get("h4", np.zeros((channels_y, 128), dtype=np.float64).tolist()),
	dtype=np.float64
	)
	h4_total += h4[0]
		
	h1_patched = np.where(h1[0] == 0, 0.1, h1[0])
	h2_patched = np.where(h2[0] == 0, 0.1, h2[0])
	h1_total_patched = np.where(h1_total == 0, 0.1, h1_total)
	h2_total_patched = np.where(h2_total == 0, 0.1, h2_total)
	h3_patched = np.where(h3[0] == 0, 0.1, h3[0])
	h4_patched = np.where(h4[0] == 0, 0.1, h4[0])
	h3_total_patched = np.where(h3_total == 0, 0.1, h3_total)
	h4_total_patched = np.where(h4_total == 0, 0.1, h4_total)
	if zaxis_type == 'log':
		h3_patched = np.log10(h3_patched)
		h4_patched = np.log10(h4_patched)
		h3_total_patched = np.log10(h3_total_patched)
		h4_total_patched = np.log10(h4_total_patched)	
	
		min_val_31 = np.floor(h3_patched.min())
		max_val_31 = np.ceil(h3_patched.max())
		min_val_32 = np.floor(h3_total_patched.min())
		max_val_32 = np.ceil(h3_total_patched.max())
		min_val_41 = np.floor(h4_patched.min())
		max_val_41 = np.ceil(h4_patched.max())
		min_val_42 = np.floor(h4_total_patched.min())
		max_val_42 = np.ceil(h4_total_patched.max())
			
		tickvals_31 = np.arange(min_val_31, max_val_31 + 1)
		ticktext_31 = [f"10^{int(v)}" for v in tickvals_31]
		tickvals_32 = np.arange(min_val_32, max_val_32 + 1)
		ticktext_32 = [f"10^{int(v)}" for v in tickvals_32]
		tickvals_41 = np.arange(min_val_41, max_val_41 + 1)
		ticktext_41 = [f"10^{int(v)}" for v in tickvals_41]
		tickvals_42 = np.arange(min_val_42, max_val_42 + 1)
		ticktext_42 = [f"10^{int(v)}" for v in tickvals_42]		
		colorbar_31 = dict(title='counts',tickvals=tickvals_31,ticktext=ticktext_31,x=0.736,y=0.79,len=0.49,thickness=20)
		colorbar_41 = dict(title='counts',tickvals=tickvals_41,ticktext=ticktext_41,x=1.0,y=0.79,len=0.49,thickness=20)
		colorbar_32 = dict(title='counts',tickvals=tickvals_32,ticktext=ticktext_32,x=0.736,y=0.24,len=0.49,thickness=20)
		colorbar_42 = dict(title='counts',tickvals=tickvals_42,ticktext=ticktext_42,x=1.0,y=0.24,len=0.49,thickness=20)

		fig.add_trace(go.Heatmap(z = np.transpose(h3_patched), y=h3[2][:-1],colorbar=colorbar_31,colorscale=heatmap_color), row = 1, col = 3)
		fig.add_trace(go.Heatmap(z = np.transpose(h4_patched), y=h4[2][:-1],colorbar=colorbar_41,colorscale=heatmap_color), row = 1, col = 4)
		fig.add_trace(go.Heatmap(z = np.transpose(h3_total_patched),y=h3[2][:-1],colorbar=colorbar_32,colorscale=heatmap_color), row = 2, col = 3)
		fig.add_trace(go.Heatmap(z = np.transpose(h4_total_patched),y=h4[2][:-1], colorbar=colorbar_42,colorscale=heatmap_color), row = 2, col = 4)
	else:
		fig.add_trace(go.Heatmap(z = np.transpose(h3[0]), y=h3[2][:-1],colorbar=dict(title="Counts",x=0.736,y=0.79,len=0.49,thickness=20),colorscale=heatmap_color), row = 1, col = 3)
		fig.add_trace(go.Heatmap(z = np.transpose(h4[0]),y=h4[2][:-1], colorbar=dict(title="Counts",x=1.0,y=0.79,len=0.49,thickness=20),colorscale=heatmap_color), row = 1, col = 4)
		fig.add_trace(go.Heatmap(z = np.transpose(h3_total),y=h3[2][:-1],colorbar=dict(title="Counts",x=0.736,y=0.24,len=0.49,thickness=20),colorscale=heatmap_color), row = 2, col = 3)
		fig.add_trace(go.Heatmap(z = np.transpose(h4_total), y=h4[2][:-1],colorbar=dict(title="Counts",x=1.0,y=0.24,len=0.49,thickness=20),colorscale=heatmap_color), row = 2, col = 4)

	fig.add_trace(go.Bar(x = h1[1][:-1], y = h1_patched, marker_color = colors[0]), row = 1, col = 1)
	fig.add_trace(go.Bar(x = h2[1][:-1], y = h2_patched, marker_color = colors[4]), row = 1, col = 2)

	fig.add_trace(go.Bar(x = h1[1][:-1], y = h1_total_patched, marker_color = colors[0]), row = 2, col = 1)
	fig.add_trace(go.Bar(x = h2[1][:-1], y = h2_total_patched, marker_color = colors[4]), row = 2, col = 2)


	fig.update_xaxes(title_text=("x [pitch 0.4 mm]"), row = 1, col = 1)
	fig.update_yaxes(type=yaxis_type,title_text=("counts"), row = 1, col = 1)
	
	fig.update_xaxes(title_text=("y [pitch 0.4 mm]"), row = 1, col = 2)
	fig.update_yaxes(type=yaxis_type,title_text=("counts"), row = 1, col = 2)
	
	fig.update_xaxes(title_text=("x [pitch 0.4 mm]"), row = 1, col = 3)
	fig.update_yaxes(title_text=("adc"), row = 1, col = 3)
	
	fig.update_xaxes(title_text=("y [pitch 0.4 mm]"), row = 1, col = 4)
	fig.update_yaxes(title_text=("adc"), row = 1, col = 4)

	fig.update_xaxes(title_text=("x [pitch 0.4 mm]"), row = 2, col = 1)
	fig.update_yaxes(type=yaxis_type,title_text=("counts"), row = 2, col = 1)
	
	fig.update_xaxes(title_text=("y [pitch 0.4 mm]"), row = 2, col = 2)
	fig.update_yaxes(type=yaxis_type,title_text=("counts"), row = 2, col = 2)
	
	fig.update_xaxes(title_text=("x [pitch 0.4 mm]"), row = 2, col = 3)
	fig.update_yaxes(title_text=("adc"), row = 2, col = 3)
	
	fig.update_xaxes(title_text=("y [pitch 0.4 mm]"), row = 2, col = 4)
	fig.update_yaxes(title_text=("adc"), row = 2, col = 4)
	fig.update_layout(height = the_height, width = the_width, showlegend = False, barmode='group',uirevision='constant',bargap=0,plot_bgcolor='white',paper_bgcolor='white')
	
	return fig, {
		"h1": h1_total.tolist(),
		"h2": h2_total.tolist(),
		"h3": h3_total.tolist(),
		"h4": h4_total.tolist(),
	}
