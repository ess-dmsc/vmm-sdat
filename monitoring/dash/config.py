import os
import glob
from dash import dcc, html
from datetime import datetime
import re
import plotly.express as px

#########################
acquire_pcapng=0
delete_root=1
path="../../build"
file_name="monitoring"
file_size=20000
number_of_files = 10
#interface="enp1s0f0"
interface="en0"
geo_file="example_monitoring.json"
calib_file="empty"
cluster_algorithm="charge2"
cs="1"
ccs="3"
dt="200"
mst="1"
spc="500"
dp="500"
crl="0.1"
cru="10"
save='[[0,1,2,3],[0,1,2,3],[0,1,2,3]]'
df="SRS"
th="0"
bc="44.444444"
tac="60"
channels_x = 256
channels_y = 256
image_x = 256
image_y = 256
max_charge = 10000
charge_scale = 0.05
bins_size=20
bins_missing_strip=4
bins_span=20
bins_max_delta_hits=20
bins_delta_plane=50
color_x0 = '#1f77b4'
color_x1 = '#17becf'
color_x2 = '#9467bd'
color_x3 = '#2ca02c'
color_y0 ='#d62728'
color_y1 ='#ff7f0e'
color_y2 ='#6D071A'
color_y3 ='#e377c2'
color_xy0 = '#000000'
color_xy1 = '#800000'
color_xy2 = '#808080'
color_xy3 = '#c0c0c0'
#color_xy0 = '#000000'
#color_xy1 = '#8f8f8f'
#color_xy2 = '#555500'
#color_xy3 = '#8c564b'
the_height=1200
the_width=3000
time_points=50
##########################

color_config = [color_x0,color_x1,color_x2, color_x3,color_y0,color_y1,color_y2,color_y3,color_xy0,color_xy1,color_xy2,color_xy3 ]

shared_data = {
	"hit_data": [],
	"plane_data": [],
	"cluster_data": [],
	"hit_updated_at": 0.0,
	"plane_updated_at": 0.0,
	"cluster_updated_at": 0.0
}


palette_map = {
	"Plotly": px.colors.qualitative.Plotly,
	"D3": px.colors.qualitative.D3,
	"Set1": px.colors.qualitative.Set1,
	"Pastel": px.colors.qualitative.Pastel,
	"Dark2": px.colors.qualitative.Dark2,
	"Bold": px.colors.qualitative.Bold,
	"Safe": px.colors.qualitative.Safe,
	"Vivid": px.colors.qualitative.Vivid,
	"Alphabet": px.colors.qualitative.Alphabet,
	"Antique": px.colors.qualitative.Antique
}

color_options=[
	{"label": "Config", "value": "Config"},
	{"label": "Plotly", "value": "Plotly"},
	{"label": "D3", "value": "D3"},
	{"label": "Set1", "value": "Set1"},
	{"label": "Pastel", "value": "Pastel"},
	{"label": "Dark2", "value": "Dark2"},
	{"label": "Bold", "value": "Bold"},
	{"label": "Safe", "value": "Safe"},
	{"label": "Vivid", "value": "Vivid"},
	{"label": "Alphabet", "value": "Alphabet"},
	{"label": "Antique", "value": "Antique"}]

def get_gradient_style(colors):
	"""Build CSS gradient style from list of colors."""
	return {
		'background': f'linear-gradient(to right, {", ".join(colors)})',
		'height': '15px',
		'borderRadius': '5px',
		'margin': '5px 0'
	}

def get_colorscale_options():
	sequential = px.colors.sequential
	palette_names = [name for name in dir(sequential) if not name.startswith("_") and isinstance(getattr(sequential, name), list)]
	
	options = []
	for name in palette_names:
		colors = getattr(sequential, name)
		options.append({
			'label': html.Div([
				html.Div(style=get_gradient_style(colors)),
				html.Div(name, style={'fontSize': '12px', 'marginTop': '2px'})
			]),
			'value': name
		})
	return options
	
def count_files(extension, name_filter):
    # Find all .txt files
	file_list = glob.glob(os.path.join(".", extension))
	file_list = [f for f in file_list if name_filter in os.path.basename(f)]
	if not file_list:
		return None
	if len(file_list) < 2:
		return None
	return file_list


def get_oldest_file(file_list):
	if not file_list:
		return None
	if len(file_list) < 2:
		return None
	return min(file_list, key=os.path.getmtime)


def cleanup_old_root_files(root_dir, max_files=10):
	# Get all .root files with full paths
	root_files = [os.path.join(root_dir, f) for f in os.listdir(root_dir) if f.endswith(".root")]
	#root_files = [os.path.join(root_dir, f) for f in os.listdir(root_dir) if f.endswith(".root") and os.path.isfile(os.path.join(root_dir, f))]
	#print(root_files)
	if len(root_files) <= max_files:
		return  # Nothing to clean up

	# Sort by filename timestamp (reuse your existing function)
	sorted_files = sorted(root_files, key=extract_timestamp_from_filename)

	# Files to delete
	files_to_delete = sorted_files[:-max_files]  # keep the newest `max_files`

	for f in files_to_delete:
		try:
			os.remove(f)
			print(f"Deleted old file: {f}")
		except Exception as e:
			print(f"Could not delete {f}: {e}")



def extract_timestamp_from_filename(file):
	match = re.search(r"_(\d{14})_", file)
	if match:
		return datetime.strptime(match.group(1), "%Y%m%d%H%M%S")
	else:
		return datetime.max  # push malformed names to the end
