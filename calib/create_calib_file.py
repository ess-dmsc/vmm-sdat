import json
import argparse
text_description =  "JSON calibration file manipulator. The tool create_calib_file.py creates a calibration file for a whole system of multiple hybrids. The user has to specify the JSON mapping file that describes the system (mapping between hybrid, FEC and VMM), the name for the new combined system calibration file, the directory with the calibration files of the individual hybrids, and the  choice of calibrations (ADC, time, time walk)."
parser = argparse.ArgumentParser(description=text_description,
								 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("geo", help="hybrid geometry description")
parser.add_argument("-d", "--directory", default=".", help="directory with hybrid calibration files")
parser.add_argument("-c", "--calib", default="calib.json", help="new calibration file")
parser.add_argument("-a", "--adc", action="store_true", help="add ADC calib")
parser.add_argument("-t", "--time", action="store_true", help="add time calib")
parser.add_argument("-w", "--timewalk", action="store_true", help="add time-walk calib")

args = parser.parse_args()
config = vars(args)
name = config["geo"]
calib = config["calib"]
directory = config["directory"]
is_adc = config["adc"]
is_time = config["time"]
is_time_walk = config["timewalk"]
lastId=0

with open(calib, 'w') as json_file:
	with open(name, "r") as geo_file:
		data  = json.load(geo_file)
		#print("\n************** File content " + name + " ******************")
		#print(json.dumps(data, sort_keys=True, indent=4))
		#print("********************************\n")
		for d in data["hybrid_mapping"]:
			if lastId != d["hybridID"]:
				vmm=0
			else:
				vmm=1
			if is_time:
				y = {"time_slopes":[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]}
				d.update(y)
				y = {"time_offsets":[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]}
				d.update(y)
			if is_adc:
				y = {"adc_slopes":[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]}
				d.update(y)
				y = {"adc_offsets":[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]}
				d.update(y)
			if is_time_walk:
				y = {"twalk_offsets":0}
				d.update(y)
				y = {"twalk_scales":0}
				d.update(y)
				y = {"twalk_means":0}
				d.update(y)	
				y = {"twalk_sigmas":0}
				d.update(y)				
			if is_adc:
				adc_file_name = "vmm_adc_calibration_" + d["hybridID"] +".json"
				with open(adc_file_name, "r") as adc_file:
					adc  = json.load(adc_file)
					for a in adc["vmm_calibration"]:
						if a["vmmID"] == vmm:
							d["adc_slopes"] = a["adc_slopes"]
							d["adc_offsets"] = a["adc_offsets"]
			if is_time:
				time_file_name = "vmm_time_calibration_" + d["hybridID"] +".json"
				with open(time_file_name, "r") as time_file:
					time  = json.load(time_file)
					for a in time["vmm_calibration"]:
						if a["vmmID"] == vmm:
							d["time_slopes"] = a["time_slopes"]
							d["time_offsets"] = a["time_offsets"]
			if is_time_walk:
				time_file_name = "vmm_twalk_calibration_" + d["hybridID"] +".json"
				with open(time_file_name, "r") as time_file:
					time  = json.load(time_file)
					for a in time["vmm_calibration"]:
						if a["vmmID"] == vmm:
							d["twalk_offsets"] = a["twalk_offsets"]
							d["twalk_scales"] = a["twalk_scales"]
							d["twalk_means"] = a["twalk_means"]
							d["twalk_sigmas"] = a["twalk_sigmas"]					
			lastId = d["hybridID"]
			#print(d["hybridID"])
		data["vmm_calibration"] = data.pop("hybrid_mapping")	
		#json.dump(data,json_file,indent=4, sort_keys=True,ensure_ascii=False)
		json.dump(data,json_file,sort_keys=True)
