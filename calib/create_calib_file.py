import json
import argparse
import os
text_description =  "JSON calibration file manipulator. The tool create_calib_file.py creates a calibration file for a whole system of multiple hybrids. The user has to specify the JSON mapping file that describes the system (mapping between hybrid, FEC and VMM), the name for the new combined system calibration file, the directory with the calibration files of the individual hybrids, and the  choice of calibrations (ADC, time, time walk)."
parser = argparse.ArgumentParser(description=text_description,
								 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("geo", help="hybrid geometry description")
parser.add_argument("-d", "--directory", default=".", help="directory with hybrid calibration files")
parser.add_argument("-c", "--calib", default="calib.json", help="new calibration file")
parser.add_argument("-a", "--adc", action="store_true", help="add ADC calib")
parser.add_argument("-t", "--time", action="store_true", help="add time calib")
parser.add_argument("-w", "--timewalk", action="store_true", help="add timewalk calib")

args = parser.parse_args()
config = vars(args)
name = config["geo"]
calib = config["calib"]
directory = config["directory"]
is_adc = config["adc"]
is_time = config["time"]
is_timewalk = config["timewalk"]
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
			if is_timewalk:
				y = {"timewalk_a":[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]}
				d.update(y)
				y = {"timewalk_b":[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]}
				d.update(y)
				y = {"timewalk_c":[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]}
				d.update(y)	
				y = {"timewalk_d":[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]}
				d.update(y)				
			if is_adc:
				for filename in os.listdir(directory):
					search_name = "vmm_adc_calibration_" + d["hybridID"]
					if filename.startswith(search_name) and filename.endswith(".json"): 
						with open(filename, "r") as file:
							the_calib  = json.load(file)
							for c in the_calib["vmm_calibration"]:
								if c["vmmID"] == vmm:
									d["adc_slopes"] = c["adc_slopes"]
									d["adc_offsets"] = c["adc_offsets"]
			if is_time:
				for filename in os.listdir(directory):
					search_name = "vmm_time_calibration_" + d["hybridID"]
					if filename.startswith(search_name) and filename.endswith(".json"): 
						with open(filename, "r") as file:
							the_calib  = json.load(file)
							for c in the_calib["vmm_calibration"]:
								if c["vmmID"] == vmm:
									d["time_slopes"] = c["time_slopes"]
									d["time_offsets"] = c["time_offsets"]
			if is_timewalk:
				for filename in os.listdir(directory):
					search_name = "vmm_timewalk_calibration_" + d["hybridID"] +".json"
					if filename.startswith(search_name) and filename.endswith(".json"): 
						with open(filename, "r") as file:						
							the_calib  = json.load(file)
							for c in the_calib["vmm_calibration"]:
								if c["vmmID"] == vmm:
									d["timewalk_a"] = c["timewalk_a"]
									d["timewalk_b"] = c["timewalk_b"]
									d["timewalk_c"] = c["timewalk_c"]
									d["timewalk_d"] = c["timewalk_d"]					
			lastId = d["hybridID"]
			#print(d["hybridID"])
		data["vmm_calibration"] = data.pop("hybrid_mapping")	
		#json.dump(data,json_file,indent=4, sort_keys=True,ensure_ascii=False)
		json.dump(data,json_file,sort_keys=True)
