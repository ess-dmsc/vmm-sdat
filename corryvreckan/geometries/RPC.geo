[W001_GEM1]
coordinates = "cartesian"
material_budget = 0.01068
number_of_pixels = 256, 256
orientation = 0deg,0deg,0deg
orientation_mode = "xyz"
pixel_pitch = 400um,400um
position = 0um,0um,0um
role = "reference"
spatial_resolution = 115um,115um
time_resolution = 15ns
type = "gem"

[TRIGGER] # use only one channel from RPC 
coordinates = "cartesian"
material_budget = 0.01068
orientation_mode = "xyz"
orientation = 0deg,0deg,0deg
position = 0um,0um,0um
role = "auxiliary"
time_resolution= 5ns
type = "trigger"

[W002_GEM2]
coordinates = "cartesian"
material_budget = 0.01068
number_of_pixels = 256, 256
orientation_mode = "xyz"
orientation = 0deg,0deg,+90deg
pixel_pitch = 400um,400um
position = 0um,0um,-547mm
spatial_resolution = 115um,115um
time_resolution = 15ns
type = "gem"


[W003_GEM3]
coordinates = "cartesian"
material_budget = 0.01068
number_of_pixels = 256, 256
orientation_mode = "xyz"
orientation = 0deg,0deg,-90deg
pixel_pitch = 400um,400um
position = 0um,0um,-958mm
spatial_resolution = 115um,115um
time_resolution = 15ns
type = "gem"

[W0014_RPC]
coordinates = "cartesian"
material_budget = 0.01068
number_of_pixels = 16,1 
orientation = 0deg,0deg,0deg
#pixel_pitch = 100000um,6200um
pixel_pitch = 620um,1000um
orientation_mode = "xyz"
position = 0um,0um,300mm
spatial_resolution = 1790um,28868um
role = "dut"
time_resolution = 15ns 
type = "rpc"

