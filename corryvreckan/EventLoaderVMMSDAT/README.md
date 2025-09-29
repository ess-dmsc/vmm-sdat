# EventLoaderVMMSDAT

**Maintainer:** [dorothea.pfeiffer@cern.ch](mailto:dorothea.pfeiffer@cern.ch)  
**Module Type:** *GLOBAL*  
**Status:** *Immature*  

---

## Overview
`EventLoaderVMMSDAT` is a Corryvreckan module for loading pre-processed VMM3a data.  
It is derived from the `EventLoaderMPGD` module (by Maryna Borysova, <maryna.borysova@weizmann.ac.il>).

The module reads ROOT files produced by the [vmm-sdat](https://github.com/ess-dmsc/vmm-sdat) analysis program.  
`vmm-sdat` generates ROOT files containing three trees:

- `hits`
- `clusters_plane`
- `clusters_detector`

This module uses the **`clusters_detector`** tree, which contains clusters already matched between the X and Y planes.

---

## Features
- Loads cluster information from ROOT trees.
- Associates clusters with detectors defined in the configuration.
- Produces histograms for monitoring:
  - Total number of clusters per event.
  - Number of clusters per detector.
  - 2D cluster maps in local and global coordinates.

---

## Configuration Parameters

Add the following section to your Corryvreckan steering file:

```toml
[EventLoaderVMMSDAT]

# General settings
log_level = "STATUS"
input_file = "/path/to/sorted_by_time0_100events.root"
tree_name = "clusters_detector"

# Number of clusters to skip/read
number_clusters_to_skip = 0
number_clusters_to_read = -1   # -1 = read all

# Event building
time_window = 2000             # Time window (ns) to group clusters into events
sort_clusters = true           # Sort clusters by time (recommended)

# Detector mapping
detector_map = [[1, "GEMXY1"], [2, "GEMXY2"], [3, "GEMXY3"], [4, "GEMXY4"]]

# Detector geometry corrections (pitch, offsets)
detector_position_scale_map   = [["GEMXY1", 0.4], ["GEMXY2", 0.4], ["GEMXY3", 0.4], ["GEMXY4", 0.4]]
detector_position_shift_x_map = [["GEMXY1", 0.0], ["GEMXY2", 0.0], ["GEMXY3", 0.0], ["GEMXY4", 0.0]]
detector_position_shift_y_map = [["GEMXY1", 0.0], ["GEMXY2", 0.0], ["GEMXY3", 0.0], ["GEMXY4", 0.0]]

# Trigger and required detectors
detector_trigger   = "GEMXY3"
detector_required  = ["GEMXY1", "GEMXY2", "GEMXY3"]

# Algorithms for position/time/charge/size
position_algorithm = "cog"     # [cog | charge2 | utpc | algo]
time_algorithm     = "cog"     # [cog | charge2 | utpc | algo]
time_choice        = "plane0"  # [plane0 | plane1 | both]
charge_choice      = "plane0"  # [plane0 | plane1 | both]
size_choice        = "plane0"  # [plane0 | plane1 | both]
```

---

## Output Histograms

- **Global**
  - Total clusters per event (`NClusterTotalEvent`)
- **Per Detector**
  - `ClusterMap`: 2D histogram of cluster hit positions in **local coordinates** (pixels).
  - `ClusterMap_global`: 2D histogram of cluster hit positions in **global coordinates** (mm).
  - `NClustersPerEvent`: 1D histogram of clusters per event.

---

## Example Usage

Minimal steering file:

```toml
[EventLoaderVMMSDAT]
input_file = "/data/vmm-sdat/sorted_by_time0_100events.root"
tree_name  = "clusters_detector"
log_level  = "STATUS"
```

Run Corryvreckan with:

```bash
corry -c your_config.conf
```

---

## Notes
- The module requires a valid `clusters_detector` tree in the ROOT input file.  
- Geometry corrections (scale, shift) should ideally come from detector descriptions, but are currently configured manually.  
- The code creates *dummy pixels* for clusters to ensure compatibility with other Corryvreckan modules.
