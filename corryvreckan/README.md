# Corryvreckan

`Corryvreckan` is a versatile, highly configurable software with a modular structure designed to reconstruct and analyse test beam and laboratory data. It is especially used for tracking.  
More information: [Corryvreckan on GitLab](https://gitlab.cern.ch/corryvreckan/corryvreckan)

---

## Use of vmm-sdat Data in Corryvreckan

With the modules provided in **EventLoaderVMMSDAT**, data clustered in [vmm-sdat](https://github.com/ess-dmsc/vmm-sdat) can be read into Corryvreckan.  
For more details, please refer to the README file in the `EventLoaderVMMSDAT` subfolder.

---

## Prerequisites

- A complete working installation of **Corryvreckan**.

---

## Example

This folder contains the following resources for testing the module:

- **`runCorry_vmmsdat.py`** – Python script to run the Corryvreckan executable (`corry`) with the configuration file `EventLoaderVMMSDAT.conf`.
- **`EventLoaderVMMSDAT.conf`** – Example configuration file.
- **`GEM_telescospe.geo`** – Example geometry file.
- **`output/`** – Empty folder that will contain histograms after running the example.
- **`example.root`** – Example vmm-sdat data file to be read into Corryvreckan.
- **`sort_tree.cpp`** – Utility script to pre-sort vmm-sdat data by timestamp.  
  By default, `EventLoaderVMMSDAT` sorts clusters during runtime. For repeated runs on large files, pre-sorting with this script can significantly improve performance.

---
