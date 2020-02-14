#!/usr/bin/python3

# VMM HDF5-to-ROOT Analysis
# --------------------------------------------------------------------
# This script is a simple example, showing how to read the data from a
# ROOT tree, generated with vmm-hdf5-to-root from Dorothea Pfeiffer.
# Furthermore a simple plot of the data, using matplotlib, is
# generated.
# --------------------------------------------------------------------
# Lucian Scharenberg
# lucian.scharenberg@cern.ch
# 18 November 2019


# --------------------------------------------------------------------
# PACKAGE HANDLING
# --------------------------------------------------------------------

import uproot
import matplotlib.pyplot as plt
import sys


# --------------------------------------------------------------------
# DATA HANDLING
# --------------------------------------------------------------------

# Get the tree in the ROOT file
tree = uproot.open(sys.argv[1])['events']

# Now get the branch of interest
adc0 = tree.array('clusters_detector.adc0').flatten()


# --------------------------------------------------------------------
# PLOT THE DATA
# --------------------------------------------------------------------

# Generate a histogram
plt.hist(adc0, bins = 2500, range = [0.0, 2500.0], histtype = 'step')

# Some labels
plt.xlabel('Energy (ADC Counts)')
plt.ylabel('Occurrence')

# Save the plot and show it
plt.savefig('tree-reader.png', dpi = 500)
plt.show()
