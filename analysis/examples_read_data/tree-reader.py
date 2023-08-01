#!/usr/bin/python3

# VMM Analysis
# --------------------------------------------------------------------
# This script is a simple example, showing how to read the data from a
# ROOT tree, generated with vmm-sdat. Furthermore a simple plot of the
# data, using matplotlib, is generated.
# --------------------------------------------------------------------
# Lucian Scharenberg
# lucian.scharenberg@cern.ch
# 18 November 2019 and 03 March 2022


# --------------------------------------------------------------------
# PACKAGE HANDLING
# --------------------------------------------------------------------

import uproot3 as uproot
import matplotlib.pyplot as plt
import sys


# --------------------------------------------------------------------
# DATA HANDLING
# --------------------------------------------------------------------

# Get the tree in the ROOT file
tree = uproot.open(sys.argv[1])['clusters_detector']

# Now get the branch of interest
adc0 = tree.array('adc0')


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
