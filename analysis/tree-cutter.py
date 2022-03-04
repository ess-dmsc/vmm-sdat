#!/usr/bin/python3

# VMM Analysis
# --------------------------------------------------------------------
# This script is a simple example, showing how to read the data from a
# ROOT tree, generated with vmm-sdat. In addition, some cuts are
# applied to the data using pandas. In this specific example, this
# means that only the position of the clusters is plotted, if the
# ADC value of the cluster is larger than a specific value.
# --------------------------------------------------------------------
# Lucian Scharenberg
# lucian.scharenberg@cern.ch
# 18 November 2019 and 03 March 2022


# --------------------------------------------------------------------
# PACKAGE HANDLING
# --------------------------------------------------------------------

import uproot3 as uproot
import pandas as pd
import matplotlib.pyplot as plt
import sys


# --------------------------------------------------------------------
# DATA HANDLING
# --------------------------------------------------------------------

# Get the tree in the ROOT file
tree = uproot.open(sys.argv[1])['clusters_detector']

# Now get the branches of interest
adc0 = tree.array('adc0')
pos0 = tree.array('pos0')
pos1 = tree.array('pos1')

# Create a pandas data frame, which is used to apply the cuts
data = {'adc0': adc0,
        'pos0': pos0,
        'pos1': pos1}
df = pd.DataFrame(data)

# Get only the events, which 'survive' the cut
events = df.query('adc0 > 500')


# --------------------------------------------------------------------
# PLOT THE DATA
# --------------------------------------------------------------------

# Create the plot of the positions
plt.scatter(events['pos0'], events['pos1'], s = 0.05)

# Some labels
plt.xlabel('Position 0 (Strip Numbers)')
plt.xlabel('Position 1 (Strip Numbers)')

# Save the plot and show it
plt.savefig('tree-cutter.png', dpi = 500)
plt.show()
