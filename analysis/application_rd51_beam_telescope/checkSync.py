#!/usr/bin/python3

# CHECK FEC SYNCHRONISATION
# --------------------------------------------------------------------
# This script allows to check the synchronisation between FECs of
# VMM3a/SRS. The method behind is the same as for time resolution
# studies: one detector functions as a time reference, while it is
# checked in the other detectors, if they have matching recorded
# interactions within a certain time window around the reference time.
# --------------------------------------------------------------------
# Lucian Scharenberg
# lucian.scharenberg@cern.ch
# 18 August 2023


# --------------------------------------------------------------------
# PACKAGE HANDLING
# --------------------------------------------------------------------

import uproot3 as uproot
import pandas as pd
import numpy as np
import plotly.express as px
import argparse


# --------------------------------------------------------------------
# DATA HANDLING
# --------------------------------------------------------------------

def time_difference(reference, dut):

    results = []
    
    i = 0
    k = 0

    while i < len(reference):
    
        if i % 10000 == 0:
    
            print(i)
    
        while k < len(dut):

            difference = dut[k] - reference[i]
                        
            if difference <= -500.0:

                k += 1
                continue

            elif (difference > -500.0 and difference < 500.0):

                results.append(difference)

                k += 1
                i += 1
                break

            elif difference >= 500.0:

                i += 1
                break

            elif i == len(reference) or k == len(dut):

                break

        if i == len(reference) or k == len(dut):

            break

    return np.asarray(results)


def get_data(filename):

    clusters = uproot.open(filename)['clusters_detector']
    time = clusters.array('time0')
    det = clusters.array('det')

    return pd.DataFrame({'time': time, 'det': det}, dtype = np.float64)


def main_process():
    
    parser = argparse.ArgumentParser(description = 'Check the FEC synchronisation')
    #parser.add_argument('-i',
    #                    help = 'Array of the last bytes of the FEC IP addresses. \
    #                            Example: In case of FECs with IPs 10.0.0.2 and 10.0.0.10, it would be \"2, 10\".',
    #                    required = True)
    parser.add_argument('-d',
                        help = 'Array of the detectors connected to the FECs (has to be matched with the FEC IPs). \
                                Example: In case det. 1 is connected to 10.0.0.2 and det. 2 is connected to 10.0.0.10, it would be \"1, 2\".',
                        required = True)
    parser.add_argument('-f',
                        help = 'ROOT file, that has to be read in.',
                        required = True)
    args = vars(parser.parse_args())
    print(args)

    filename = args['f']
    #ip_addresses = [int(i) for i in args['i'].split(' ')]
    detectors = [int(i) for i in args['d'].split(' ')]

    df = get_data(filename)
    
    reference = np.sort(df.query('det == @detectors[0]')['time'].values)

    results_time = []
    results_info = []

    for index, item in enumerate(detectors[1:]):

        print(item)
        temp = np.sort(df.query('det == @item')['time'].values)
        results_time.append(time_difference(reference, temp))
        #results_info.append(['Det = ' + str(item) + ', IP = ' + str(ip_addresses[index])] * len(results_time[index]))
        results_info.append([str(item)] * len(results_time[index]))

    return pd.DataFrame(dict(series = np.concatenate(results_info), data = np.concatenate(results_time)))


# --------------------------------------------------------------------
# PLOT THE DATA
# --------------------------------------------------------------------

fig = px.histogram(main_process(), x = 'data',
                   color = 'series', barmode = 'overlay',
                   labels = {'data': 'Time Difference / ns',
                             #'series': 'Identifier'})
                             'series': 'Detector'})
fig.show()
