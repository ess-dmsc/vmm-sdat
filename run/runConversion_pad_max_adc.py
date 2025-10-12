#!/usr/bin/python
import os
import subprocess
import re
import sys

try:
        args = ['../build/convertFile', '-f', './example_pad.pcapng', '-geo', 'example_geometry_pad.json', '-bc', '44.444', '-tac', '60', '-th','0', '-cs','1', '-dt', '200', '-mp0', '1', '-mp1','1', '-spc', '1000', '-dp', '100', '-save', '[[1],[0],[1]]', '-algo', '10', '-info', 'max_adc_pads', '-df','SRS', '-log', 'INFO']
        subprocess.call(args)

except OSError:
	pass