#!/usr/bin/python
import os
import subprocess
import re
import sys

try:
        args = ['../build/convertFile', '-f', './example_xyu.pcapng', '-geo', './example_geometry_xyu.json', '-bc', '40', '-tac', '60', '-th','0', '-cs','1', '-ccs', '3', '-dt', '200', '-mst', '1', '-spc', '500', '-dp', '200', '-coin', 'center-of-masss', '-crl', '0.2', '-cru', '10', '-save', '[[0],[0],[0]]', '-algo', '4', '-info', 'xyu', '-df','SRS','-log', 'INFO']	
        subprocess.call(args)

except OSError:
	pass