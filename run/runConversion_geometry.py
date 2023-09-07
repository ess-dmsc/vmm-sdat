#!/usr/bin/python
import os
import subprocess
import re
import sys
	
try:	
	args = ['../build/convertFile', '-f', 'example.pcapng', '-geo', 'example_geometry.json',
'-axis', '[[1,0],0],[[1,1],0]', '-bc', '44.444', '-tac', '60', '-th','0', '-cs','1', '-ccs', '3', '-dt', '100', '-mst', '1', '-spc', '200', '-dp', '200', '-coin', 'center-of-masss', '-crl', '0.2', '-cru', '10', '-save', '[[1],[1],[1]]', '-json','0', '-algo', '0', '-info', 'CALIB', '-cal', './example_calib.json', '-df', 'SRS']	
	subprocess.call(args)


except OSError:
	pass
