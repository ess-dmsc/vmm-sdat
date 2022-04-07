#!/usr/bin/python
import os
import subprocess
import re
import sys

files = ["./example_LET.pcapng"]


for file in files:	
	try:	
		args = ['../build/convertFile', '-f', file, '-geo', '../run/LET_geometry.json',
	'-bc', '44.44', '-tac', '60', '-th','0', '-cs','1', '-ccs', '2', '-dt', '200', '-mst', '0', '-spc', '500', '-dp', '500', '-coin', 'center-of-masss', '-crl', '0.1', '-cru', '10', '-save', '[[0,1],[0,1],[0,1]]', '-json','0', '-algo', '0', '-info', 'LET', '-df', 'ESS']	
		subprocess.call(args)


	except OSError:
		pass
