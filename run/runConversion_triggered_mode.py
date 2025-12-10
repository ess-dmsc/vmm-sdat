#!/usr/bin/python
import os
import subprocess
import re
import sys
	
try:	
	args = ['../build/convertFile', '-f', 'example_triggered_mode.pcapng','-geo', 'example_geometry_trigger.json', '-bc', '44.444', '-tac', '60', '-th','0', '-cs','1', '-ccs', '2', '-dt', '100', '-mst', '1', '-spc', '200', '-dp', '200', '-coin', 'center-of-masss', '-crl', '0.2', '-cru', '10', '-save', '[[1],[1],[1]]', '-json','0', '-algo', '0','-df', 'TRG','-log', 'INFO']		
	subprocess.call(args)


except OSError:
	pass