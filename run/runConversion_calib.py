#!/usr/bin/python
import os
import subprocess
import re
import sys
	
try:	
	args = ['../build/convertFile', '-f', 'example.pcapng', '-vmm', '[[1,0,1,0],[1,0,1,1],[1,1,2,0],[1,1,2,1]]',
'-axis', '[[1,0],0],[[1,1],0]', '-bc', '44.444', '-tac', '60', '-th','0', '-cs','1', '-ccs', '3', '-dt', '100', '-mst', '1', '-spc', '200', '-dp', '200', '-coin', 'center-of-masss', '-crl', '0.2', '-cru', '10', '-save', '[[1],[1],[1]]', '-json','0', '-algo', '0', '-info', 'CALIB', '-cal', './example_calib.json', '-df', 'SRS', '-cahi','1','-log', 'INFO']		
	subprocess.call(args)


except OSError:
	pass

