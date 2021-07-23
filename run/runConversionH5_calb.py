#!/usr/bin/python
import os
import subprocess
import re
import sys
	
try:	
	args = ['../build/convertFile', '-f', 'example.h5', '-vmm', '[[1,0,2,2],[1,0,2,3],[1,0,2,0],[1,0,2,1],[1,1,2,8],[1,1,2,9],[1,1,2,6],[1,1,2,7]]',
'-axis', '[[1,0],0],[[1,1],0]', '-bc', '40', '-tac', '60', '-th','0', '-cs','2', '-ccs', '3', '-dt', '200', '-mst', '3', '-spc', '500', '-dp', '200', '-coin', 'center-of-masss', '-crl', '0.5', '-cru', '10', '-save', '[[],[],[1]]', '-json','0', '-algo', '0', '-info', 'CALIB', '-cal', './example_calib.json']	
	subprocess.call(args)


except OSError:
	pass




