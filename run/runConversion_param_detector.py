#!/usr/bin/python
import os
import subprocess
import re
import sys
	
try:	
	args = ['../build/convertFile', '-f', 'example.pcapng', '-geo', 'example_geometry.json', '-bc', '44.444', '-tac', '60', '-th','[0,1]', '-cs','[1,2]', '-ccs', '[3,4]', '-dt', '[100,200]', '-mst', '[1,2]', '-spc', '[200,300]', '-dp', '[200,500]', '-coin', 'center-of-masss', '-crl', '[0.2,0.4]', '-cru', '[10,5]','-save', '[[1],[1],[1]]', '-json','0', '-algo', '0', '-info', 'EXAMPLE', '-df', 'SRS','-log', 'INFO']	
	subprocess.call(args)


except OSError:
	pass




