#!/usr/bin/python
import os
import subprocess
import re
import sys
	
try:	
	args = ['../build/convertFile', '-f', 'example_maxiroc.pcapng', '-geo', 'example_maxiroc.json', '-bc', '44.444', '-tac', '60', '-th','0', '-cs','1', '-ccs', '3', '-dt', '100', '-mst', '1', '-spc', '200', '-dp', '200', '-coin', 'center-of-masss', '-crl', '0.2', '-cru', '10', '-save', '[[1],[1],[1]]', '-json','0', '-algo', '0', '-df', 'MAX','-log', 'INFO']		
	subprocess.call(args)


except OSError:
	pass
