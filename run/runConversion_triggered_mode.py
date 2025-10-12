#!/usr/bin/python
import os
import subprocess
import re
import sys
	
try:	
	args = ['../build/convertFile', '-f', 'example_triggered_mode.pcapng', '-vmm', '[[1,0,2,2],[1,0,2,3],[1,0,2,4],[1,0,2,5],[1,1,2,10],[1,1,2,11],[1,1,2,12],[1,1,2,13]]',
'-axis', '[[1,0],0],[[1,1],0]', '-bc', '44.444', '-tac', '60', '-th','0', '-cs','1', '-ccs', '2', '-dt', '100', '-mst', '1', '-spc', '200', '-dp', '200', '-coin', 'center-of-masss', '-crl', '0.2', '-cru', '10', '-save', '[[1],[1],[1]]', '-json','0', '-algo', '0','-df', 'TRG','-log', 'INFO']		
	subprocess.call(args)


except OSError:
	pass

