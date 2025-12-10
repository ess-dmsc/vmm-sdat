#!/usr/bin/python
import os
import subprocess
import re
import sys

try:
	args = ['../build/convertFile', '-f', 'example_rpc.pcapng', '-geo', 'example_rpc.json', '-bc', '44.444', '-tac', '60', '-th','0', '-cs','1', '-ccs', '3', '-dt', '100', '-mst', '[1,1,1,1,1]', '-spc', '200', '-dp', '200', '-coin', 'center-of-masss', '-crl', '0.2', '-cru', '10', '-save', '[[0,1,2,3,10],[0,1,2,3,10],[0,1,2,3,10]]', '-algo', '4', '-info', 'rpc', '-df', 'SRS','-log', 'INFO']	       
	subprocess.call(args)

except OSError:
	pass