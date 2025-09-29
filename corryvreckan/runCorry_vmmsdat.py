#!/usr/bin/python
import os
import subprocess
import re
import sys
	
try:	
	args = ['/Users/dpfeiffe/programming/vmm-data/corryvreckan/bin/corry', '-c', './EventLoaderVMMSDAT.conf']	
	
	subprocess.call(args)

except OSError:
	pass




