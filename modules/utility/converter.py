from field2hdf5 import field2hdf5
import numpy as np
import glob,os
import os.path
import time

#MAIN PROGRAM
#-------------------------------------------------------------------
#Define path to FOLDER 
#workdir  = "/data/toni/"
workdir  = "/share/drive/toni/Re160s10/case2/"
#workdir  = "/share/drive/toni/Re160s80/case1/SS/"
#workdir  = "/data2/toni/turb/"
#Define main name of file
fname_main = "s10b_" 
#fname_main = "s80my1101resx2_" 
#Define variable to convert
#Options available:
#Tfxz,upxz,vpxz,wpxz,o1xz,o2xz,o3xz
#variable = ".Tfxz"
#variable = ".o2xz"
#variable = ".o3xz"
variable_v = ['.upxz','.vpxz','.wpxz','.Tfxz','.o1xz','.o2xz','.o3xz']
ifirst = 1
ilast = 17
#Show file list before proceeding
print "File list to convert: \n"
for ifile in range(ifirst,ilast+1):
	for variable in variable_v:
		filename = workdir + fname_main + '%03d' %ifile + variable
		print "%s \n" %(filename)

#Let check list before proceeding
print "Waiting 10 s before to proceed, Ctrl+C to cancel program."
time.sleep(10.0)

for ifile in range(ifirst,ilast+1):
	for variable in variable_v:
		#filemame is the complete path to file
		filename = workdir + fname_main + '%03d' %ifile + variable
		print "Working with file = %s" %(filename)
		field2hdf5(filename)
		#remove binary field after conversion
		if os.path.isfile(filename) == True:
			os.remove(filename)

print "Conversion of '%s' files finished!" % (variable)


