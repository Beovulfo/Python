"""
Module for converting binary file for field XZ generated by TOFIS_LOMA
into HDF5 compressed file
"""
def field2hdf5(filename):
	""" 
	This function reads the XZ field generated by TOFIS fortran program and 
	converts it to same filename + .h5, using gzip compression. Furthermore
	this new file includes: xvec,zvec,y,jspecy,Re,time and the field itself.
	{y,xvec,zvec,time} = readfieldxz(filename)
	"""
	import numpy as np
  	import h5py
	import os.path
	
	if os.path.isfile(filename+'.h5')==True:
		print "This file already exists man!Nothing done."
		return
	f = open(filename,'rb')
	#Create dtypes for proper reading from Fortran unformatted
	# binary file
	#Declaring types
	yfmap = np.dtype([('y','float64'),('fmap','float64')])
	#uw00 = np.dtype([('u00','float32'),('w00','float32')])
	rec1 = np.dtype([('dummy1','uint32'), \
                ('time','float32'),('Re','float32'), \
		('alp','float32'),('bet','float32'), \
		('mgalx','uint32'),('my','uint32'),('mgalz','uint32'),\
		('nspec','uint32'),('nacum','uint32'),\
		('dummy2','uint32')])
	#Read first record
	RECORD1=np.fromfile(f,rec1,1)

	#Check if first record is ok...
	if RECORD1['dummy1'] != RECORD1['dummy2']:
		print "File read not good...!"
		return
	else:
    		print "File RECORD1 read correctly :)"

	mgalx=RECORD1['mgalx'][0]
	my=RECORD1['my'][0]
	mgalz=RECORD1['mgalz'][0]
	nspec=RECORD1['nspec'][0]


	print "nspec= %s" % nspec

	rec2 = np.dtype([('dummy1','uint32'), \
 		('jspecy','uint32',nspec),  \
		('yfmap',yfmap,my),	\
		('dummy2','uint32')])
	#READ record2
	RECORD2=np.fromfile(f,rec2,1)

	#Check if record2 is ok...
	if RECORD2['dummy1'] != RECORD2['dummy2']:
		print "File read not good...!"
		return
	else:
    		print "File RECORD2 read correctly :)"

	#Save y vector amd jspecy
	y = RECORD2['yfmap']['y'][0,]
        jspecy = RECORD2['jspecy'][0,]
        #Define plane info
	planey = np.ndarray(shape=(mgalx,mgalz),order='F',dtype=float)

	#Create type "recplane"
	recplaney = np.dtype([('dummy1','uint32'), \
                 ('data','float32',[mgalz,mgalx]), \
                 ('dummy2','uint32')])


	#Read all planes Y info
	FIELD1=np.ndarray(shape=(nspec,mgalz,mgalx),\
	dtype=float)
	for j in range(nspec):
    		readrec = np.fromfile(f,recplaney,1)
    		planeydata = readrec['data']
		FIELD1[j,:,:] = planeydata[:,:]
	f.close()
        #FIELD1.shape=(nspec,mgalz,mgalx)
        #Create vector X and Z
	Lx = 2*3.1415/RECORD1['alp']
	Lz = 2*3.1415/RECORD1['bet']
	x = np.linspace(-Lx/2.,Lx/2.,mgalx)
	z = np.linspace(-Lz/2.,Lz/2.,mgalz)


	hf5 = h5py.File(filename+'.h5','w')
	hf5.create_dataset('field',data=FIELD1,compression='gzip')
	hf5.create_dataset('xvec',data=x)
	hf5.create_dataset('zvec',data=z)
	hf5.create_dataset('y',data=y)
	hf5.create_dataset('time',data=RECORD1['time'])
	hf5.create_dataset('Re',data=RECORD1['Re'])
	hf5.create_dataset('jspecy',data=jspecy)
	hf5.close()
	
        del FIELD1
        print 'Data from time = %s' % RECORD1['time']
        print 'mgalx = %s, my = %s, mgalz = %s' % (RECORD1['mgalx'], \
		RECORD1['my'],RECORD1['mgalz'])
	#Return y, FIELD
	print "File conversion finished. Congratulations"
	return  

