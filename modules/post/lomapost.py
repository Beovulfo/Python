import glob
import os
from readfiles import *
import matplotlib.pyplot as plt
import pandas as pd
import h5py
#import matplotlib as mpl
import numpy as np
from IPython.display import clear_output
from IPython.display import Image

def set_fig_props():
	plt.rcParams['figure.figsize'] = (10.0, 8.0)
	plt.rcParams['savefig.dpi'] = 300
	plt.rcParams['figure.facecolor'] = 'white'
	sizefont= 30
	font = {'family' : 'normal','weight' : 'normal','size'   : sizefont}
	plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern'],'size':sizefont})
	plt.rc('text', usetex=True)
	plt.rc('xtick', labelsize=sizefont) 
	plt.rc('ytick', labelsize=sizefont) 
	plt.rcParams['xtick.major.pad']='8'
	plt.rcParams['ytick.major.pad']='8'

def set_fig_props2():
        plt.rcParams['figure.figsize'] = (10.0, 8.0)
        plt.rcParams['savefig.dpi'] = 300
        plt.rcParams['figure.facecolor'] = 'white'
        sizefont= 12
        font = {'family' : 'normal','weight' : 'normal','size'   : sizefont}
        plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern'],'size':sizefont})
        plt.rc('text', usetex=True)
        plt.rc('xtick', labelsize=sizefont)
        plt.rc('ytick', labelsize=sizefont)
        plt.rcParams['xtick.major.pad']='8'
        plt.rcParams['ytick.major.pad']='8'


def create_l_files(p_folders):
	print p_folders                    #print path folders
	njobs = len(p_folders)          #save number of jobs in njobs
	l_sta_files = list(range(njobs))                   #initialize l_sta_files
	l_spe_files = list(range(njobs))  
	#Save all sta files list on l_sta_files for every job
	for ijob in range(njobs):
	    l_sta_files[ijob] = sorted(glob.glob(p_folders[ijob]+'/*.sta'))
	    l_spe_files[ijob] = sorted(glob.glob(p_folders[ijob]+'/*.spe'))
	    #Checking:
	    #print "job(%s): 1st file = %s, last file = %s" %(ijob+1,l_sta_files[ijob][0],l_sta_files[ijob][-1])   
	    #print "job(%s): 1st file = %s, last file = %s" %(ijob+1,l_spe_files[ijob][0],l_spe_files[ijob][-1])  
#Total number of files for every job
	nfiles=list(range(njobs))
	nfiles2=list(range(njobs))
	for ijob in range(njobs):
    		nfiles[ijob] = len(l_sta_files[ijob])
    		nfiles2[ijob] = len(l_spe_files[ijob])
#Same files STA-SPE!!
	return l_sta_files,l_spe_files,nfiles,nfiles2

def read_l_files(l_sta_files,l_sta_opt):
	njobs = np.shape(l_sta_files)[0]
	stats = []
	for ijob in range(njobs):
    		for j in range(len(l_sta_files[ijob])):
        		stats.append([])
#for every job in l_folders
	for ijob in range(njobs):
    #for every stafiles in l_sta_files list:
    		for fsta,ista in zip(l_sta_files[ijob],range(len(l_sta_files[ijob]))):
			temp = worksta(fsta,l_sta_opt[ijob]) #save the stats on the right place
			stats[ijob].append(temp)
	return stats

def read_l_files_sp(l_spe_files,plane=-1,nvar=12):
	njobs = np.shape(l_spe_files)[0]
	spe = []
	for ijob in range(njobs):
    		for j in range(len(l_spe_files[ijob])):
        		spe.append([])
	#for every job in l_folders
	for ijob in range(njobs):
	 #for every spefiles in l_spe_files list:
    		for fspe,ispe in zip(l_spe_files[ijob],range(len(l_spe_files[ijob]))):
			temp = workspe(fspe,plane,nvar) #save the stats on the right place
			spe[ijob].append(temp)
	return spe

def read_l_files_sp_all(l_spe_files,nvar=12):
	spe = []
    	for fspe,ispe in zip(l_spe_files,range(len(l_spe_files))):
		temp = workspeall(fspe,nvar) #save the stats on the right place
		spe.append(temp)
	return spe




def load_sta_files(p_folders,l_sta_opt):
	#path to results folder (p_results), p_folders=subfolders to study,l_sta_opt=sta file version
	#0 == incompressible, 1== lomacte, without mum, 2==loma with compressible parts of dissipation
	njobs = len(p_folders)
	[l_sta_files,l_spe_files,nfiles,nfiles2]=create_l_files(p_folders)	
	stats = read_l_files(l_sta_files,l_sta_opt)
	return stats

def load_spe_files(p_folders,plane=-1,nvar=12):
	#path to results folder (p_results), p_folders=subfolders to study,l_sta_opt=sta file version
	#0 == incompressible, 1== lomacte, without mum, 2==loma with compressible parts of dissipation
	njobs = len(p_folders)
	[l_sta_files,l_spe_files,nfiles,nfiles2]=create_l_files(p_folders)	
	spe = read_l_files_sp(l_spe_files,plane,nvar)
	return spe

def load_spe_files_all(p_folder,ifirst,ilast,nvar=12):
	#path to results folder (p_results), p_folders=subfolders to study,l_sta_opt=sta file version
	#0 == incompressible, 1== lomacte, without mum, 2==loma with compressible parts of dissipation
	njobs = len(p_folder)
	l_spe_files = sorted(glob.glob(p_folder+'/*.spe'))
	spe = read_l_files_sp_all(l_spe_files[ifirst:ilast],nvar)
	return spe





def plot1Dts(x,y,styles,xlabel,ylabel,xlen,ylen,xlim=[0.0,0.0],ylim=[0.0,0.0],saveopt=[False,'/home','kk']):
#Preparing plot
	PRINT = saveopt[0]
	save_folder = saveopt[1]
	namefig = saveopt[2]
	fig = plt.figure(figsize=(xlen,ylen)) #Define figure with size
	ax=fig.add_subplot(1,1,1)
#Plot all jobs in the same graph
	for ijob in range(len(styles)):
    		plt.plot(x[ijob],y[ijob],styles[ijob],lw=2.0)
	ax.set_xlabel(xlabel)
	ax.set_ylabel(ylabel)
	if not xlim[0]==xlim[1]:
		ax.set_xlim(xlim)
	if not ylim[0]==ylim[1]:
		ax.set_ylim(ylim)
	plt.grid('on')
	if PRINT == 'True':
    		pathname = save_folder + namefig
    		fig.savefig(pathname,bbox_inches='tight')
		print "Figure has been saved"
def plot1Dts2(x,y,y2,styles,xlabel,ylabel,xlen,ylen,xlim=[0.0,0.0],ylim=[0.0,0.0],saveopt=[False,'/home','kk']):
#Preparing plot
	PRINT = saveopt[0]
	save_folder = saveopt[1]
	namefig = saveopt[2]
	fig = plt.figure(figsize=(xlen,ylen)) #Define figure with size
	ax=fig.add_subplot(1,1,1)
#Plot all jobs in the same graph
	for ijob in range(len(styles)):
    		plt.plot(x[ijob],y[ijob],styles[ijob],lw=2.0)
    		plt.plot(x[ijob],y2[ijob],styles[ijob],lw=2.0,ls='dashed')
		#ax.fill_between(x[ijob], y[ijob], y2[ijob], where=y>=y2, facecolor='green', interpolate=True)
	ax.set_xlabel(xlabel)
	ax.set_ylabel(ylabel)
	if not xlim[0]==xlim[1]:
		ax.set_xlim(xlim)
	if not ylim[0]==ylim[1]:
		ax.set_ylim(ylim)
	plt.grid('on')
	if PRINT == 'True':
    		pathname = save_folder + namefig
    		fig.savefig(pathname,bbox_inches='tight')
		print "Figure has been saved"





def pcolorts(time,stats,variable,frame_max,xlabel,ylabel,xlen,ylen,xlim=[0.0,0.0],ylim=[0.0,0.0],saveopt=[False,'/home','kk']):
#Preparing plot
        PRINT = saveopt[0]
        save_folder = saveopt[1]
        namefig = saveopt[2]
        fig = plt.figure(figsize=(xlen,ylen)) #Define figure with size
#Plot all jobs in the same graph
        for ijob in range(len(time)):
                ax=fig.add_subplot(len(time),1,ijob+1)
                y = stats[ijob][0]['y'][:]
                kk = np.zeros([frame_max[ijob],len(y)])
                for j in range(0,frame_max[ijob]):
                        kk[j,:] = stats[ijob][j][variable][:]
                plt.pcolor(time[ijob][0:frame_max[ijob]],y,kk.transpose(),edgecolor='None',cmap='Accent')
                #plt.contour(time[ijob][0:frame_max[ijob]],y,kk.transpose(),20,lw=4.0,color='b-')
                ax.set_xlabel(xlabel)
                ax.set_ylabel(ylabel)
                if not xlim[0]==xlim[1]:
                        ax.set_xlim(xlim)
                if not ylim[0]==ylim[1]:
                        ax.set_ylim(ylim)
                plt.grid('on')
        if PRINT == 'True':
                pathname = save_folder + namefig
                fig.savefig(pathname,bbox_inches='tight')
                print "Figure has been saved"


def contoursts(time,stats,variable,frame_max,xlabel,ylabel,xlen,ylen,xlim=[0.0,0.0],ylim=[0.0,0.0],ncontours=10,saveopt=[False,'/home','kk']):
#Preparing plot
    PRINT = saveopt[0]
    save_folder = saveopt[1]
    namefig = saveopt[2]
    fig = plt.figure(figsize=(xlen,ylen)) #Define figure with size
    ax=fig.add_subplot(1,1,1)
    collist = ['b','g','r']
    CS =[1 for x in range(len(time))]
#Plot all jobs in the same graph
    for ijob in range(len(time)):
        y = stats[ijob][0]['y'][:]
        kk = np.zeros([frame_max[ijob],len(y)])
        for j in range(0,frame_max[ijob]):
            kk[j,:] = stats[ijob][j][variable][:]
	print collist[ijob]
        CS[ijob] =plt.contour(time[ijob][0:frame_max[ijob]],y,kk.transpose(),ncontours,colors=collist[ijob])
    cbar = plt.colorbar(CS[0])
    for ijob in range(1,len(time)):
        cbar.add_lines(CS[ijob])
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if not xlim[0]==xlim[1]:
        ax.set_xlim(xlim)
    if not ylim[0]==ylim[1]:
        ax.set_ylim(ylim)
    plt.grid('on')
    if PRINT == 'True':
        pathname = save_folder + namefig
        fig.savefig(pathname,bbox_inches='tight')
        print "Figure has been saved"



#import lomapost
def averageprofile(stats,variable,scaling,xilim,njobs,ifirst,ilast):
    """
        Average stats from ifirst to ilast and return one only curve
    """
    Ny = 1000
    xi = np.linspace(-xilim,xilim,Ny)
    fun_ave = [[] for x in range(njobs)]
    faverage = np.zeros([Ny])
    for ijob in range(njobs):
        y = stats[ijob][0]['y'] #need vertical axis
    #new sta files list from minimum time:
        faverage = np.zeros([Ny])
        for ista in range(ifirst[ijob],ilast[ijob]+1):
            fun   = stats[ijob][ista][variable]
            delta = stats[ijob][ista][scaling]
            #Calculate scaled function
            yxi = y/delta
            scaled = np.interp(xi,yxi,fun)
            #accumulate data
            faverage = faverage + scaled
        fun_ave[ijob]=faverage/(ilast[ijob]-ifirst[ijob]+1)
        faverage = np.zeros([Ny])
    return xi,fun_ave

def averagestats(yxi,fun_ave,ncases,jobs):
    """
        Define with jobs the way you want to average stats
        jobs=[[0,1,2],[3,4,5],[6,7,8]], jobs[icase][ijob]
    """
    fun = [[] for x in range(ncases)]
    fsum = np.zeros(len(yxi))
    for icase in range(ncases):
        for ijob in jobs[icase]:
            fsum = fsum + fun_ave[ijob]
        fun[icase]=fsum/(len(jobs[icase]))
        fsum=np.zeros(len(yxi))
    return yxi,fun


def plotprofiles(xi,fun,styles,xlabel,ylabel,xlen,ylen,xlim=[0.0,0.0],ylim=[0.0,0.0],saveopt=[False,'/home','kk'],refx=[],refy=[],refcolor=[]):
    """
        Make plots of all functions listed in "fun" with same horizontal axis "xi"
    """
#Preparing plot
    PRINT = saveopt[0]
    save_folder = saveopt[1]
    namefig = saveopt[2]
    fig = plt.figure(figsize=(xlen,ylen)) #Define figure with size
    ax=fig.add_subplot(1,1,1)
#Plot all jobs in the same graph
    for f,style in zip(fun,styles):
        plt.plot(xi,f,style,lw=2.0)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if not xlim[0]==xlim[1]:
        ax.set_xlim(xlim)
    if not ylim[0]==ylim[1]:
        ax.set_ylim(ylim)
    if not refx ==[]:
	i=0
	for x,y in zip(refx,refy):
		plt.plot(x,y,'o',color=refcolor[i],ms = 7,mfc=refcolor[i],alpha=0.8)
		i+=1
    plt.grid('on')
    if PRINT == True:
        pathname = save_folder + namefig
        fig.savefig(pathname,bbox_inches='tight')
        print "Figure has been saved in %s" %(pathname)


def plotprofiles2(xi,fun,styles,xlabel,ylabel,xlen,ylen,xlim=[0.0,0.0],ylim=[0.0,0.0],saveopt=[False,'/home','kk']):
    """
        Make plots of all functions listed in "fun" with horizontal axis listed in "xi"
    """
#Preparing plot
    PRINT = saveopt[0]
    save_folder = saveopt[1]
    namefig = saveopt[2]
    fig = plt.figure(figsize=(xlen,ylen)) #Define figure with size
    ax=fig.add_subplot(1,1,1)
#Plot all jobs in the same graph
    i = 0
    for x,f in zip(xi,fun):
        plt.plot(x,f,styles[i],lw=2.0)
        i += 1
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if not xlim[0]==xlim[1]:
        ax.set_xlim(xlim)
    if not ylim[0]==ylim[1]:
        ax.set_ylim(ylim)
    plt.grid('on')
    if PRINT == True:
        pathname = save_folder + namefig
        fig.savefig(pathname,bbox_inches='tight')
        print "Figure has been saved"

def getSScase3(p_folders,lcases,variable,scaling,SS1,SS2,Npoints=500,derivative=False,ferror='std',fdmean= False):
    """
    returns one for case wanted to compute, using its own scaling for each run!
    Distribution of runs is given in lcases[icase](ex. lcases = [0,1,2] o [3,5]]
    only for one case at a time
    return only one profile per case and its Error of the mean.
    http://teacher.nsrl.rochester.edu/phy_labs/AppendixB/AppendixB.html
    """
    ncases = len(lcases)
    fycases = [[] for i in range(ncases)]
    ycases = [[] for i in range(ncases)]
    ferrorcases = [[] for i in range(ncases)]
    funint = np.zeros(Npoints);
    dfunint = np.zeros(Npoints);
    lim = 15.0;
    if scaling == 'dw':
        lim=2.0
    ydelta = np.linspace(-lim,lim,Npoints)
    nacum = 0
    fcum  = np.zeros(Npoints)
    fcum2 = np.zeros(Npoints)
    dfcum = np.zeros(Npoints)
    icount = 1
    damax = 0.0
    #Obtain dispersion on scaling
    Ntimes= 50
    [tint,deltaint,deltastd] = gettimeaverage(p_folders,lcases,scaling,SS1,SS2,Ntimes)
    #---------------------------------------
    for ijob in lcases:
        path = p_folders[ijob]
        print path
        with h5py.File(path,"r") as f:
            time = np.array(f['time']);
            if variable =='vmd':
                fun = np.array(f['vm'])
                delta = np.array(f[scaling])
                fun = fun*delta
            else:
                fun = np.array(f[variable]);
            y = np.array(f['y'])
            if variable=='rum':
                rhom = np.array(f['rhom'])
                fun = np.divide(fun,rhom)
            elif variable=='R12':
                ruv = np.array(f['ruv'])
                rum = np.array(f['rum'])
                rvm = np.array(f['rvm'])
                rhom = np.array(f['rhom'])
                fun = np.divide(np.subtract(ruv,np.divide(np.multiply(rum,rvm),rhom)),rhom)
            elif variable=='R11':
                ruu = np.array(f['ruu'])
                rum = np.array(f['rum'])
                #rvm = np.array(f['rvm'])
                rhom = np.array(f['rhom'])
                fun = np.divide(np.subtract(ruu,np.divide(np.multiply(rum,rum),rhom)),rhom)
            elif variable=='R22':
                rvv = np.array(f['rvv'])
                rvm = np.array(f['rvm'])
                #rvm = np.array(f['rvm'])
                rhom = np.array(f['rhom'])
                fun = np.divide(np.subtract(rvv,np.divide(np.multiply(rvm,rvm),rhom)),rhom)
            elif variable=='R33':
                rww = np.array(f['rww'])
                rwm = np.array(f['rwm'])
                #rvm = np.array(f['rvm'])
                rhom = np.array(f['rhom'])
                fun = np.divide(np.subtract(rww,np.divide(np.multiply(rwm,rwm),rhom)),rhom)

                        
            #print "computed rum/rhom"
            delta = np.array(f[scaling])
            
            it1 = where(time>=SS1[ijob])[0][0]
            it2 = where(time>=SS2[ijob])[0][0]
            #print "SS1[%s]  = %s" %(ijob,SS1[ijob])
            #print "SS2[%s] = %s" %(ijob,SS2[ijob])
            for itime in range(it1,it2):
                dd = delta[itime]
                if fdmean == True:
                    dd = np.interp(time[itime],tint,deltaint)
                #edelta = 0.01*dd
                edelta = np.interp(time[itime],tint,deltastd)
                if derivative == True:
                    #ed = der1(y,fun[:,itime])*y/dd**2
                    funint = np.interp(ydelta,y/dd,dd*der1(y,fun[:,itime]))
                    dfunint = 0.0*funint#funint**2/dd**2*edelta**2
                    if itime==it1:
                        print "Taking derivative and multiplying by scaling"
                else:
                    if variable=='w3m' or variable[0:2] =='ep':
                        #ed = fun[:,itime]*y/dd**2
                        funint = np.interp(ydelta,y/dd,dd*fun[:,itime])
                        #edint   = np.interp(ydelta,y/dd,ed)
                        #dfunint = edint**2*edelta**2
                        dfunint = 0.0*funint
                        if itime==it1: #First time
                            print "Multiplying by scaling"

                    else:
                        funint = np.interp(ydelta,y/dd,fun[:,itime])
                        dfunint = 0.0*funint
                nacum+= 1
                fcum  = fcum  + funint
                fcum2 = fcum2 + funint**2
                dfcum = dfcum + dfunint
                #damax = max(damax,np.max(np.abs(funint)))
    funmean = fcum/nacum
    e1square = np.abs(fcum2-nacum*funmean**2)/(nacum-1)
    #e1square = np.abs(fcum2/nacum-funmean**2)
    #e2square = np.mean(dfcum)
    e2square = dfcum/nacum
    #funstd  = np.sqrt((e1square+e2square)/nacum)
    funstd = np.sqrt((e1square+e2square))
    #calculate ERROR of the mean= std/sqrt(N)
    print "Max of standard deviation = %s o/o" %(np.nanmax(np.abs(funstd))/np.nanmax(np.abs(funmean))*100)
    print "Number of profiles = %s" %(nacum)
    if ferror == 'std':
        print "Giving std-deviation error"
        return ydelta,funmean,funstd
    #funstd = np.sqrt(e2square+)
    else: 
        print "Giving error of the mean"
        return ydelta,funmean,funstd/np.sqrt(nacum)


def gettimeaverage(p_folders,lcases,variable,SS1,SS2,Npoints=500,ferror='std'):
    """
    returns one for case wanted to compute.
    Distribution of runs is given in lcases[icase](ex. lcases = [[0,1,2],[3,5]])
    only for one case at a time
    return only one profile per case and its std
    """
    nruns = len(lcases)
    #Obtain dispersion on scaling
    Ntimes= Npoints
    funcum  = np.zeros(Ntimes) 
    funcum2 = np.zeros(Ntimes)
    tmin=10000;
    tmax=0;
    ltimes =[]
    lfuns =[]
    for ijob in lcases:
        path = p_folders[ijob]
        with h5py.File(path,"r") as f:
            time = np.array(f['time']);
            if variable =='dyeta':
                fun  = np.nanmax(np.array(f[variable]),axis=0)
                Lx=2*np.pi/np.float32(f['alp'])
                print "Lx=%s"%Lx
            elif variable=='eta':
                fun  = np.nanmin(np.array(f[variable]),axis=0)
            else:
                fun  = np.array(f[variable]);
            it1 = where(time>=SS1[ijob])[0][0]
            it2 = where(time>=SS2[ijob])[0][0]
            tmin = min(tmin,time[it1])
            tmax = max(tmax,time[it2])
            ltimes.append(time)
            lfuns.append(fun)
    tint = np.linspace(tmin,tmax,Ntimes)
    for ij in range(nruns):
        funint   = np.interp(tint,ltimes[ij],lfuns[ij])
        funcum  += funint
        funcum2 += funint**2
    funmean = funcum/nruns
    funstd = np.sqrt(np.abs(funcum2-nruns*funmean**2)/(nruns-1))
    if ferror=='std':
        return tint,funmean,funstd
    else:
    #funstd = np.sqrt(e2square+)
        return tint,funmean,funstd/np.sqrt(nruns)


def getdeltavstime(p_folders,lcases,variable,SS1,SS2,Npoints=500,ferror='std',scaling='dw'):
    """
    returns one for case wanted to compute.
    Distribution of runs is given in lcases[icase](ex. lcases = [[0,1,2],[3,5]])
    only for one case at a time
    return only one profile per case and its std
    """
    nruns = len(lcases)
    #Obtain dispersion on scaling
    Ntimes= Npoints
    funcum  = np.zeros(Ntimes) 
    funcum2 = np.zeros(Ntimes)
    tmin=10000;
    tmax=0;
    ltimes =[]
    lfuns =[]
    for ijob in lcases:
        path = p_folders[ijob]
        with h5py.File(path,"r") as f:
            time = np.array(f['time']);
            fun  = np.array(f[variable]);
            it1 = where(time>=SS1[ijob])[0][0]
            it2 = where(time>=SS2[ijob])[0][0]
            tmin = min(tmin,time[it1])
            tmax = max(tmax,time[it2])
            ltimes.append(time)
            lfuns.append(fun)
    tint = np.linspace(tmin,tmax,Ntimes)
    for ij in range(nruns):
        funint   = np.interp(tint,ltimes[ij],lfuns[ij])
        funcum  += funint
        funcum2 += funint**2
    funmean = funcum/nruns
    funstd = np.sqrt(np.abs(funcum2-nruns*funmean**2)/(nruns-1))
    #---------------------------------------
    if ferror=='std':
        return tint,funmean,funstd
    else:
    #funstd = np.sqrt(e2square+)
        return tint,funmean,funstd/np.sqrt(nruns)


def createmesh3(path,Ly,npoints,lim1,k,mright,fac2,Dy0,Dyf,write=False):
    """ 
    Create a mesh in 3 pices.
    Path : where to save the mesh
    Ly : size of the mesh
    npoints: number of points of half of the domain.
    lim1:length of the first piece (piece-1)
    k: Jump of Dy from mid to end of piece-1 (kref=2.5)
    mright: value of the second derivative of y (between 0.001 and 0.03)
    fac2: lim2=fac2*Ly (fac2=0.9)
    Dy0: mesh size in the middle (0.025dw0 is a good number)
    Dyf: mesh size in the end (same order of magnitud than middle)
    """
    #Define the linear system according to the asked requirements of the mesh (by hand)
    A = np.array([[lim1**3, lim1**4],[3*lim1**2,4*lim1**3]])
    b = np.array([(k-1)*Dy0,mright])
    x = np.linalg.solve(A, b)
    L = Ly/2;
    y=np.linspace(0,L,npoints)
    Dy1 = Dy0 + x[0]*y**3+x[1]*y**4
    D = x[0]
    E = x[1]
    Dylim =  Dy0 + x[0]*lim1**3+x[1]*lim1**4
    ilim = np.where(y<=lim1)[0][-1]
    #Second piece
    #--------------
    lim2=fac2*Ly/2;
    #Dy(y) = A + B*y
    #Dy(lim=Dy1[-1])
    #Dy'(lim) = mright
    #1)
    #Dy (lim) = A + B*lim
    #dotDy(lim)=B
    B = mright;
    A = Dylim-mright*lim1
    #y2=np.linspace(0,Ly/2,npoints)
#Dy2=0
    Dy2 = A+B*y
    Dy2lim2=A+B*lim2
    ilim2 = np.where(y<=lim2)[0][-1]
    #Third limit
    #--------------
    #Dy3 = A + B y + C y**2 + Dy**3
    #Dy2(lim2)= A + B*lim2 + C*lim2**2 + D*lim2**3
    #B = B' + 2*C*lim2+3*D*lim2**2
    #Dy3(L)= A + B*y(L) + C*y(L)**2 + D*y(L)**3
    #0 = B' + 2*C*y(L)+ 3*D*y(L)**2
    A = [[1,lim2,lim2**2,lim2**3],[0.0,1,2*lim2,3*lim2**2],[1,L,L**2,L**3],[0,1,2*L,3*L**2]]
    b = [Dy2lim2,B,Dyf,0]
    x = np.linalg.solve(A, b)
    #y3 = np.linspace(lim2,L,npoints3)
    #y3 = np.linspace(lim2,L,80)
    Dy3 = x[0]+x[1]*y+x[2]*y**2+x[3]*y**3
    #Put All together
    Dy = np.r_[Dy1[0:ilim],Dy2[ilim:ilim2],Dy3[ilim2:]]
    y2= np.array(cumtrapz(Dy))
    y2=y2/y2[-1]*(L)
    
    ynew = np.zeros(2*len(y2)+1)
    mynew = len(ynew)
    ynew[0:mynew/2] = -y2[::-1]
    ynew[mynew/2+1:mynew] = y2[0:len(y2)+1]
    #Write file
    if write==True:
        np.savetxt(path, ynew)
	print "mesh saved at %s"%(path)
    return ynew
