import glob
import os
from readfiles import *
import matplotlib.pyplot as plt
import pandas as pd
#import matplotlib as mpl
import numpy as np
from IPython.display import clear_output
from IPython.display import Image

def set_fig_props():
	plt.rcParams['figure.figsize'] = (10.0, 8.0)
	plt.rcParams['savefig.dpi'] = 300
	plt.rcParams['figure.facecolor'] = 'white'
	sizefont= 24
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
			temp = workstaHZ(fsta) #save the stats on the right place
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




def load_sta_files(p_folders,l_sta_opt=0):
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
