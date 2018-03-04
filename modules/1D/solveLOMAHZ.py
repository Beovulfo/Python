from scipy import *
from scipy.integrate import cumtrapz
from scipy.integrate import simps
from scipy.integrate import romb
from math  import *
from pylab import *
import numpy as np
import compact




def isNaN(num):
    return num != num

def myinteg(y,u):
	result = np.zeros(len(y))
	for j in range(0,len(y)-2):
		result[j+1]=result[j]+0.5*(y[j+1]-y[j])*(u[j+1]+u[j])
	#result = result - 0.5* result[-1]
	return result

def der1(y,u):
	"""computes first derivative using Central Finite Differences
        """
	result = array(u)
        for j in range(len(y)):
 		if j==0:
			result[j]=(u[j+1]-u[j])/(y[j+1]-y[j])
		elif j==(len(y)-1):
			result[j]=(u[j]-u[j-1])/(y[j]-y[j-1])
		else:
			result[j]=(u[j+1]-u[j-1])/(y[j+1]-y[j-1])
	
	return result

def der2(y,u):
	"""Computes second derivative using central differences
	"""
	result = zeros(len(y))
	for j in range(len(y)):
		if j==0:
			result[j] = (u[j+2]-2*u[j+1]+u[j])/((y[j+2]-y[j])/2.)**2
		elif j==(len(y)-1):
			result[j] = (u[j-2]-2*u[j-1]+u[j])/((y[j-2]-y[j])/2.)**2
		else:
			result[j]=(u[j+1]-2*u[j]+u[j-1])/((y[j+1]-y[j-1])/2.)**2
        return result

def integ(y,u):
        result = cumtrapz(u,y,initial=0.0)
	#return append(result,result[-1])
	return result

def geninit(y,rho0,s,option=None):
	""" Generate from y and s the initial arrays
		(ru00,rv00,t00) = geninit(y,rho0,s,option=None)
	"""
        my = len(y)
	u00 = zeros(my)
	T00 = zeros(my) 
	rhou00 = zeros(my)
	rhov00 = zeros(my)
        dm0 = 1.0
	u0 = 0.5
	rho00 = rho0*(1.0 + (s-1.0)/(s+1.0)*np.tanh(-y/(2*dm0)))
	u00 = u0*np.tanh(-y/(2*dm0))
	T00 = 1.0/rho00
        rhou00 = rho00*u00 
	return rhou00,rhov00,T00

def calcdw(y,um):
	DU =abs(um.max()-um.min())
	m = (abs(der1(y,um))).max()
	if m ==0:
		return 0.0
	else:
		return DU/m 


def dtvisc(y,re):
	""" Estiamtes Dtvisc and Dtconv for umax=1
	"""
	#Estimate Dt
	Dtvisc = re/2.*(diff(y)[0])**2
	Dtconv = diff(y)[0]/1.0
	return Dtvisc, Dtconv

def calcdmcomp(y,um,rum,rhom):
	""" Pantano definition
	"""
	from numpy import trapz
	DU = um.max()-um.min()
	if um.any() != rum.any(): #Variable density definition
		r0 = (rhom[0]+rhom[-1])/2.0
		favum = rum/rhom
		result = 1./(r0*DU**2.0)*trapz(rhom*(0.5*DU-favum)*(0.5*DU+favum),y)
	else: #constant density case
		r0 = 1.0
		favum = um
		rhom = um/um #ones
		result = 1./(r0*DU**2.0)*trapz(rhom*(0.5*DU-favum)*(0.5*DU+favum),y)
	return result



def calcdm(y,u00):
	integ1 = trapz(u00**2.0,y)/(u00[-1]-u00[0])**2.0
	return 0.25*(y[-1]-y[0])-integ1	

#def calcdm(y,u00):
#	integ1 = trapz(u00**2.0,y)/(u00[-1]-u00[0])**2.0
#	return 0.25*(y[-1]-y[0])-integ1	


def lomacte(y,CFL,nstep,rho0,s,RKopt, param=0.4,param0 =-0.2,verb = False,time0=None, ru00=[], rv00=[], t00=[]):
	""" This functions evolve rhou00,T00 and updates rhov00 using drho,
	    reproduces 00 mode evolution of loma
		(time,y,rhou00,rhov00,T00,timev,dmv) = evolve(my,CFL,nstep,deltaT,RKopt,time0,y,ru00,rv00,t00)
	"""
	if verb == True:
		print "Running Runge-Kutta evolution of LOMAHZ_1D..."
	#RKNEW
	xi    = [0   ,-17.0/60.0,  -5.0/12.0]
	gamma = [8./15.,    5./12.,  3./4.    ]
	beta  = [gamma[0]-param0,    param,  (param*(gamma[0]*gamma[2]+gamma[1]*gamma[2]+xi[1]*gamma[2]+gamma[0]*xi[2])-gamma[0]*gamma[1]*gamma[2]+gamma[0]*gamma[1]*xi[2])/(param-gamma[0]*gamma[1])]
	alpha = [param0,   gamma[1]+xi[1]-param, gamma[2]+xi[2]-beta[2]   ]

	#------------Physics -------------------------#
        Reynolds = 160.0;Pr=0.7
	ire = 1./Reynolds #inverse of Reynolds
	ipe = 1./(Reynolds*Pr) #inverse of Peclet
        my = len(y)
        L  = y[-1]
	#---------------------------------------------#
	#Initialization of variables (arrays)
	time=0.0
	dm=zeros(nstep)
	timev = zeros(nstep)
	T00last = zeros(my) 
	T00last1 = zeros(my) 
	drho = zeros(my) 
	drhodt = zeros(my) 
	drholast = zeros(my) 
	intdrho = zeros(my) 
	DT00 = zeros(my) 

	if time0 is None:
		(rhou00,rhov00,T00) = geninit(y,rho0,s)
	else: #restarting
		print "Restarting from previous run..."
		rhou00 = ru00
		rhov00 = rv00
		T00    = t00
		time   = time0

	[prem1,prem3,dt11,dt12,dt21,dt22,fmap]= compact.derivadas(y,my)
	rhotop = 1.0/T00[-1]
	rhobot = 1.0/T00[0]
	#
	#density ratio
	sratio = s
	if verb == True:
		print "density ratio = %s" % sratio
	#Start RHS
	rhs1=zeros([my,3]) #RHS rhou00
	rhs4=zeros([my,3]) #RHS T00
	rhomin = min(rhobot,rhotop)
        sigma = 0.0
        Tmax = 1.0/rhomin
        #Calculating viscous timestep
	Dmu =Tmax**(sigma+1.0)/(Reynolds*Pr)*(1.0/np.min(diff(y))**2)
	#Start istep 
	for istep in range(nstep):
		Dc  =np.max(rhov00*T00)*(1.0/np.min(diff(y)))
	#Obtain timestep
		Dt = CFL/(Dc+Dmu)
		if istep==0 and verb == True:
			#	Dt = 0.01*Dtvisc #first step is important
			print "Initial Dt = %s" % Dt
		for rkstep in range(3):
		#For each rkstep
		#Define the RHS's
			rhs1[:,rkstep] = -compact.deryr(rhou00*rhov00*T00,dt12,prem1,fmap,my) +ire *compact.deryyr(rhou00*T00,dt22,prem3,my)
			rhs4[:,rkstep] = -rhov00*T00*compact.deryr(T00,dt12,prem1,fmap,my) + ipe*T00*compact.deryyr(T00,dt22,prem3,my)  
			#Define Dt*RKparameters
			dtxi = Dt* xi[rkstep]
			dtgamma = Dt* gamma[rkstep]
			dtbeta = Dt* beta[rkstep]
                        #Calculations
			T00last = T00
			#rhou00=evolverkstep(rhou00,rhs1[:,rkstep],rhs1[:,rkstep-1],dtgamma,dtxi]
			rhou00  = rhou00[:] + dtgamma*rhs1[:,rkstep] + dtxi*rhs1[:,rkstep-1]
			DT00   = dtgamma*rhs4[:,rkstep] + dtxi*rhs4[:,rkstep-1]
			T00 = T00 + DT00	
                        #drhodt=-Drho/betaDt
                        drhodt = DT00/(beta[rkstep]*Dt*T00last*T00)
                        kk = compact.inty80(list(drhodt),fmap,dt11,dt12,prem1)
                        M = kk[-1] #flow outcome
                        #M = trapz(drhodt,y) #flow outcome
                        #drhodt = cumtrapz(drhodt,y,initial=0)-M/(1.0+sratio**(-0.5))
                        drhodt = kk-M/(1.0+sratio**(-0.5))
                        rhov00 = drhodt-alpha[rkstep]/beta[rkstep]*rhov00
		#istep loop
		time += Dt
		dm[istep]=calcdmcomp(y,rhou00*T00,rhou00,1.0/T00)
		timev[istep]=time
		if isNaN(dm[istep]):
			print 'NaN found in this simulation!'
			return time,rhou00,rhov00, T00,timev,dm
	if verb == True:
		print "last Dt = %s" % Dt
		print "final time = %s, dm = %s " % (timev[-1],dm[-1])
	return time, rhou00,rhov00, T00,timev,dm

		#=======================END OF evolve function ==================================#

def evolverkstep(var,rhs,rhslast,dtgamma,dtxi):
	return var + dtgamma*rhs+dtxi*rhslast

#Manual cumtrapz
def mycumtrapz(y,f):
	"""make cumtrapz same way used on Loma Code
	  returns cumintegral
	"""
	cumintegral = zeros(len(f))
	cumintegral[0] = 0.0
	for j in range(len(y)-1):
    		cumintegral[j+1]=cumintegral[j]+0.5*(y[j+1]-y[j])*(f[j+1]+f[j])
	return cumintegral
	

#TEST FUNCTIONS
def test():
	#TEST derivatives
	my = 1025
	y = linspace(-2*pi,2*pi,my)
	f = sin(y)
	df_anal = cos(y)
	d2f_anal = -sin(y)
	df = der1(y,f)
	d2f = der2(y,f)
	print "Testing derivatives with f = sin(y)"
	figure()
	plot(y,df_anal,'r',y,df,'b--')
	suptitle('First derivative comparition')
	xlabel('y')
	figure()
	plot(y,d2f_anal,'r',y,d2f,'b--')
	suptitle('Second derivative comparition')
	xlabel('y')
	e1 = sum(abs(df-df_anal))/my
	e2 = sum(abs(d2f-d2f_anal))/my
	assert  e1 < 1e-3, \
		"First derivative not good enough. Error = %20.15e" % e1
	assert  e2 < 1e-3, \
		"Second derivative not good enough.Error = %20.15e" % e2

	#TEST integ1
	intf_anal = -cos(y)
	intf = integ(y,f) - 1.0 #-1 is the primitive at j=0	
	intfcumtrapz = cumtrapz(f,y,initial = 0)-1
	e3 = sum(abs(intf-intf_anal))/my
	figure()
	plot(y,intf_anal,'r',y,intf,'b--',y,intfcumtrapz,'r--')
	suptitle('Integration')
	xlabel('y')
	assert e3 < 1e-3, \
		"Itegration looks like SHIT. Error = %20.15e" %e3
	

def loma(y,CFL,nstep,rho0,s,RKopt, param=0.4,param0 =-0.2,verb = False,time0=None, ru00=[], rv00=[], t00=[]):
	""" This functions evolve rhou00,T00 and updates rhov00 using drho,
	    reproduces 00 mode evolution of loma
		(time,y,rhou00,rhov00,T00,timev,dmv) = evolve(my,CFL,nstep,deltaT,RKopt,time0,y,ru00,rv00,t00)
	"""
	if verb == True:
		print "Running Runge-Kutta evolution of LOMAHZ_1D..."
	#RKNEW
	xi    = [0   ,-17.0/60.0,  -5.0/12.0]
	gamma = [8./15.,    5./12.,  3./4.    ]
	beta  = [gamma[0]-param0,    param,  (param*(gamma[0]*gamma[2]+gamma[1]*gamma[2]+xi[1]*gamma[2]+gamma[0]*xi[2])-gamma[0]*gamma[1]*gamma[2]+gamma[0]*gamma[1]*xi[2])/(param-gamma[0]*gamma[1])]
	alpha = [param0,   gamma[1]+xi[1]-param, gamma[2]+xi[2]-beta[2]   ]

	#------------Physics -------------------------#
        Reynolds = 160.0;Pr=0.7
	ire = 1./Reynolds #inverse of Reynolds
	ipe = 1./(Reynolds*Pr) #inverse of Peclet
        my = len(y)
        L  = y[-1]
	#---------------------------------------------#
	#Initialization of variables (arrays)
	time=0.0
	dm=zeros(nstep)
	timev = zeros(nstep)
	T00last = zeros(my) 
	T00last1 = zeros(my) 
	drho = zeros(my) 
	drhodt = zeros(my) 
	drholast = zeros(my) 
	intdrho = zeros(my) 
	DT00 = zeros(my) 

	if time0 is None:
		(rhou00,rhov00,T00) = geninit(y,rho0,s)
	else: #restarting
		print "Restarting from previous run..."
		rhou00 = ru00
		rhov00 = rv00
		T00    = t00
		time   = time0

	[prem1,prem3,dt11,dt12,dt21,dt22,fmap]= compact.derivadas(y,my)
	rhotop = 1.0/T00[-1]
	rhobot = 1.0/T00[0]
	#
	#density ratio
	sratio = s
	if verb == True:
		print "density ratio = %s" % sratio
	#Start RHS
	rhs1=zeros([my,3]) #RHS rhou00
	rhs4=zeros([my,3]) #RHS T00
	rhomin = min(rhobot,rhotop)
        sigma = 0.0
        Tmax = 1.0/rhomin
        #Calculating viscous timestep
	Dmu =Tmax**(sigma+1.0)/(Reynolds*Pr)*(1.0/np.min(diff(y))**2)
	#Start istep 
	for istep in range(nstep):
		Dc  =np.max(rhov00*T00)*(1.0/np.min(diff(y)))
	#Obtain timestep
		Dt = CFL/(Dc+Dmu)
		if istep==0 and verb == True:
			#	Dt = 0.01*Dtvisc #first step is important
			print "Initial Dt = %s" % Dt
		for rkstep in range(3):
		#For each rkstep
		#Define the RHS's
			rhs1[:,rkstep] = -compact.deryr(rhou00*rhov00*T00,dt12,prem1,fmap,my) +ire *compact.deryyr(rhou00*T00,dt22,prem3,my)
			rhs4[:,rkstep] = -rhov00*T00*compact.deryr(T00,dt12,prem1,fmap,my) + ipe*T00*compact.deryyr(T00,dt22,prem3,my)  
			#Define Dt*RKparameters
			dtxi = Dt* xi[rkstep]
			dtgamma = Dt* gamma[rkstep]
			dtbeta = Dt* beta[rkstep]
                        #Calculations
			T00last = T00
			rhou00  = rhou00[:] + dtgamma*rhs1[:,rkstep] + dtxi*rhs1[:,rkstep-1]
			DT00   = dtgamma*rhs4[:,rkstep] + dtxi*rhs4[:,rkstep-1]
			T00 = T00 + DT00	
                        #drhodt=-Drho/betaDt
                        #drhodt = DT00/(beta[rkstep]*Dt*T00last*T00)
                        drhodt =- DT00/(beta[rkstep]*Dt*T00last*T00)-\
	                          drhodt*alpha[rkstep]/beta[rkstep]
                        kk = compact.inty80(list(drhodt),fmap,dt11,dt12,prem1)
                        M = kk[-1] #flow outcome
                        #M = trapz(drhodt,y) #flow outcome
                        #drhodt = cumtrapz(drhodt,y,initial=0)-M/(1.0+sratio**(-0.5))
                        rhov00 = -kk+M/(1.0+sratio**(-0.5))
		#istep loop
		time += Dt
		dm[istep]=calcdmcomp(y,rhou00*T00,rhou00,1.0/T00)
		timev[istep]=time
		if isNaN(dm[istep]):
			print 'NaN found in this simulation!'
			return time,rhou00,rhov00, T00,timev,dm
	if verb == True:
		print "last Dt = %s" % Dt
		print "final time = %s, dm = %s " % (timev[-1],dm[-1])
	return time, rhou00,rhov00, T00,timev,dm

		#=======================END OF evolve function ==================================#



def lomahz(y,CFL,nstep,rho0,s,RKopt, param=0.4,param0 =-0.2,verb = False,time0=None, ru00=[], rv00=[], t00=[]):
	""" This functions evolve rhou00,T00 and updates rhov00 using drho,
	    reproduces 00 mode evolution of loma
		(time,y,rhou00,rhov00,T00,timev,dmv) = evolve(my,CFL,nstep,deltaT,RKopt,time0,y,ru00,rv00,t00)
	"""
	if verb == True:
		print "Running Runge-Kutta evolution of LOMAHZ_1D..."
	#RKNEW
	xi    = [0   ,-17.0/60.0,  -5.0/12.0]
	gamma = [8./15.,    5./12.,  3./4.    ]
	beta  = [gamma[0]-param0,    param,  (param*(gamma[0]*gamma[2]+gamma[1]*gamma[2]+xi[1]*gamma[2]+gamma[0]*xi[2])-gamma[0]*gamma[1]*gamma[2]+gamma[0]*gamma[1]*xi[2])/(param-gamma[0]*gamma[1])]
	alpha = [param0,   gamma[1]+xi[1]-param, gamma[2]+xi[2]-beta[2]   ]

	#------------Physics -------------------------#
        Reynolds = 160.0;Pr=0.7
	ire = 1./Reynolds #inverse of Reynolds
	ipe = 1./(Reynolds*Pr) #inverse of Peclet
        my = len(y)
        L  = y[-1]
	#---------------------------------------------#
	#Initialization of variables (arrays)
	time=0.0
	dm=zeros(nstep)
	timev = zeros(nstep)
	T00last = zeros(my) 
	T00last1 = zeros(my) 
	drho = zeros(my) 
	drhodt = zeros(my) 
	drholast = zeros(my) 
	intdrho = zeros(my) 
	DT00 = zeros(my) 

	if time0 is None:
		(rhou00,rhov00,T00) = geninit(y,rho0,s)
	else: #restarting
		print "Restarting from previous run..."
		rhou00 = ru00
		rhov00 = rv00
		T00    = t00
		time   = time0

	[prem1,prem3,dt11,dt12,dt21,dt22,fmap]= compact.derivadas(y,my)
	rhotop = 1.0/T00[-1]
	rhobot = 1.0/T00[0]
	#
	#density ratio
	sratio = s
	if verb == True:
		print "density ratio = %s" % sratio
	#Start RHS
	rhs1=zeros([my,3]) #RHS rhou00
	rhs4=zeros([my,3]) #RHS T00
	rhomin = min(rhobot,rhotop)
        sigma = 0.0
        Tmax = 1.0/rhomin
        #Calculating viscous timestep
	Dmu =Tmax**(sigma+1.0)/(Reynolds*Pr)*(1.0/np.min(diff(y))**2)
	#Start istep 
	for istep in range(nstep):
		Dc  =np.max(rhov00*T00)*(1.0/np.min(diff(y)))
	#Obtain timestep
		Dt = CFL/(Dc+Dmu)
		if istep==0 and verb == True:
			#	Dt = 0.01*Dtvisc #first step is important
			print "Initial Dt = %s" % Dt
		for rkstep in range(3):
		#For each rkstep
		#Define the RHS's
			rhs1[:,rkstep] = -compact.deryr(rhou00*rhov00*T00,dt12,prem1,fmap,my) +ire *compact.deryyr(rhou00*T00,dt22,prem3,my)
			rhs4[:,rkstep] = -rhov00*T00*compact.deryr(T00,dt12,prem1,fmap,my) + ipe*T00*compact.deryyr(T00,dt22,prem3,my)  
			#Define Dt*RKparameters
			dtxi = Dt* xi[rkstep]
			dtgamma = Dt* gamma[rkstep]
			dtbeta = Dt* beta[rkstep]
                        #Calculations
			T00last = T00
			#rhou00=evolverkstep(rhou00,rhs1[:,rkstep],rhs1[:,rkstep-1],dtgamma,dtxi]
			rhou00  = rhou00[:] + dtgamma*rhs1[:,rkstep] + dtxi*rhs1[:,rkstep-1]
			DT00   = dtgamma*rhs4[:,rkstep] + dtxi*rhs4[:,rkstep-1]
			T00 = T00 + DT00	
                        #drhodt=-Drho/betaDt
                        drhodt = DT00/(beta[rkstep]*Dt*T00last*T00)
                        kk = compact.inty80(list(drhodt),fmap,dt11,dt12,prem1)
                        M = kk[-1] #flow outcome
                        #M = trapz(drhodt,y) #flow outcome
                        #drhodt = cumtrapz(drhodt,y,initial=0)-M/(1.0+sratio**(-0.5))
                        drhodt = kk-M/(1.0+sratio**(-0.5))
                        rhov00 = drhodt-alpha[rkstep]/beta[rkstep]*rhov00
		#istep loop
		time += Dt
		dm[istep]=calcdmcomp(y,rhou00*T00,rhou00,1.0/T00)
		timev[istep]=time
		if isNaN(dm[istep]):
			print 'NaN found in this simulation!'
			return time,rhou00,rhov00, T00,timev,dm
	if verb == True:
		print "last Dt = %s" % Dt
		print "final time = %s, dm = %s " % (timev[-1],dm[-1])
	return time, rhou00,rhov00, T00,timev,dm

		#=======================END OF evolve function ==================================#


