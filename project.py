import numpy as np
import matplotlib.pyplot as plt
import time
import sys
from mpl_toolkits.mplot3d import Aces3D
import easygui as gui
import scipy as sp

global k = 3.86e5
global D = 10.
global alpha = 0.105
global gam = 2.211
global Bext = np.array(0., 0., 0.,)
global beta = 0.05
global Mp = np.array([0., 0., 0.])
global x0 = np.array([1., 0., 0.])
global x1 = np.array([-1., 0., 0.])
global z0 = np.array([0., 0., 1.])
global z1 = np.array([1., 0., -1.])

#first order differential equation for magnetization
def LLG(M):
	global alpha
	global gamma
	part1 = np.cross(M, Heff(M))
	part2 = alpha * (np.cross(M, np.cross(M, Heff(M))))/Norm(M)
	result = -1.* gamma * (part1 + part2)
	return result

#returns normalization factor
def Norm(M):
	return np.sqrt((M[0]**2) + (M[1]**2) + (M[2]**2))

def Heff(M):
	global k
	global D
	global Bext
	result = np.array([M[0]+Bext[0], 0.+Bext[1], -1.*D*M[2])
	return result

def STT(M):
	global beta
	global Mp
	result = beta * np.cross(M, (np.cross(M, Mp)))
	return result

def Full(M):
	return LLG(M) + STT(M)

def Dist(M1, M2):
	M = M1 - M2
	result = np.sqrt((M[0]**2)+(M[1]**2)+(M[2]**2))

def RK4(f, M, h):
	k1 = h * f(M)
	k2 = h * f(M + (k1/2.))
	k3 = h * f(M + (k2/2.))
	k4 = h * f(M + k3)
	result = M + (k1+k4)/6. + (k2+k3)/3.
	return result

def RK4A(f, m0, v0, prec = 10e-10, lim = 10e5):
	N = int(lim)
	h = 1e-5
	t = np.empty(N)
	t[0] = 0.
	M = np.empty([N, 3])
	M[0] = np.array([m0[0], m0[1], m0[2])
	delta = prec
	i = 0
	done = False
	do while(done==False):
		global x0
		global x1
		global z0
		global z1

		M0 = RK4(f, M[i], h)
		M1 = RK4(f, M0, h)
		M2 = RK4(f, M[i], 2*h)
		dis = Dist(M1, M2)	
	
		if((dis-dis)==dis):
			rho = 10.
		else:
			rho = 30. * h * delta/ dis
		if(rho>=1):
			i++
			M[i] = M1
			h = h * np.power(rho, .25)
		elif(rho<1)
			h = h * np.power(rho, .25)
		
		if(i==(N-1)):
			endpoint = M[i]
			print("Not Enough Points!")
			done == True
	
		dx0 = Dist(M[i], x0)
		dx1 = Dist(M[i], x1)
		dz0 = Dist(M[i], z0)
		dz1 = Dist(M[i], z1)

		if(dx0<prec):
			M[i] = x0
			np.resize(M,[i+1,3])
			done = True
		elif(dx1<prec):
			M[i] = x1
			np.resize(M,[i+1,3])
			done = True
		elif(dz0<prec):
			M[i] = z0
			np.resize(M,[i+1,3])
			done = True	
		elif(dz1<prec):
			M[i] = z1
			np.resize(M,[i+1,3])
			done = True

		#naisu naisu over raisu include the limiting cycle
	return M
	#derp this M is a trajectory of magnetization

#input an array or trajectories
def plot(M):
	for k in range(M):
		Mag = M[k] #1 trajectory
		x = np.empty(len(Mag))
		y = np.empty(len(Mag))
		z = np.empty(len(Mag))
		phi = np.empty(len(Mag))
		theta = np.empty(len(Mag))
		for i in range(Mag):
			x[i] = Mag[i][0]
			y[i] = Mag[i][1]
			z[i] = Mag[i][2]
			
			if(x[i]>0):
				phi[i] = np.arctan(y[i]/x[i])
			elif(x[i]<0):
				phi[i] = np.arctan(y[i]/x[i]) + np.pi
			elif(x[i]==0):
				phi[i] = np.sign(y[i])*np.pi/2.
			r = np.sqrt((x[i]**2) + (y[i]**2))
			
			if(z[i]>0):
				theta[i] = np.arctan(r/z[i])
			elif(z[i]<0):
				theta[i] = (np.pi/2.) - np.arctan(r/z[i])
			elif((z[i]+z[i])==z[i]):
				theta[i] = np.pi/2.

		



















