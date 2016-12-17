#Change notes

#2
    #Improved point generation: 3d rather than 2d slice
    #Input parameters rather than calculated fields
    #Efficiency improvements
    #Hardcoded real-numbers

#2-1
    #Hardcoded 8 initial poisn
    #Refined simpleplot interface

#2-2
    #working on an array plot that solves for a range of polarization angles and magnitudes

import numpy as np
import matplotlib.pyplot as plt
import time
import sys
from mpl_toolkits.mplot3d import Axes3D
import easygui as gui
import scipy as sp

#Reads functions (f), startpoints(r0), restriction(D), effective field (heff), plate magnetization(Mp), ratio array, time step(h) and number of points(n)
def rk4(f, r0, D, hext, mp, ratio, h, n, path):
    t = (np.arange(n+1)*h).tolist() #creates a list of n+1 times
    r = [[0. for dim in range(3)] for num in range(n+1)] # creates a 3d array for rotation at each point
    r[0][0], r[0][1], r[0][2] = r0 # sets initial position
    prog = 0.
    totalprog= 0.
    direction = ""
    if(path%2==0):
        direction = "forwards"
    else:
        direction = "backwards"
    for i in range(n):
        k1 = 0.
        k2 = 0.
        k3 = 0.
        k4 = 0.
        heff = heffective(r[i], hext, D)
        for j in range(len(f)):
            k1 = k1 + f[j](t[i], r[i], heff, ratio[j], mp)
            k2 = k2 + f[j](t[i]+(h/2.), r[i]+((h/2.)*k1), heff, ratio[j], mp) 
            k3 = k3 + f[j](t[i]+(h/2.), r[i]+((h/2.)*k2), heff, ratio[j], mp)
            k4 = k4 + f[j](t[i]+h, r[i] + (h*k3), heff, ratio[j], mp)
           
            r[i] = np.array(r[i]/normalize(r[i])) #normalizes final position
            r[i+1] = r[i] + (h/6.)*(k1 + 2*k2 + 2*k3 + k4) #sets next position
        prog = 100*(float(i)/float(n))
        totalprog = (100*(path/36.)) + (prog*(1./32.))

        sys.stdout.write("Calculating path {0:1.0f}/18 moving ".format((int(path)/2)+1) + direction + " in time:\" %{0:6.2f}".format(prog) + "\tTotal Completion: %{0:6.2f} \r".format(totalprog))
        sys.stdout.flush()

    r[n] = np.array(r[n]/normalize(r[n])) #normalizes final position 
    return r

def plot(n, h, f, D=10., lam = 2.211, ad = 0.05, astt = 0.05, hext = np.array([0.,0.,0.]), mp = np.array([0.,0.,0.]), startpoints=False, endpoints=False, savefig=False, name=False): #takes initial position INPUT INITIAL CONDITION AS FLOATS   
    #Clock Start
    start = time.clock()
    num=int(n)
    #Generates array if starting points
    r0=[(1,0,0),(0.9961,0,.0871),(.9848, 0, .1736),(.9659, 0, .2588),(.9397, 0, .3420),(.9063, 0, .4226),(.8660, 0, .5),(.8192, 0, .5736),(.7660, 0, .6428),(.7071,0,.7071),(.6428, 0, .7660),(.5736, 0, .8192),(.5, 0 ,.8660),(.4226, 0, .9063),(.3420, 0, .9397),(.2588, 0, .9659),(.1736, 0, .9848),(.0871, 0, .9961)]
    for i in range(len(r0)):
        r0[i]=np.array(r0[i])
    r0 = np.array(r0)
    #Establishing ratios of different forces, magnitude of magnitization
    ratio = [lam, -1.*lam*ad, -1*lam*astt]
    
    #Generates arrays for receiving data from rk4 function(r[individual particle][time][dimension])
    rf = [[[0. for k in range(3)] for j in range(n+1)] for i in range(len(r0))]
    rr = [[[0. for k in range(3)] for j in range(n+1)] for i in range(len(r0))]
    r  = [[[0. for k in range(3)] for j in range(2*(n+1))] for i in range(len(r0))]
    #Splits output into 3-d, cylindrical, and scaled cylindrical
    x = [[0] for points in range(len(r0))]
    y = [[0] for points in range(len(r0))] 
    z = [[0] for points in range(len(r0))]
    theta = [[0] for points in range(len(r0))]
    phi = [[0] for points in range(len(r0))]
    
    #Sets up completion meter
    perc = -1
    #Begins array or each generated point (i)
    for i in range(len(r0)):
        n = num
        #convert all values to floats
        for j in range(3):
            r0[i][j] = float(r0[i][j])

        #normalizes starting point by magnitization magnitude
        r0[i] = np.array(r0[i]/normalize(r0[i]))
       
        #Generates all points in time for particle [i] 
        rf[i] = rk4(f, r0[i], D, hext, mp, ratio, h, n, 2*i)
        rr[i] = rk4(f, r0[i], D, hext, mp, ratio, (-1*h), n, (2*i)+1)
        rr[i] = np.flipud(rr[i])   

        for j in range(0,n+1,1):
            r[i][j] = rr[i][j]
            r[i][j+(n+1)] = rf[i][j]
        #Sets position at time zero (should be one of the earlier generated points)
        x[i][0] = r[i][0][0]
        y[i][0] = r[i][0][1]
        z[i][0] = r[i][0][2]
        
        n = 2*(n+1)
        #Adds a point for each point in time
        for j in range(1,n,1): #plots all points
            x[i].append(r[i][j][0])
            y[i].append(r[i][j][1])
            z[i].append(r[i][j][2])
            
            #Calculates the angle in cylindrical coordinates
            if(r[i][j][0]<0)and(r[i][j][1]>0):
                theta[i].append(np.arctan(r[i][j][1]/r[i][j][0])+np.pi)
            elif(r[i][j][0]<0)and(r[i][j][1]<0):
                theta[i].append(np.arctan(r[i][j][1]/r[i][j][0])-np.pi)
            elif(r[i][j][0]==0)and(r[i][j][1]>0):
                    theta[i].append(np.pi*.5)
            elif(r[i][j][0]==0)and(r[i][j][1]<0):
                    theta[i].append(np.pi*-.5)    
            elif((r[i][j][0]==0)and(r[i][j][1]==0)):
                    theta[i].append(0)
            else:
                theta[i].append(np.arctan(r[i][j][1]/r[i][j][0]))
            
            #Rotates the entire system by pi halves
            theta[i][j] = theta[i][j]+np.pi/2.
            
            if theta[i][j]>np.pi:
                theta[i][j] = theta[i][j] - (2. * np.pi)
            
            #Appends the scaled 2d image
            phi[i].append(theta[i][j]*np.sqrt(1-((z[i][j])**2)))            
        
        #Sets up masked terms for 2d projections
        zm = z
        zm[i] = np.ma.array(z[i])    
        theta[i] = np.ma.array(theta[i])
        phi[i] = np.ma.array(phi[i])

        #MASKING
        for j in range(len(theta[i])-1):
            if ((theta[i][j+1]-theta[i][j])>3) or ((theta[i][j+1]-theta[i][j])<-3):
                zm[i][j], zm[i][j+1]=np.ma.masked, np.ma.masked
                theta[i][j],theta[i][j+1]=np.ma.masked, np.ma.masked
                phi[i][j],phi[i][j+1]=np.ma.masked, np.ma.masked

    print "\nPlotting..."

    #Setup figures
    fig = plt.figure(figsize=(10,10))
    ax3d = fig.add_subplot(221, projection = '3d')
    plt.title('3D Image')
    ax2d = fig.add_subplot(222)
    plt.title('2D Image')
    ax2dp = fig.add_subplot(223)
    plt.title('Scaled 2D Image')
    
    #Plot each point
    for i in range(len(r0)):    #plots all points, starting points, and terminii
        #Plot color based on endpoints
        if (theta[i][n-1]>((np.pi/2.)-0.1)) and (theta[i][n-1]<(np.pi/2.)+0.1) and (z[i][n-1]<.1) and (z[i][n-1]>-.1):        
            ax3d.plot(x[i],y[i],z[i], '-b', label='Spin Magnetization')       
            ax2d.plot(theta[i],zm[i], '-b', label='Spin Magnetization')
            ax2dp.plot(phi[i],zm[i], '-b', label='Spin Magnetization')
        elif (theta[i][n-1]>(-1*np.pi/2)-.1) and (theta[i][n-1]<(-1*np.pi/2)+.1) and (z[i][n-1]<.1) and (z[i][n-1]>-.1):
            ax3d.plot(x[i],y[i],z[i], '-r', label='Spin Magnetization')
            ax2d.plot(theta[i],zm[i], '-r', label='Spin Magnetization')
            ax2dp.plot(phi[i],zm[i], '-r', label='Spin Magnetization')
        
        #If no endpoint, plot black
        else:
            ax3d.plot(x[i],y[i],z[i], '-k', label='Spin Magnetization')
            ax2d.plot(theta[i],zm[i], '-k', label='Spin Magnetization')
            ax2dp.plot(phi[i],zm[i], '-k', label='Spin Magnetization')
       
       #Optional plot start/endpoints
        if startpoints==True:
            ax3d.plot([x[i][0]], [y[i][0]], [z[i][0]], 'og', label='start points')
            ax2d.plot([theta[i][0]], [zm[i][0]], 'og', label='start points')
            ax2dp.plot([phi[i][0]], [zm[i][0]], 'og', label='start points')
        if endpoints==True:
            ax3d.plot([x[i][n-1]], [y[i][n-1]], [z[i][n-1]], 'or', label='end points')
            ax2d.plot([theta[i][n-1]], [zm[i][n-1]], 'or', label='end points')
            ax2dp.plot([phi[i][n-1]], [zm[i][n-1]], 'or', label='end points')
        
    ax3d.set_xlabel("x-axis")
    ax3d.set_ylabel("y-axis")
    ax3d.set_zlabel("z-axis")
    ax3d.set_xlim3d(-1.2, 1.2)
    ax3d.set_ylim3d(-1.2, 1.2)
    ax3d.set_zlim3d(-1.2, 1.2)
    
    ax2d.set_xlabel("Phi-axis")
    ax2d.set_ylabel("Z-axis")
    ax2d.set_xlim(-np.pi,np.pi)
    ax2d.set_ylim(-1.2, 1.2)
    
    ax2dp.set_xlabel("Scaled Phi-axis")
    ax2dp.set_ylabel("Z-axis")
    ax2dp.set_xlim(-np.pi,np.pi)
    ax2dp.set_ylim(-1.2, 1.2)
    
    elapsed = time.clock()-start 

    Stats = "Elapsed time: " + str(int((elapsed-(elapsed%60))/60))+ " minutes and {0:6.2f}".format(elapsed%60)+" seconds\nPoints: " +str(len(r0)) + "\nNumber of Steps: " + str(num) + "\nStep Size: " + str(h) + "s\nD value: " + str(D) + "\nDamping ratio = " + str(ad) + "\nSpin Torque Current Ratio: " + str(astt) + "\nMagnitization Polarization: " + str(mp) + "\nExternal Field: " + str(hext)

   
    print Stats
    
    if (savefig==False):
        plt.show()
        while(1):
            save = str(raw_input("Save image? (y/n): "))
            if (save == "y"):
                if (name==False):
                    fig.savefig("C:/Users/Chris/Documents/Physics/Research/Spin Analysis/"+time.strftime("%d%m%y")+"at"+time.strftime("%H%M%S")+".png")
                    print("Saved image!")
                    Text = open("C:/Users/Chris/Documents/Physics/Research/Spin Analysis/"+time.strftime("%d%m%y")+"at"+time.strftime("%H%M%S")+".txt", "w")
                    Text.write(Stats)
                    Text.close()
                else:
                    fig.savefig("C:/Users/Chris/Documents/Physics/Research/Spin Analysis/"+name+".png")
                    print("Saved image!")
                    Text = open("C:/Users/Chris/Documents/Physics/Research/Spin Analysis/"+name+".txt", "w")
                    Text.write(Stats)
                    Text.close()
                break
            elif (save == "n"):
                print("Image not saved")
                break
            else:
                print "Not Valid input"
    elif (savefig==True):
        if (name==False):
            fig.savefig("C:/Users/Chris/Documents/Physics/Research/Spin Analysis/"+time.strftime("%d%m%y")+"at"+time.strftime("%H%M%S")+".png")
            print("Saved image!")
            Text = open("C:/Users/Chris/Documents/Physics/Research/Spin Analysis/"+time.strftime("%d%m%y")+"at"+time.strftime("%H%M%S")+".txt", "w")
            Text.write(Stats)
            Text.close()
        else:
            fig.savefig("C:/Users/Chris/Documents/Physics/Research/Spin Analysis/"+name+".png")
            print("Saved image!")
            Text = open("C:/Users/Chris/Documents/Physics/Research/Spin Analysis/"+name+".txt", "w")
            Text.write(Stats)
            Text.close()

    plt.cla()
    plt.close

def heffective(M, hext, D):
    k = 3.86e5
    m0 = 4e-7*np.pi
    #e3
    #((2*k)/(m0*Ms))
    heff = np.array([M[0]+hext[0], 0., -1.*D*M[2]+hext[2]])
    return heff

def normalize(r3): #returns normalization factor
    norm = np.sqrt(r3[0]**2+r3[1]**2+r3[2]**2)  
    if norm == 0:
        norm = 1
    return norm

def Procession(t, r, heff, rat, mp): # calculates the movement due to Procession energy
    rot = rat*np.cross(r,heff) #the rotation vector is the cross of the magnetization and the gradient of the energy field
    return rot

def Damping(t, r, heff, rat, mp): # calculates resitance due to damping
    rot = rat * np.cross(r, (np.cross(r,heff))) /normalize(r)
    return rot
    
def STT(t, r, heff, rat, mp):
    rot = rat * np.cross(r, (np.cross(r, mp))) /normalize(mp)
    return rot

def SimplePlot():
    msg = "Enter your desired values"
    title = "Spin Analysis Plot"
    FieldNames = ["Number of timesteps", "Timestep size (s)"]
    FieldValues = []
    FieldValues = gui.multenterbox(msg, title, FieldNames)
    while True:
        if FieldValues == None: break
        error = ""
        for i in range(len(FieldNames)):
            if FieldValues[i].strip()=="":
                error = error + (FieldNames[i]+" is a required value.\n")
        if error == "": break
        FieldValues = gui.multenterbox(msg, title, FieldNames)
    
    n = int(float(FieldValues[0]))
    h = float(FieldValues[1])
    f=[Procession, Damping, STT]

    msg = "Enter your desired values (all values optional)"
    title = "Spin Analysis Plot"
    FieldNames = ["D value (restrictive field, default 10)", "Damping ratio (default .05)", "Spin Torque current ratio", "Magnetic Polarization angle (0-90)", "External field (default 0,0,0, enter in form 'hx,hy,hz')"]
    FieldValues = []
    FieldValues = gui.multenterbox(msg, title, FieldNames)

    msg ="What key points do you want to view?"
    title = "Key points"
    choices = ["1. None", "2. Both", "3. Startpoints", "4. Endpoints"]
    choice = gui.multchoicebox(msg, title, choices)
    choice = int((choice.split("."))[0])

    spoints = False 
    epoints = False

    if (choice == 2):
        spoints = True
        epoints = True

    elif (choice == 3):
        spoints = True

    elif (choice == 4):
        epoints = True

    D = 10.
    da = 0.05
    sa = 0.05
    hext = '0,0,0'
    mp =  0.
    hextf = [0.,0.,0.]
    mpf = [0., 0., 0.]


    if (FieldValues[0]!=""):
        D = float(FieldValues[0])
    if (FieldValues[1]!=""):
        da = float(FieldValues[1])
    if (FieldValues[2]!=""):
        sa = float(FieldValues[2])
    if (FieldValues[3]!=""):
        mp = (FieldValues[3])
    if (FieldValues[4]!=""):
        hext = (FieldValues[4])
    hext = hext.split(',')
    mp = (float(mp)/360)*2*np.pi
    mpf = np.array([np.cos(mp), 0., np.sin(mp)])
    mpf = mpf/normalize(mpf)
    for i in range(3):
        hextf[i] = float(hext[i])

    plot(n, h, f, D=D, ad=da, astt=sa, hext = np.array(hextf), mp = mpf, startpoints = spoints, endpoints = epoints)

def ManyPlot(minangle, maxangle, numangle, minrat, maxrat, numrat, ist=1, jst=1):
    
    if (minangle<0 or maxangle>90 or minrat<0 or numangle<=0 or numrat<=0):
        print "Invalid values angle (Angle must be 0-90, Ratio must be >0, Numbers of each must be at least 1)"
        sys.exit
    minangle = minangle*(2*np.pi/360)
    maxangle = maxangle*(2*np.pi/360)
    angles = np.linspace(minangle, maxangle, numangle)
    ratios = np.linspace(minrat, maxrat, numrat)
    for i in range(int(ist-1), numangle,1):
        for j in range(int(jst-1), numrat,1):
            plot(int(5e3), .001, [Procession, Damping, STT], astt = ratios[j], mp = np.array([np.cos(angles[i]), 0., np.sin(angles[i])]), savefig=True, name = (str(i) + "-" + str(j)))
            print("Completed " + str(1+j+(i*numrat)) + " of " + str(numrat*numangle))
            print("")
    
