# -*- coding: utf-8 -*-
"""
Created on Wed Apr  2 15:25:07 2025

@author: ellen
"""

import numpy as np
import astropy.units as u
from astropy.constants import G

# import plotting modules
import matplotlib.pyplot as plt
import matplotlib

import OrbitCOM
from GalaxyMass import ComponentMass


'''
The goal of this code is to do two things:
    1) Determine the fraction of the total stellar mass that is contributed
    by the bulge and the disk of both the Milky Way and M31
    2) Determine the total angular momentum that is contributed by the bulge
    and the disk of both the Milky Way and M31
  where both of these will be determined as a function of the distance from
  the center of the merger. This will be done as follows:
    1) The time and snapshot of the merger will be measured using
    code from previous homeworks
    2) At that time and snapshot, the center of mass will be determined
    3) Using the center of mass codes, I will determine the mass of each
    component about this center of mass
    4) These masses will be plotted as a function of distance from this
    center of mass
    5) The velocity of each component about the center of mass will then
    be determined
    6) Finally, using the velocity and the mass, we will be able to determine
    the angular momentum of each component
'''

#first, we are determining the time at which they would merge
#previous work has implied it would be around 5 Gyr from now
#I used this in mind to create the files (from snapshot 445-800) for both
#the disk and the bulge

MW_disk = np.genfromtxt('Orbit_MWMerge.txt', dtype = None, names = True, skip_header = 0)
M31_disk = np.genfromtxt('Orbit_M31Merge.txt', dtype = None, names = True, skip_header = 0)
MW_bulge = np.genfromtxt('Orbit_MWMerge_Bulge.txt', dtype = None, names = True, skip_header = 0)
M31_bulge = np.genfromtxt('Orbit_M31Merge_Bulge.txt', dtype = None, names = True, skip_header = 0)

#To do this, I would use the OrbitCOM files to create a merge .txt file from
#the first merged snapshot until the end of the simulation (about 800)
#This will be the first plot that I have to show the start of the merge
#and these will be the files used later
#These files will be made twice since I can choose the particle type
#within the function

def vectors(v1, v2):
    #the vectors can be for position or velocity
    difference = v1 - v2
    #now need to sum the squares of all of these vectors to find the magnitudes
    magsum = 0 #initial to add on to
    for i in difference:
        magsum = np.sum(difference**2)
    mag = np.sqrt(magsum) #now sqrt so that we have the sqrt of all the sum of the 
    #squared differences
    return mag #this is what we can use now to find the differences between
    #the different galaxies
    
MWpos_disk = np.zeros([len(MW_disk['x']), 3]) #initial position array at x pos
for i in range(0, len(MW_disk['x'])):
    MWpos_disk[i] += [MW_disk['x'][i], MW_disk['y'][i], MW_disk['z'][i]]
    
MWpos_bulge = np.zeros([len(MW_bulge['x']), 3]) #initial position array at x pos
for i in range(0, len(MW_bulge['x'])):
    MWpos_bulge[i] += [MW_bulge['x'][i], MW_bulge['y'][i], MW_bulge['z'][i]]

MWvel_disk = np.zeros([len(MW_disk['vx']), 3]) #literally exact same idea but velocities
for i in range(0, len(MW_disk['vx'])):
    MWvel_disk[i] += [MW_disk['vx'][i], MW_disk['vy'][i], MW_disk['vz'][i]]
    
MWvel_bulge = np.zeros([len(MW_bulge['vx']), 3]) #literally exact same idea but velocities
for i in range(0, len(MW_bulge['vx'])):
    MWvel_bulge[i] += [MW_bulge['vx'][i], MW_bulge['vy'][i], MW_bulge['vz'][i]]

M31pos_disk = np.zeros([len(M31_disk['x']),3])
for i in range(0, len(M31_disk['x'])):
    M31pos_disk[i] += [M31_disk['x'][i], M31_disk['y'][i], M31_disk['z'][i]]
    
M31pos_bulge = np.zeros([len(M31_bulge['x']),3])
for i in range(0, len(M31_bulge['x'])):
    M31pos_bulge[i] += [M31_bulge['x'][i], M31_bulge['y'][i], M31_bulge['z'][i]]
    
M31vel_disk = np.zeros([len(M31_disk['vx']), 3]) 
for i in range(0, len(M31_disk['vx'])):
    M31vel_disk[i] += [M31_disk['vx'][i], M31_disk['vy'][i], M31_disk['vz'][i]]
    
M31vel_bulge = np.zeros([len(M31_bulge['vx']), 3]) 
for i in range(0, len(M31_bulge['vx'])):
    M31vel_bulge[i] += [M31_bulge['vx'][i], M31_bulge['vy'][i], M31_bulge['vz'][i]]
    
MWM31pos_disk = np.zeros(len(MW_disk['x'])) #empty array to store
for i in range(len(MW_disk['x'])):
    MWM31pos_disk[i] = vectors(MWpos_disk[i], M31pos_disk[i])
    
MWM31vel_disk = np.zeros(len(MW_disk['vx']))#empty array again
for i in range(len(MW_disk['vx'])):
    MWM31vel_disk[i] = vectors(MWvel_disk[i], M31vel_disk[i])
    
MWM31pos_bulge = np.zeros(len(MW_bulge['x'])) #empty array to store
for i in range(len(MW_bulge['x'])):
    MWM31pos_bulge[i] = vectors(MWpos_bulge[i], M31pos_bulge[i])
    
MWM31vel_bulge = np.zeros(len(MW_bulge['vx']))#empty array again
for i in range(len(MW_bulge['vx'])):
    MWM31vel_bulge[i] = vectors(MWvel_bulge[i], M31vel_bulge[i])
    
plt.plot(MW_disk['t'], MWM31pos_disk, color = 'red', lw=3)
plt.xlabel("time (Gyr)")
plt.ylabel("Position (kpc)")
plt.title("MW and M31 disk position")
plt.savefig("MWM31diskpositions.png")
plt.figure()

plt.plot(MW_disk['t'], MWM31vel_disk, color = 'blue', lw=3)
plt.xlabel("time (Gyr)")
plt.ylabel("Velocity (km/s)")
plt.title("MW and M31 disk velocity")
plt.savefig("MWM31diskvelocities.png")
plt.figure()

plt.plot(MW_bulge['t'], MWM31pos_bulge, color = 'goldenrod', lw=3)
plt.xlabel("time (Gyr)")
plt.ylabel("Position (kpc)")
plt.title("MW and M31 bulge position")
plt.savefig("MWM31bulgepositions.png")
plt.figure()

plt.plot(MW_bulge['t'], MWM31vel_bulge, color = 'green', lw=3)
plt.xlabel("time (Gyr)")
plt.ylabel("Position (km/s)")
plt.title("MW and M31 bulge velocity")
plt.savefig("MWM31bulgevelocities.png")

#Now that we will have the orbits of the disk and the bulge for both the
#MW and M31, we will plot each of these components as a function of their
#distance from the center of mass. 
#Since the files contain their x, y, and z data as a function of time,
#we can create a radial component that we analyze and compare it
#to the mass that would be determined through GalaxyMass & mass component (?)

#General outline for getting mass:
MW_bulge = ComponentMass('MW_000.txt', 3)
MW_disk = ComponentMass('MW_000.txt', 2)
M31_bulge = ComponentMass('M31_000.txt', 3)
M31_disk = ComponentMass('M31_000.txt', 2)
#not exactly sure if this is what I should be doing, which is what I plan
#to check in class tomorrow  

#We also have the velocity components about the center of mass from the
#OrbitCOM files that we received, and we can then further plot the velocity
#against the mass that we got earlier

#example:

#plt.plot(MWpos_disk, MW_disk)

#From here we can create some sort of iterative equation to account for the 
#angular momentum at all of the radial components and all of the velocities

#this would be my unique function
# distance = summationofxyz like in the homeworks where its the magnitude
#radial velocity is what we've accounted for (is my understanding)

def angularmomentum(distance, velocity):
    angular_momentum = velocity/distance
    return angular_momentum
#then from here I would iterate using a for loop for every point of the radius
#and this would then be measured against the function of time










