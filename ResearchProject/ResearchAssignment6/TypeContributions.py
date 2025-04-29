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

import MassProfile


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

MW_disk = np.genfromtxt('Orbit_MWMerge.txt', dtype = None, names = True, skip_header = 0, comments = '#')
M31_disk = np.genfromtxt('Orbit_M31Merge.txt', dtype = None, names = True, skip_header = 0, comments = '#')
MW_bulge = np.genfromtxt('Orbit_MWMerge_Bulge.txt', dtype = None, names = True, skip_header = 0, comments = '#')
M31_bulge = np.genfromtxt('Orbit_M31Merge_Bulge.txt', dtype = None, names = True, skip_header = 0, comments = '#')

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
    
#the following functions gather the positions and velocities of the MW and M31's
#bulge and disk so that way we can find their difference and confirm
#that this is the correct snapshot and to understand the context we're 
#starting in
    
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
    
#now that we have the info, we can plot them to show that the two galaxies
#have merged for the bulge and disk positions and velocities
    
figure, axis = plt.subplots(2,2)
    
axis[0,0].set_xlabel("time from present (Gyr)")
axis[0,0].set_ylabel("Position (kpc)")
axis[0,0].set_title("MW and M31 disk position separation")
axis[0,0].plot(MW_disk['t'], MWM31pos_disk, color = 'red', 
               label = 'Milky Way and M31 disk separation')

axis[0,1].set_xlabel("time from present (Gyr)")
axis[0,1].set_ylabel("Velocity (km/s)")
axis[0,1].set_title("MW and M31 disk velocity")
axis[0,1].plot(MW_disk['t'], MWM31vel_disk, color = 'blue', 
               label = 'Milky Way and M31 velocity difference')

axis[1,0].set_xlabel("time from present (Gyr)")
axis[1,0].set_ylabel("Position (kpc)")
axis[1,0].set_title("MW and M31 bulge position")
axis[1,0].plot(MW_bulge['t'], MWM31pos_bulge, color = 'goldenrod',
               label = 'Milky Way and M31 bulge separation')

axis[1,1].set_xlabel("time from present (Gyr)")
axis[1,1].set_ylabel("Position (km/s)")
axis[1,1].set_title("MW and M31 bulge velocity")
axis[1,1].plot(MW_bulge['t'], MWM31vel_bulge, color = 'green', 
               label = 'Milky Way and M31 velocity difference')

plt.tight_layout()
plt.savefig('MWM31differences.png')

#These plots show that our snapshot start is correct and that this is
#where the merger begins. This is evidence for the snapshot we are using
#moving forward for the mass, radius, and velocity of the disk & bulge

#here's the radii that we want to look over to plug into the
#mass enclosed function, since this should encompass a majority of the merger
#remnant, at least for the parts we are interested in looking at. This was
#tweaked with until I felt like the mass profiles had flattened out

radius_array = np.arange(0, 150)

#now we need the masses using the mass profile code and the mass enclosed function
#we are using the low res snapshot 445 that was accessed from nimoy
#this snapshot contains that mass, position, and velocity information for all
#the particles in each galaxy at that point in time. We can use this to 
#get the mass profiles for the particles as a function of their distance from the
#center of mass

MW_mass = MassProfile.MassProfile('MW', 445)
MW_mass_disk = MW_mass.MassEnclosed(2, radius_array).value
MW_mass_bulge = MW_mass.MassEnclosed(3, radius_array).value

M31_mass = MassProfile.MassProfile('M31', 445)
M31_mass_disk = M31_mass.MassEnclosed(2, radius_array).value
M31_mass_bulge = M31_mass.MassEnclosed(3, radius_array).value

#now that I have the masses, I can plot them against the radii to get the
#mass profiles. I used fig, ax so that I can put them all in the same on for the
#appropriate masses so we can appropriately compare the contributions to the
#total mass

figure, axis = plt.subplots(1,2, figsize = (8, 6))

axis[0].set_xlabel('Distance from center of mass (kpc)')
axis[0].set_ylabel('Cumulative Mass (M$_{\odot}$)')
axis[0].set_title('Milky Way and M31 cumulative \n bulge mass from center of mass')
axis[0].plot(radius_array, MW_mass_bulge, lw = 2, color = 'indigo', label =
             'Milky Way bulge mass')
axis[0].plot(radius_array, M31_mass_bulge, lw = 2, color = 'forestgreen', 
             label = 'M31 bulge mass')
axis[0].legend()

axis[1].set_xlabel('Distance from center of mass(kpc)')
axis[1].set_ylabel('Culumative Mass (M$_{\odot}$)')
axis[1].set_title('Milky Way and M31 cumulative \n disk mass from center of mass')
axis[1].plot(radius_array, M31_mass_disk, lw = 2, color = 'forestgreen',
             label = 'M31 disk mass')
axis[1].plot(radius_array, MW_mass_disk, lw = 2, color = 'indigo', label = 
             'Milky Way disk mass')
axis[1].legend()
plt.tight_layout()
plt.savefig('MassProfilesMWM31.png')

#now that we have the mass profiles, we can make the velocity profiles
#we are doing this the exact same way, except that we are using
#the circular velocity from the mass profile code that we have imported

MW_velocity_disk = MW_mass.CircularVelocity(radius_array, 2).value
MW_velocity_bulge = MW_mass.CircularVelocity(radius_array, 3).value

M31_velocity_disk = M31_mass.CircularVelocity(radius_array, 2).value
M31_velocity_bulge = M31_mass.CircularVelocity(radius_array, 3).value

#same plotting situation as above and I used the same colors as before
#so that way the comparisons are more obvious

figure, axis = plt.subplots(1,2, figsize = (8,6))
axis[0].set_xlabel('Distance from center of mass (kpc)')
axis[0].set_ylabel('Circular Velocity (km/s)')
axis[0].set_title('Milky Way and M31 disk circular \n velocities from center of mass')
axis[0].plot(radius_array, MW_velocity_disk, lw = 2, color = 'indigo', label = 
             'Milky Way disk velocity')
axis[0].plot(radius_array, M31_velocity_disk, lw = 2, color = 'forestgreen', 
             label = 'M31 disk velocity')
axis[0].legend()

axis[1].set_xlabel('Distance from center of mass (kpc)')
axis[1].set_ylabel('Circular Velocity (km/s)')
axis[1].set_title('Milky Way and M31 bulge circular \n velocities from center of mass')
axis[1].plot(radius_array, MW_velocity_bulge, lw = 2, color = 'indigo', label =
             'Milky Way bulge velocity')
axis[1].plot(radius_array, M31_velocity_bulge, lw = 2, color = 'forestgreen', 
             label = 'M31 bulge velocity')
axis[1].legend()
plt.tight_layout()
plt.savefig('MWM31velocityprofiles.png')
plt.figure()


#From here we can create an iterative function to account for the 
#angular momentum at all of the radial components and all of the velocities

#the unique function that I have created takes the crucial components
#of angular momentum adn multiplies them together as is standard for the
#angular momentum equation. This makes it difficult to exactly cite a paper
#as this is more of a standard and known equation. 
#The goal is to iterate through each point in the file at the same i
#so that at the first point we multiply the first radius, the first velocity,
#and the first mass. Then it all gets added to an angular momentum array
#so that we can plot it against whatever we want to understand what
#the angular momentum depends on. Then we can also repeat it for the disk
#and bulge of each galaxy within the merger remnant
new_radius_array = np.arange(0, 150, step = 1)

angular_momentum = np.zeros(len(new_radius_array))
def angularmomentum(distance, velocity, mass):
    for i in range(len(distance)):
        angular_momentum[i] += 0.5*velocity[i]*((distance[i])**2)*mass[i]
    return angular_momentum

#angular momentum is plotted against the distance from
#the center since this is something both galaxies would have in common
#The angular momentum function is unique to me and to this project

MWdisk_L = angularmomentum(new_radius_array, MW_velocity_disk, MW_mass_disk)
plt.semilogy(radius_array, MWdisk_L, label = 'MW disk L', lw = 2, color = 'indigo')

angular_momentum=np.zeros(len(new_radius_array))
M31_disk_L = angularmomentum(new_radius_array, M31_velocity_disk, M31_mass_disk)
plt.semilogy(radius_array, M31_disk_L, label = 'M31 disk L', lw = 2, color = 'forestgreen')
plt.xlabel('Distance from center of mass (kpc)')
plt.ylabel('Angular momentum (M$_{\odot}$kpc$^{2}$s$^{-1}$)')
plt.title('Angular Momentum vs distance from COM: Disks')
plt.legend()
plt.savefig('Ang_mom_disk.png')
plt.figure()

angular_momentum = np.zeros(len(new_radius_array))
MWbulge_L = angularmomentum(new_radius_array, MW_velocity_bulge, MW_mass_bulge)
plt.semilogy(radius_array, MWbulge_L, label = 'MW bulge L', lw = 2, color = 'indigo')

angular_momentum = np.zeros(len(new_radius_array))
M31_bulge_L = angularmomentum(new_radius_array, M31_velocity_bulge, M31_mass_bulge)
plt.semilogy(radius_array, M31_bulge_L, label = 'M31 bulge L', color = 'forestgreen')
plt.xlabel('Distance from center of mass (kpc)')
plt.ylabel('Angular momentum (M$_{\odot}$kpc$^{2}$s$^{-1}$)')
plt.title('Angular momentum vs distance from COM: Bulges')
plt.legend()
plt.savefig('Ang_mom_bulge.png')





