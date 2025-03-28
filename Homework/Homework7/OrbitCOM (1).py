

# Homework 6 Template
# G. Besla & R. Li

# import modules
import numpy as np
import astropy.units as u
from astropy.constants import G

# import plotting modules
import matplotlib.pyplot as plt
import matplotlib

# my modules
from readme import Read
# Step 1: modify CenterOfMass so that COM_P now takes a parameter specifying 
# by how much to decrease RMAX instead of a factor of 2
from CenterOfMass2 import CenterOfMass


def OrbitCOM(galaxyname, start, end, n):
    """function that loops over all the desired snapshots to compute the COM 
    pos and vel as a function of time.
    inputs:
        galaxy name: the name of the specific galaxy we're looking at (MW, M31, M33)
        start: the first snapshot we want to look at
        end: the last snapshot we want to look at
        n: integer to define the intervals we want to return the COM
          
    outputs: 
        we will save the information to a file
    """
    
    # compose the filename for output
    fileout = "Orbit_"+str(galaxyname)+'.txt'
    
    #  set tolerance and VolDec for calculating COM_P in CenterOfMass
    # for M33 that is stripped more, use different values for VolDec
    delta, volDec = 0.1, 2
    M33delta, M33volDec = 0.1, 4
    
    # generate the snapshot id sequence 
    # it is always a good idea to also check if the input is eligible (not required)
    snaps_ids = np.arange(start, end+1, n)
    #I had to add +1 otherwise it went to 795, not 800
    
    
    # initialize the array for orbital info: t, x, y, z, vx, vy, vz of COM
    orbit_array = np.zeros([len(snaps_ids), 7])
    
    
    
    # a for loop 
    for i, snap_id in enumerate(snaps_ids): # loop over files
        
        # compose the data filename (be careful about the folder)
        ilbl = '000' + str(snap_id) 
        #remove all but last 3 numbers so we know what snapid it is
        ilbl = ilbl[-3:]
        
        #compiling it so it matches the format of other files and is 
        #easily discernible
        filename = './'+galaxyname+'/'+ galaxyname + '_' + str(ilbl) + '.txt'
        
        # Initialize an instance of CenterOfMass class, using disk particles
        CenterofMass = CenterOfMass(filename, 2)


        # Store the COM pos and vel. Remember that now COM_P required VolDec
        #have to specify if it's m33 or not because of the different values
        if galaxyname == 'M33':
            COM_P = CenterofMass.COM_P(M33delta, M33volDec)
            COM_V = CenterofMass.COM_V(COM_P[0], COM_P[1], COM_P[2])
        else:
            COM_P = CenterofMass.COM_P(delta, volDec)
            COM_V = CenterofMass.COM_V(COM_P[0], COM_P[1], COM_P[2])
    
        # store the time, pos, vel in ith element of the orbit array,  without units (.value) 
        # note that you can store 
        # a[i] = var1, *tuple(array1)
        #there should be 7 that we use since we made an array w 7 columns
        time = CenterofMass.time / 1000
        orbit_array[i][0] = time.value
        orbit_array[i][1] = COM_P[0].value #x value
        orbit_array[i][2] = COM_P[1].value #y value
        orbit_array[i][3] = COM_P[2].value #z value
        orbit_array[i][4] = COM_V[0].value #vx
        orbit_array[i][5] = COM_V[1].value #vy
        orbit_array[i][6] = COM_V[2].value #vz
        
        # print snap_id to see the progress
        #print(snap_id) commenting it out bcuz I don't wanna see 800 three times lol
        
    # write the data to a file
    # we do this because we don't want to have to repeat this process 
    # this code should only have to be called once per galaxy.
    np.savetxt(fileout, orbit_array, fmt = "%11.3f"*7, comments='#',
               header="{:>10s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}"\
                      .format('t', 'x', 'y', 'z', 'vx', 'vy', 'vz'))


# Recover the orbits and generate the COM files for each galaxy
# read in 800 snapshots in intervals of n=5
# Note: This might take a little while - test your code with a smaller number of snapshots first! 
#OrbitCOM('MW', 0, 800, 5)
#OrbitCOM('M31', 0, 800, 5)
#OrbitCOM('M33', 0, 800, 5)
#I commented these out bcuz they already uploaded data and if I don't my code
#takes literally 5 years to run

#checked them all and they worked yay now I can go to next step

# Read in the data files for the orbits of each galaxy that you just created
# headers:  t, x, y, z, vx, vy, vz
# using np.genfromtxt
#plug in the txt files I saved using the functions I ran a few lines earlier
MWdata = np.genfromtxt("Orbit_MW.txt", dtype = None, names = True, skip_header = 0)
M31data = np.genfromtxt("Orbit_M31.txt", dtype = None, names = True, skip_header = 0)
M33data = np.genfromtxt("Orbit_M33.txt", dtype = None, names = True, skip_header = 0)

# function to compute the magnitude of the difference between two vectors 
# You can use this function to return both the relative position and relative velocity for two 
# galaxies over the entire orbit  
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

# Determine the magnitude of the relative position and velocities 

# of MW and M31
#first we need to find the positions and the velocities of both of these
MWpos = np.zeros([len(MWdata['x']), 3]) #initial position array at x pos
for i in range(0, len(MWdata['x'])):
    MWpos[i] += [MWdata['x'][i], MWdata['y'][i], MWdata['z'][i]]

MWvel = np.zeros([len(MWdata['vx']), 3]) #literally exact same idea but velocities
for i in range(0, len(MWdata['vx'])):
    MWvel[i] += [MWdata['vx'][i], MWdata['vy'][i], MWdata['vz'][i]]

M31pos = np.zeros([len(M31data['x']),3])
for i in range(0, len(M31data['x'])):
    M31pos[i] += [M31data['x'][i], M31data['y'][i], M31data['z'][i]]
    
M31vel = np.zeros([len(M31data['vx']), 3]) 
for i in range(0, len(M31data['vx'])):
    M31vel[i] += [M31data['vx'][i], M31data['vy'][i], M31data['vz'][i]]
    
MWM31pos = np.zeros(len(MWdata['x'])) #empty array to store
for i in range(len(MWdata['x'])):
    MWM31pos[i] = vectors(MWpos[i], M31pos[i])
    
MWM31vel = np.zeros(len(MWdata['vx']))#empty array again
for i in range(len(MWdata['vx'])):
    MWM31vel[i] = vectors(MWvel[i], M31vel[i])


# of M33 and M31
#don't need to reload M31 so just putting M33 here
M33pos = np.zeros([len(M33data['x']),3])
for i in range(0, len(M33data['x'])):
    M33pos[i] += [M33data['x'][i], M33data['y'][i], M33data['z'][i]]
    
M33vel = np.zeros([len(M33data['vx']), 3]) 
for i in range(0, len(M33data['vx'])):
    M33vel[i] += [M33data['vx'][i], M33data['vy'][i], M33data['vz'][i]]

M31M33pos = np.zeros(len(M31data['x'])) #empty array again
for i in range(len(M33data['x'])):
    M31M33pos[i] = vectors(M33pos[i], M31pos[i])
    
M33M31vel = np.zeros(len(M31data['vx']))#empty array again
for i in range(len(M33data['vx'])):
    M33M31vel[i] = vectors(M33vel[i], M31vel[i])

# Plot the Orbit of the galaxies 
#################################
#MW and M31
plt.plot(MWdata['t'], MWM31pos, color = 'turquoise', lw=3)
plt.xlabel("time (Gyr)")
plt.ylabel("Position (kpc)")
plt.title("MW and M31 position")
plt.savefig("MWM31positions.png")
plt.figure()

#M31 M33
plt.plot(M31data['t'], M31M33pos, color = 'rebeccapurple', lw=3)
plt.xlabel("time (Gyr)")
plt.ylabel("Position (kpc)")
plt.title("M31 and M33 positions")
plt.savefig("M31M33positions.png")
plt.figure()

# Plot the orbital velocities of the galaxies 
#################################
#MW and M31
plt.plot(MWdata['t'], MWM31vel, color = 'goldenrod', lw=3)
plt.xlabel("time (Gyr)")
plt.ylabel("Velocity (km/s)")
plt.title("MW and M31 velocity")
plt.savefig("MWM31velocities.png")
plt.figure()

#M31 and M33
plt.plot(M31data['t'], M33M31vel, color = 'maroon', lw=3)
plt.xlabel("time (Gyr)")
plt.ylabel("Velocity (km/s)")
plt.title("M33 and M31 velocity")
plt.savefig("M31M33velocities.png")

'''
Answers to the questions:
    1) They will have two close encounters before eventually merging on the
    third as that is when the position reaches zero, implying no difference
    between them and they are merged
    2) As the two galaxies get closer together, their velocities increase
    Since velocity increases as distance decreases, they are inversely related
    3) M31 and MW merge after about 6 billion years. At this time, the orbits
    of M31 and M33 continue to have close encounters, though they get closer over
    time, too, meaning they will likely merge soon. With the start and end we've
    chosen, we are unable to exactly see their merging, though we can certainly
    see traces in their positions and even their velocities as the velocity approaches
    zero and then greatly shoots up before somewhat evening back out again in
    resonance with their encounters. 
    
    
'''


