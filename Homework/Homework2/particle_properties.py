# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 22:00:42 2025

@author: ellen
"""

import numpy as np
import astropy.units as u
from readme import Read
def particleproperties(file, type, number):
    
    #the goal of this code is to use the data we have been given in a file
    #to deduce different properties about the astronomical object we have
    #the data for. This way, we can understand the overall dynamics,
    #motions, and physics of the object we're looking at. For my commenting,
    #my comments will typically be above the line of code they are describing.
    
    #first step is taking in all of the information from the readme file.
    #I had called it a different name than the ReadFile since it made more
    #sense to me to do it that way.
    #We must take the data like this since this is how we returned it in
    #its function, and it allows us to do many things with the information
    #in the file than what we just do in this code so far. 
    
    time, total_particles, data = Read(file)
    
    #now we find the x, y, and z positions for all of the particles
    #we can do this by having python look through the file (now called data)
    #and searching for the column titled x, y, and z accordingly.
    #It will then store all of these values in the appropriate arrays, which
    #we can use later to index and find magnitudes for. These positions
    #are then put into units of kpc using astropy. 
    
    x = data['x']*u.kpc
    y = data['y']*u.kpc
    z = data['z']*u.kpc
    
    #This equation finds the magnitude of the distance by taking the 
    #magnitude of the function through the square root of all of the distances
    #squared. It is already in units of km/s. 
    
    distance_magnitude = (np.sqrt((x**2)+(y**2)+(z**2)))
    
    #This does the same fundamental thing as the previous x, y, and z functions
    #so it can store all of the velocities in arrays for later use. They are
    #then put into units of km/s. 
    
    vx = data['vx']*u.km/u.s
    vy = data['vy']*u.km/u.s
    vz = data['vz']*u.km/u.s
    
    #Similarly to the distance magnitude, this equation takes the magnitude
    #of the velocities so it finds the net value. It is already in units of km/s.
    
    velocity_magnitude = (np.sqrt((vx**2)+(vy**2)+(vz**2)))
    
    #This creates an index of all of the masses under the column 'm' inside
    #the data file and converts them to the mass of the sun. We must multiply
    #it by 10**9 as they are 10**10 solar masses and we must account for this
    
    mass = data['m']*(10**9)*u.M_sun
    
    #this creates a new index to store specific values of the data to a new
    #index we want to make. In this instance, the index is dependent on the
    #type that is a condition to the function. There are 3 possible options
    #and 3 possible values in the file, so it will look for those types we
    #input. It then creates the new x, y, and z arrays, similar to before, that
    #are further dependent on the number. This number is the particle number
    #that we are looking for (in the case of the homework, it is #100).
    #that way it searches the index to confirm that it is the correct type and
    #it takes it at that specific number. It also confirms they are in units of
    #kpc.
    
    #np.around(value, 3)
    index = np.where(data['type'] == type)
    xnew = data['x'][index][number]*u.kpc
    ynew = data['y'][index][number]*u.kpc
    znew = data['z'][index][number]*u.kpc
    
    #Now we are making a new distance magnitude that is dependent on the 
    #number we inputted. I kept the previous value as it is a general case,
    #so that in future iterations of this code it can be adjusted to meet 
    #the various coding needs. It is already in units of kpc.
    
    new_distance_magnitude = (np.sqrt((xnew**2)+(ynew**2)+(znew**2)))
    
    #this will read out the distance that we have just calculated so we know
    #the magnitude of the particles at this specific number and rounds it to
    #3 decimal places
    
    print("The distance is", np.round(new_distance_magnitude, 3))
    
    #Just like the distances, we are taking the velocities according to the
    #index at the specific number that we have inputted and making sure
    #it is in units of km/s
    
    vxnew = data['vx'][index][number]*u.km/u.s
    vynew = data['vy'][index][number]*u.km/u.s
    vznew = data['vz'][index][number]*u.km/u.s
    
    #Taking the velocity magnitude similarly to how we did before, but we 
    #are using the new velocities that we found that are already in units of km/s. 
    
    new_velocity_magnitude = (np.sqrt((vxnew**2)+(vynew**2)+(vznew**2)))
    
    #Prints the velocity magnitude out so we know the velocity of the particles
    #at that specific number and rounds it to 3 decimal places. 
    
    print("The velocity is", np.round(new_velocity_magnitude,3))
    
    #Now, we have to find the mass at the index and number as we have before
    #and multiply it by 10**9 for the same reasons as before and make sure it
    #is in units of solar masses. 
    
    newmass = data['m'][index][number]*(10**9)*u.M_sun
    
    #This print statement tells us what the mass value is for the number at
    #the specific index. It also sets the precedent of rounding to 3 decimal
    #places despite it not being necessary in our given homework assignment. 
    
    print("The mass is", np.round(newmass, 3))
    
    #Finally, we must convert the distance into units of light years using the
    #to function. We then print the distance and put it into 3 decimal places
    #so that it is rounded. 
    
    lyrdistance = new_distance_magnitude.to(u.lyr)
    print("The distance in light years is", np.round(lyrdistance,3))
    
particleproperties('MW_000.txt', 2, 99)

#I found it important to keep the type so that way we can confirm that the
#line we are wanting to look at (the number) is that type, otherwise there
#will be some sort of error. I think it can help so that I can better 
#understand the dynamics if I am able to confirm the type of object it is. 
