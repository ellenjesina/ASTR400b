# -*- coding: utf-8 -*-
"""
Created on Tue Feb 18 12:42:20 2025

@author: ellen
"""

import numpy as np
import astropy.units as u
from readme import Read
from CenterOfMass import CenterOfMass
import matplotlib.pyplot as plt

class MassProfile:
    
    def __init__(self, galaxyname, snapshotnumber):
            #now we need to reconstruct the file name
            #add a string of the file number to the value '000'
            ilbl = '000' + str(snapshotnumber)
            #now remove all but the last 3 digits
            ilbl = ilbl[-3:]
            self.filename = "%s_"%(galaxyname) + ilbl + '.txt'
            
            #as always must initialize
            self.time, self.total, self.data = Read(self.filename)                                                                                             
            
            # store the mass & positions
            self.m = self.data['m']*u.Msun
            self.x = self.data['x']*u.kpc
            self.y = self.data['y']*u.kpc
            self.z = self.data['z']*u.kpc
            
            #store the name of the galaxy for later (ahem, m33)
            self.gname = galaxyname
            
    def MassEnclosed(self, ptype, radius):
        
        #using previous homework's classes
        centerofmass = CenterOfMass(self.filename, ptype)
        COM = centerofmass.COM_P(0.1)
        
        #find where this is the particle type so we can take the appropriate mass
        index = np.where(self.data['type'] == ptype)
        
        #initializing mass as zeros
        #I went through this a million times with the units and decided to add
        #them here, otherwise I ran into too many issues where it kept squaring
        #the solar masses. So I decided the masses would be initialized in the
        #array to uncomplicate things later
        massarray = np.zeros(len(radius))*u.Msun
        #calling the earlier mass
        mass = self.m[index]
        
        #determining center of mass position with magnitude
        center_radius = np.sqrt((self.x[index] - COM[0])**2 +
                                (self.y[index] - COM[1])**2 +
                                (self.z[index] - COM[2])**2).value
        
        #now we can go through the for loop, assuming radius is an array
        for i in range(len(radius)):
            #this makes it so it's within the radial distance given
            mass_index = np.where(center_radius < radius[i]) 
            #now we are finding the masses in that radius so we can summate
            massinradius = mass[mass_index] 
            #summating the masses
            mass_enclosed = np.sum(massinradius)
            #now we have to add it to the mass array we had earlier
            #and put it in the correct order of magnitude
            massarray[i] = mass_enclosed*1e10
        #gonna need this later   
        return massarray
    
    def TotalMassEnclosed(self, radii):
        #need to make a blank array to add it to by particle type
        #it's gonna be zeroes so we have to do it in the length of the array
        #of radii that we are given
        massparticle = np.zeros(len(radii))*u.Msun
        
        #adding in the caveat of M33 (so no bulge)
        if self.gname == 'M33':
            #since the bulge is particle type 4, I can only go through 1-3
            for i in range(1,3):
                #adding the masses to the mass array I created earlier
                #then adding it to the particle mass
                massarray = self.MassEnclosed(i, radii)
                #now we're adding this mass to the particle mass
                massparticle += massarray
                #print(massparticle.shape)
                
        #for the other galaxies (MW and M31), we don't need a caveat so it's 
        #in the else function, where the range CAN include the bulge (1-4)
        else:
             for i in range(1,4):   
                #then we add its masses
                massarray = self.MassEnclosed(i, radii)
                #now we have to add the mass to the particle mass
                massparticle += massarray
                #this should ensure that the mass we're adding is dependent
                #on the particle so that we don't try to sum over a mass
                #that doesn't exist for M33
                
        return massparticle

    def HernquistMass(self, radius, a, Mhalo):
        """ Function that defines the Hernquist 1990 dark matter mass profile 
    Inputs:
        r: astropy quantity
            Galactocentric distance in kpc
        a: astropy quantity
            scale radius of the Hernquist profile in kpc
        m_halo: float
            total halo mass in units of 1e12 Msun 
        
    Ouputs:
        mass:  astropy quantity
            total mass within the input radius r in Msun
            """
        #this is the same idea as what we did in Lab 4 so I just copied and
        #pasted what I had done. However, I didn't need to correct for the
        #Mhalo units since I do that elsewhere and don't need the redundancy
        b = Mhalo #constants
        c = radius**2/((a + radius)**2)
    
        mass = b*c
        
        return mass
    
    def CircularVelocity(self, radii, ptype):
        
        #in general, v^2 = GM/r, so gotta find G and use the enclosed mass
        #from earlier I beleive
        
        #~conversion~
        from astropy.constants import G
        G = G.to(u.kpc*u.km**2/u.s**2/u.Msun)
        
        #now doing the mass based on particle type from earlier
        mass = self.MassEnclosed(ptype, radii)
        velocity = np.sqrt((G*mass)/radii)
        velocity = np.round(velocity, 3)
        
        return velocity
    
    def TotalCircularVelocity(self, radii):
        #same equation as with circular velocity I just have to use
        #the entire mass enclosed function
        #it's basically like there's a mass - circular velocity relation
        #and a total mass to total circular velocity relation
        
        #~more conversion~
        from astropy.constants import G
        G = G.to(u.kpc*u.km**2/u.s**2/u.Msun)
        
        #using the total mass function from earlier to find the total velocity
        total_mass = self.TotalMassEnclosed(radii)
        
        #equation for velocity but taking all of the mass so it's total velocity
        totalvelocity = np.sqrt((G*total_mass)/radii)
        totalvelocity = np.round(totalvelocity, 3)
        
        #aaaand gonna need this later
        return totalvelocity
    
    def HernquistCircularSpeed(self, radius, a, Mhalo):
        #same equations as before
        
        #~even more conversion why didn't I do this earlier I don't know~
        from astropy.constants import G
        G = G.to(u.kpc*u.km**2/u.s**2/u.Msun)
        
        #using the hernquist mass for hernquist velocity (duh lol)
        mass = self.HernquistMass(radius, a, Mhalo)
        
        #and velocity equation yay
        hernvelocity = np.sqrt((G*mass)/radius)
        hernvelocity = np.round(hernvelocity, 3)
        return hernvelocity
    

#need to put the plotting in a function so that I can call the function to plot
#for every radius I have in the array
#just seemed a little less complicated than calling everything a million times

def PlottheGalaxy(radius_array):
    
    #let's start off with the milky way which I will call milky naturally
    
    milky = MassProfile('MW', 0) #snapshot is 0
    
    #now finding the masses of the various components using their particle types
    halo_mass = milky.MassEnclosed(1, radius_array)
    disk_mass = milky.MassEnclosed(2, radius_array)
    bulge_mass = milky.MassEnclosed(3, radius_array)
    total_mass = milky.TotalMassEnclosed(radius_array)
    
    #gonna plot these before I figure out the hernquist profile
    plt.plot(radius_array, halo_mass, color = 'red', label = 'Halo Mass')
    plt.plot(radius_array, disk_mass, color = 'orange', label = 'Disk Mass')
    plt.plot(radius_array, bulge_mass, color = 'darkgoldenrod', label = 'Bulge Mass')
    plt.plot(radius_array, total_mass, color = 'darkgreen', label = 'Total Mass')
    
    #finally for the masses, we need the hernquist mass to plot in this
    #we will have to sum over all of the mass tho
    #this is also where I added in the mass because when I put it in the array
    #it got way too complicated so I just initialized it here
    hernquist_array = np.zeros(len(halo_mass)) *u.Msun #this is an array we will
    #add all of the masses to so that it's all the masses on the graph
    #we have to find the masses in the dark matter part, which is where the
    #particle type is 1
    darkmassindex = np.where(milky.data['type'] == 1) #dark matter
    #now we go through and add to the hernquist array we initialized
    #to be able to plot the function overall
    for i in range(len(halo_mass)):
        hernquist_array[i] += milky.HernquistMass(radius_array[i], 60*u.kpc, 
                                            np.sum(milky.m[darkmassindex])*1e10)
        #we used a = 60 kpc in class so I'm assuming it's the same, and it works
    #now I can plot this mass too!
    plt.plot(radius_array, hernquist_array, color = 'blue')
        
    #formatting it now :) 
    #and gotta make it logged
    plt.semilogy()
    plt.legend()
    plt.xlabel("Radius in kpc")
    plt.ylabel("Mass enclosed in Msun")
    plt.title("Milky Way enclosed mass by particle type")
    plt.savefig('Masses for Milky Way')
    #resetting the image for the velocities
    plt.figure()
    
    #repeating for M33! pretty easy to just copy and paste, and since I had
    #already explained it a lot earlier I am a little less specific with my
    #commenting throughout
    
    M33 = MassProfile('M33', 0) #snapshot is 0
    
    #now finding the masses of the various components
    halo_mass = M33.MassEnclosed(1, radius_array)
    disk_mass = M33.MassEnclosed(2, radius_array)
    #no bulge mass so it's empty
    total_mass = M33.TotalMassEnclosed(radius_array)
    
    #gonna plot these before I figure out the hernquist profile
    #and gotta make it logged
    plt.plot(radius_array, halo_mass, color = 'red', label = 'Halo Mass')
    plt.plot(radius_array, disk_mass, color = 'orange', label = 'Disk Mass')
    
    plt.plot(radius_array, total_mass, color = 'darkgreen', label = 'Total Mass')
    
    #finally for the masses, we need the hernquist mass to plot in this
    #we will have to sum over all of the mass tho
    hernquist_array = np.zeros(len(halo_mass))*u.Msun #this is an array we will
    #add all of the masses to so that it's all the masses on the graph
    #we have to find the masses in the dark matter part, which is where the
    #particle type is 1
    darkmassindex = np.where(M33.data['type'] == 1) #dark matter
    for i in range(len(halo_mass)):
        hernquist_array[i] += M33.HernquistMass(radius_array[i], 60*u.kpc, 
                                            np.sum(M33.m[darkmassindex])*1e10)
        #we used a = 60 kpc in class so I'm assuming it's the same
    plt.plot(radius_array, hernquist_array, color = 'blue', label = 'Hernquist Mass')
        
    #formatting it now :) 
    plt.semilogy()
    plt.legend()
    plt.xlabel("Radius in kpc")
    plt.ylabel("Mass enclosed in Msun")
    plt.title("M33 enclosed mass by particle type")
    plt.savefig('Masses for M33')
    #resetting the image for the velocities
    plt.figure()
    
    #repeating for M31! also can just copy and paste very nice
    
    M31 = MassProfile('M31', 0) #snapshot is 0
    
    #now finding the masses of the various components
    halo_mass = M31.MassEnclosed(1, radius_array)
    disk_mass = M31.MassEnclosed(2, radius_array)
    bulge_mass = M31.MassEnclosed(3, radius_array)
    total_mass = M31.TotalMassEnclosed(radius_array)
    
    #gonna plot these before I figure out the hernquist profile
    #and gotta make it logged
    plt.plot(radius_array, halo_mass, color = 'red', label = 'Halo Mass')
    plt.plot(radius_array, disk_mass, color = 'orange', label = 'Disk Mass')
    plt.plot(radius_array, bulge_mass, color = 'darkgoldenrod', label = 'Bulge Mass')
    plt.plot(radius_array, total_mass, color = 'darkgreen', label = 'Total Mass')
    
    #finally for the masses, we need the hernquist mass to plot in this
    #we will have to sum over all of the mass tho
    hernquist_array = np.zeros(len(halo_mass))*u.Msun #this is an array we will
    #add all of the masses to so that it's all the masses on the graph
    #we have to find the masses in the dark matter part, which is where the
    #particle type is 1
    darkmassindex = np.where(M31.data['type'] == 1) #dark matter
    for i in range(len(halo_mass)):
        hernquist_array[i] += M31.HernquistMass(radius_array[i], 60*u.kpc, 
                                            np.sum(M31.m[darkmassindex])*1e10)
        #we used a = 60 kpc in class so I'm assuming it's the same
    plt.plot(radius_array, hernquist_array, color = 'blue', label = 'Hernquist Mass')
        
    #formatting it now :) 
    plt.semilogy()
    plt.legend()
    plt.xlabel("Radius in kpc")
    plt.ylabel("Mass enclosed in Msun")
    plt.title("M31 enclosed mass by particle type")
    plt.savefig('Masses for M31')
    #resetting the image for the velocities
    plt.figure()
    
    #and now it's time to plot the rotation curves!
    #I feel like we can use a lot of similar components of the previous code
    #starting with the milky way
    #since it does end up being so similar, I also limited comments since
    #it's basically just the same ideas
    #at a radius of 30 kpc
    
    milky = MassProfile('MW', 0) #snapshot is 0
    
    #particle types
    halo_velocity = milky.CircularVelocity(radius_array, 1)
    disk_velocity = milky.CircularVelocity(radius_array, 2)
    bulge_velocity = milky.CircularVelocity(radius_array, 3)
    total_velocity = milky.TotalCircularVelocity(radius_array)
    
    #now gonna plot again
    plt.plot(radius_array, halo_velocity, color = 'red', label = 'Halo Velocity')
    plt.plot(radius_array, disk_velocity, color = 'orange', label = 'Disk Velocity')
    plt.plot(radius_array, bulge_velocity, color = 'darkgoldenrod', label = 'Bulge Velocity')
    plt.plot(radius_array, total_velocity, color = 'darkgreen', label = 'Total Velocity')
    
    #and finally, we gotta do the hernquist velocity (pretty much same as before)
    hernquist_array = np.zeros(len(halo_velocity))*u.km/u.s #this is an array we will
    #add all of the velocities so that it's all the velocities on the graph
    #we have to find the velocities in the dark matter part, which is where the
    #particle type is 1
    darkvelocityindex = np.where(milky.data['type'] == 1) #dark matter
    for i in range(len(halo_velocity)):
        hernquist_array[i] += milky.HernquistCircularSpeed(radius_array[i], 60*u.kpc, 
                                            np.sum(milky.m[darkvelocityindex])*1e10)
        #we used a = 60 kpc in class so I'm assuming it's the same, also still
        #used the dark matter mass for the halo
    plt.plot(radius_array, hernquist_array, color = 'blue', label = 'Hernquist Velocity')
    
    plt.semilogy()
    plt.legend()
    plt.xlabel("Radius in kpc")
    plt.ylabel("Velocity in km/s")
    plt.title("Milky Way Rotation Curve")
    plt.savefig('Milky Way Rotation Curve')
    #resetting the image for the other velocities
    plt.figure()
    
    #and the best part is being able to copy and paste all of this
    #for the other two galaxies! yayyyy
    
    #M33 is up first
    M33 = MassProfile('M33', 0) #snapshot is 0
    
    halo_velocity = M33.CircularVelocity(radius_array, 1)
    disk_velocity = M33.CircularVelocity(radius_array, 2)
    #no bulge again yahoo
    total_velocity = M33.TotalCircularVelocity(radius_array)
    
    #now gonna plot again
    plt.plot(radius_array, halo_velocity, color = 'red', label = 'Halo Velocity')
    plt.plot(radius_array, disk_velocity, color = 'orange', label = 'Disk Velocity')
    
    plt.plot(radius_array, total_velocity, color = 'darkgreen', label = 'Total Velocity')
    
    #and finally, we gotta do the hernquist velocity (pretty much same as before)
    hernquist_array = np.zeros(len(halo_velocity))*u.km/u.s #this is an array we will
    #add all of the velocities so that it's all the velocities on the graph
    #we have to find the velocities in the dark matter part, which is where the
    #particle type is 1
    darkvelocityindex = np.where(M33.data['type'] == 1) #dark matter
    for i in range(len(halo_velocity)):
        hernquist_array[i] += M33.HernquistCircularSpeed(radius_array[i], 60*u.kpc, 
                                            np.sum(M33.m[darkvelocityindex])*1e10)
        #we used a = 60 kpc in class so I'm assuming it's the same, also still
        #used the dark matter mass for the halo
    plt.plot(radius_array, hernquist_array, color = 'blue', label = 'Hernquist Velocity')
    
    plt.semilogy()
    plt.legend()
    plt.xlabel("Radius in kpc")
    plt.ylabel("Velocity in km/s")
    plt.title("M33 Rotation Curve")
    plt.savefig('M33 Rotation Curve')
    #resetting the image for the other velocities
    plt.figure()
    
    #and once again for M31 :)
    M31 = MassProfile('M31', 0) #snapshot is 0
    
    halo_velocity = M31.CircularVelocity(radius_array, 1)
    disk_velocity = M31.CircularVelocity(radius_array, 2)
    bulge_velocity = M31.CircularVelocity(radius_array, 3)
    total_velocity = M31.TotalCircularVelocity(radius_array)
    
    #now gonna plot again
    plt.plot(radius_array, halo_velocity, color = 'red', label = 'Halo Velocity')
    plt.plot(radius_array, disk_velocity, color = 'orange', label = 'Disk Velocity')
    plt.plot(radius_array, bulge_velocity, color = 'darkgoldenrod', label = 'Bulge Velocity')
    plt.plot(radius_array, total_velocity, color = 'darkgreen', label = 'Total Velocity')
    
    #and finally, we gotta do the hernquist velocity (pretty much same as before)
    hernquist_array = np.zeros(len(halo_velocity))*u.km/u.s #this is an array we will
    #add all of the velocities so that it's all the velocities on the graph
    #we have to find the velocities in the dark matter part, which is where the
    #particle type is 1
    darkvelocityindex = np.where(M31.data['type'] == 1) #dark matter
    for i in range(len(halo_velocity)):
        hernquist_array[i] += M31.HernquistCircularSpeed(radius_array[i], 60*u.kpc, 
                                            np.sum(M31.m[darkvelocityindex])*1e10)
        #we used a = 60 kpc in class so I'm assuming it's the same, also still
        #used the dark matter mass for the halo
    plt.plot(radius_array, hernquist_array, color = 'blue', label = 'Hernquist Velocity')
    
    plt.semilogy()
    plt.legend()
    plt.xlabel("Radius in kpc")
    plt.ylabel("Velocity in km/s")
    plt.title("M31 Rotation Curve")
    plt.savefig('M31 Rotation Curve')
    #resetting the image for the other velocities
    #plt.figure() don't need it since it's the last one
    
#now to use the function with the radii that I have picked

if __name__ == '__main__':
    radius_array =  np.arange(0.1, 30, 2)*u.kpc
    PlottheGalaxy(radius_array)
    
#YAY IT WORKEDDDDD
#sorry I'm so casual in my comments, I tend to do the stream of consciousness
#way of typing when it comes to commenting, but if it's an issue I can adjust
#for the future