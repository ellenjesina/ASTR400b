# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 13:06:07 2025

@author: ellen
"""

import numpy as np
import astropy.units as u
from readme import Read
import pandas as pd
import matplotlib.pyplot as plt

def ComponentMass(filename, type):
    #The goal of this function is to return the total mass
    #of any desired galaxy component
    #Inputs: 
        #file name, which will tell us which file we are reading
        #type, which will tell us the galaxy component that we are measuring
    #Outputs: 
        #the total mass of the galaxy component that we are measuring,
        #in units of Solar Mass which we will have divided by 1e2 to
        #put it in units of 10^12 solar masses
    
    time, total_particles, data = Read(filename)
    
    #this stores the new data when the data is equivalent to the
    #type that we are looking for to measure
    
    index = np.where(data['type'] == type)
    
    #we are now creating new data that is dependent on the type we are looking
    #for so that we can only take account of the masses for this type
    #instead of all of the mass of the galaxy
    
    newdata = data[index]
    newmass = newdata['m']
    
    #now we are summing over the entire mass of the specific type
    #that we are looking at to account for the total mass
    #which we then divide by 1e2 to account of 10^12 solar masses
    
    totalmass = (np.sum(newmass)*u.Msun) / 1e2
    
    #now we are rounding the total mass to 3 decimal points so that
    #it stays clean (sig figs)
    
    totalnewmass = np.round(totalmass, 3)
    
    return totalnewmass

#Now using the function to find the total galaxy's masses and sum them
#together. This will fill in the table, and we can use it to get fbar,
#and it will tell us about the total stellar mass, dark matter mass, and
#the overall total mass

#The following code does as follows:
    #calculates the total mass of the halo of a galaxy,
    #calculates the total mass of the disk of a galaxy,
    #calculates the total mass of the bulge of a galaxy,
    #calculates the total mass of the galaxy by combining these components,
    #and calculates fbar (stellar mass/total mass) and rounds it to 3 
    #decimal places. This is done through utilizing the function I have 
    #previously written and using the correct text files to get the
    #proper data and mass measurements. 

#as stated above, for the milky way

MW1 = ComponentMass('MW_000.txt', 1)
MW2 = ComponentMass('MW_000.txt', 2)
MW3 = ComponentMass('MW_000.txt', 3)
MWtotal = MW1 + MW2 + MW3
MWfbar = np.round(((MW2 + MW3) / MWtotal), 3)

#as stated above, for M31

M311 = ComponentMass('M31_000.txt', 1)
M312 = ComponentMass('M31_000.txt', 2)
M313 = ComponentMass('M31_000.txt', 3)
M31total = M311 + M312 + M313
M31fbar = np.round(((M312 + M313) / M31total), 3)

#as stated above, for M33

M331 = ComponentMass('M33_000.txt', 1)
M332 = ComponentMass('M33_000.txt', 2)
M333 = ComponentMass('M33_000.txt', 3)
M33total = M331 + M332 + M333
M33fbar = np.round(((M332 + M333) / M33total), 3)

#as stated above, but for the local group, meaning it is the 
#summation of specific components from each of the galaxies
#and then rounded to three decimal places.

local1 = np.round((MW1 + M311 + M331),3)
local2 = np.round((MW2 + M312 + M332), 3)
local3 = np.round((MW3 + M313 + M333), 3)
total_local_group = MWtotal + M31total + M33total
localfbar = np.round(((local2 + local3) / total_local_group), 3)

#compiling the table through putting in all of the data into individual
#lines, and then using pandas dataframe capabilities to put the data
#into a table format
tablenames = ['Galaxy Name', r'Halo Mass 10$^{12}$ M$_{\odot}$', 
              r'Disk Mass 10$^{12}$ M$_{\odot}$',
              r'Bulge Mass 10$^{12}$ M$_{\odot}$', 
              r'Total Mass 10$^{12}$ M$_{\odot}$', r'f$_{bar}$']
line1 = ['Milky Way', MW1.value, MW2.value, MW3.value, 
         MWtotal.value, MWfbar.value]
line2 = ['M31', M311.value, M312.value, M313.value, 
         M31total.value, M31fbar.value]
line3 = ['M33', M331.value, M332.value, M333.value, 
         M33total.value, M33fbar.value]
line4 = ['Local Group', local1.value, local2.value, local3.value, 
         total_local_group.value, localfbar]

#using the pandas function pd.DataFrame to compile it into a table
#where we have the lines first, then we assign the column values as the
#previously stated table names

dataframe = pd.DataFrame([line1, line2, line3, line4], columns = tablenames)
#print(dataframe), used this to double check the table data before
#compiling it into a plot. Kept in case I needed to double check it again. 

#now we are going to need to make a figure of the table, so to do this
#we must utilize matplotlib. Must make subplots and turn off the axes 
#so that matplotlib doesn't try to make a plot, but instead lets us use
#the ax.table function so that we can create a table. We already have
#all of the data stored in the dataframe variable, so we just have to
#assign the values to this table. 

fig, ax = plt.subplots()
ax.axis('off') #turns off the axes so we can just see the table
table = ax.table(cellText = dataframe.values, colLabels = dataframe.columns,
                 loc = 'center')
table.auto_set_font_size(False) #we want to be able to set the value, not ax
table.set_fontsize(7) #this seemed to fit the best
table.scale(1.25, 1.25) #played around with scaling until it fit the text
#and the table in the screen well

#now we just need to show the data and save it as a figure and voila!
plt.show()
plt.savefig('GalaxyMassTable')

