# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 13:00:19 2025

@author: ellen
"""

import numpy as np
import astropy.units as u
#

def Read(filename):
    #opening the file with the intention of reading it
    file = open(filename, 'r')
    #reading the first line
    line1 = file.readline()
    #assigning variables to the information in line 1 and then splitting it
    label, value = line1.split()
    #assigning a variable to the value, and putting it in Myr
    time = float(value)*u.Myr
    #now we are reading the second line, which we will then split
    line2 = file.readline()
    #assigning the second label and values to this line
    label2, value2 = line2.split()
    #assigning a variable to the value that we found
    total_particles = float(value2)
    #we are done using this file for now, so we can close it
    file.close()
    data = np.genfromtxt(filename, dtype = None, names = True, skip_header = 3)
    #we are now returning this specific information so that when we call
    #this function in the future, the first three lines are alrady 
    #accounted for and can be used however we need them to. It also
    #starts us off on line 4 which should be the start of the data we
    #generall want to work with
    return time, total_particles, data

#commented out for now, but can replace this with various text files so that
#we can make sure they read in properly. 
#Read("M31_000.txt")