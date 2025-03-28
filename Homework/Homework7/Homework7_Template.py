
# # Homework 7 Template
# 
# Rixin Li & G . Besla
# 
# Make edits where instructed - look for "****", which indicates where you need to 
# add code. 




# import necessary modules
# numpy provides powerful multi-dimensional arrays to hold and manipulate data
import numpy as np
# matplotlib provides powerful functions for plotting figures
import matplotlib.pyplot as plt
# astropy provides unit system and constants for astronomical calculations
import astropy.units as u
import astropy.constants as const
# import Latex module so we can display the results with symbols
from IPython.display import Latex

# **** import CenterOfMass to determine the COM pos/vel of M33
from CenterOfMass2 import CenterOfMass

# **** import the GalaxyMass to determine the mass of M31 for each component
from GalaxyMass import ComponentMass

# # M33AnalyticOrbit




class M33AnalyticOrbit:
    """ Calculate the analytical orbit of M33 around M31 """
    
    def __init__(self, filename): # **** add inputs
        '''
        inputs:
            filename: the name that we want the file to be
        outputs:
            none; this is setting the constants and such througout when we
            use a file to run this class
        '''
        #get the gravitational constant (the value is 4.498502151575286e-06)
        self.G = const.G.to(u.kpc**3/u.Msun/u.Gyr**2).value
        
        ### **** store the output file name
        self.output = filename
        
        ### get the current pos/vel of M33 
        # **** create an instance of the  CenterOfMass class for M33 
        M33_COM = CenterOfMass('M33_000.txt', 2) #looking for disk
        # **** store the position VECTOR of the M33 COM (.value to get rid of units)
        delta = 0.1
        voldec = 4
        M33_position = M33_COM.COM_P(delta, voldec).value
        # **** store the velocity VECTOR of the M33 COM (.value to get rid of units)
        M33_velocity = M33_COM.COM_V(M33_position[0], M33_position[1], M33_position[2]).value
        
        ### get the current pos/vel of M31 
        # **** create an instance of the  CenterOfMass class for M31 
        M31_COM = CenterOfMass('M31_000.txt', 2) #looking for disk
        # **** store the position VECTOR of the M31 COM (.value to get rid of units)
        delta = 0.1
        voldec = 2
        M31_position = M31_COM.COM_P(delta, voldec).value
        # **** store the velocity VECTOR of the M31 COM (.value to get rid of units)
        M31_velocity = M31_COM.COM_V(M31_position[0], M31_position[1], M31_position[2]).value
        
        ### store the DIFFERENCE between the vectors posM33 - posM31
        # **** create two VECTORs self.r0 and self.v0 and have them be the
        # relative position and velocity VECTORS of M33
        self.r0 = (M33_position - M31_position)
        self.v0 = (M33_velocity - M31_velocity)
        
        ### get the mass of each component in M31 !!!!!!!!!
        ### disk
        # **** self.rdisk = scale length (no units) but it is 5 kpc?
        self.rdisk = 5 
        # **** self.Mdisk set with ComponentMass function. Remember to *1e12 to get the right units. Use the right ptype
        self.Mdisk = ComponentMass('M31_000.txt', 2) * 1e12
        ### bulge
        # **** self.rbulge = set scale length (no units) but it is still 1 kpc ...?
        self.rbulge = 1
        # **** self.Mbulge  set with ComponentMass function. Remember to *1e12 to get the right units Use the right ptype
        self.Mbulge = ComponentMass('M31_000.txt', 3)  * 1e12
        # Halo
        # **** self.rhalo = set scale length from HW5 (no units)
        self.rhalo = 60 #from hw 5
        # **** self.Mhalo set with ComponentMass function. Remember to *1e12 to get the right units. Use the right ptype
        self.Mhalo = ComponentMass('M31_000.txt', 1) * 1e12
    
    
    def HernquistAccel(self, M, ra, r): # it is easiest if you take as an input the position VECTOR 
        """
        inputs:
            M: total halo or bulge mass
            ra: corresponding scale length
            r: a vector of the positions in x, y, and z that we are evaluating
            
        outputs:
            Hern: Hernquist acceleration vector (x,y,z)
        """
        
        ### **** Store the magnitude of the position vector
        rmag = np.sqrt(r[0]**2 + r[1]**2 + r[2]**2)
        
        ### *** Store the Acceleration
        a = -self.G * M #constants, using the self defined G
        b = rmag * ((ra + rmag)**2)
        Hern =   (a/b) * r #follow the formula in the HW instructions
        # NOTE: we want an acceleration VECTOR so you need to make sure that in the Hernquist equation you 
        # use  -G*M/(rmag *(ra + rmag)**2) * r --> where the last r is a VECTOR 
        return Hern
    
    def MiyamotoNagaiAccel(self, M, r_d, r):# it is easiest if you take as an input a position VECTOR  r 
        """ 
        inputs:
            M: Mass of the disk
            r_d: previously defined self.rdisk, or the radius of the disk
            r: vectory of the positions we're evaluating (x,y,z)
        outputs:
            acc: an array of the accelerations for x,y, and z of the disk
        """

        ### Acceleration **** follow the formula in the HW instructions
        z_d = r_d / 5.0
        R = np.sqrt(r[0]**2 + r[1]**2)
        B = r_d + np.sqrt(r[2]**2 + (z_d)**2)
        a = -(self.G * M)
        b = (R**2 + B**2)**(1.5)
        c = a/b
        zstuff = b / (np.sqrt(r[2]**2 + (z_d)**2))
        z_array = np.array([1, 1, zstuff])
        acc = c * r * z_array
        
        
        # AGAIN note that we want a VECTOR to be returned  (see Hernquist instructions)
        # this can be tricky given that the z component is different than in the x or y directions. 
        # we can deal with this by multiplying the whole thing by an extra array that accounts for the 
        # differences in the z direction:
        #  multiply the whle thing by :   np.array([1,1,ZSTUFF]) 
        # where ZSTUFF are the terms associated with the z direction
        #I did all the stuff before I saw the Zstuff buit this seems to work
        #I think so if it doesn't I'll come back up and fix it
        
       
        return acc
        # the np.array allows for a different value for the z component of the acceleration
     
    
    def M31Accel(self, r_array): # input should include the position vector, r
        """
        sums all the acceleration vectors from each of the galaxy components
        inputs:
            r_array: the 3D position vector of where we want to analyze
        outputs:
            a_sum: sum of the accelerations of the halo, bulge, and disk at
            the evaluated point
        """

        ### Call the previous functions for the halo, bulge and disk
        # **** these functions will take as inputs variable we defined in the initialization of the class like 
        # self.rdisk etc. 
        acc_M31halo = self.HernquistAccel(self.Mhalo, self.rhalo, r_array)
        acc_M31disk = self.HernquistAccel(self.Mdisk, self.rdisk, r_array)
        acc_M31bulge = self.HernquistAccel(self.Mbulge, self.rbulge, r_array)
            
            # return the SUM of the output of the acceleration functions - this will return a VECTOR 
        a_sum = acc_M31bulge + acc_M31disk + acc_M31halo
        return a_sum
    
    
    
    def LeapFrog(self, dt, r, v): # take as input r and v, which are VECTORS. Assume it is ONE vector at a time
        """ 
        goal is to integrate forward in time to understand the orbit of M33 specifically
        inputs:
            dt: time interval for integration
            r: starting position vector, relative to M31
            v: velocity vector v of M33 relative to M31
        outputs:
            rnew: the new radius after iterating through time
            vnew: the new velocity after iterating through time
        """
        
        # predict the position at the next half timestep
        rhalf = r + (v*(dt/2))
        
        # predict the final velocity at the next timestep using the acceleration field at the rhalf position 
        a = self.M31Accel(rhalf).value
        vnew = v + (a*dt)
        
        # predict the final position using the average of the current velocity and the final velocity
        # this accounts for the fact that we don't know how the speed changes from the current timestep to the 
        # next, so we approximate it using the average expected speed over the time interval dt. 
        rnew = rhalf + (vnew*(dt/2))
        
        return rnew, vnew # return the new position and velcoity vectors
    
    def OrbitIntegration(self, t0, dt, tmax):
        """ 
        Using the previous function (LeapFrog) to iterate through time into
        the future to predict orbits
        inputs:
            t0 = the initial time we're starting at
            dt = the time interval for integration
            tmax = the time we are wanting to end at, or at what point in the
            future we're measuring to
        outputs:
            none; saves the text to a file for us to use and plot
        """

        # initialize the time to the input starting time
        t = t0
        
        # initialize an empty array of size :  rows int(tmax/dt)+2  , columns 7
        orbit = np.zeros((int(tmax/dt)+2, 7))
        
        # initialize the first row of the orbit
        orbit[0] = t0, *tuple(self.r0), *tuple(self.v0)
        r = self.r0
        v = self.v0
        # this above is equivalent to 
        # orbit[0] = t0, self.r0[0], self.r0[1], self.r0[2], self.v0[0], self.v0[1], self.v0[2]

        # initialize a counter for the orbit.  
        i = 1 # since we already set the 0th values, we start the counter at 1
        
        # start the integration (advancing in time steps and computing LeapFrog at each step)
        while (t<tmax):  # as long as t has not exceeded the maximal time 
            
            # **** advance the time by one timestep, dt
            t += dt
            # **** store the new time in the first column of the ith row
            orbit[i,0] = t
            
            # ***** advance the position and velocity using the LeapFrog scheme
            # remember that LeapFrog returns a position vector and a velocity vector  
            # as an example, if a function returns three vectors you would call the function and store 
            # the variable like:     a,b,c = function(input)
            rnew , vnew = self.LeapFrog(dt, r, v)
            
    
            # ****  store the new position vector into the columns with indexes 1,2,3 of the ith row of orbit
            # TIP:  if you want columns 5-7 of the Nth row of an array called A, you would write : 
            # A[n, 5:8] 
            # where the syntax is row n, start at column 5 and end BEFORE column 8
            
            
            # ****  store the new position vector into the columns with indexes 1,2,3 of the ith row of orbit
            orbit[i, 1:4] = rnew[0], rnew[1], rnew[2]
            orbit[i, 4:7] = vnew[0], vnew[1], vnew[2]
            
            # **** update counter i , where i is keeping track of the number of rows (i.e. the number of time steps)
            i += 1
            r = rnew
            v = vnew
        
        print(self.output)
        # write the data to a file
        np.savetxt('C:/Users/ellen/OneDrive/Desktop/U of A/spring 25/astr400b/' + self.output, orbit, fmt = "%11.3f"*7, comments='#', 
                   header="{:>10s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}"\
                   .format('t', 'x', 'y', 'z', 'vx', 'vy', 'vz'))
        
        # there is no return function ok :)
#gotta make the orbits first using all the text files that I got


#let's do the M33 orbit and plot it

M33orbit = M33AnalyticOrbit('M33prediction.txt') #this is the name we want
#and we had to pick the file we wanna save it in

M33orbit.OrbitIntegration(0, 0.1, 10) # ~integrating with what we were given ~

M33prediction = np.genfromtxt('M33prediction.txt',comments = '#', names = True)
#setting it to be the prediction now

#now we gotta make the r and v and THEN we can plot lol

r33 = np.sqrt((M33prediction['x'])**2 + (M33prediction['y'])**2 + (M33prediction['z'])**2)
v33 = np.sqrt((M33prediction['vx'])**2 + (M33prediction['vy'])**2 + (M33prediction['vz'])**2)

M33_HW6orbit = np.genfromtxt('Orbit_M33.txt', comments = '#', names = True)
M31_HW6orbit = np.genfromtxt('Orbit_M31.txt', comments = '#', names = True)

#now we have to calculate the magnitude between m31 and m33, but it could
#be any galaxy so i'll call em galone and galtwo
def magnitude(galone, galtwo):
    #seems easier to separate then calculate so I'm separating here
    x = (galone['x'] - galtwo['x'])**2
    y = (galone['y'] - galtwo['y'])**2
    z = (galone['z'] - galtwo['z'])**2
    vx = (galone['vx'] - galtwo['vx'])**2
    vy = (galone['vy'] - galtwo['vy'])**2
    vz = (galone['vz'] - galtwo['vz'])**2
    position = np.sqrt(x + y + z)
    velocity = np.sqrt(vx + vy + vz)
    
    return position, velocity

#taking these positions and velocities
position, velocity = magnitude(M33_HW6orbit, M31_HW6orbit)

#answering question one by plotting positions over each other
plt.figure()
plt.plot(M33prediction['t'], r33, label = 'M33-M31 prediction')
plt.plot(M33_HW6orbit['t'], position, color = 'red', label = 'M33 - M31 hw6 simulation')
plt.xlabel('time (Gyr)')
plt.ylabel('avg separation (kpc)')
plt.legend(loc = 'upper left')
plt.savefig('HW7_position_projections.png')

plt.figure()

#answering question one by plotting velocities over each other
plt.plot(M33prediction['t'], v33, color = 'goldenrod', label = 'M33-M31 prediction')
plt.plot(M33_HW6orbit['t'], velocity, color = 'navy', label = 'M33 - M31 hw6 simulation')
plt.xlabel('Time (Gyr)')
plt.ylabel('Avg velocity difference (kpc/s)')
plt.legend()
plt.savefig('HW7_velocity_projections.png')

#the questions :)
'''
Question 2: These plots are not similar beyond about 1.5-2 Gyr. After that,
they really don't have much similarity and if anything are almost complete
opposites. The predictions are much simpler than the simulations from HW6
    
Question 3: I think that the reason there is this discrepancy is because 
we can consider more than just the interactions with one another (M33 and M31)
in HW6, for example we have some consideration with the Milky Way. For this
homework, though, we don't have that, specifically in the hernquist acceleration.
    
Question 4: I guess I sort of answered this in question 3 but I would add it
into the hernquist acceleration as well as the Miyamoto acceleration since
these are the ones that take mass and distance into account, but only for the
relation between M31 and M33, not with the MW. When there is this taken more
into account, like in HW6, we see a more accurate plot. 
    
     
'''



