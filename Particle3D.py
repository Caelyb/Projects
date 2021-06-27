"""
CMod Project A - Solar System: Particle3D
This class simulates the interactions of individual particles in the system
"""
import math
import numpy as np

class Particle3D(object):


    """
    Class to describe 3D particles.

    Methods:
    * formatted output
    * kinetic energy
    * velocity update
    * first order position updates
    * second order position updates
    * static method to read and sort input files data
    * static method to find relative position
    * static method to find the force
    * static method to find the potential energy
    """

    def __init__(self, label, mass, pos, vel):

        """
        Initialise a Particle3D instance

        * particle type
        * mass as a float
        * position as a numpy array of x, y and z coordinates
        * velocity as a numpy array of x, y and z coordinates
        * system specific parameter De as a float
        * system specific parameter re as a float
        * system specific parameter alpha as a float
        """

        self.label = label
        self.mass = mass
        self.position = pos
        self.velocity = vel




    def __str__(self):

        """
        Define output format.
        For particle H2, p=(2.0, 0.5, 1.0) this will print as
        " h2 2.0 0.5 1.0 "
        """

        return str(self.label),(float(self.position[0])),(float(self.position[1])), (float(self.position[2]))


    def kinetic_energy(self):

        """
        Return kinetic energy as
        1/2*mass*velocity^2
        """

        return (0.5*self.mass*(np.linalg.norm(self.velocity))**2)


    def velocity_update(self,force, delta_t):

        """
        velocity update,
        v(t+dt) = v(t) + dt*F(t)/m

        :dt: timestep as a float
        :force: force on particle as a float
        """

        self.velocity = self.velocity+delta_t*(force/self.mass)


    def position_update(self,force, delta_t):

        """
        Second-order position update,
        x(t+dt) = x(t) + dt*v(t) + 1/2*dt^2*F(t)

        :force: current force as float
        :delta_t: timestep as float
        """

        self.position = self.position + (delta_t*self.velocity)+(((delta_t**2)*force)/(2*self.mass))


    def from_file(file_handle):

        """
        Read files and convert data to floats

        :file handle 1: particle data and constants
        :file handle 2: initial position and velocity data
        :return Particle3D instance
        """
        line = file_handle.readline()
        items = line.split()


        position = np.array([float(items[2]),float(items[3]),float(items[4])])
        velocity = np.array([float(items[5]),float(items[6]),float(items[7])])
        return Particle3D(items[0],float(items[1]),position,velocity)


    def relative_position(particle1,particle2):

        """
        Find and return the distance between the two particles

        :particle1,particle2: Particle3D instance
        :return relative position as a numpy array
        """

        rel_pos=particle1.position-particle2.position

        return rel_pos
