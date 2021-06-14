"""
CMod Project A - Solar System: Solar3D_nbody
This class simulates the interactions between the particles in the system
"""

import math
import numpy as np
from Particle3D import Particle3D


class  n_body(object):

    """
    Class to describe 3D particles interactions.

    Methods:
    * formatted output
    * total kinetic energy of system
    * potential energy between two bodies
    * total potential energy
    * force between two particles in system
    * total force of system
    * centre of mass velocity correction
    * magnitude of seperation from first particle - written to file
    * vector seperation from first particle
    * magniude of seperation from first particle
    * initilsing initial inputs from file and Particle3D class
    """

    def __init__(self, particles, label):
        """
        Initialise a n_body instance
        * system name
        * system of paricles
        """
        self.particles = particles
        self.label = label


    def __str__(self):
        """
        Define output format of system of particles name
        """

        return str(self.label)


    def total_kinetic(self):
        """
        Return kinetic energy as the sum of kinetic energs of bodies
        :uses kinetic_energy from Particle3D class

        """
        ke = 0
        for i in range(len(self.particles)):
            ke += self.particles[i].kinetic_energy()

        return ke


    def pot_energy_method(p1,p2,G):
        """
        Return potential energy between two particles in system where
        potential energy = -G*m1*m2/|r|
        """
        r = np.linalg.norm(Particle3D.relative_position(p1,p2))
        pe = -(G*p1.mass*p2.mass)/(r)

        return pe


    def pot_energy(self,G):
        """
        Return sum of all potential energy between particles in system
        :uses potential energy method method
        """
        pe=0

        for i in range(len(self.particles)-1):
            for j in range(i+1,len(self.particles)):
                pe += n_body.pot_energy_method(self.particles[i],self.particles[j],G)

        return pe


    def force_method(p1,p2, G):
        """
        Return force between two particles in system where
        F=-G*m1*m2/|r|^3
        """
        r = Particle3D.relative_position(p1,p2)
        f = -(G*(p1.mass*p2.mass))/(np.linalg.norm(r)**3)

        return f*r


    def forces (self,G):
        """
        Return list of forceinteraction, with each force interaction the sum of
        the interactions on that particle using newtons third law.
        :uses force method method
        """
        f_t=[]

        for i in range(len(self.particles)):
            f_t.append(np.zeros(3))

        for i in range(len(self.particles)-1):
            for j in range(i+1,len(self.particles)):
                f = n_body.force_method(self.particles[i],self.particles[j],G)
                f_t[i]+=f
                f_t[j]-=f

        return f_t


    def initial_velocity_correction(self):
        """
        Correct for initial centre of mass of system where
        momentum = the sum of mi*vi, and
        Velocity of centre of mass = momentum * 1/sum of mi
        This is subtracted from the initial velocity
        """
        tot_momentum=np.zeros(3)
        tot_mass=0

        for i in range(len(self.particles)):
            tot_momentum +=self.particles[i].mass*self.particles[i].velocity
            tot_mass +=self.particles[i].mass

        vel_com=(tot_momentum/tot_mass)

        for i in range(len(self.particles)):
            self.particles[i].velocity=self.particles[i].velocity-vel_com


    def apo_perapises(self,outfile_seperation):
        """
        Calculate magnitude of seperation between the first particles (the sun)
        and the other particles and output this to a file
        """
        n=len(self.particles)
        seperation=[]

        for i in range (n-2):
            seperation.append(np.linalg.norm(Particle3D.relative_position(self.particles[0],self.particles[i+1])))

        #This calculates the magnetude of the seperation between the moon and earth
        seperation.append((np.linalg.norm(Particle3D.relative_position(self.particles[3],self.particles[n-1]))))
        outfile_seperation.write(str(seperation)+"\n")


    def period_seperation(self):
        """
        Return the seperation between the first particles (the sun)
        and the other particles and return it as a list
        """
        sep=[]

        for i in range (len(self.particles)-2):

            sep.append((Particle3D.relative_position(self.particles[0],self.particles[i+1])))

        #This calculates the seperation between the moon and earth
        sep.append(((Particle3D.relative_position(self.particles[3],self.particles[(len(self.particles))-1]))))

        return(sep)


    def seperation_apo_pera(self):
        """
        Return the magnitude between seperation between the first particles
        (the sun) and the other particles and return it as a list
        """
        n= len(self.particles)
        diff=np.zeros(n-1)

        for i in range (n-2):
            diff[i]=(np.linalg.norm(Particle3D.relative_position(self.particles[0],self.particles[i+1])))

        #This calculates the seperation between the moon and earth
        diff[n-2]=((np.linalg.norm(Particle3D.relative_position(self.particles[3],self.particles[n-1]))))

        return (diff)


    def file_input(file_handle):
        """
        Read files, send to particle 3D from file method to create its Particle
        3D instance class, using this to create the Solar3D_nbody instance class
        containing Particle 3D instance classes.

        :file handle: initial position and velocity data
        :send file to Particle 3D from file method
        :return Solar3d_nbody instance containing Particle3D instances
        """
        system = []
        n = len(file_handle.readlines())
        file_handle.seek(0)

        for i in range(n):
            k = Particle3D.from_file(file_handle)
            system.append(k)

        return n_body(system, "solar system")
