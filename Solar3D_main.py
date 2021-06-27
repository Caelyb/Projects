"""
CMod Project A - Solar System: Solar3D_main
This contains two clasess, one to run the code, the other to implement the user
interface
"""
import numpy as np
import math as m
from Solar3D_nbody import n_body
from Particle3D import Particle3D
import matplotlib.pyplot as plt

def main(Number, dt, conversion, c, answer):

    """
    Extracting timestep, Gravitational constant and metre to (chosen distance unit) conversion factor
    from Solar_constants.txt file
    """

    time = 0

    G = 6.674e-11
    N = Number
    delta_t = dt
    conv = conversion
    com = c
    choice = answer
    delta_t_day = 86400/dt

    """
    extract celestial bodies; label, mass, position(x,y,z) and velocity(x,y,z)
    from Solar_Set_up.txt file
    then export data to nbody class which creates system object which contains a list
    of particle objects, one for each body
    """

    file_handle = open("Solar_Set_up.txt","r")
    # number of bodies
    n = len(file_handle.readlines())

     # move reader back to line 1 of the file so future methods read from start of the file
    file_handle.seek(0)

    #initate initial conditions from file in particle3D and Solar3D_nbody classes
    system = n_body.file_input(file_handle)

    # calculates intial forces
    force = system.forces(G)

    # opens file for particle seperations to be stored
    outfile_seperation = open("seperation.txt","w")
    # opens vmd data file for particle positions to be stored
    outfile_vmd = open("vmd_file.xyz" , "w")

    #correct for centre of mass if option selected in UI
    if com == 1:
        system.initial_velocity_correction()

    #set up empty list
    positions = []
    energy_list = []
    time_list = []
    time_period=[]
    new_time=[]
    end_period=[]
    diff_list=[]
    ay_total=[]

    #if selected to calculate orbital data in UI set up lists with values of zero
    if choice == 3 or choice == 4:
        for i in range (n-1):
            ay_total.append(0)
            time_period.append(0)
            new_time.append(0)
            end_period.append(0)

        diff_initial = system.period_seperation()

        diff_list.append(system.seperation_apo_pera())

    """
    Main iteration sequence (velocity verlet), for each timestep;
    1. The number of objects and time is recorded in the vmd file(every 5th step)
    2. Positions of all bodies updated
    3. The positions of each body is recorded in the vmd file(every 20th step)
    4. Forces on all bodies are calculated
    5. Velocity of each body is updated based on the force and previous velocity
    6. (if selected)Total energy of the system is calculated and appended to a list
    7. Time is recorded and force and time are updated before the iteration starts again
    """

    for i in range(N):

        #add partilce count and timestep every 20th iteration to vmd file
        if i%20 == 0:
            outfile_vmd.write(str(n)+"\n" + "point="+str(i/20)+"\n")

        #update positions and add positions to position list
        for j in range (n):
            positions.append(system.particles[j].position)
            system.particles[j].position_update(force[j],delta_t)

            #add particles positins evry 20th step to vmd outfile
            if i%20 == 0:
                x,y,z,a = system.particles[j].__str__()
                outfile_vmd.write("{} {} {} {}\n".format(x,y*conv,z*conv,a*conv))

        force_new = system.forces(G)

        #update velocity
        for j in range (n):
            system.particles[j].velocity_update(0.5*(force[j]+force_new[j]),delta_t)

        if choice == 3 or choice == 4:
            system.apo_perapises(outfile_seperation)
            diff = system.seperation_apo_pera()
            diff_list.append(diff)


        if choice == 2 or choice == 4:
            energy = float(system.total_kinetic() + system.pot_energy(G))
            energy_list.append(energy)
            time_list.append(float(time))

        force = force_new
        time = time + delta_t

        #if selected to calculate orbital data in UI, calculate orbital period
        if choice == 3 or choice == 4:
            diff_part = system.period_seperation()
            for j in range (n-1):
                """
                using scalar product to calculate the angle between the current
                and previous position, this is summed up and when it complete one
                rotaion (2*pi) it calculates the period
                it uses its own time for each particle, which is reset after each
                period calculation
                """
                scalar_product = np.inner(diff_initial[j],diff_part[j])
                dot = np.linalg.norm(diff_initial[j])*np.linalg.norm(diff_part[j])
                ay = m.acos(scalar_product/dot)

                ay_total[j] = ay_total[j]+ay
                time_period[j] = time_period[j]+1
                diff_initial[j] = diff_part[j]

                if ay_total[j] >= (2*m.pi):
                    #calculate period and reset time and angle to 0
                    period_t = (time_period[j]/delta_t_day)
                    ay_total[j] = 0
                    time_period[j] = 0
                    end_period[j] = period_t

    #if selected to calculate orbital data in UI
    if choice == 3 or choice == 4:
        print()
        print("If the period appears as 0 the number of iterations needs to be increased")
        print()
        for j in range (n-1):
            #find name of paricle
            name,x,y,z = system.particles[j+1].__str__()
            #find max and min seperation for particles and print out orbital data
            print (name)
            print ("Apoapsis",max(map(lambda xx: xx[j], diff_list)), "km")
            print ("Periapsis",min(map(lambda xx: xx[j], diff_list)), "km")
            print ("Period", end_period[j], "earth days")

            print()

    """
    Creates lits of energy/average_energy and plots them against time
    in days to show energy fluctuations of the system
    """
    # if selected to calculate energy fluctuations plot in UI
    if choice == 2 or choice == 4:

        #calculate average energy
        norm_energy = []
        energy_average = sum(energy_list)/len(energy_list)
        for i in energy_list:
            norm_energy.append(i/energy_average)

        time_convert = []
        for i in range(N):
            time_convert.append(time_list[i]*(1/delta_t))

        print("Average system energy: "+str(energy_average))

        plt.suptitle('Average Energy vs time')
        plt.xlabel('time (days)')
        plt.ylabel('Percentage of average energy)')
        plt.plot(time_convert, norm_energy)
        plt.show()

    print()
    print("To view an animation of the solar system please type: " +str("vmd vmd_file.xyz\ninto the terminal and be sure to check the recomended display\nsettings in the readme.txt file"))
    print()

def UserInterface():
    """
    The User interface allows the user to choose the Number and size of the timesteps
    then asks which distance units to output the data in
    It also gives the option to correct for centre of mass motion
    finally it gives the user the option to decide which parts of the program they
    want to run, which may save them time if they are not interested in the
    total sytem energy for example
    """
    print()
    print("Please provide only integer inputs")
    print()
    input1 = input("Please choose number of iterations,\n (10000 is recommended) \n answer: ")
    N = int(input1)
    print()

    input2 = input("Please select a timestep in earth hours: \n (12 is recommended) \n answer: ")
    dt = 3600*float(input2)
    print()

    input3 = input("Which distance units would you like your output (vmd file) in? \n 1) Kilometeres \n 2) Astronomical units(recommended) \n 3) Solar radii \n 4) Lightyears \n answer: ")
# inputs need to be converted to floats so they work with if statements
    input_3 = float(input3)

# This checks answer is an option and returns error message and quits program if answer is not an availible choice
    if input_3 != 1:
        if input_3 != 2:
            if input_3 != 3:
                if input_3 != 4:
                    print("User input is not recognised")
                    quit()
# conversion factor chosen based on user choice
    conv = 0
    answer3 = float(input3)
    if answer3 == 1:
        conv = 1e-3
    if answer3 == 2:
        conv = 6.68459e-12
    if answer3 == 3:
        conv = 1.4368e-9
    if answer3 == 4:
        conv = 1.05702341e-16
    print()

    input4 = input("Would you like to correct for initial centre of mass velocity? \n 1) yes \n 2) no \n answer: ")
    com = float(input4)
    print()

    input5 = input("Availible simulation options: \n 1) Run basic simulation \n 2) Run simulation energy fluctuations plot \n 3) Run simulation with Orbital data \n 4) Run simulation with both orbital and energy data\n 5) exit program \n answer: ")

    input_5 = float(input5)

    if input_5 != 1:
        if input_5 != 2:
            if input_5 != 3:
                if input_5 != 4:
                    if input_5 != 5:
                        print("User input is not recognised")
                        quit()

    answer5 = float(input5)

    if answer5 == 5:
        quit()

    main(N, dt, conv, com, answer5)

UserInterface()
