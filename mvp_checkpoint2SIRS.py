#Checkpoint2 - SIRS model

import numpy as np
import math
import matplotlib as plt
from matplotlib import pyplot as plt
import time


# number of recovered
n1 = 0.0
# number of susceptable
n2 = 0.99
# number of infected
n3 = 0.01

def init(i):
# randomly initialise with 1% infected and no recovered
    array = np.random.choice(a = (2,1,0),size = (i,i), p = (n1,n2,n3))
    return array

def init2(i):
# alternate intialisation the has coloumb of infected for better looking wavep1s,s
    array = np.ones((i,i))
    for j in range(i):
        array[j][0] = 0
    return array

def init_v(i,n4):
    n2 = 0.99-n4
    array = np.random.choice(a = (2,1,0,3),size = (i,i), p = (n1,n2,n3,n4))
    return array

def neighbours(grid,i,j,k):
# checks if any of the four direct neighbours are infected
    inf = grid[(j+1)%i][k]*grid[(j-1)%i][k]*grid[j][(k+1)%i]*grid[j][(k-1)%i]
    return inf

def rules(grid,i,p1,p2,p3):
# takes copy of array to be updated
    update = np.copy(grid)
    rando = np.random.rand()
# randomly selects coordinates
    j = int(np.random.rand()*i)
    k = int(np.random.rand()*i)
# if infected - chance to recover
    if grid[j][k] == 0:
        if p2 > rando:
            update[j][k] = 2
# if susceptable - chance to get infected
    if grid[j][k] == 1:
        i = neighbours(grid,i,j,k)
        if i == 0:
            if p1 > rando:
                update[j][k] = 0
# if recovered - chance to loose immunity
    if grid[j][k] == 2:
        if p3 > rando:
            update[j][k] = 1
    return update

def pandemic(grid,i):
# goes through array and counts infected sites
    infected = 0
    for a in range(i):
        for b in range(i):
            if grid[a][b] == 0:
                infected += 1
    return infected

def main(array,i,N,p1,p2,p3):
    infected = []
    var = 0
    grid = array
# n is sweep number , p is changing random array element
    for n in range(N):
        for p in range(i**2):
            new = rules(grid,i,p1,p2,p3)
            grid = new
# animation
        plt.cla()
        im = plt.imshow(grid, animated = True)
        plt.draw()
        plt.pause(0.0001)

        infected.append(pandemic(grid,i))
        #print("step : "+str(n))
    #print("Average Number of Infected : "+str(infected/N))
    avgI = np.average(infected)
# calculates variance
    for a in range(N):
        var += (infected[a]**2 - avgI**2)/(i**2)
# average fraction of infected
    frac = avgI/float(i**2)
    return var,frac,infected

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def User_Interface():
# asks user what scenario they require
    input1 = int(input("\n Menu :\n 1) Custom state \n 2) Absorbing state \n 3) Equilibrium state \n 4) Cyclic wave state \n 5) Heat map of infected fraction + varience map \n 6) Precise Variance plot \n 7) Herd Immunity simulation \n : "))
    if input1 != 1:
        if input1 != 2:
            if input1 != 3:
                if input1 != 4:
                    if input1 != 5:
                        if input1 != 6:
                            if input1 != 7:
                                print("Error input not recognised")
# i : array dimensions  -  N : number of sweeps (after equilibrium)
    i = int(input(" \n Please input array dimensions : \n : "))
    N = int(input(" \n Please input number of sweeps : \n : "))

# Custom probability set up
    if input1 == 1:
        array = init(i)
        p1 = float(input(" \n Please input infection probability : \n : "))
        p2 = float(input(" \n Please input recovery probability : \n : "))
        p3 = float(input(" \n Please input immunity loss probability : \n : "))
        main(array,i,N,p1,p2,p3)

# absorbing state
    if input1 == 2:
        p1 = 0.9
        p2 = 0.05
        p3 = 0.
        array = init(i)
        main(array,i,N,p1,p2,p3)

# Equilibrium state
    if input1 == 3:
        p1 = 0.75
        p2 = 0.3
        p3 = 0.2
        array = init(i)
        main(array,i,N,p1,p2,p3)

# cyclic wave state
    if input1 == 4:
        p1 = 0.8
        p2 = 0.1
        p3 = 0.01
        array = init2(i)
        main(array,i,N,p1,p2,p3)

# runs over a range of p1 and p3 values to create heat map of infected No at those values (+ varience)
    if input1 == 5:
# opens files to store data - optional
        file1 = open('infected.dat','w')
        file2 = open('var.dat','w')
        size = 20
# charts to store data for heat maps
        chart1 = np.zeros((size,size))
        chart2 = np.zeros((size,size))
        p2 = 0.5
# r keeps track of run
        r = 0
# goes over p1 and p3 range
        for x in range(size):
            p1 = x/float(size)
            for y in range(size):
                p3 = y/float(size)
                array = init(i)
# 100 sweeps to equilibriate
                main(array,i,100,p1,p2,p3)
# runs for 1000 sweeps to get data
                var,frac,infected = main(array,i,N,p1,p2,p3)
                r += 1
# put data into arrays for ploting and files to keep data
                chart1[x][y] = frac
                chart2[x][y] = var
                file1.write('%d\n'%(frac*(i**2)))
                file2.write('%f\n'%(var))
                print("run : "+str(r))
                print("fraction of infected : "+str(frac))
                #print(chart)
        file1.close()
        file2.close()

# plots
        plt.imshow(chart1,extent=(0,1,0,1),interpolation='none',cmap='hot',origin='lower')
        plt.colorbar()
        plt.xlabel('p1')
        plt.ylabel('p3')
        plt.show()

        max = 1
        plt.imshow(chart2,extent=(0,1,0,1),interpolation='none',cmap='hot',origin='lower')
        plt.colorbar()
        plt.xlabel('p1')
        plt.ylabel('p3')
        plt.show()

    if input1 == 6:
        file1 = open('infected_p.dat','w')
        file2 = open('var_p.dat','w')

        p2 = 0.5
        p3 = 0.5
        r = 0

        p1s = []
        variences = []
        inf_lists = []

        for v in range(20,50):
            array = init(i)
            p1 = v/100.
            p1s.append(p1)
# equlibriation
            #main(array,i,100,p1,p2,p3)
# run main for data
            var,frac,infected = main(array,i,N,p1,p2,p3)
            variences.append(var)
            inf_lists.append(infected)
# count runs
            r += 1
            print("run : "+str(r))
# write into files
            file1.write('%d\n'%(frac*(i**2)))
            file2.write('%f\n'%(var))
        file1.close()
        file2.close()

# Bootstrap method
        errors = np.empty(30)
        for y in range(30):
            fakevar = np.empty(20)
            for z in range(20):
                resample = sklearn.utils.resample(inf_lists[y],n_samples=N)
                fakevar[z] = np.var(resample)
            errors[y] = np.var(fakevar)


        plt.errorbar(p1s,variences,errors)
        plt.title('infection probablilty vs Varience')
        plt.xlabel('Infection probability')
        plt.ylabel('Varience')
        plt.savefig("P1_Vs_Var")
        plt.show()

# vaccination scenario
    if input1 == 7:
        file1 = open('Fraction_infected_herd.dat','w')
# lists of lists (5 lists of 100)
        fracs = []
        vacs = []
        for y in range(5):
            p1 = 0.5
            p2 = 0.5
            p3 = 0.5
# lists of 100 values : fraction = infected fraction
            fraction = []
            vaccinated = []
            r = 0

            for n in range(100):
# n4 = vaccinated percentage
                n4 = n/100.
                array = init_v(i,n4)
# equilibration time 100 sweeps
                main(array,i,100,p1,p2,p3)
                var,frac,infected = main(array,i,N,p1,p2,p3)
# append data to lists
                fraction.append(frac)
                vaccinated.append(n4)
                r += 1
                print("Run "+str(r)+" completed")
            fracs.append(fraction)
            vacs.append(vaccinated)
# reset lists to empty for each of the 5 runs
            fraction = []
            vaccinated = []

# lists for plotting (length of 100 each)
        errors = []
        avgV = []
        avgF = []
        for t in range(100):
# lists 1 and 2 takes the corresponding data points from each of the 5 runs
            list1 = [fracs[0][t],fracs[1][t],fracs[2][t],fracs[3][t],fracs[4][t]]
            list2 = [vacs[0][t],vacs[1][t],vacs[2][t],vacs[3][t],vacs[4][t]]
            file1.write('%f\n'%(np.average(list1)))
            avgV.append(np.average(list2))
            avgF.append(np.average(list1))
# standard deviation used for error
            error = np.std(list1)
            errors.append(error)
        file1.close()

        plt.errorbar(avgV,avgF,yerr=errors)
        plt.title('Herd imunity Test')
        plt.xlabel('Vaccinated percentage')
        plt.ylabel('Fraction infected')
        plt.savefig("Herd_immunity")
        plt.show()

User_Interface()
