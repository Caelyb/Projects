# checkpoint 2 - MVP - Game of life
import numpy as np
import math
import matplotlib as plt
from matplotlib import pyplot as plt
import pylab

# array dimensions
i = 50
# number of steps
N = 0

# initialise random array
def init1(i):
    array = np.random.choice(a = (1,0),size = (i,i), p = (0.5,0.5))
    return array

# initialise oscillator
def init2(i):
    array = np.zeros((i,i))
    array[(i+24)%i][(i+25)%i] = 1
    array[(i+25)%i][(i+25)%i] = 1
    array[(i+26)%i][(i+25)%i] = 1
    return array

# initialise glider
def init3(i):
    array = np.zeros((i,i))
    array[4][4] = 1.
    array[4][3] = 1.
    array[4][5] = 1.
    array[3][5] = 1.
    array[2][4] = 1.
    return array

# method to count live cell neighbours
def neighbours(array,i,j,k):
    livecount = 0
    if array[j][(k+1)%i] == 1:
        livecount += 1
    if array[j][(k-1)%i] == 1:
        livecount += 1
    if array[(j-1)%i][k] == 1:
        livecount += 1
    if array[(j-1)%i][(k-1)%i] == 1:
        livecount += 1
    if array[(j-1)%i][(k+1)%i] == 1:
        livecount += 1
    if array[(j+1)%i][(k+1)%i] == 1:
        livecount += 1
    if array[(j+1)%i][(k-1)%i] == 1:
        livecount += 1
    if array[(j+1)%i][k] == 1:
        livecount += 1
    return livecount

# method to apply rules based on live neighbours
def rules(array,i):
    update = np.copy(array)
    for j in range(i):
        for k in range(i):
            livecount = neighbours(array,i,j,k)
# rule 1
            if livecount < 2:
                update[j][k] = 0
# rule 2
            if array[j][k] == 1:
                if livecount == 2 or livecount == 3:
                    update[j][k] = 1
# rule 3
            if livecount > 3:
                update[j][k] = 0
# rule 4
            if livecount == 3:
                update[j][k] = 1
    return update

def sitecount(array,i):
# counts alive elements in the array
    count = 0
    for a in range(i):
        for b in range(i):
            if array[a][b] == 1:
                count += 1
    return count

def glide(grid,i):
# finds alive elements positions (j,K) and uses pythagoras to find r magnitude
    posj = 0
    posk = 0
    posr = 0
    for j in range(i):
        for k in range(i):
            if grid[j][k] == 1.:
                if j > 2 and k > 2:
                    if j < i-3 and k < i-3:
                        posj = j
                        posk = k
                        posr += (math.sqrt(posj**2+posk**2))
    r = posr/5
    print(r)
    return r

def time(list,N):
# checks over 5 sweeps, if same number of alive elements then prints equilibrium time
    eq = 0
    for n in range(N):
        if list[(n+4)%N] == list[(n+3)%N] == list[(n+2)%N] == list[(n+1)%N] == list[n]:
# will get error if equilibrium is reached in last few n updates
            eq = n
            break
    print("iterations until equilibrium : "+str(eq))
    return eq

# main used for random initialisation
def main1(array,i,N):
# sets up empty lists for data and grid so array can be updated properly
    list = []
    grid = array
    positions = []
    for n in range(N):
# applies game of life rules to array
        new = rules(grid,i)
        grid = new
# animation
        plt.cla()
        im = plt.imshow(grid, animated = True)
        plt.draw()
        plt.pause(0.0001)
# active sites - equilibrium tracker
        active = sitecount(grid,i)
        list.append(active)
        #print("time step : "+str(n+1))
# prints equilibrium time
    t = time(list,N)
    return t

# main used for glider measurments
def main2(array,i,N):
    list = []
    grid = array
    positions = []
    for n in range(N):
        new = rules(grid,i)
        grid = new
# animation
        plt.cla()
        im = plt.imshow(grid, animated = True)
        plt.draw()
        plt.pause(0.0001)

        print("time step : "+str(n+1))
# glider measurements
        r = glide(grid,i)
        positions.append(r)
    #print(positions)
    dist = positions[N-1]-positions[0]
    speed = dist/N
# speed in r - diagonal speed
    print("Glider Speed(r/n^-1) : "+str(speed))

def User_Interface():
    input1 = int(input("\n Menu :\n 1) Game of life \n 2) Oscilator Demonstration \n 3) Glider Demonstration \n 4) Produce Histogram : "))
    if input1 != 1:
        if input1 != 2:
            if input1 != 3:
                if input1 != 4:
                    print("Error input not recognised")
# i is array dimensions , N is number of sweeps
    i = int(input(" \n Please input array dimensions : \n : "))
    N = int(input(" \n Please input Number of steps : \n : "))

# normal game of life - random intialisation
    if input1 == 1:
        array = init1(i)
        main1(array,i,N)

# oscilator demonstration
    if input1 == 2:
        array = init2(i)
        main1(array,i,N)

# glider demonstration + measurements
    if input1 == 3:
        array = init3(i)
        main2(array,i,N)

    if input1 == 4:
        file1 = open('Equilibration_times.dat','w')
        time = []
        r = 0
# run simulation 100 times and put eq times into histogram and datafile
        for a in range(100):
            array = init1(i)
            t = main1(array,i,N)
            time.append(t)
            file1.write('%d\n'%(t))
            r += 1
            print("Run "+str(r)+" Completed")
        file1.close()

# make histogram
        hist = pylab.hist(time, bins = 'auto')

        pylab.xlabel("Time")
        pylab.ylabel("Number of candidates")
        pylab.title("Game of life - Histogram")
        pylab.grid(True)
        pylab.savefig("GoL_Histogram_.pdf")
        pylab.show()

User_Interface()
