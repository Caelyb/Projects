#MVP - Checkpoint1 Ising model, Magnet at ultra low temperatures
import numpy as np
import math
import matplotlib
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation

d = 5000
# number of steps

def init(i):
    assembly = np.random.choice(a = (1,-1),size = (i,i), p = (0.5,0.5))
    return assembly
# generate array with random mixture of +1 or -1
def init2(i):
    assembly = np.ones((i,i))
    return assembly
# approximate equlibrium state for glauber at low T

def doublecount(s1,l1,s2,l2,assembly,i):
# used to adjust energy for nearest neighbours in kawasaki
    adj = 0
    if (s1 == s2):
        if (l1 == (l2+1)%i) or (l1 == (l2-1)%i):
            #print("adjustment")
            adj += 4

    if (l1 == l2):
        if (s1 == (s2+1)%i) or (s1 == (s2-1)%i):
            #print("adjustment")
            adj += 4

    return adj

def DeltaE(s,l,assembly,i):
    E = 2.0*assembly[s][l]*(assembly[(s+1)%i][l] + assembly[(s-1)%i][l] + assembly[s][(l+1)%i] + assembly[s][(l-1)%i])
#  works out number of nearest neighbours and uses multiplication laws to ensure correct + or - signs
# % used so that wraps around periodic boundaries
    return E

def Glauber(T,i,assembly):
    #print("")
    l = int(np.random.rand()*i)
    s = int(np.random.rand()*i)
    #print("shell : "+str(s)+"\n"+"item : "+str(l))
# Randomly select item from array of dimesions length = i
    dE = DeltaE(s,l,assembly,i)
    #print("energy change: "+str(E))
    if dE <= 0:
        assembly[s][l] = (-1)*assembly[s][l]
        #print("flip")
# condition for auto flip
    else:
        p = math.exp(-dE/T)
        r = np.random.rand()
        if p > r:
            assembly[s][l] = (-1)*assembly[s][l]
            #print("flip")
# p < r not needed as there is no change to the array

#assembly = init(5)
#print(assembly)
#Glauber(2,5,assembly)
#print(assembly)
# potential flip based on Boltmann weight

def Kawasaki(T,i,assembly):
    l1 = int(np.random.rand()*i)
    s1 = int(np.random.rand()*i)
    #print("")
    #print("shell : "+str(s1)+"\n"+"item : "+str(l1))
    l2 = int(np.random.rand()*i)
    s2 = int(np.random.rand()*i)
    #print("shell : "+str(s2)+"\n"+"item : "+str(l2))
    if (assembly[s1][l1] == assembly[s2][l2]):
        pass
        #print("unchanged")
    else:
# Randomly select items from array of dimesions length = i
        E1 = DeltaE(s1,l1,assembly,i)
        E2 = DeltaE(s2,l2,assembly,i)
        adj = doublecount(s1,l1,s2,l2,assembly,i)
        E = E1 + E2 + adj
        #print("energy change: "+str(E))
        if E <= 0:
            x = assembly[s1][l1]
            assembly[s1][l1] = assembly[s2][l2]
            assembly[s2][l2] = x
            #print("flip")
# condition for auto flip
        else:
            p = math.exp(-E/T)
            r = np.random.rand()
            if p > r:
                x = assembly[s1][l1]
                assembly[s1][l1] = assembly[s2][l2]
                assembly[s2][l2] = x
                #print("flip")

def mag(assembly,i,T):
# magnetisation method
    m = 0
    for j in range(i):
        for k in range(i):
            m += assembly[j][k]
    #print("magnitisation : " + str(m))
    return m

def suspect(assembly,i,Mlist,T):
# susceptability method
    sus = []
#reads mag data from file and appends to list
    m_avg = np.average(Mlist)
# average magnetism value needed for susceptability data
    #print("average magnetisation : " + str(m_avg))
    for g in range(int(d/10.)):
        s = (float(np.average((np.asarray(Mlist))**2))-m_avg**2)/(i*T)
        #print(s)
        sus.append(s)
    s_avg = np.average(sus)
    return m_avg, s_avg

def energy(assembly,i,T):
    energy = 0
    for j in range(i):
        for k in range(i):
            if assembly[j][k] == assembly[(j+1)%i][k]:
                energy -= 1
            else:
                energy += 1
            if assembly[j][k] == assembly[j][(k+1)%i]:
                energy -= 1
            else:
                energy +=1
# only record nearest neighbours above and to the right to avoid double counting energy
    return energy

def heatCapacity(assembly,i,T,Eng):
    c = 0
    E_avg = np.average(Eng)
    for j in range(int(d/10.)):
        c += (np.average(np.asarray(Eng)**2)-E_avg**2)/((i**2)*(T**2))
    C_avg = (c/(d/10))
    return C_avg , E_avg




# see main_G(below) for comments, main_K is the same except uses Kawasaki dynamics
def main_K(T,i,assembly):
    Eng = []
    for n in range(d):
        for k in range(i**2):
            Kawasaki(T,i,assembly)
        if(n%10 == 0):
            #mag(assembly,i,T)
            Eng.append(energy(assembly,i,T))
            plt.cla()
            im = plt.imshow(assembly, animated = True)
            plt.draw()
            plt.pause(0.0001)

    C, E_avg = heatCapacity(assembly,i,T,Eng)
    return C, E_avg
# magnetisation and susceptability are not relevent for kawasaki dynamics (both remain constant)

def main_G(T,i,assembly):
# main method using Glauber dynamics
    Mlist = []
    E = []
    for n in range(d):
        for k in range(i**2):
# n = number of iterations, k = number of flips per sweep
            Glauber(T,i,assembly)
        if(n%10 == 0):
# update measurements every 10 steps
            Mlist.append(mag(assembly,i,T))
            E.append(energy(assembly,i,T))
            plt.cla()
# clear axes - makes code run faster
            im = plt.imshow(assembly, animated = True)
            plt.draw()
            plt.pause(0.0001)
# show animation
    m_avg, s_avg = suspect(assembly,i,Mlist,T)
# get average magnetisation and susceptability for given temperature
    C, E_avg = heatCapacity(assembly,i,T,E)
# close files so they can be opened in read mode instead of write
    return m_avg, s_avg, C, E_avg

def User_Interface():

    input1 = int(input("\nAvailible options :\n 1) Glauber method \n 2) Kawasaki method \n  : "))
    if input1 != 1:
        if input1 != 2:
            print("error input not recognised")
# asks user for dynamic method

    input2 = int(input("\nMenu :\n 1) 1-3 Temperature range (0.1K intervals) \n 2) set custom Temperature \n  : "))
    if input2 != 1:
        if input2 != 2:
            print("error input not recognised")
# asks user for varied temperature run (with graphs) or shorter custom temp run
    if input2 == 2:
        input3 = float(input("\nMenu :\n Please select a temperature (Kelvin): "))
        #if not input3 in [1,2,3,4,5,6,7,8,9,10]:
            #print("error input not recognised")
        T = input3
# user supplies custom T
    input4 = int(input("\nMenu :\n Please select system dimensions : "))
    i = input4
# sets up assembly of + and - spins using predefined function ^^line 13from matplotlib import pyplot as plt
    if input2 == 1:
        m_averages = []
        s_averages = []
        c = []
        e_avg = []
        Temps = []
# empty lists for graphs

        if input1 == 1:
            assembly = init2(i)
            E_G = open('Energy_G.dat','w')
            Cap_G = open('HeatCapacity_G.dat','w')
            M_G = open('mag_G.dat','w')
            S_G = open('susceptability_G.dat','w')
            Temperatures = open('Temps_G.dat','w')
            for o in range(10,31,1):
# temp range (o = t*10) as easier to use integers in loop
                for b in range((i**2)*100):
                    Glauber(((o)/10.0),i,assembly)
# runs for (i**2)*10 flips (100 sweeps) so run starts in equilibrium state
                m_avg, s_avg, C, E_avg = main_G((o/10.0),i, assembly)

                m_averages.append(np.absolute(m_avg))
                s_averages.append(s_avg)
                c.append(float(C))
                e_avg.append(E_avg)
                Temps.append((o)/10.0)
# append data to lists for graphs
                Temperatures.write('%f\n'%((o)/10.0))
                M_G.write('%d\n'%(m_avg))
                S_G.write('%f\n'%(s_avg))
                E_G.write('%f\n'%(E_avg))
                Cap_G.write('%f\n'%(C))
                print("temp : "+str((o)/10.0)+"K")
# write in data files
            S_G.close()
            M_G.close()
            E_G.close()
            Cap_G.close()
            Temperatures.close()

            plt.close()
            plt.plot(Temps, m_averages)
            plt.title('Magnetisation Vs Temperature')
            plt.xlabel('Temperature')
            plt.ylabel('Magnetisation')
            plt.savefig("Magnetisation graph")
            plt.show()
            # print and save graph of magnetisation vs temp, must change name by hand or will be over writen each run
            plt.plot(Temps, s_averages)
            plt.title('Susceptability Vs Temperature')
            plt.xlabel('Temperature')
            plt.ylabel('Magnetic Susceptability')
            plt.savefig("susceptability graph")
            plt.show()

# print statement so user can see how far through run they are
        if input1 == 2:
            E = open('Energy_K.dat','w')
            Cap = open('HeatCapacity_K.dat','w')
            Temperatures = open('Temps_k.dat','w')
            assembly = init(i)
            for o in range(10,31,1):
                for b in range((i**2)*200):
                    Kawasaki((o/10.0),i,assembly)
                C, E_avg = main_K(((o)/10.0),i,assembly)
                #m_averages.append(np.absolute(m_avg))
                #s_averages.append(s_avg)
                c.append(float(C))
                e_avg.append(E_avg)
                Temps.append(float(o)/10.0)

                Temperatures.write('%f\n'%((o)/10.0))
                E.write('%f\n'%(E_avg))
                Cap.write('%f\n'%(C))
                print("temp : "+str(float(o)/10.0)+"K")
            plt.close()
            E.close()
            Cap.close()
            Temperatures.close()
# the same but using kawasaki instead of glauber


# print and save graph of susceptability vs temp, must change name by hand or will be over writen each run
        plt.plot(Temps, e_avg)
        plt.title('Energy Vs Temperature_K')
        plt.xlabel('Temperature')
        plt.ylabel('Energy')
        plt.savefig("Energy graph")
        plt.show()
# energy graph
        plt.plot(Temps, c)
        plt.title('heat capacity Vs Temperature_K')
        plt.xlabel('Temperature')
        plt.ylabel('heat capacity')
        plt.savefig("Heat capacity graph")
        plt.show()
# heat capacity graph

    if input2 == 2:
        if input1 == 1:
            assembly = init2(i)
            #for b in range((i**2)*100):
                #Glauber(T,i,assembly)
            main_G(T,i,assembly)
        if input1 == 2:
            #for p in range((i**2)*100):
                #Kawasaki(T,i,assembly)
            assembly = init(i)
            main_K(T,i,assembly)
# runs method based on dynmics used for a custom temperature and assembly dimensions
User_Interface()
