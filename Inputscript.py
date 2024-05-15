import numpy
import gsd.hoomd
import itertools
import math
import os
import random


N_star = 8000 # no. of stars in the system
f = 4 # functionality
n_m = 16 # no. of monomers in each arm
N = N_star * (f * n_m + 1) # total no. of monomer
N_m = f * n_m + 1 # no. of monomer in each star
N_b = 20 #
ni = N
sticky_flag = 1 # just a flag to indicate whether to put the sticky bead or not 
# if (sticky_flag == 0):
#     N += N_b
# elif (sticky_flag == 1):
#     N += int(2 * N_b)
# ni = N
while (ni <= N):

    N_particles = ni # for fcc 4*m**3, for bcc 3*m**3
    r0 = 100 # lattice constant
    K = math.ceil(N_particles ** (1 / 3)) # no. of particle/dimension
    L =790 # box length
    sigma1 = 8
    sigma2 = 100 # for LJ class of potentials  
    nbnp = int(N_particles / N_b) # to store the num of particles per chain  
    #==============================================run the C scripts for the initial positions============================================
    seed = random.getrandbits(19)
    command1 = "gcc -std=c99 init_pos_random.c -o init -lm"
    command2 = "./init "+str(n_m)+" "+str(f)+" "+str(N_star)+" "+str(int(10*L))+" "+str(int(10*L))+" "+str(int(10*L))+" "+str(int(r0))+" "+str(sigma1)+" "+str(sigma2)+" "+str(sticky_flag)+" "+str(seed)
    rc = os.system(command1)
    rc2 = os.system(command2)
    #=====================================================XXXXXXXXXX======================================================================


    # Intialise systemb
    #N_b = 50
    #m = 1200 # typically m is number which is the measure of number of unit cell/dimension
    #x = numpy.linspace(L , L , K, endpoint=False)
    posfile = f"Init_pos_N_{N_particles}_C_{N_star}.dat"
    initpos = numpy.loadtxt(posfile, unpack=False, delimiter='\t')
    position = list(initpos)
    #position = list(zip(x,x,x))
    frame = gsd.hoomd.Frame()
    frame.particles.N = N_particles
    frame.particles.position = position[0:N_particles]
    frame.particles.types = ['A','B']
    #frame.particles.typeid = list(numpy.random.randint(2, size=N_particles))
    frame.particles.typeid = [0] * int(N_particles)
    frame.configuration.box = [L, L, L, 0, 0, 0]
    frame.particles.diameter = [1.0] * N_particles
    
    for pi in range(0,N_particles,N_m):
        for arms in range(1,f+1,1):
            if (sticky_flag == 1):
                frame.particles.typeid[pi + arms * n_m] = 1
                frame.particles.diameter[pi + arms * n_m] = 0.1 # Line not needed for the RevCross potential

    frame.bonds.N = N_particles - N_star
    frame.bonds.types = ['A-A','A-B']
    frame.bonds.typeid = [0] * (N_particles - N_star)
    bi = 0 #bonds array index
    if (N_particles % N_star == 0):
        ndgrp = numpy.ones(((N_particles - N_star),2))
        for ci in range(0,N_particles, N_m):
            for ai in range(0,f,1): #ai is arm index : 0 1 2 3 ......(f-1)
                ndgrp[bi][0] = ci
                ndgrp[bi][1] = ci + ai * n_m + 1
                bi += 1
                for li in range(ci + ai * n_m + 1,ci + ai * n_m + n_m,1):
                    ndgrp[bi][0] = li
                    ndgrp[bi][1] = li + 1
                    if(li == ci + ai * n_m + n_m - 1):
                       frame.bonds.typeid[bi] = 1 
                    bi += 1

    
    frame.bonds.group = ndgrp


    

    frame.angles.N = (N_particles - N_star)
    frame.angles.types = ['A-A-A','A-O-A']
    frame.angles.typeid = [0] * (N_particles - N_star)

    angi = 0 #angles array index
    if (N_particles % N_star == 0):
        anggrp = numpy.ones(((N_particles - N_star),3))
        for ci in range(0,N_particles, N_m):
            for ai in range(0,f,1): #ai is arm index : 0 1 2 3 ......(f-1)
                anggrp[angi][0] = ci + ai * n_m + 1
                anggrp[angi][1] = ci
                anggrp[angi][2] = ci + ((ai + 1) % f) * n_m + 1
                frame.angles.typeid[angi] = 1
                angi += 1
                for li in range(ci + ai * n_m + 1,ci + ai * n_m + n_m ,1):
                    if(li == ci + ai * n_m + 1):
                        anggrp[angi][0] = ci
                        anggrp[angi][1] = li
                        anggrp[angi][2] = li + 1
                        angi += 1
                    else:
                        anggrp[angi][0] = li - 1
                        anggrp[angi][1] = li
                        anggrp[angi][2] = li + 1
                        angi += 1
    

    
    frame.angles.group = anggrp

    frame.validate()
    #Write the frame into a file
    filename = 'Initial_config_N_' + str(N_particles)+'_L_'+str(L)+'.gsd'
    with gsd.hoomd.open(name=filename, mode='x') as f:
        f.append(frame)

    ni += 50

#============================================Initial config ends here======================================================================

