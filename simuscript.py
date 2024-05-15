
import hoomd
import matplotlib
import numpy
import gsd.hoomd
import itertools
import math
import os
import random

sigma=1.0
N = 520000
ni = 520000
N_b = 20
L = 650

def runLangevin(device, L, N_star, N_arm, f, del_t, e_gauss, init_file, N_steps):

    N_particles = N_star * (f * N_arm + 1) # 1 for the central binding monomer

    if (device == 'gpu' or device == 'GPU') :
        simu = hoomd.Simulation(device=hoomd.device.GPU(),seed=3)
    else:
        simu = hoomd.Simulation(device=hoomd.device.CPU(),seed=3)
    
    integrator = hoomd.md.Integrator(dt=del_t,integrate_rotational_dof=False)
    
    # Open the initial frame

    simu.create_state_from_gsd(filename=init_file, frame=-1)

    # Get the sigma from the file itself

    fileobj = gsd.hoomd.open(init_file, mode = 'r')
    
    for a_frame in fileobj:
        sig_B, sig_A = numpy.unique(a_frame.particles.diameter)



    # Define forces 

    cell = hoomd.md.nlist.Tree(buffer=0.1)
    lj = hoomd.md.pair.LJ(nlist=cell, mode = 'shift')
    lj.params[('A','A')] = dict(epsilon = 1, sigma = sig_A)
    lj.params[('B','B')] = dict(epsilon = 0,sigma = sig_B)
    lj.params[('A','B')] = dict(epsilon = 0, sigma = (sig_A + sig_B)/2)
    lj.r_cut[('A','A')] = math.pow(2.,1./6.) * sig_A  
    lj.r_cut[('B','B')] = False
    lj.r_cut[('A','B')] = False

    gauss = hoomd.md.pair.Gaussian(nlist = cell, default_r_cut = None, mode = 'none')
    gauss.params[('B', 'B')] = dict(epsilon = -e_gauss, sigma = 0.05)
    gauss.r_cut[('B', 'B')] = 3.0
    gauss.params[('A', 'B')] = dict(epsilon = 0.0, sigma = 1.0)
    gauss.r_cut[('A', 'B')] = False
    gauss.params[('A', 'A')] = dict(epsilon = 0.0, sigma = 1.0)
    gauss.r_cut[('A', 'A')] = False

    #bonds-----------------------------------------------------------------------  
    fene = hoomd.md.bond.FENEWCA()
    fene.params['A-A'] = dict(k = 30, r0 = 1.5, epsilon = 1.0, sigma = 1.0, delta = 0.0)
    fene.params['A-B'] = dict(k = 0.0, r0 = 0.765, epsilon = 0.0, sigma = 0.55, delta = 0.0)
    #----------------------------------------------------------------------------- 
    harmonic2 = hoomd.md.angle.Harmonic()
    harmonic2.params['A-A-A'] = dict(k = 20, t0 = 3.141592) 
    harmonic2.params['A-O-A'] = dict(k = 50, t0 = 1.91055375)

    harm = hoomd.md.bond.Harmonic()
    harm.params['A-A'] = dict(k = 0.0, r0 = 1.5)
    harm.params['A-B'] = dict(k = 30000, r0 = 0.50)

    # Set-up the integrator and the thermostat

    integrator.forces.append(lj)
    integrator.forces.append(gauss)
    integrator.forces.append(harmonic2)
    integrator.forces.append(harm)
    integrator.forces.append(fene)

    langevin = hoomd.md.methods.Langevin(filter=hoomd.filter.All(), kT=1.0)
    integrator.methods.append(langevin)
    simu.operations.integrator = integrator

    # Test run to randomize the velocity(which does not really makes much of a difference)

    simu.run(200000)

    #===============================================

    thermodynamic_properties = hoomd.md.compute.ThermodynamicQuantities(filter=hoomd.filter.All()) # Computer for thermodynamic values
    zero_momentum = hoomd.md.update.ZeroMomentum(hoomd.trigger.Periodic(1000)) # To rescale the momentum

    simu.operations.updaters.append(zero_momentum)
    simu.operations.computes.append(thermodynamic_properties)

    # I/O operations with the ongoing simulation 

    #-----------Some parameters for the file naming and compression of the system-----------------------------
    rho_i = N_particles / (L**3)
    #rho_f = N_particles / ((N_arm  * 10)**3)
    rho_f = 0.40
    
    comp_file_name = 'Comp_traj_N_{}_rho_i_{:.3f}_rho_f_{:.3f}.gsd'.format(N_particles,rho_i,rho_f)
    gsd_writer2 = hoomd.write.GSD(trigger=hoomd.trigger.Periodic(50000),filename=comp_file_name,mode='xb',filter=hoomd.filter.All(),dynamic=['property','momentum','topology','attribute'])
    simu.operations.writers.append(gsd_writer2)
    gsd_writer2.write_diameter = True
    simu.run(100000)
    gsd_writer2.flush()
    simu.operations.writers.remove(gsd_writer2)
    #Starting to compress the system

    ramp = hoomd.variant.Ramp(A=0, B=1, t_start=simu.timestep, t_ramp=100000)
    rho = simu.state.N_particles / simu.state.box.volume
    initial_box = simu.state.box
    final_box = hoomd.Box.from_box(initial_box)  # make a copy of initial_box
    final_box.volume = simu.state.N_particles / rho_f
    box_resize_trigger = hoomd.trigger.Periodic(2)
    box_resize = hoomd.update.BoxResize(box1=initial_box, box2=final_box, variant=ramp, trigger=box_resize_trigger)
    simu.operations.updaters.append(box_resize)
    
    simu.run(100001)

    simu.operations.updaters.remove(box_resize)
    
    #simu.run(N_steps * 20)

    equi_file_name = 'Equil_traj_N_{}_rho_{:.3f}.gsd'.format(N_particles,rho_f)
    gsd_writer3 = hoomd.write.GSD(trigger=hoomd.trigger.Periodic(int((N_steps)/1000)),filename=equi_file_name,mode='xb',filter=hoomd.filter.All(),dynamic=['property','momentum','topology','attribute'])
    simu.operations.writers.append(gsd_writer3)
    gsd_writer3.write_diameter = True
    simu.run(N_steps)
    gsd_writer3.flush()

    #equi_file_name = 'Stat_traj_N_{}_rho_{:.3f}.gsd'.format(N_particles,rho_f)
    #gsd_writer4 = hoomd.write.GSD(trigger=hoomd.trigger.Periodic(100),filename=equi_file_name,mode='xb',filter=hoomd.filter.All(),dynamic=['property','momentum','topology','attribute'])
    #simu.operations.writers.append(gsd_writer4)
    #gsd_writer4.write_diameter = True
    #simu.run(100000)
    #gsd_writer4.flush()

# This function is to change the monomer density on an existing simulation trajectory file

def changeRho_NVT(device, L, N_star, N_arm, f, del_t, e_gauss, init_file, N_steps, final_rho):
    N_particles = N_star * (f * N_arm + 1) # 1 for the central binding monomer

    if (device == 'gpu' or device == 'GPU') :
        simu = hoomd.Simulation(device=hoomd.device.GPU(),seed=3)
    else:
        simu = hoomd.Simulation(device=hoomd.device.CPU(),seed=3)
    
    integrator = hoomd.md.Integrator(dt=del_t,integrate_rotational_dof=False)
    
    # Open the initial frame

    simu.create_state_from_gsd(filename=init_file, frame=-1)

    # Get the sigma from the file itself

    fileobj = gsd.hoomd.open(init_file, mode = 'r')
    
    for a_frame in fileobj:
        sig_B, sig_A = numpy.unique(a_frame.particles.diameter)



    # Define forces 

    cell = hoomd.md.nlist.Tree(buffer=0.1,exclusions=('bond',))
    lj = hoomd.md.pair.ExpandedLJ(nlist=cell, mode = 'shift')
    lj.params[('A','A')] = dict(epsilon = 10, sigma = sig_A, delta = 0.0)
    lj.params[('B','B')] = dict(epsilon = 0,sigma = sig_B, delta = 0.0)
    lj.params[('A','B')] = dict(epsilon = 0, sigma = (sig_A + sig_B)/2, delta = 0.0)
    lj.r_cut[('A','A')] = math.pow(2.,1./6.) * sig_A  
    lj.r_cut[('B','B')] = sig_B
    lj.r_cut[('A','B')] = math.pow(2, 1./6.) * (sig_B + sig_A)/2

    gauss = hoomd.md.pair.ExpandedGaussian(nlist = cell, default_r_cut = None, mode = 'none')
    gauss.params[('B', 'B')] = dict(epsilon = -e_gauss, sigma = (sig_B / 10.) * 3.0, delta = sig_B)
    gauss.r_cut[('B', 'B')] = 3.0
    gauss.params[('A', 'B')] = dict(epsilon = 0.0, sigma = 1.0, delta = 0.0)
    gauss.r_cut[('A', 'B')] = 0.0
    gauss.params[('A', 'A')] = dict(epsilon = 0.0, sigma = 1.0, delta = 0.0)
    gauss.r_cut[('A', 'A')] = 0.0

    #bonds-----------------------------------------------------------------------  
    fene = hoomd.md.bond.FENEWCA()
    fene.params['A-A'] = dict(k = 30, r0 = 1.2, epsilon = 0.000, sigma = 1.0, delta = 0.97)
    fene.params['A-B'] = dict(k = 300000, r0 = 0.7, epsilon = 0.000, sigma = 0.1, delta = 0.50)
    #----------------------------------------------------------------------------- 
    harmonic2 = hoomd.md.angle.Harmonic()
    harmonic2.params['A-A-A'] = dict(k = 20, t0 = 3.141592) 
    harmonic2.params['A-O-A'] = dict(k = 50, t0 = 1.91055375)

    # Set-up the integrator and the thermostat

    integrator.forces.append(lj)
    integrator.forces.append(gauss)
    integrator.forces.append(harmonic2)
    integrator.forces.append(fene)

    langevin = hoomd.md.methods.Langevin(filter=hoomd.filter.All(), kT=1.0)
    integrator.methods.append(langevin)
    simu.operations.integrator = integrator

    # Test run to initialize the system before the final monomer density

    simu.run(400)
    if (simu.state.N_particles / simu.state.box.volume > final_rho):
        prename = 'Expansion'
    else :
        prename = 'Compression'

    equi_file = prename + '_traj_N_{}_rho_{:.3f}.gsd'.format(N_particles,final_rho)
    gsd_writer2 = hoomd.write.GSD(trigger=hoomd.trigger.Periodic(20000),filename=equi_file,mode='xb',filter=hoomd.filter.All(),dynamic=['property','momentum','topology','attribute'])
    simu.operations.writers.append(gsd_writer2)
    gsd_writer2.write_diameter = True
    simu.run(4000)
 

    #Starting to compress/expand the system

    ramp = hoomd.variant.Ramp(A=0, B=1, t_start=simu.timestep, t_ramp=N_steps)
    rho = simu.state.N_particles / simu.state.box.volume
    initial_box = simu.state.box
    final_box = hoomd.Box.from_box(initial_box)  # make a copy of initial_box
    final_box.volume = simu.state.N_particles / final_rho
    box_resize_trigger = hoomd.trigger.Periodic(2)
    box_resize = hoomd.update.BoxResize(box1=initial_box, box2=final_box, variant=ramp, trigger=box_resize_trigger)
    simu.operations.updaters.append(box_resize)
    
    simu.run(N_steps+1)

    simu.operations.updaters.remove(box_resize)
    gsd_writer2.flush()
    simu.operations.writers.remove(gsd_writer2)


def main():
    device = 'GPU'
    L = 790
    total_N = 520000
    N_star = 8000
    N_arm = 16
    f = 4
    dt = 0.002
    k_angle = 10
    k_bond = 10
    eps_gauss = 15.0
    init_file_name = "Initial_config_N_{0:d}_L_{1:d}.gsd".format(total_N,L)
    #init_file_name = "Equil_traj_N_65000_rho_0.400.gsd"
    N_steps = 300000000
    if (os.path.exists("equili.txt")):
        print("Opps...The code is not running")
    else:
        print("The run parameters are :")
        print("========================")
        print("\nN : {0:d} \nk : {1:d} \nrho_f : {2:.5f}\nepsilon_gauss : {3:.3f}\n".format((N_star * (N_arm * f + 1)),10, 0.20, eps_gauss))
        runLangevin(device,L,N_star,N_arm,f,dt,eps_gauss,init_file_name,N_steps)

    #changeRho_NVT('GPU',200,N_star,N_arm,f,dt,eps_gauss,'Equil_traj_N_65000_rho_0.400.gsd',N_steps,0.008)

if __name__ == "__main__":
    main()

