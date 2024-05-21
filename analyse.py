import freud
import gsd.hoomd
import numpy as np


#freud.parallel.set_num_threads(nthreads=8)
#num_threads = freud.parallel.get_num_threads()
#print("Number of threads:", num_threads)

freud.parallel.NumThreads(N=8) # Parallelization by the Freud package



#================ Parameters for the analysis============================

filename = 'Equil_traj_N_13000_rho_0.250.gsd'
num_of_stars = 1000
star_monomer = 13
eps = 25
init_frame_rdf = 2000
init_frame_msd = 1200
final_frame = 3000
r_max = 7.0
num_bins = 1000
trajectory = gsd.hoomd.open(name=filename, mode='r')
#=========================================================================


# Function for the RDF and MSD calculation of the star polymers

#def calc_rdf_gsd(trajectory, init_frame, final_frame, r_max, n_bins):
    # Open the GSD trajectory file
    
#freud.parallel.set_num_threads(nthreads=None)
#freud.parallel.NumThreads(N=32)
# Stores the box dimensions
box = trajectory[-1].configuration.box
position = []
for frame in range(init_frame_rdf,final_frame):
    part_type = trajectory[frame].particles.typeid
    position.append(trajectory[frame].particles.position[part_type == 1])
print(np.shape(position))
pos_arr = np.concatenate(position)
rdf = freud.density.RDF(num_bins, r_max)
rdf.compute(system=(box, pos_arr), reset=False)
    #return rdf

#ef calc_msd_com_gsd(trajectory, star_n, N_star, init_frame, fin_frame):
#reud.parallel.set_num_threads(nthreads=None)
#freud.parallel.NumThreads(N=32)
#trajectory = gsd.hoomd.open(name=file_name, mode='r')
com_pos = []
box = trajectory[-1].configuration.box
for frame_i in range(init_frame_msd,final_frame):
    _pos = []
    for star_i in range(1,num_of_stars+1):
        _pos.append(np.sum(trajectory[frame_i].particles.position[(star_i - 1) * star_monomer : star_i * star_monomer], axis = 0) / star_monomer)
    com_pos.append(_pos)   
msd = freud.msd.MSD(box, mode = 'window')
msd.compute(com_pos)
    #return msd





#msqd = calc_msd_com_gsd(trajectory, star_monomer, num_of_stars, init_frame_msd, final_frame)

#rdf_data = calc_rdf_gsd(trajectory, init_frame_rdf, final_frame, r_max, num_bins)



# Saving Data for the calculated RDF amd MSD data

msd_filename = "MSD_com_N_s_{}_S_n_{}_eps_{}.dat".format(num_of_stars,star_monomer,eps)
rdf_filename = "RDF_com_N_s_{}_S_n_{}_eps_{}.dat".format(num_of_stars,star_monomer,eps)

with open(msd_filename, 'w') as f:
    for i in range(len(msd.msd)):
        f.write(f"{i} {msd.msd[i]}\n")

print("MSD data saved to rdf_data.txt")

with open(rdf_filename, 'w') as f:
    for i in range(len(rdf.rdf)):
        f.write(f"{rdf.bin_centers[i]} {rdf.rdf[i]}\n")

print("RDF data saved to rdf_data.txt")

