import gsd.hoomd
import numpy as np
import argparse


parser = argparse.ArgumentParser(description='Process a GSD file.')

parser.add_argument('filename', type=str, help='The path to the GSD file')
parser.add_argument('Nm', type=int, help='An integer parameter Nm')

args = parser.parse_args()



trajectory1 = gsd.hoomd.open(args.filename)

Nm = args.Nm # The number of monomers in a polymer
N = trajectory1[-1].particles.N # Total number of particles in the sytstem


Rg = np.empty(len(trajectory1)) # To store the Rg of the individual star polymers
MSD = np.empty(len(trajectory1)) # To store the MSD of the star polymers
box = trajectory1[0].configuration.box[0:3]
for frame in range(len(trajectory1)):
    types =  trajectory1[frame].particles.typeid
    positions = trajectory1[frame].particles.position[types==0]
    images = trajectory1[frame].particles.image[types==0]

    unwrap_pos = positions + images * box
    star_pos = unwrap_pos.reshape(len(unwrap_pos)//(Nm-4),Nm-4,3)
    rg_stars = np.std(star_pos, axis=1)
    frame_rg2 = np.mean(np.sum(rg_stars ** 2,axis=1))
    

    Rg[frame] = frame_rg2


monomer_pos = np.empty((len(trajectory1),N,3))
monomer_image = np.empty((len(trajectory1),N,3))
for frame in range(len(trajectory1)):
    monomer_pos[frame] = trajectory1[frame].particles.position
    monomer_image[frame] = trajectory1[frame].particles.image
unwrapped_pos = monomer_pos + monomer_image * box # Unwrapped positions of all the particles in all the frames 
unwrapped_pos_star = unwrapped_pos.reshape((len(trajectory1),N//Nm,Nm,3)) # Reshaped for each star (N_frames, N_stars, Nm, 3)
unwrapped_star_com = np.mean(unwrapped_pos_star, axis = -2) # Calculate the COM for each stars
for del_t in range(len(trajectory1)):
    count = 0
    dr = 0
    for frame_i in range(len(trajectory1)//2):
        if (frame_i + del_t < len(trajectory1)):
            count += 1
            star_dr = unwrapped_star_com[frame_i] - unwrapped_star_com[frame_i + del_t]
            dr += np.mean(star_dr * star_dr)
    dr = dr / count
    MSD[del_t] = dr


if args.filename.endswith('.gsd'):
    msd_filename = args.filename[:-4] + '_MSD.dat'
    Rg_filename = args.filename[:-4] + '_Rg.dat'

with open(msd_filename, 'w') as f:
    for i in range(len(trajectory1)):
        f.write(f"{i} {MSD[i]}\n")
with open(Rg_filename, 'w') as f:
    for i in range(len(trajectory1)):
        f.write(f"{i} {np.sqrt(Rg[i])}\n")

print("data saved!")
