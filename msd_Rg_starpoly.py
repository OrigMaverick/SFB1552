import gsd.hoomd
import numpy as np 
trajectory1 = gsd.hoomd.open('../gauss_pot_eps_25/Equil_traj_N_65000_rho_0.400.gsd')
Rg4 = np.empty(len(trajectory1))
MSD_star = np.empty(len(trajectory1))
N_star = 1000
N = 65000
N_m = 65
box = trajectory[0].configuration.box[0:3]

star_cm = np.empty((len(trajectory1),N_star,3))
inv_N_m = 1 / N_m
inv_N_star = 1 / N_star

for frame in range(len(trajectory1)):
    positions = trajectory1[frame].particles.position
    images = trajectory1[frame].particles.image
    frame_rg2 = 0
    
    for start_id in range(0, N, N_m):
        segment_positions = positions[start_id:start_id + N_m] + box * images[start_id:start_id + N_m]
        R_cm = np.mean(segment_positions, axis=0)
        star_cm[frame][start_id//N_m] = R_cm
        R_g = np.sum(np.mean((segment_positions - R_cm) ** 2, axis=0))
        frame_rg2 += R_g * inv_N_star
    

    Rg4[frame] = frame_rg2



for del_t in range(len(trajectory1)):
    count = 0
    dr = 0
    for frame_i in range(len(trajectory1)//2):
        if(frame_i + del_t < len(trajectory1)):
            count += 1
            dr += np.mean(np.mean((star_cm[frame_i] - star_cm[frame_i + del_t])**2, axis = 0))
    dr = dr / count
    MSD_star[del_t] = dr


with open('msd.dat', 'w') as f:
    for i in range(len(trajectory1)):
        f.write(f"{i} {MSD_star[i]}\n")

with open('Rg_t.dat', 'w') as f:
    for i in range(len(trajectory1)):
        f.write(f"{i} {Rg4[i]}\n")

print("data saved!")
