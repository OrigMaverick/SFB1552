import numpy as np
import os
import matplotlib.pyplot as plt
import freud
import gsd.hoomd


traj_file = 'Equil_traj_N_13000_rho_0.250.gsd'
trajectory = gsd.hoomd.open(traj_file, mode = 'r')
freud.parallel.set_num_threads(8)
box = trajectory[-1].configuration.box
N_star = 1000
star_n = 13
com_pos = []

for frame_i in range(len(trajectory)):
    
    typeid = trajectory[frame_i].particles.typeid  # Assuming 3D positions
    aq = freud.locality.AABBQuery(box, trajectory[frame_i].particles.position[typeid == 1])
    distances = []
    
    if (frame_i % 100 == 0):
        filename_frame = "NList_B_frame_"+str(frame_i)+".dat"
        f = open(filename_frame, 'a')
        for bonds in aq.query(trajectory[frame_i].particles.position[typeid == 1], dict(num_neighbors = 4, exclude_ii = True)):
            distances.append(bonds[2])
            f.write(f"{bonds[0]}\t{bonds[1]}\t{bonds[2]}\n")
    
        f.close()




def calc_Coord_num(blist, arm_n):
    c_num = 0
    for i in range(0, len(blist), arm_n):
        for j in range(i, i + arm_n):
            if blist[j] < 0.2:
                c_num += 1
    return c_num

def coord_num_distrib(blist, arm_n, frame_num):
    filename = f"Frame_{frame_num}_CN_dist.dat"
    with open(filename, "w") as outputfile:
        for i in range(0, len(blist), arm_n):
            c_num = sum(1 for j in range(i, i + arm_n) if blist[j] < 0.2)
            outputfile.write(f"{c_num}\n")




def load_and_average_frames(start_frame=0, end_frame=3000, step=100):
    # Initialize an empty list to store arrays
    arrays = []

    # Iterate over each frame
    for i in range(start_frame, end_frame+1, step):
        filename = f"Frame_{i}_CN_dist.dat"
        #filepath = os.path.join(directory, filename)
        if os.path.isfile(filename):
            # Load the data from the file
            array = np.loadtxt(filename)
            # Append the array to the list
            arrays.append(array)
    print(np.shape(arrays))
    # Calculate the average of all arrays
    if arrays:
        average_array = np.mean(arrays, axis=0)
        return average_array
    else:
        return None  # Return None if no files were found

def main():
    total_frame = 3000
    frame_skip = 100
    frame_num = 0
    N_nn = 4
    N_B_part = 4000
    star_n = 13
    arm_n = 16
    N_star = 1000

    B_list = np.zeros(N_nn * N_B_part, dtype=float)
    with open("CoordinationNum_t_.dat", "w") as output:
        while frame_num < total_frame:
            filename = f"NList_B_frame_{frame_num}.dat"
            with open(filename, "r") as nlist_filename:
                for line in range(N_B_part * N_nn):
                    _, _, B_list[line] = map(float, nlist_filename.readline().split())
            coord_num = calc_Coord_num(B_list, arm_n)
            cnum = coord_num / N_star
            output.write(f"{frame_num}\t{cnum}\n")
            coord_num_distrib(B_list, arm_n, frame_num)
            frame_num += frame_skip


    average_result = load_and_average_frames()
    hist, bin_edges = np.histogram(average_result, bins = 500, density=False)
    with open('Time_avg_CN_dist.dat', 'w') as f:
        for i in range(len(hist)):
            f.write(f"{bin_edges[i]} {hist[i]}\n")

    print("data saved!")


main()

