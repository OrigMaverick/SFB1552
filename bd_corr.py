import sys
import freud
import numpy as np
import gsd.hoomd
import math

# Check if the correct number of arguments is provided
if len(sys.argv) != 2:
    print("Usage: python3 snippet.py 'string_input'")
    sys.exit(1)

# Get the filename from the command line arguments
equi_file = sys.argv[1]

# Remove the last four characters
outfile = equi_file[:-4]+'_Bd_corr.dat'

traj = gsd.hoomd.open(equi_file, mode = 'r')
box = traj[-1].configuration.box[0:3]



lifetime = []
t_init = 0
t = 50
t_final = len(traj) - t
window = t_final - t_init

distances = []
for t2 in range(len(traj)):
    distances1 = []
    aq = freud.locality.AABBQuery(box, traj[t2].particles.position)
    for bonds in aq.query(traj[t2].particles.position, dict(r_max = 0.2, exclude_ii = True)):
        distances1.append(bonds)
    distances.append(distances1)


for t1 in range(0,t):
    count = 0
    for t2 in range(t_init,t_final):
        count2 = 0
        #if t1 == 0:
        norm_fact = len(distances[t2])
        set1 = set((x[0], x[1]) for x in distances[t2 + t1])
        set2 = set((x[0], x[1]) for x in distances[t2])
        count2 = len(set1.intersection(set2))
        count += count2 / norm_fact
    
    lifetime.append(count/window)

with open(outfile, 'w') as f:
    for i in range(len(lifetime)):
        f.write(f"{i} {lifetime[i]}\n")

print("Data saved!")

