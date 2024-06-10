import sys
import freud
import numpy as np
import gsd.hoomd
import math


if len(sys.argv) != 2:
    print("Usage: python3 snippet.py 'string_input'")
    sys.exit(1)


equi_file = sys.argv[1]


outfile1 = equi_file[:-4]+'_Bd_corr.dat'
outfile2 = equi_file[:-4]+'_star_Bd_corr.dat'

traj = gsd.hoomd.open(equi_file, mode = 'r')
box = traj[-1].configuration.box[0:3]

b_index = list(traj[-1].particles.typeid)
b_index = [i for i, n in enumerate(b_index) if n == 1]
b_index = np.reshape(b_index,(int(traj[-1].particles.N)//65,4))


lifetime = []
lifetime1 = []
t_init = 0
t = 50
dt = traj[1].configuration.step - traj[0].configuration.step
t_final = len(traj) - t
window = t_final - t_init

distances = []
for t2 in range(len(traj)):
    distances1 = []
    aq = freud.locality.AABBQuery(box, traj[t2].particles.position)
    for bonds in aq.query(traj[t2].particles.position, dict(r_max = 0.2, exclude_ii = True)):
        distances1.append(bonds)
    distances.append(distances1)

norm_factor = []

for t2 in range(t_init,t_final):
        norm_fact2 = 0
        set2 = set((x[0], x[1]) for x in distances[t2])
        a1 = [i[0] for i in set2]
        for i in range(len(b_index)):
                setb = set(x for x in b_index[i])
                if len(setb.intersection(a1)) == 4 :
                    norm_fact2 += 1
        norm_factor.append(norm_fact2)

for t1 in range(0,t):
    count = 0
    count_st = 0
    for t2 in range(t_init,t_final):
        count2 = 0
        count3 = 0
        norm_fact2 = 0
        #if t1 == 0:
        norm_fact = len(distances[t2])
        set1 = set((x[0], x[1]) for x in distances[t2 + t1])
        set2 = set((x[0], x[1]) for x in distances[t2])
        a = set1.intersection(set2)
        count2 = len(a)
        a = [i[0] for i in a]
        a = set(a)
        for i in range(len(b_index)):
                setb = set(x for x in b_index[i])
                if len(setb.intersection(a)) == 4 :
                    count3 += 1
        

        count += count2 / norm_fact
        count_st += count3 / norm_factor[t2]
    
    lifetime.append(count/window)
    lifetime1.append(count_st/window)

with open(outfile1, 'w') as f:
    for i in range(len(lifetime)):
        f.write(f"{i*dt} {lifetime[i]}\n")

with open(outfile2, 'w') as f:
    for i in range(len(lifetime1)):
        f.write(f"{i*dt} {lifetime1[i]}\n")

print("Data saved!")

