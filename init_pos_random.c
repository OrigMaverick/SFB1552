#include <stdlib.h>
#include <stdio.h>
#include "RandomGenerator.h"
#include <math.h>

#define twoPI 6.28318530718
#define PIby4 0.78539816339
#define PIby2 1.57079632679
typedef struct
{
    int pid;
    float px,py,pz,sigma;
}particles;

/**
 * This function places particles in a cubic box, some points are generic to all the particles
 * 1. The neighbouring particles lies in a spherical surface of radius r_0
 * 2. If the particle coincided with the excluded volume of any other particle, it redifines its coordinates
 * 3. The particle cannot go out of the box, it redefines its coordinates again
 * 4. The polymer chain always lies in a reduced cubic box than the actual confinement, just to reduce the occurance of (3)
 * 
*/
void linearChain(particles *part, float xi, float xf, float yi, float yf, float zi, float zf, float r_0, int chain_N, int part_id_i, float sigma1 , float sigma2, int stck_flag, long *seed)
{
    float dx,dy,dz,dr,sigma;
    // Reduce the end points by 2*r_0
    {
        xf -= r_0;
        yf -= r_0;
        zf -= r_0;
        xi += r_0;
        yi += r_0;
        zi += r_0;
    }
    float lx = fabsf(xf - xi);
    float ly = fabsf(yf - yi);
    float lz = fabsf(zf - zi);
    float rand_theta, rand_u,x,y,z;
    int error = 0;
    for (int ni = part_id_i; ni < (part_id_i + chain_N); ni++)
    {
        // first bead in the chain
        if (ni % chain_N == 0)
        {
            x = ran2(seed) * (xf - xi) + xi;
            y = ran2(seed) * (yf - yi) + yi;
            z = ran2(seed) * (zf - zi) + zi;
            if (stck_flag == 0)
                sigma = sigma1;
            else
                sigma = sigma2;
            
        }
        // all other beads in the chain
        else
        {
            rerun :

            rand_u = 2. * ran2(seed) - 1.;
            rand_theta = twoPI * ran2(seed);
            if (stck_flag == 0 && ni == part_id_i + 1)
            {
                x = part[ni - 1].px + sqrtf(1. - rand_u * rand_u) * cosf(rand_theta) * 0.5;
                y = part[ni - 1].py + sqrtf(1. - rand_u * rand_u) * sinf(rand_theta) * 0.5;
                z = part[ni - 1].pz + rand_u * 0.5;
                sigma = sigma2;
            }
            else if(stck_flag == 1 && ni == (part_id_i + chain_N - 1))
            {
                x = part[ni - 1].px + sqrtf(1. - rand_u * rand_u) * cosf(rand_theta) * 0.5;
                y = part[ni - 1].py + sqrtf(1. - rand_u * rand_u) * sinf(rand_theta) * 0.5;
                z = part[ni - 1].pz + rand_u * 0.5;
                sigma = sigma1;
            }
            else
            {
                x = part[ni - 1].px + sqrtf(1. - rand_u * rand_u) * cosf(rand_theta) * r_0;
                y = part[ni - 1].py + sqrtf(1. - rand_u * rand_u) * sinf(rand_theta) * r_0;
                z = part[ni - 1].pz + rand_u * r_0;
                sigma = sigma2;
            }
            
            // check for the particles out of the box
            x = ((x > xf) ? (2. * xf - x) : ((x < xi) ? (2. * xi - x) : x + 0.));
            y = ((y > yf) ? (2. * yf - y) : ((y < yi) ? (2. * yi - y) : y + 0.));
            z = ((z > zf) ? (2. * zf - z) : ((z < zi) ? (2. * zi - z) : z + 0.));
        
            for (int j = part_id_i; j < (part_id_i + chain_N); j++)
            {
                if ( j == ni )
                    continue;
                
                dx = x - part[j].px;
                dy = y - part[j].py;
                dz = z - part[j].pz;

                dr = sqrtf(dx * dx + dy * dy + dz * dz);
                if(dr < 0.5 * (part[ni].sigma + part[j].sigma))
                {
                    goto rerun;
                }
            }
            
        }
        part[ni].sigma = sigma;
        part[ni].pid = ni;
        part[ni].px = x;
        part[ni].py = y;
        part[ni].pz = z;
    }
    
}
/**
void starPoly(particles *part, float xi, float xf, float yi, float yf, float zi, float zf, float r_0, int N_m, int f, int part_id_i, float sigma1 , float sigma2, int stck_flag, long *seed)
{
    float dx, dy, dz, dr, sigma, rand_theta, rand_u, x, y, z, lx, ly, lz;
    int n_m;
    // Reduce the end points by 2*r_0
    // {
    //     xf -= r_0;
    //     yf -= r_0;
    //     zf -= r_0;
    //     xi += r_0;
    //     yi += r_0;
    //     zi += r_0;
    // }
    lx = fabsf(xf - xi);
    ly = fabsf(yf - yi);
    lz = fabsf(zf - zi);
    n_m = (int)(N_m - 1) / f;
    int error = 0;
    int arms = 1, ni = part_id_i + 1;

    // Place the first particle in the centre of the box
    part[part_id_i].pid = part_id_i;
    part[part_id_i].px = 0.5 * (xi + xf);
    part[part_id_i].py = 0.5 * (yi + yf);
    part[part_id_i].pz = 0.5 * (zi + zf);
    part[part_id_i].sigma = sigma2;

    while (arms <= f && ni < (part_id_i + N_m))
    {
        for (int i = ni; i < (ni + n_m); i++)
        {
            rerun :

            rand_u = 2. * ran2(seed) - 1.;
            rand_theta = ran2(seed) * twoPI;
            if(stck_flag == 1 && i == part_id_i + (arms * n_m))
            {
                x = part[i - 1].px + sqrtf(1. - rand_u * rand_u) * cosf(rand_theta) * 0.50;
                y = part[i - 1].py + sqrtf(1. - rand_u * rand_u) * sinf(rand_theta) * 0.50;
                z = part[i - 1].pz + rand_u * 0.50;
                sigma = sigma1;
            }
            else if (i == ni)
            {
                //rand_theta = arms * PIby4;
                x = part[part_id_i].px + sqrtf(1. - rand_u * rand_u) * cosf(rand_theta) * r_0;
                y = part[part_id_i].py + sqrtf(1. - rand_u * rand_u) * sinf(rand_theta) * r_0;
                z = part[part_id_i].pz + rand_u * r_0;
                sigma = sigma2;
            }
            else
            {
                //rand_theta = arms * PIby4
                x = part[i - 1].px + sqrtf(1. - rand_u * rand_u) * cosf(rand_theta) * r_0;
                y = part[i - 1].py + sqrtf(1. - rand_u * rand_u) * sinf(rand_theta) * r_0;
                z = part[i - 1].pz + rand_u * r_0;
                sigma = sigma2;
            }
            
            // check for the particles out of the box
            x = ((x > xf) ? (2. * xf - x) : ((x < xi) ? (2. * xi - x) : x + 0.));
            y = ((y > yf) ? (2. * yf - y) : ((y < yi) ? (2. * yi - y) : y + 0.));
            z = ((z > zf) ? (2. * zf - z) : ((z < zi) ? (2. * zi - z) : z + 0.));
        
            for (int j = part_id_i; j < i-1; j++)
            {
                dx = x - part[j].px;
                dy = y - part[j].py;
                dz = z - part[j].pz;

                dr = sqrtf(dx * dx + dy * dy + dz * dz);
                if(dr < 0.5 * (part[j].sigma + sigma))
                {
                    goto rerun;
                }
            }
            part[i].sigma = sigma;
            part[i].pid = i;
            part[i].px = x;
            part[i].py = y;
            part[i].pz = z;
        }
        
        arms += 1;
        ni += n_m;
 
    }
    
}
*/

void starPoly(particles *part, float xi, float xf, float yi, float yf, float zi, float zf, float r_0, int N_m, int f, int part_id_i, float sigma1 , float sigma2, int stck_flag, long *seed)
{
    float dx, dy, dz, dr, sigma, rand_theta, rand_phi, rand_u, rand_v, rand_r, x, y, z, lx, ly, lz, X, Y, Z;
    int n_m;
    // Reduce the end points by 2*r_0
    // {
    //     xf -= r_0;
    //     yf -= r_0;
    //     zf -= r_0;
    //     xi += r_0;
    //     yi += r_0;
    //     zi += r_0;
    // }
    lx = fabsf(xf - xi);
    ly = fabsf(yf - yi);
    lz = fabsf(zf - zi);
    n_m = (int)(N_m - 1) / f;
    int error = 0;
    int arms = 1, ni = part_id_i + 1;

    // Place the first particle in the centre of the box
    part[part_id_i].pid = part_id_i;
    part[part_id_i].px = 0.5 * (xi + xf);
    part[part_id_i].py = 0.5 * (yi + yf);
    part[part_id_i].pz = 0.5 * (zi + zf);
    part[part_id_i].sigma = sigma2;

    while (arms <= f && ni < (part_id_i + N_m))
    {
        for (int i = ni; i < (ni + n_m); i++)
        {
            rerun :

            rand_u = ran2(seed);
            rand_v = ran2(seed);
            rand_phi = acosf(2. * rand_v - 1.);
            //rand_r = cbrtf(ran2(seed));
            rand_theta = (rand_u * (arms-1) + 0.5) * PIby2;
            X = sinf(rand_phi) * cosf(rand_theta);
            Y = sinf(rand_phi) * sinf(rand_theta);
            Z = cosf(rand_phi);
            if(stck_flag == 1 && i == part_id_i + (arms * n_m))
            {
                x = part[i - 1].px + X * 0.50;
                y = part[i - 1].py + Y * 0.50;
                z = part[i - 1].pz + Z * 0.50;
                sigma = sigma1;
            }
            else if (i == ni)
            {
                //rand_theta = arms * PIby4;
                x = part[part_id_i].px + X * r_0;
                y = part[part_id_i].py + Y * r_0;
                z = part[part_id_i].pz + Z * r_0;
                sigma = sigma2;
            }
            else
            {
                //rand_theta = arms * PIby4
                x = part[i - 1].px + X * r_0;
                y = part[i - 1].py + Y * r_0;
                z = part[i - 1].pz + Z * r_0;
                sigma = sigma2;
            }
            
            // check for the particles out of the box
            x = ((x > xf) ? (2. * xf - x) : ((x < xi) ? (2. * xi - x) : x + 0.));
            y = ((y > yf) ? (2. * yf - y) : ((y < yi) ? (2. * yi - y) : y + 0.));
            z = ((z > zf) ? (2. * zf - z) : ((z < zi) ? (2. * zi - z) : z + 0.));
        
            for (int j = part_id_i; j < i-1; j++)
            {
                dx = x - part[j].px;
                dy = y - part[j].py;
                dz = z - part[j].pz;

                dr = sqrtf(dx * dx + dy * dy + dz * dz);
                if(dr < 0.5 * (part[j].sigma + sigma))
                {
                    goto rerun;
                }
            }
            part[i].sigma = sigma;
            part[i].pid = i;
            part[i].px = x;
            part[i].py = y;
            part[i].pz = z;
        }
        
        arms += 1;
        ni += n_m;
 
    }
}
int main(int argc, char* argv[])
{
    // int total_N, chain_N, num_chain, nindex = 0, sticker_flag;
    // float Lx,Ly,Lz,r_0,div_per_dim,sigma1,sigma2;
    // long seed;
    // char filename1[1024];

    int total_N, N_m, N_star, n_m, nindex = 0, sticker_flag, f;
    float Lx,Ly,Lz,r_0,div_per_dim,sigma1,sigma2;
    long seed;
    char filename1[1024];

    /**
     * Take the inputs from the command line 
    */
    /**
     * total_N = The total number of particles in the simulation box
     * num_chain = the number of chains 
     * Lx,Ly,Lz = actual dimension of the simu box
     * r_0 = mean distance between particles
     * sigma = parameter from LJ force field(diameter of the particles)
     * chain_N = num. of particles in each chain
     * sticker_flag = 0: the sticker is only at one of the ends
     *                1: sticker present at both the ends
     *                2: sticker absent
     * star_N = num
     * 
    */
    // total_N = atoi(argv[1]); 
    // num_chain = atoi(argv[2]);
    // Lx = atof(argv[3]) * 0.1;
    // Ly = atof(argv[4]) * 0.1;
    // Lz = atof(argv[5]) * 0.1;
    // r_0 = atof(argv[6]) * 0.1;
    // sigma1 = atof(argv[7]) * 0.1;
    // sigma2 = atof(argv[8]) * 0.1;
    // sticker_flag = atoi(argv[9]);
    // seed = -1 * atol(argv[10]);
    // div_per_dim = ceilf(powf(num_chain,0.3333));
    // chain_N = total_N / num_chain;

    //total_N = atoi(argv[1]); 
    n_m = atoi(argv[1]);
    f = atoi(argv[2]);
    N_star = atoi(argv[3]);
    Lx = atof(argv[4]) * 0.1;
    Ly = atof(argv[5]) * 0.1;
    Lz = atof(argv[6]) * 0.1;
    r_0 = atof(argv[7]) * 0.01;
    sigma1 = atof(argv[8]) * 0.01;
    sigma2 = atof(argv[9]) * 0.01;
    sticker_flag = atoi(argv[10]);
    seed = -1 * atol(argv[11]);
    div_per_dim = ceilf(powf(N_star,0.3333));
    N_m = f * n_m + 1;

    total_N = N_m * N_star;
    //sigma1 += 0.1;
    //sigma2 += 0.1;
    
    //total_N = (sticker_flag == 0) ? (total_N + 1 * num_chain) : ((sticker_flag == 1) ? (total_N + 2 * num_chain): total_N + 0); // adjusting the number of particles with 
                                                                                                                                // with no. of sticky beads

    float random = ran2(&seed);


    particles *part1 = (particles*)malloc(total_N*sizeof(particles));

    nindex = 0;
    while (nindex < total_N)
    {
        for (int z = 0; z < (int)div_per_dim; z++)
        {
            for (int y = 0; y < (int)div_per_dim; y++)
            {
                for (int x = 0; x < (int)div_per_dim; x++)
                {
                    if(nindex >= total_N)
                        goto out;
                    //linearChain(part1,(-Lx / 2.) + x * (Lx / div_per_dim),(-Lx / 2.) + (x + 1) * (Lx / div_per_dim), \
                                     (-Ly / 2.) + y * (Ly / div_per_dim),(-Ly / 2.) + (y + 1) * (Ly / div_per_dim),  \
                                     (-Lz / 2.) + z * (Lz / div_per_dim),(-Lz / 2.) + (z + 1) * (Lz / div_per_dim), \
                                     r_0, chain_N, nindex, sigma1, sigma2, sticker_flag, &seed);
                    starPoly(part1,(-Lx / 2.) + x * (Lx / div_per_dim),(-Lx / 2.) + (x + 1) * (Lx / div_per_dim), \
                                     (-Ly / 2.) + y * (Ly / div_per_dim),(-Ly / 2.) + (y + 1) * (Ly / div_per_dim),  \
                                     (-Lz / 2.) + z * (Lz / div_per_dim),(-Lz / 2.) + (z + 1) * (Lz / div_per_dim), \
                                     r_0, N_m, f, nindex, sigma1, sigma2, sticker_flag, &seed);
                    
                    nindex += N_m;
                }
                
            }
            
        }
        
    }
    out :
    sprintf(filename1,"Init_pos_N_%d_C_%d.dat",total_N,N_star);
    FILE* posfile = fopen(filename1,"w");
    for (int i = 0; i < total_N; i++)
    {
        fprintf(posfile,"%f\t%f\t%f\n",part1[i].px,part1[i].py,part1[i].pz);
    }
    fclose(posfile);
    free(part1);
    return 0;

}

