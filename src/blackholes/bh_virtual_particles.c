#include <gsl/gsl_math.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#include "../main/allvars.h"
#include "../main/proto.h"

#define DEG_TO_RAD(deg) ((deg) * M_PI / 180.0)

// Robert Jenkins' 96 bit Mix Function
unsigned long mix(unsigned long a, unsigned long b, unsigned long c)
{
    a=a-b;  a=a-c;  a=a^(c >> 13);
    b=b-c;  b=b-a;  b=b^(a << 8);
    c=c-a;  c=c-b;  c=c^(b >> 13);
    a=a-b;  a=a-c;  a=a^(c >> 12);
    b=b-c;  b=b-a;  b=b^(a << 16);
    c=c-a;  c=c-b;  c=c^(b >> 5);
    a=a-b;  a=a-c;  a=a^(c >> 3);
    b=b-c;  b=b-a;  b=b^(a << 10);
    c=c-a;  c=c-b;  c=c^(b >> 15);
    return c;
}

void create_particles(void)
{
  double Pj, Mj, Vj; 
  
  Pj = All.PJet / (All.UnitEnergy_in_cgs / All.UnitTime_in_s);
  Mj = All.MJet; /// (All.UnitMass_in_g);
  Vj = All.VJet; /// (All.UnitVelocity_in_cm_per_s);

#ifdef BURST_MODE
  if((Pj * (All.Time) >= All.FeedbackCount * Mj * Vj*Vj) && (All.FeedbackFlag < 0))  
    All.FeedbackFlag = 1;
#endif 
  
  if(All.FeedbackFlag > 0)
    {
      int i;
      double theta, phi, x, y, z;
      
      int particles_spawned = 0;
      int tot_particles_spawned = 2;
      if(ThisTask == 0)
        particles_spawned = 2;
      mpi_printf("\nJETS: Kicking Particles -> FeedbackCount:%d\n\n", All.FeedbackCount);      
      
      int *list;

      if(All.MaxID == 0) 
        calculate_maxid();

      list = mymalloc("list", NTask * sizeof(int));

      MPI_Allgather(&particles_spawned, 1, MPI_INT, list, 1, MPI_INT, MPI_COMM_WORLD);

      MyIDType newid = All.MaxID + 1;

      for(i = 0; i < ThisTask; i++)
        newid += list[i];

      myfree(list);

      if(particles_spawned > 0)
        {
          unsigned long seed = mix(clock(), time(NULL), getpid());
      
          srand(seed);

          // Generate random angles for phi (azimuthal angle) and theta (polar angle)
          phi = DEG_TO_RAD(rand() % 360);
          if(rand() % 2 == 0)
            theta = DEG_TO_RAD(rand() % 21 - 10);
          else
            theta = DEG_TO_RAD(rand() % 21 + 170);
          // Calculate x, y, and z components
          x = cos(phi) * sin(theta);
          y = sin(phi) * sin(theta);
          z = cos(theta);
        }

      for(i = 0; i < particles_spawned; i++)
        {
/* Assign properties to the spawned particles */

/* Assign new unique IDs to the spawned particles */
          P[NumPart + i].ID = newid;

          newid++;

          //assign mass
          P[NumPart + i].Mass = Mj;
          //assign pos
          P[NumPart + i].Pos[0] = 150;
          P[NumPart + i].Pos[1] = 150;
          P[NumPart + i].Pos[2] = 300;
          //assign vel 
          if(i%2 == 0)
            {
              P[NumPart + i].Vel[0] = x*Vj;
              P[NumPart + i].Vel[1] = y*Vj;
              P[NumPart + i].Vel[2] = z*Vj;
            }
          if(i%2 == 1)
            {
              P[NumPart + i].Vel[0] = -x*Vj;
              P[NumPart + i].Vel[1] = -y*Vj;
              P[NumPart + i].Vel[2] = -z*Vj;
            }
          //assign acc
          P[NumPart + i].GravAccel[0] = 0;
          P[NumPart + i].GravAccel[1] = 0;
          P[NumPart + i].GravAccel[2] = 0;
          //assign type
          P[NumPart + i].Type = 5;
          //assign ti_current
          P[NumPart + i].Ti_Current = All.Ti_Current;
          //assign bh_ids
          P[NumPart + i].BhID = NumBh + i;
          BhP[NumBh + i].PID  = NumPart + i;
          //assing density loop properties
          BhP[NumBh + i].Hsml = 5;
          BhP[NumBh + i].DestroyFlag = -1;
          BhP[NumBh + i].DensityFlag = 1;
        }

      All.MaxID      += tot_particles_spawned;
      All.TotNumPart += tot_particles_spawned;
      All.TotNumBh   += tot_particles_spawned;
      NumPart        += particles_spawned;
      NumBh          += particles_spawned;

#ifdef BURST_MODE
        All.FeedbackFlag = -1;
        All.FeedbackCount++;
#endif
    }  
}

void destroy_particles(void)
{
 int i;
 double x, y, z, r2;
 double Radius = 10;

 for(i=0; i<NumBh; i++)
   { 
     if(BhP[i].DestroyFlag < 0)
       {
         x = PPB(i).Pos[0] - 150;
         y = PPB(i).Pos[1] - 150;
         z = PPB(i).Pos[2] - 300;

         r2 = x*x + y*y + z*z;

         if(r2 > Radius*Radius)
           {
             BhP[i].DestroyFlag = 1;
             printf("\nJETS: Destroying Particles\n\n");
           }
       }
     else
       BhP[i].DestroyFlag = 2;
   }
}
