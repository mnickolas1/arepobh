#include <gsl/gsl_math.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../main/allvars.h"
#include "../main/proto.h"

#define DEG_TO_RAD(deg) ((deg) * M_PI / 180.0)

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
      double x,y,z;

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
          double phi = ((double)rand() / RAND_MAX) * 2 * M_PI;

          // Generate a random value for cosine(theta) to ensure uniform distribution in the cone
          double cos_theta = ((double)rand() / RAND_MAX) * (cos(DEG_TO_RAD(15)) - 1) + 1;
          double theta = acos(cos_theta); 

          // Convert spherical coordinates to Cartesian coordinates
          x = sin(theta) * cos(phi); 
          y = sin(theta) * sin(phi);  
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
          P[NumPart + i].Pos[0] = boxHalf_X;
          P[NumPart + i].Pos[1] = boxHalf_Y;
          P[NumPart + i].Pos[2] = boxHalf_Z;
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
          BhP[NumBh + i].Hsml = 10;
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

 int local_destroyed = 0, global_destroyed = 0;

 for(i=0; i<NumBh; i++)
   { 
     if(BhP[i].DestroyFlag < 0)
       {
         x = PPB(i).Pos[0] - boxHalf_X;
         y = PPB(i).Pos[1] - boxHalf_Y;
         z = PPB(i).Pos[2] - boxHalf_Z;

         r2 = x*x + y*y + z*z;

         if(r2 > Radius*Radius)
           {
             BhP[i].DestroyFlag = 1;
             local_destroyed++;
           }
       }
     else
       BhP[i].DestroyFlag = 2;
   }

  MPI_Reduce(&local_destroyed, &global_destroyed, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  
  if(global_destroyed > 0) 
    mpi_printf("\nJETS: Destroying %d particles across all tasks\n\n", global_destroyed);
}
