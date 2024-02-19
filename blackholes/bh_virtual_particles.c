#include <gsl/gsl_math.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../main/allvars.h"
#include "../main/proto.h"
#include <time.h>

#define DEG_TO_RAD(deg) ((deg) * M_PI / 180.0)


void create_particles(void)
{
  double Pj, Mj, Vj; 
  
  Pj = All.PJet / (All.UnitEnergy_in_cgs / All.UnitTime_in_s);
  Mj = All.MJet; /// (All.UnitMass_in_g);
  Vj = All.VJet; /// (All.UnitVelocity_in_cm_per_s);

#ifdef BURST_MODE
  if((Pj * (All.Time - All.LastFeedbackTime) >= Mj * Vj*Vj) && (All.FeedbackFlag < 0))  
    All.FeedbackFlag = 1;
#endif 
  
  if(All.FeedbackFlag > 0)
    {

      int i;
      int particles_spawned = 0;
      int tot_particles_spawned = 2;
      if(ThisTask < tot_particles_spawned)
        particles_spawned = 1;
      
      int *list;

      if(All.MaxID == 0) 
        calculate_maxid();

      list = mymalloc("list", NTask * sizeof(int));

      MPI_Allgather(&particles_spawned, 1, MPI_INT, list, 1, MPI_INT, MPI_COMM_WORLD);

      MyIDType newid = All.MaxID + 1;

      for(i = 0; i < ThisTask; i++)
        newid += list[i];

      myfree(list);

      for(i = 0; i < particles_spawned; i++)
        {
/* Assign new unique IDs to the spawned particles */
          P[NumPart + i].ID = newid;

          newid++;

/* Assign properties to the spawned particles */       
          srand(time(NULL));

          double phi, theta;

          // Generate random angles for phi (azimuthal angle) and theta (polar angle)
          phi = DEG_TO_RAD(rand() % 360);
          if(rand() % 2 == 0)
            theta = DEG_TO_RAD(rand() % 21 - 10);
          else
            theta = DEG_TO_RAD(rand() % 21 + 170);
          // Calculate x, y, and z components
          double x = cos(phi) * sin(theta);
          double y = sin(phi) * sin(theta);
          double z = cos(theta);

          //assign mass
          P[NumPart + i].Mass = Mj;
          //assign pos
          P[NumPart + i].Pos[0] = 450;
          P[NumPart + i].Pos[1] = 150;
          P[NumPart + i].Pos[2] = 150;
          //assign vel 
          if(ThisTask == 0)
            {
              P[NumPart + i].Vel[0] = x*Vj;
              P[NumPart + i].Vel[1] = y*Vj;
              P[NumPart + i].Vel[2] = z*Vj;
            }
          else
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
        All.LastFeedbackTime = All.Time;
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
         z = PPB(i).Pos[2] - 150;

         r2 = x*x + y*y + z*z;

         if(r2 > Radius*Radius && BhP[i].Hsml < 5)
           BhP[i].DestroyFlag = 1;
       }
     else
       BhP[i].DestroyFlag = 2;
   }
}
