#include <gsl/gsl_math.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../main/allvars.h"
#include "../main/proto.h"


void create_particles(void)
{
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
          P[NumPart + i].Mass = 0.01;
          P[NumPart + i].Pos[0] = 0.1;
          P[NumPart + i].Pos[1] = 0;
          P[NumPart + i].Pos[2] = 0;
          P[NumPart + i].Vel[0] = 1;
          P[NumPart + i].Vel[1] = 0;
          P[NumPart + i].Vel[2] = 0;
          P[NumPart + i].GravAccel[0] = 0;
          P[NumPart + i].GravAccel[1] = 0;
          P[NumPart + i].GravAccel[2] = 0;

          P[NumPart + i].Type = 5;
          
          P[NumPart + i].BhID = NumBh + i;
          BhP[NumBh + i].PID  = NumPart + i;
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

#ifdef BURST_MODE
  if(All.PJet * (All.Time - All.LastFeedbackTime) >= All.MJet * All.VJet*All.VJet)  
    All.FeedbackFlag = 1;
#endif   
}

void destroy_particles(void)
{
 int i;
 int Radius = 10;

 for(i=0; i<NumBh; i++)
   { 
     if(BhP[i].DestroyFlag < 0)
       {
         if(PPB(i).Pos[0] > Radius)
           BhP[i].DestroyFlag = 1;
       }
     
     else
       BhP[i].DestroyFlag = 2;
   }
}