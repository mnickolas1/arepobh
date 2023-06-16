#include <gsl/gsl_math.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../main/allvars.h"
#include "../main/proto.h"

void perform_end_of_step_bh_physics(void)
{
  int i;
  double pj;
  double kick_vector[3];

/*find cone particles to kick*/
    if(All.Time >= All.FeedbackTime)
      {   
        if(All.FeedbackFlag > 0)
          {
            int queue = 0;
            int queue_all = 0;
            for(i = 0; i < NumGas; i++)
              {
                if(SphP[i].JetQueue > queue) //JetQueue gives the priority list for cone particles
                  queue = SphP[i].JetQueue;
              }
            MPI_Allreduce(&queue, &queue_all, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
            MPI_Barrier(MPI_COMM_WORLD); // synchronize all tasks
              
            struct pv_update_data pvd;
            if(All.ComovingIntegrationOn)
              {
                pvd.atime    = All.Time;
                pvd.hubble_a = hubble_function(All.Time);
                pvd.a3inv    = 1 / (All.Time * All.Time * All.Time);
              }
            else
              pvd.atime = pvd.hubble_a = pvd.a3inv = 1.0;  
/*kick the particles*/            
            for(i = 0; i < NumGas; i++)
              {
                if(SphP[i].JetQueue == queue_all)
                  {  
                    kick_vector[0] = SphP[i].BhKickVector[0];
                    kick_vector[1] = SphP[i].BhKickVector[1];
                    kick_vector[2] = SphP[i].BhKickVector[2];

                    pj = P[i].Mass * All.VJet; 

                    /*update momentum*/
                    SphP[i].Momentum[0] = kick_vector[0] * pj / sqrt(pow(kick_vector[0], 2) + pow(kick_vector[1], 2) + pow(kick_vector[2], 2));
                    SphP[i].Momentum[1] = kick_vector[1] * pj / sqrt(pow(kick_vector[0], 2) + pow(kick_vector[1], 2) + pow(kick_vector[2], 2));
                    SphP[i].Momentum[2] = kick_vector[2] * pj / sqrt(pow(kick_vector[0], 2) + pow(kick_vector[1], 2) + pow(kick_vector[2], 2));   
                 
                    /*update velocities*/
                    update_primitive_variables_single(P, SphP, i, &pvd);  

                    /*update total energy*/
                    SphP[i].Energy = SphP[i].Utherm * P[i].Mass + 0.5 * P[i].Mass * (pow(P[i].Vel[0], 2) + pow(P[i].Vel[1], 2) + pow(P[i].Vel[2], 2));                 
                    /*update internal energy*/
                    update_internal_energy(P, SphP, i, &pvd);
                    /*update pressure*/
                    set_pressure_of_cell_internal(P, SphP, i);
#ifdef PASSIVE_SCALARS                 
                    /*tracer field advected passively*/
                    SphP[i].PScalars[0] = 1;
                    SphP[i].PConservedScalars[0] = P[i].Mass;
                  }
              }     
#endif
#ifdef BURST_MODE
        All.FeedbackFlag = -1;
        All.LastFeedbackTime = All.Time;
#endif
          }
      }
#ifdef BURST_MODE
  if(All.PJet * (All.Time - All.LastFeedbackTime) >= All.MJet * All.VJet*All.VJet)  
    All.FeedbackFlag = 1;
#endif   
}
