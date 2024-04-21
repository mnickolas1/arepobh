#include <gsl/gsl_math.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../main/allvars.h"
#include "../main/proto.h"


/*THIS PART ADAPTED FROM GADGET4*/
/* fall back to cubic spline kernel */
#if !defined(CUBIC_SPLINE_KERNEL) && !defined(WENDLAND_C2_KERNEL) && !defined(WENDLAND_C4_KERNEL) && !defined(WENDLAND_C6_KERNEL)
#define CUBIC_SPLINE_KERNEL 
#endif

/* fall back to three dimensions */
#if !defined(TWODIMS) && !defined(ONEDIMS)
#define THREEDIMS
#endif

/* Norms */
#ifdef CUBIC_SPLINE_KERNEL

#ifdef THREEDIMS
#define NORM (8.0 / M_PI) /*!< For 3D-normalized kernel */
#endif

#ifdef TWODIMS
#define NORM (40.0 / (7.0 * M_PI)) /*!< For 2D-normalized kernel */
#endif

#ifdef ONEDIMS
#define NORM (4.0 / 3.0) /*!< For 1D-normalized kernel */
#endif

#endif /* CUBIC_SPLINE_KERNEL */

#ifdef WENDLAND_C2_KERNEL

#ifdef THREEDIMS
#define NORM (21.0 / (2.0 * M_PI)) /*!< For 3D-normalized kernel */
#endif

#ifdef TWODIMS
#define NORM (7.0 / M_PI) /*!< For 2D-normalized kernel */
#endif

#ifdef ONEDIMS
#define NORM (5.0 / 4.0) /*!< For 1D-normalized kernel */
#endif

#endif /* WENDLAND_C2_KERNEL */

#ifdef WENDLAND_C4_KERNEL

#ifdef THREEDIMS
#define NORM (495.0 / (32.0 * M_PI)) /*!< For 3D-normalized kernel */
#endif

#ifdef TWODIMS
#define NORM (9.0 / M_PI) /*!< For 2D-normalized kernel */
#endif

#ifdef ONEDIMS
#define NORM (3.0 / 2.0) /*!< For 1D-normalized kernel */
#endif

#endif /* WENDLAND_C4_KERNEL */

#ifdef WENDLAND_C6_KERNEL

#ifdef THREEDIMS
#define NORM (1365.0 / (64.0 * M_PI)) /*!< For 3D-normalized kernel */
#endif

#ifdef TWODIMS
#define NORM (78.0 / (7.0 * M_PI)) /*!< For 2D-normalized kernel */
#endif

#ifdef ONEDIMS
#define NORM (55.0 / 32.0) /*!< For 1D-normalized kernel */
#endif

#endif /* WENDLAND_C6_KERNEL */

static int int_compare(const void *a, const void *b);

/*sph loop kernel function -> u < 1 */
void kernel(double u, double hinv3, double hinv4, double *wk, double *dwk)
{
#ifdef CUBIC_SPLINE_KERNEL
#if defined(WENDLAND_C2_KERNEL) || defined(WENDLAND_C4_KERNEL) || defined(WENDLAND_C6_KERNEL)
#error "Only one SPH kernel can be used"
#endif
  if(u < 0.5)
    {
      *dwk = u * (18.0 * u - 12.0);
      
      *wk = (1.0 + 6.0 * (u - 1.0) * u * u);
    }
  else
    {
      double t1 = (1.0 - u);
      double t2 = t1 * t1;
      
      *dwk = -6.0 * t2;
      
      *wk = 2.0 * t2 * t1;
    }
#endif

#ifdef WENDLAND_C2_KERNEL /* Dehnen & Aly 2012 */
#ifdef ONEDIMS
  double t1 = (1.0 - u);
  double t2 = (t1 * t1);

  
  *dwk = -12.0 * u * t2;
  
  *wk = t2 * t1 * (1.0 + u * 3.0);

#else /* 2d or 3d */
  double t1 = (1.0 - u);
  double t2 = (t1 * t1);
  double t4 = t2 * t2;
  
  *dwk = -20.0 * u * t2 * t1;
  
  *wk = t4 * (1.0 + u * 4.0);

#endif
#endif /* WENDLAND_C2_KERNEL */

#ifdef WENDLAND_C4_KERNEL /* Dehnen & Aly 2012 */
#ifdef ONEDIMS
  double t1 = (1.0 - u);
  double t2 = t1 * t1;
  double t4 = t2 * t2;
  double t5 = t4 * t1;

  
  *dwk = -14.0 * t4 * (4.0 * u + 1) * u;
  
  *wk = t5 * (1.0 + u * (5.0 + 8.0 * u));

#else /* 2d or 3d */
  double t1 = (1.0 - u);
  double t2 = (t1 * t1);
  double t4 = t2 * t2;
  double t6 = t2 * t2 * t2;
  
  *dwk = -56.0 / 3.0 * u * t4 * t1 * (5.0 * u + 1);
  
  *wk = t6 * (1.0 + u * (6.0 + 35.0 / 3.0 * u));

#endif
#endif /* WENDLAND_C4_KERNEL */

#ifdef WENDLAND_C6_KERNEL /* Dehnen & Aly 2012 */
#ifdef ONEDIMS
  double t1 = (1.0 - u);
  double t2 = (t1 * t1);
  double t4 = t2 * t2;
  double t6 = t4 * t2;
  double t7 = t4 * t2 * t1;
  
  *dwk = -6.0 * u * t6 * (3.0 + u * (18.0 + 35.0 * u));
  
  *wk = t7 * (1.0 + u * (7.0 + u * (19.0 + 21.0 * u)));

#else /* 2d or 3d */
  double t1 = (1.0 - u);
  double t2 = (t1 * t1);
  double t4 = t2 * t2;
  double t7 = t4 * t2 * t1;
  double t8 = t4 * t4;
  
  *dwk = -22.0 * u * (1.0 + u * (7.0 + 16.0 * u)) * t7;
  
  *wk = t8 * (1.0 + u * (8.0 + u * (25.0 + 32.0 * u)));

#endif
#endif /* WENDLAND_C6_KERNEL */
  
  *dwk *= NORM * hinv4;
  
  *wk *= NORM * hinv3;
}
/*THIS PART ADAPTED FROM GADGET4*/

void active_virtual_part_init_alloc(struct ActiveVirtualPart *AVP, const char *name, int *MaxPart)
{
  AVP->NActiveParticles   = 0;
  AVP->ActiveParticleList = 0;
  AVP->MaxPart  = MaxPart;

  char Identifier[200];
  Identifier[199] = 0;

  snprintf(Identifier, 199, "NextActiveParticle%s", AVP->Name);
  AVP->ActiveParticleList = (int *)mymalloc_movable(&AVP->ActiveParticleList, Identifier, *(AVP->MaxPart) * sizeof(int));
}

void active_virtual_part_set(struct ActiveVirtualPart *AVP)
{
  AVP->NActiveParticles = 0;
  
  for(int i=0; i<NumBh; i++)
    if(BhP[i].DestroyFlag == 1)
      {
        AVP->ActiveParticleList[AVP->NActiveParticles] = i;
        AVP->NActiveParticles++;
      }
}

void virtual_part_feedback(void)
{
  int idx, i;
  double Ftherm;

  Ftherm  = pow(10,7)*BOLTZMANN / GAMMA_MINUS1 / PROTONMASS / 0.6; //add thermal energy ~ 10^7 K
  Ftherm /= (All.UnitEnergy_in_cgs / All.UnitMass_in_g);  

  struct pv_update_data pvd;
  if(All.ComovingIntegrationOn)
    {
      pvd.atime    = All.Time;
      pvd.hubble_a = hubble_function(All.Time);
      pvd.a3inv    = 1 / (All.Time * All.Time * All.Time);
    }
  else
    pvd.atime = pvd.hubble_a = pvd.a3inv = 1.0;

  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
      continue;
      if(SphP[i].F < 0)
      continue;
  
      SphP[i].Momentum[0] += SphP[i].FMomentum[0];
      SphP[i].Momentum[1] += SphP[i].FMomentum[1];
      SphP[i].Momentum[2] += SphP[i].FMomentum[2];
      P[i].Mass           += SphP[i].FMass;
      
      /*update velocities*/
      update_primitive_variables_single(P, SphP, i, &pvd);
      /*update total energy*/
      SphP[i].Energy = (SphP[i].Utherm + Ftherm) * P[i].Mass + 0.5 * P[i].Mass * (pow(P[i].Vel[0], 2) + pow(P[i].Vel[1], 2) + pow(P[i].Vel[2], 2)); 
      /*update internal energy*/
      update_internal_energy(P, SphP, i, &pvd);
      /*update pressure*/
      set_pressure_of_cell_internal(P, SphP, i);

#ifdef PASSIVE_SCALARS                 
      /*tracer field advected passively*/
      SphP[i].PScalars[0] = 1;
      SphP[i].PConservedScalars[0] = P[i].Mass;
#endif
      /*set feed flags to zero*/
      SphP[i].FMomentum[0] = SphP[i].FMomentum[1] = SphP[i].FMomentum[2] = SphP[i].FMass = 0;
      SphP[i].F = -1;
    }
}
