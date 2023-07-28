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

/*static int int_compare(const void *a, const void *b);*/

/*update bh-timestep at prior_mesh_construction*/
/*void update_bh_timesteps(void)
{
  int i;

  for(i = 0; i < NumBh; i++)
    { 
      if(BhP[i].DestroyFlag == 1)
        BhP[i].TimeBinBh = 0;
      else
        BhP[i].TimeBinBh = 29;
    }
  reconstruct_bh_timebins();
  update_list_of_active_bh_particles();
}*/

/*call this function as the reconstruct_timebins() bh version*/
/*void reconstruct_bh_timebins(void)
{
  int i, bin;

  for(bin = 0; bin < TIMEBINS; bin++)
    {
      TimeBinsBh.TimeBinCount[bin]   = 0;
      TimeBinsBh.FirstInTimeBin[bin] = -1;
      TimeBinsBh.LastInTimeBin[bin]  = -1;
    }
  
  for(i = 0; i < NumBh; i++)
    {
      
      bin = BhP[i].TimeBinBh;

      if(TimeBinsBh.TimeBinCount[bin] > 0)
        {
          TimeBinsBh.PrevInTimeBin[i]                                  = TimeBinsBh.LastInTimeBin[bin];
          TimeBinsBh.NextInTimeBin[i]                                  = -1;
          TimeBinsBh.NextInTimeBin[TimeBinsBh.LastInTimeBin[bin]]      = i;
          TimeBinsBh.LastInTimeBin[bin]                                = i;
        }
      else
        {
          TimeBinsBh.FirstInTimeBin[bin] = TimeBinsBh.LastInTimeBin[bin] = i;
          TimeBinsBh.PrevInTimeBin[i] = TimeBinsBh.NextInTimeBin[i] = -1;
        }
      TimeBinsBh.TimeBinCount[bin]++;
    }
}*/

/*call this function after updating the bh-timebin to the ngb condition*/
/*void update_list_of_active_bh_particles(void)
{
  int i, n;
  TimeBinsBh.NActiveParticles = 0;
  for(n = 0; n < TIMEBINS; n++)
    {
      if(TimeBinSynchronized[n]) 
        {
          for(i = TimeBinsBh.FirstInTimeBin[n]; i >= 0; i = TimeBinsBh.NextInTimeBin[i])
            {
              TimeBinsBh.ActiveParticleList[TimeBinsBh.NActiveParticles] = i;
              TimeBinsBh.NActiveParticles++;  
            }
        }
    }

    mysort(TimeBinsBh.ActiveParticleList, TimeBinsBh.NActiveParticles, sizeof(int), int_compare);*/

  /*n = 1;
  int in;
  long long out;

  in = TimeBinsBh.NActiveParticles;

  sumup_large_ints(n, &in, &out);

  TimeBinsBh.GlobalNActiveParticles = out;*/
}

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
  for(int i=0; i<NumBh; i++)
    if(BhP[i].DestroyFlag == 1)
      {
        AVP[AVP->NActiveParticles] = i;
        AVP->NActiveParticles++;
      }
}

/*
static int int_compare(const void *a, const void *b)
{
  if(*((int *)a) < *((int *)b))
    return -1;

  if(*((int *)a) > *((int *)b))
    return +1;

  return 0;
}*/
