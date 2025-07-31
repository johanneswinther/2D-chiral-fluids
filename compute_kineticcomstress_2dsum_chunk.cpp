/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.

   This file is part of custom LAMMPS code for simulations of chiral matter.
    author: Johannes Winther
    e-mail: johannesiwk@gmail.com
    github: https://github.com/mandadapu-group/2D-chiral-fluids
------------------------------------------------------------------------- */

#include "domain.h"
#include <cstdlib>
#include <cstring>
#include <cmath>
#include "compute_kineticcomstress_2dsum_chunk.h"
#include "neighbor.h"
#include "atom.h"
#include "update.h"
#include "comm.h"
#include "force.h"
#include "pair.h"
#include "bond.h"
#include "modify.h"
#include "fix.h"
#include "memory.h"
#include "error.h"
#include "compute_chunk_atom.h"
#include "compute_property_atom.h"
#include <array>  // Required for std::array

using namespace LAMMPS_NS;

enum { ONCE, NFREQ, EVERY };

/* ----------------------------------------------------------------------

Compute the per-molecule (per-chunk) contributions to the 3x3 kinetic center-of-mass coarse
grained stress tensor for a 2D system of molecules.

---------------------------------------------------------------------- */

ComputeKineticcomstress2DSumChunk::ComputeKineticcomstress2DSumChunk(LAMMPS *lmp, int narg, char **arg) :
  ComputeChunk(lmp, narg, arg),
  id_temp(NULL), stress(nullptr),stress_all(nullptr), massproc(nullptr), masstotal(nullptr), vcm(nullptr), vcmall(nullptr)
{
  if (narg < 4) error->all(FLERR,"Illegal compute kineticcomstress/chunk command");

  vector_flag = 1;
  extvector = 0;
  size_vector = 4;
  
  ComputeKineticcomstress2DSumChunk::init();
  ComputeKineticcomstress2DSumChunk::allocate();
}

/* ---------------------------------------------------------------------- */

ComputeKineticcomstress2DSumChunk::~ComputeKineticcomstress2DSumChunk()
{
  delete [] id_temp;
  memory->destroy(stress);
  memory->destroy(stress_all);
  memory->destroy(massproc);
  memory->destroy(masstotal);
  memory->destroy(vcm);
  memory->destroy(vcmall);
}

/* ---------------------------------------------------------------------- */

void ComputeKineticcomstress2DSumChunk::allocate()
{
  ComputeChunk::allocate();

  memory->destroy(stress);
  memory->destroy(stress_all);
  memory->destroy(massproc);
  memory->destroy(masstotal);
  memory->destroy(vcm);
  memory->destroy(vcmall);

  maxchunk = nchunk;
  if (maxchunk == 0) error->all(FLERR, "No chunks found for kineticcomstress2dsum/chunk");

  memory->create(stress,4,"kineticcomstress2dsum/chunk:stress");
  memory->create(stress_all,4,"kineticcomstress2dsum/chunk:stress_all");
  memory->create(massproc, maxchunk, "kineticcomstress2dsum/chunk:massproc");
  memory->create(masstotal, maxchunk, "kineticcomstress2dsum/chunk:masstotal");
  memory->create(vcm, maxchunk, 2, "kineticcomstress2dsum/chunk:vcm");
  memory->create(vcmall, maxchunk, 2, "kineticcomstress2dsum/chunk:vcmall");
  memory->create(stress_comp, 4, "kineticcomstress2dsum/chunk:stress_comp");
  memory->create(vcm_comp, maxchunk, 2, "kineticcomstress2dsum/chunk:vcm_comp");
  vector=stress_all;
}


/* ---------------------------------------------------------------------- */

void ComputeKineticcomstress2DSumChunk::setup()
{
  // one-time calculation of per-chunk mass
  // done in setup, so that ComputeChunkAtom::setup() is already called

  if (firstflag && cchunk->idsflag == ONCE) {
    compute_vector();
    firstflag = massneed = 0;
  }
}

/* ---------------------------------------------------------------------- */

void ComputeKineticcomstress2DSumChunk::compute_vector()
{
  int index;
  double massone;

  ComputeChunk::compute_vector();
  int *ichunk = cchunk->ichunk;

  for (int i=0; i<4; i++) {
    stress_all[i]=0.0;
    stress_comp[i]=0.0;
  }

  for (int i = 0; i < nchunk; i++) {
    for (int d = 0; d < 2; d++) {
      vcm[i][d] = 0.0;
      vcmall[i][d] = 0.0;
      vcm_comp[i][d] = 0.0;
    }
  }
  double **v = atom->v;
  int *mask = atom->mask;
  int *type = atom->type;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int nlocal = atom->nlocal;


  // For printing time (debugging)
  // bigint timestep = update->ntimestep;
  // double dt = update->dt;  // if available
  // double simtime = timestep * dt;
  // if (comm->me == 0) {
  //   printf("Current timestep: %lld, simulation time: %g\n", timestep, simtime);
  //   fflush(stdout);
  // }

  // Compute VCM for each chunk on this processor
  // This piece of code below is adapted from "compute_vcm_chunk.cpp"
  if (massneed)
    for (int i = 0; i < nchunk; i++) massproc[i] = 0.0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      index = ichunk[i] - 1;
      if (index < 0) continue;
      if (rmass)
        massone = rmass[i];
      else
        massone = mass[type[i]];
      // temperature->remove_bias(i, v[i]); //remove macroscopic velocity bias

      //normal summation
      // vcm[index][0] += v[i][0] * massone;
      // vcm[index][1] += v[i][1] * massone;
      // vcm[index][2] += v[i][2] * massone;
      
      //Kahan summation
      for (int d = 0; d < 2; d++) {
        // Kahan summation of vcm[index][d] += v[i][d] * massone
        double term = v[i][d] * massone;
        double y    = term - vcm_comp[index][d];
        double t    = vcm[index][d] + y;
        vcm_comp[index][d] = (t - vcm[index][d]) - y;
        vcm[index][d]      = t;
      }
      // temperature->restore_bias(i, v[i]); //restore macroscopic velocity bias
      if (massneed) massproc[index] += massone;
    }
  //This piece of code above is adapted from "compute_vcm_chunk.cpp"

  
  // Sum mass and VCM over all processors
  MPI_Allreduce(&vcm[0][0], &vcmall[0][0], 2 * nchunk, MPI_DOUBLE, MPI_SUM, world);
  if (massneed) MPI_Allreduce(massproc, masstotal, nchunk, MPI_DOUBLE, MPI_SUM, world);

  // Normalize VCM by mass
  for (int i = 0; i < nchunk; i++) {
    if (masstotal[i] > 0.0) {
      vcmall[i][0] /= masstotal[i];
      vcmall[i][1] /= masstotal[i];

    } else
      vcmall[i][0] = vcmall[i][1] = 0.0;
  }

  double Lx = domain->boxhi[0] - domain->boxlo[0];
  double Ly = domain->boxhi[1] - domain->boxlo[1];

  inv_volume = 1.0 / (Lx*Ly);


  //normal summation
  // for (int i = 0; i < nchunk; i++) {
  //     stress_all[0] += masstotal[i]*vcmall[i][0]*vcmall[i][0]; // T_xx
  //     stress_all[1] += masstotal[i]*vcmall[i][0]*vcmall[i][1]; // T_xy
  //     stress_all[2] += masstotal[i]*vcmall[i][0]*vcmall[i][2]; // T_xz
  //     stress_all[3] += masstotal[i]*vcmall[i][1]*vcmall[i][0]; // T_yx
  //     stress_all[4] += masstotal[i]*vcmall[i][1]*vcmall[i][1]; // T_yy
  //     stress_all[5] += masstotal[i]*vcmall[i][1]*vcmall[i][2]; // T_yz
  //     stress_all[6] += masstotal[i]*vcmall[i][2]*vcmall[i][0]; // T_zx
  //     stress_all[7] += masstotal[i]*vcmall[i][2]*vcmall[i][1]; // T_zy 
  //     stress_all[8] += masstotal[i]*vcmall[i][2]*vcmall[i][2]; // T_zz
  // }

  //Kahan summation
  for (int i = 0; i < nchunk; i++) {
    // precompute the 9 stress contributions for chunk i
    double m   = masstotal[i];
    double vx  = vcmall[i][0];
    double vy  = vcmall[i][1];
    double terms[4] = {
      m*vx*vx,   // T_xx
      m*vx*vy,   // T_xy
      m*vy*vx,   // T_yx
      m*vy*vy,   // T_yy
    };
  
    // Kahan summation over all 4 components
    for (int c = 0; c < 4; c++) {
      double y = terms[c] - stress_comp[c];
      double t = stress_all[c] + y;
      stress_comp[c] = (t - stress_all[c]) - y;
      stress_all[c] = t;
    }
  }

  for (int i=0; i<4; i++) {
    stress_all[i] *= inv_volume; // Normalize by volume
  } 
}

double ComputeKineticcomstress2DSumChunk::memory_usage()
{
  double bytes = (double) (4*2+maxchunk*2+maxchunk*2*2)* sizeof(double);
  return bytes;
}
