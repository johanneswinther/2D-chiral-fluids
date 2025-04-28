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
#include "compute_kineticcomstress_chunk.h"
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
grained stress tensor for a 2D/3D system of molecules.

---------------------------------------------------------------------- */

ComputeKineticcomstressChunk::ComputeKineticcomstressChunk(LAMMPS *lmp, int narg, char **arg) :
  ComputeChunk(lmp, narg, arg),
  id_temp(NULL), stress(nullptr), massproc(nullptr), masstotal(nullptr), vcm(nullptr), vcmall(nullptr)
{
  if (narg < 4) error->all(FLERR,"Illegal compute kineticcomstress/chunk command");

  array_flag = 1;
  size_array_cols = 9;
  size_array_rows = 0;
  size_array_rows_variable = 1;
  extarray = 0;

  ComputeKineticcomstressChunk::init();
  ComputeKineticcomstressChunk::allocate();
}

/* ---------------------------------------------------------------------- */

ComputeKineticcomstressChunk::~ComputeKineticcomstressChunk()
{
  delete [] id_temp;
  memory->destroy(stress);
  memory->destroy(massproc);
  memory->destroy(masstotal);
  memory->destroy(vcm);
  memory->destroy(vcmall);
}

/* ---------------------------------------------------------------------- */

void ComputeKineticcomstressChunk::allocate()
{
  ComputeChunk::allocate();

  memory->destroy(stress);
  memory->destroy(massproc);
  memory->destroy(masstotal);
  memory->destroy(vcm);
  memory->destroy(vcmall);

  maxchunk = nchunk;
  if (maxchunk == 0) error->all(FLERR, "No chunks found for kineticcomstress/chunk");

  memory->create(stress,maxchunk,9,"kineticcomstress/chunk:stress");
  memory->create(massproc, maxchunk, "kineticcomstress/chunk:massproc");
  memory->create(masstotal, maxchunk, "kineticcomstress/chunk:masstotal");
  memory->create(vcm, maxchunk, 3, "kineticcomstress/chunk:vcm");
  memory->create(vcmall, maxchunk, 3, "kineticcomstress/chunk:vcmall");
  array=stress;
}


/* ---------------------------------------------------------------------- */

void ComputeKineticcomstressChunk::setup()
{
  // one-time calculation of per-chunk mass
  // done in setup, so that ComputeChunkAtom::setup() is already called

  if (firstflag && cchunk->idsflag == ONCE) {
    compute_array();
    firstflag = massneed = 0;
  }
}

/* ---------------------------------------------------------------------- */

void ComputeKineticcomstressChunk::compute_array()
{
  int index;
  double massone;

  ComputeChunk::compute_array();
  int *ichunk = cchunk->ichunk;

  for (int i = 0; i < nchunk; i++) {
    stress[i][0]=0.0;
    stress[i][1]=0.0;
    stress[i][2]=0.0;
    stress[i][3]=0.0;
    stress[i][4]=0.0;
    stress[i][5]=0.0;
    stress[i][6]=0.0;
    stress[i][7]=0.0;
    stress[i][8]=0.0;
    vcm[i][0] = 0.0;
    vcm[i][1] = 0.0;
    vcm[i][2] = 0.0;
    vcmall[i][0] = 0.0;
    vcmall[i][1] = 0.0;
    vcmall[i][2] = 0.0;
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
      vcm[index][0] += v[i][0] * massone;
      vcm[index][1] += v[i][1] * massone;
      if (massneed) massproc[index] += massone;
    }
  //This piece of code above is adapted from "compute_vcm_chunk.cpp"

  
  // Sum mass and VCM over all processors
  MPI_Allreduce(&vcm[0][0], &vcmall[0][0], 3 * nchunk, MPI_DOUBLE, MPI_SUM, world);
  if (massneed) MPI_Allreduce(massproc, masstotal, nchunk, MPI_DOUBLE, MPI_SUM, world);

  // Normalize VCM by mass
  for (int i = 0; i < nchunk; i++) {
    if (masstotal[i] > 0.0) {
      vcmall[i][0] /= masstotal[i];
      vcmall[i][1] /= masstotal[i];
      vcmall[i][2] /= masstotal[i];

    } else
      vcmall[i][0] = vcmall[i][1] = vcmall[i][2] = 0.0;
  }

  double Lx = domain->boxhi[0] - domain->boxlo[0];
  double Ly = domain->boxhi[1] - domain->boxlo[1];
  double Lz = domain->boxhi[2] - domain->boxlo[2];

  // 2D or 3D
  if (domain->dimension == 2){
    inv_volume = 1.0 / (Lx*Ly);
  } else {
    inv_volume = 1.0 / (Lx*Ly*Lz);
  }

  // inv_volume =1.0;

  for (int i = 0; i < nchunk; i++) {
    if (masstotal[i] > 0.0) {
      stress[i][0] = masstotal[i]*vcmall[i][0]*vcmall[i][0]*inv_volume; // T_xx
      stress[i][1] = masstotal[i]*vcmall[i][0]*vcmall[i][1]*inv_volume; // T_xy
      stress[i][2] = masstotal[i]*vcmall[i][0]*vcmall[i][2]*inv_volume; // T_xz
      stress[i][3] = masstotal[i]*vcmall[i][1]*vcmall[i][0]*inv_volume; // T_yx
      stress[i][4] = masstotal[i]*vcmall[i][1]*vcmall[i][1]*inv_volume; // T_yy 
      stress[i][5] = masstotal[i]*vcmall[i][1]*vcmall[i][2]*inv_volume; // T_yz
      stress[i][6] = masstotal[i]*vcmall[i][2]*vcmall[i][0]*inv_volume; // T_zx
      stress[i][7] = masstotal[i]*vcmall[i][2]*vcmall[i][1]*inv_volume; // T_zy 
      stress[i][8] = masstotal[i]*vcmall[i][2]*vcmall[i][2]*inv_volume; // T_zz
    } else
      stress[i][0] = stress[i][1] = stress[i][2] = stress[i][3] = stress[i][4] = stress[i][5] = stress[i][6] = stress[i][7] = stress[i][8]= 0.0;
  }
  size_array_rows = nchunk;
}

double ComputeKineticcomstressChunk::memory_usage()
{
  double bytes = sizeof(double)*(maxchunk*9+maxchunk*2+2*3*maxchunk);
  return bytes;
}
