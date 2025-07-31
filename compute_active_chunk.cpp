/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.

   This file is part of custom LAMMPS code for simulations of active matter.
    author: Johannes Winther
    e-mail: johannesiwk@gmail.com
    github: https://github.com/johanneswinther/2D-chiral-fluids
------------------------------------------------------------------------- */

#include "compute_active_chunk.h"

#include "atom.h"
#include "compute_chunk_atom.h"
#include "domain.h"
#include "error.h"
#include "memory.h"
#include "force.h"
#include "update.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeActiveChunk::ComputeActiveChunk(LAMMPS *lmp, int narg, char **arg) :
    ComputeChunk(lmp, narg, arg), massproc(nullptr), masstotal(nullptr), com(nullptr),
    comall(nullptr), stress(nullptr), stressall(nullptr), summedstress(nullptr), k_torque(0.0)
{
  if (narg != 5) error->all(FLERR, "Illegal compute active/chunk command");

  vector_flag = 1;
  size_vector = 3;
  size_array_rows = 0;
  size_array_rows_variable = 1;
  extvector = 0;

  k_torque = utils::numeric(FLERR,arg[4],false,lmp);

  ComputeActiveChunk::init();
  ComputeActiveChunk::allocate();
}

/* ---------------------------------------------------------------------- */

ComputeActiveChunk::~ComputeActiveChunk()
{
  memory->destroy(massproc);
  memory->destroy(masstotal);
  memory->destroy(com);
  memory->destroy(comall);
  memory->destroy(stress);
  memory->destroy(stressall);
  memory->destroy(summedstress);
}

/* ---------------------------------------------------------------------- */

void ComputeActiveChunk::compute_vector()
{
  int i, index;
  double dx, dy, dz, massone;
  double unwrap[3];

  ComputeChunk::compute_vector();
  int *ichunk = cchunk->ichunk;

  invoked_vector = update->ntimestep;

  // zero local per-chunk values

  for (i = 0; i < nchunk; i++) {
    massproc[i] = 0.0;
    com[i][0] = com[i][1] = com[i][2] = 0.0;
    stress[i][0] = stress[i][1] = stress[i][2] = stress[i][3] = 0.0;
  }

  // compute COM for each chunk

  double **x = atom->x;
  int *mask = atom->mask;
  int *type = atom->type;
  imageint *image = atom->image;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int nlocal = atom->nlocal;

  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      index = ichunk[i] - 1;
      if (index < 0) continue;
      if (rmass)
        massone = rmass[i];
      else
        massone = mass[type[i]];
      domain->unmap(x[i], image[i], unwrap);
      massproc[index] += massone;
      com[index][0] += unwrap[0] * massone;
      com[index][1] += unwrap[1] * massone;
      com[index][2] += unwrap[2] * massone;
    }

  MPI_Allreduce(massproc, masstotal, nchunk, MPI_DOUBLE, MPI_SUM, world);
  MPI_Allreduce(&com[0][0], &comall[0][0], 3 * nchunk, MPI_DOUBLE, MPI_SUM, world);

  for (i = 0; i < nchunk; i++) {
    if (masstotal[i] > 0.0) {
      comall[i][0] /= masstotal[i];
      comall[i][1] /= masstotal[i];
      comall[i][2] /= masstotal[i];
    }
  }

  // compute torque on each chunk

  double **f = atom->f;

  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      index = ichunk[i] - 1;
      if (index < 0) continue;
      domain->unmap(x[i], image[i], unwrap);
      dx = unwrap[0] - comall[index][0];
      dy = unwrap[1] - comall[index][1];
      fx = -k_torque*dy;
      fy = k_torque*dx;
      stress[index][0] += fx*dx;
      stress[index][1] += fx*dy;
      stress[index][2] += fy*dx;
      stress[index][3] += fy*dy;
    }

  MPI_Allreduce(&stress[0][0], &stressall[0][0], 4 * nchunk, MPI_DOUBLE, MPI_SUM, world);


  summedstress[0] = summedstress[1] = summedstress[2] = summedstress[3]=0;

  for (int i = 0; i < nchunk; i++) {
    for (int j = 0; j < 4; j++) {
      summedstress[j] += stressall[i][j];
    }
  }

  nktv2p = force->nktv2p;
  double Lx = domain->boxhi[0] - domain->boxlo[0];
  double Ly = domain->boxhi[1] - domain->boxlo[1];
  inv_vol = 1.0 / (Lx*Ly);

  for (int i = 0; i < 4; i++)
    summedstress[i] *= inv_vol * nktv2p;
}
/* ----------------------------------------------------------------------
   free and reallocate per-chunk arrays
------------------------------------------------------------------------- */

void ComputeActiveChunk::allocate()
{
  ComputeChunk::allocate();
  memory->destroy(massproc);
  memory->destroy(masstotal);
  memory->destroy(com);
  memory->destroy(comall);
  memory->destroy(stress);
  memory->destroy(stressall);
  memory->destroy(summedstress);
  maxchunk = nchunk;
  memory->create(massproc, maxchunk, "active/chunk:massproc");
  memory->create(masstotal, maxchunk, "active/chunk:masstotal");
  memory->create(com, maxchunk, 3, "active/chunk:com");
  memory->create(comall, maxchunk, 3, "active/chunk:comall");
  memory->create(stress, maxchunk, 4, "active/chunk:stress");
  memory->create(stressall, maxchunk, 4, "active/chunk:stressall");
  memory->create(summedstress, 4, "active/chunk:summedstress");
  vector = summedstress;
}

/* ----------------------------------------------------------------------
   memory usage of local data
------------------------------------------------------------------------- */

double ComputeActiveChunk::memory_usage()
{
  double bytes = (double) maxchunk * 2 * sizeof(double) + ComputeChunk::memory_usage();
  bytes += (double) maxchunk * 2 * 3 * sizeof(double);
  bytes += (double) maxchunk * 2 * 4 * sizeof(double);
  return bytes;
}
