/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.

   This file is part of custom LAMMPS code for simulations of chiral matter and is based on compute_angmom_chunk.cpp.
    author: Johannes Winther
    e-mail: johannesiwk@gmail.com
    github: https://github.com/mandadapu-group/2D-chiral-fluids
------------------------------------------------------------------------- */


#include "compute_spintemp_sum_chunk.h"

#include "atom.h"
#include "compute_chunk_atom.h"
#include "domain.h"
#include "error.h"
#include "memory.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeSpintempSumChunk::ComputeSpintempSumChunk(LAMMPS *lmp, int narg, char **arg) :
    ComputeChunk(lmp, narg, arg), massproc(nullptr), masstotal(nullptr), com(nullptr),
    comall(nullptr), angmom(nullptr), angmomall(nullptr), spintempall(nullptr)
{
  if (narg != 4) error->all(FLERR, "Illegal compute angvel/chunk command");

  vector_flag = 1;
  size_vector = 3;
  size_array_rows = 0;
  size_array_rows_variable = 1;
  extvector = 0;

  ComputeSpintempSumChunk::init();
  ComputeSpintempSumChunk::allocate();
}

/* ---------------------------------------------------------------------- */

ComputeSpintempSumChunk::~ComputeSpintempSumChunk()
{
  memory->destroy(massproc);
  memory->destroy(masstotal);
  memory->destroy(com);
  memory->destroy(comall);
  memory->destroy(angmom);
  memory->destroy(angmomall);
  memory->destroy(spintempall);
}

/* ---------------------------------------------------------------------- */

void ComputeSpintempSumChunk::compute_vector()
{
  ComputeChunk::compute_vector();

  int i, index;
  double dx, dy, dz, massone;
  double unwrap[3];
  int *ichunk = cchunk->ichunk;

  // zero local per-chunk values

  for (i = 0; i < nchunk; i++) {
    massproc[i] = 0.0;
    com[i][0] = com[i][1] = com[i][2] = 0.0;
    angmom[i][0] = angmom[i][1] = angmom[i][2] = 0.0;
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

  // compute angmom for each chunk

  double **v = atom->v;

  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      index = ichunk[i] - 1;
      if (index < 0) continue;
      domain->unmap(x[i], image[i], unwrap);
      dx = unwrap[0] - comall[index][0];
      dy = unwrap[1] - comall[index][1];
      dz = unwrap[2] - comall[index][2];

      if (rmass)
        massone = rmass[i];
      else
        massone = mass[type[i]];
      
      angmom[index][0] += (dy * v[i][2] - dz * v[i][1])*massone;
      angmom[index][1] += (dz * v[i][0] - dx * v[i][2])*massone;
      angmom[index][2] += (dx * v[i][1] - dy * v[i][0])*massone;

    }
  
  MPI_Allreduce(&angmom[0][0], &angmomall[0][0], 3 * nchunk, MPI_DOUBLE, MPI_SUM, world);

  //Now that I have the angular velocity, I can sum it

  avg_mass = 0.0;
  spintempall[0] = spintempall[1] = spintempall[2] = 0.0;

  for (int i = 0; i < nchunk; i++) {
    spintempall[0] += angmomall[i][0]*angmomall[i][0];
    spintempall[1] += angmomall[i][1]*angmomall[i][1];
    spintempall[2] += angmomall[i][2]*angmomall[i][2];
    avg_mass += masstotal[i];
  }
  avg_mass /= nchunk;

  double Lx = domain->boxhi[0] - domain->boxlo[0];
  double Ly = domain->boxhi[1] - domain->boxlo[1];
  double Lz = domain->boxhi[2] - domain->boxlo[2];

  // 2D or 3D

  if (domain->dimension == 2){
    inv_vol = 1.0 / (Lx*Ly);
  } else {
    inv_vol = 1.0 / (Lx*Ly*Lz);
  }

  rho_bead = avg_mass*nchunk*inv_vol;

  scale = avg_mass/(rho_bead*rho_bead*rho_bead)*inv_vol;

  //Multiply by inv vol^2

  for (i = 0; i < 3; i++) {
    spintempall[i] *= scale;
  }

}

/* ----------------------------------------------------------------------
   free and reallocate per-chunk arrays
------------------------------------------------------------------------- */

void ComputeSpintempSumChunk::allocate()
{
  memory->destroy(massproc);
  memory->destroy(masstotal);
  memory->destroy(com);
  memory->destroy(comall);
  memory->destroy(angmom);
  memory->destroy(angmomall);
  maxchunk = nchunk;
  memory->create(massproc, maxchunk, "spintemp/chunk:massproc");
  memory->create(masstotal, maxchunk, "spintemp/chunk:masstotal");
  memory->create(com, maxchunk, 3, "spintemp/chunk:com");
  memory->create(comall, maxchunk, 3, "spintemp/chunk:comall");
  memory->create(angmom, maxchunk, 3, "spintemp/chunk:angmom");
  memory->create(angmomall, maxchunk, 3, "spintemp/chunk:angmomall");
  memory->create(spintempall, 3, "spintemp/chunk:spintempall");
  vector = spintempall;
}

/* ----------------------------------------------------------------------
   memory usage of local data
------------------------------------------------------------------------- */

double ComputeSpintempSumChunk::memory_usage()
{
  double bytes = ComputeChunk::memory_usage();
  bytes += (bigint) maxchunk * 2 * sizeof(double);
  bytes += (double) maxchunk * 2 * 3 * sizeof(double);
  bytes += (double) maxchunk * 2 * 3 * sizeof(double);
  return bytes;
}
