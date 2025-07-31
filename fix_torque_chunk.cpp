// clang-format off
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


#include "fix_torque_chunk.h"

#include "atom.h"
#include "comm.h"
#include "compute_chunk_atom.h"
#include "compute_com_chunk.h"
#include "error.h"
#include "memory.h"
#include "modify.h"
#include "respa.h"
#include "update.h"
#include "domain.h"

#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

// ----------------------------------------------------------------------
// This fix is used to torque individual chunks in 2d, around their center of mass.
// ----------------------------------------------------------------------

/* ---------------------------------------------------------------------- */

FixTorqueChunk::FixTorqueChunk(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  idchunk(nullptr), idcom(nullptr)
{
  if (narg != 6) utils::missing_cmd_args(FLERR, "fix torque/chunk", error);

  restart_global = 1;
  scalar_flag = 1;
  global_freq = 1;
  extscalar = 1;
  energy_global_flag = 1;
  respa_level_support = 1;
  ilevel_respa = 0;

  k_torque = utils::numeric(FLERR,arg[3],false,lmp);

  idchunk = utils::strdup(arg[4]);
  idcom = utils::strdup(arg[5]);

  nchunk = 0;
}

/* ---------------------------------------------------------------------- */

FixTorqueChunk::~FixTorqueChunk()
{

  // decrement lock counter in compute chunk/atom, it if still exists

  cchunk = dynamic_cast<ComputeChunkAtom *>(modify->get_compute_by_id(idchunk));
  if (cchunk) {
    cchunk->unlock(this);
    cchunk->lockcount--;
  }
  delete[] idchunk;
  delete[] idcom;
}

/* ---------------------------------------------------------------------- */

int FixTorqueChunk::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixTorqueChunk::init()
{
  // current indices for idchunk and idcom

  cchunk = dynamic_cast<ComputeChunkAtom *>(modify->get_compute_by_id(idchunk));
  if (!cchunk)
    error->all(FLERR,"Chunk/atom compute {} does not exist or is not chunk/atom style", idchunk);

  ccom = dynamic_cast<ComputeCOMChunk *>(modify->get_compute_by_id(idcom));
  if (!ccom)
    error->all(FLERR,"Com/chunk compute {} does not exist or is not com/chunk style", idcom);

  // check that idchunk is consistent with ccom->idchunk

  if (ccom && (strcmp(idchunk,ccom->idchunk) != 0))
    error->all(FLERR,"Fix torque/chunk chunk ID {} not the same as compute com/chunk chunk ID {}",
               idchunk, ccom->idchunk);

  if (utils::strmatch(update->integrate_style,"^respa")) {
    ilevel_respa = (dynamic_cast<Respa *>(update->integrate))->nlevels-1;
    if (respa_level >= 0) ilevel_respa = MIN(respa_level,ilevel_respa);
  }
}

/* ---------------------------------------------------------------------- */

void FixTorqueChunk::setup(int vflag)
{
  if (utils::strmatch(update->integrate_style,"^verlet"))
    post_force(vflag);
  else {
    (dynamic_cast<Respa *>(update->integrate))->copy_flevel_f(ilevel_respa);
    post_force_respa(vflag,ilevel_respa,0);
    (dynamic_cast<Respa *>(update->integrate))->copy_f_flevel(ilevel_respa);
  }
}

/* ---------------------------------------------------------------------- */

void FixTorqueChunk::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixTorqueChunk::post_force(int /*vflag*/)
{
  // calculate current centers of mass for each chunk
  // extract pointers from idchunk and idcom
  // compute COM for each chunk
  ccom->compute_array();

  nchunk = cchunk->nchunk;
  int *ichunk = cchunk->ichunk;
  double **com = ccom->array;
  double **f = atom->f;
  int nlocal = atom->nlocal;

  //Calculate vectors from COM to each atom in the chunk
  double **x = atom->x;
  int *mask = atom->mask;
  imageint *image = atom->image;
  int i, index;
  double dx,dy; 
  double unwrap[3];

  for (i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      index = ichunk[i] - 1;
      if (index < 0) continue;
      domain->unmap(x[i], image[i], unwrap);
      dx = unwrap[0] - com[index][0];
      dy = unwrap[1] - com[index][1];
      
      f[i][0] -= k_torque*dy;
      f[i][1] += k_torque*dx;
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixTorqueChunk::post_force_respa(int vflag, int ilevel, int /*iloop*/)
{
  if (ilevel == ilevel_respa) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixTorqueChunk::min_post_force(int vflag)
{
  post_force(vflag);
}

