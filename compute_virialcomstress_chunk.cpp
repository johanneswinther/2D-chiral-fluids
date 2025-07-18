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
#include "compute_virialcomstress_chunk.h"
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
#include <cstring>
#include "math_const.h"
#include "memory.h"
#include "neigh_list.h"
#include "respa.h"
#include "pair_lj_cut.h"
#include <iostream>
#include <inttypes.h>


using namespace LAMMPS_NS;

enum { ONCE, NFREQ, EVERY };

/* ----------------------------------------------------------------------

Compute the per-molecule (per-chunk) contributions to the 3x3 virial center-of-mass coarse
grained stress tensor for a 2D/3D system of molecules.

---------------------------------------------------------------------- */

ComputeVirialcomstressChunk::ComputeVirialcomstressChunk(LAMMPS *lmp, int narg, char **arg) :
  ComputeChunk(lmp, narg, arg),
  id_temp(NULL), stress(nullptr),stress_all(nullptr), massproc(nullptr), masstotal(nullptr),xcm_all(nullptr),xcm(nullptr)
{
  if (narg < 4) error->all(FLERR,"Illegal compute virialcomstress/chunk command");

  array_flag = 1;
  size_array_cols = 9;
  size_array_rows = nchunk;
  size_array_rows_variable = 1;
  extarray = 0;

  perchunk_flag = 1;
  // presschunkflag = 1;
  // peratom_flag = 1;
  timeflag = 1;
  comm_reverse = 4;
  nmax = 0;
  keflag = 0; /*I changed something here*/
  pairflag = 0;
  bondflag = 0;

  ComputeVirialcomstressChunk::init();
  ComputeVirialcomstressChunk::allocate();
}

/* ---------------------------------------------------------------------- */

ComputeVirialcomstressChunk::~ComputeVirialcomstressChunk()
{
  delete [] id_temp;
  memory->destroy(xcm);
  memory->destroy(xcm_all);
  memory->destroy(stress);
  memory->destroy(stress_all);
  memory->destroy(massproc);
  memory->destroy(masstotal);
}

/* ---------------------------------------------------------------------- */

void ComputeVirialcomstressChunk::allocate()
{
  ComputeChunk::allocate();

  memory->destroy(xcm);
  memory->destroy(xcm_all);
  memory->destroy(stress);
  memory->destroy(stress_all);
  memory->destroy(massproc);
  memory->destroy(masstotal);

  int maxatom=atom->nmax;
  maxchunk = nchunk;
  if (maxchunk == 0) error->all(FLERR, "No chunks found for virialcomstress/chunk");

  memory->create(stress,maxchunk,9,"virialcomstress/chunk:stress");
  memory->create(stress_all,maxchunk,9,"virialcomstress/chunk:stress_all");
  memory->create(massproc,maxchunk,"virialcomstress/chunk:massproc");
  memory->create(masstotal,maxchunk,"virialcomstress/chunk:masstotal");
  memory->create(xcm,maxchunk,3,"virialcomstress/chunk:xcm");
  memory->create(xcm_all,maxchunk,3,"virialcomstress/chunk:xcm_all");
  array=stress_all;
}

void ComputeVirialcomstressChunk::init_list(int /*id*/, NeighList *ptr) {
  list = ptr;
}

void ComputeVirialcomstressChunk::init() {
  //neighbor_request_id = modify->request_neighbor(this);
  if (force->pair == nullptr)
    error->all(FLERR, "No pair style is defined for compute stress/spherical");
  if (force->pair->single_enable == 0)
    error->all(FLERR, "Pair style does not support compute stress/spherical");
  neighbor->add_request(this, NeighConst::REQ_OCCASIONAL);
}

/* ---------------------------------------------------------------------- */

void ComputeVirialcomstressChunk::setup()
{
  // one-time calculation of per-chunk mass
  // done in setup, so that ComputeChunkAtom::setup() is already called

  if (firstflag && cchunk->idsflag == ONCE) {
    compute_array();
    firstflag = massneed = 0;
  }
}

/* ---------------------------------------------------------------------- */

void ComputeVirialcomstressChunk::compute_array()
{
  ComputeChunk::compute_array();
  size_array_rows = nchunk;

  //invoked_vector = update->ntimestep;
    
  Pair *pair = force->pair;
  double **cutsq = pair->cutsq;  // Get cutoff distance squared

  int l, m;
  int index;
  double massone;

  int i, j, ii, jj, inum, jnum, itype, jtype;
  double xtmp, ytmp, ztmp, delx, dely, delz, fpair;
  double rsq, r2inv, r6inv, forcelj, factor_lj,factor_coul,eng;
  int *ilist, *jlist, *numneigh, **firstneigh;

  double **x = atom->x;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_lj = force->special_lj;
  double *special_coul = force->special_coul;
  int newton_pair = force->newton_pair;

  neighbor->build_one(list);

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  int *mask = atom->mask;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int npair = nlocal;

  int *ichunk = cchunk->ichunk;

  //Zero out elements
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

    stress_all[i][0]=0.0;
    stress_all[i][1]=0.0;
    stress_all[i][2]=0.0;
    stress_all[i][3]=0.0;
    stress_all[i][4]=0.0;
    stress_all[i][5]=0.0;
    stress_all[i][6]=0.0;
    stress_all[i][7]=0.0;
    stress_all[i][8]=0.0;

    xcm[i][0]=0.0;
    xcm[i][1]=0.0;
    xcm[i][2]=0.0;

    xcm_all[i][0]=0.0;
    xcm_all[i][1]=0.0;
    xcm_all[i][2]=0.0;
  }
  
  // zero local per-chunk values
  if (massneed)
    for (int i = 0; i < nchunk; i++) massproc[i] = 0.0;

  // compute COM for each chunk and do this by unwrapping the coordinates
  imageint *image = atom->image;
  double unwrap[3];

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      index = ichunk[i] - 1;
      if (index < 0) continue;
      if (rmass)
        massone = rmass[i];
      else
        massone = mass[type[i]];
      domain->unmap(x[i], image[i], unwrap);
      xcm[index][0] += unwrap[0] * massone;
      xcm[index][1] += unwrap[1] * massone;
      xcm[index][2] += unwrap[2] * massone;
      if (massneed) massproc[index] += massone;
    }

  MPI_Allreduce(&xcm[0][0], &xcm_all[0][0], 3 * nchunk, MPI_DOUBLE, MPI_SUM, world);
  if (massneed) MPI_Allreduce(massproc, masstotal, nchunk, MPI_DOUBLE, MPI_SUM, world);

  double Lx = domain->boxhi[0] - domain->boxlo[0];
  double Ly = domain->boxhi[1] - domain->boxlo[1];
  double Lz = domain->boxhi[2] - domain->boxlo[2];

    // // 2D or 3D
  nktv2p = force->nktv2p;
  if (domain->dimension == 2){
    inv_volume = 1.0 / (Lx*Ly);
  } else {
    inv_volume = 1.0 / (Lx*Ly*Lz);
  }

  scale=0.5*nktv2p*inv_volume;

  // Get forces for each atom on the processor and bin them into their corresponding intermolecular interactions
  for (ii = 0; ii < inum; ii++) {

    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];

    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    // Check if atom i is in the group
    if (!(mask[i] & groupbit)) continue;

    // Molecule index (is shared between processors when using atom style "full")
    int i_index;
    i_index= atom->molecule[i] - 1;
    double com_i[3];

    //Get COM for molecule i
    for (int d = 0; d < 3; ++d) {
      com_i[d] = xcm_all[i_index][d] / masstotal[i_index];
    }

    for (jj = 0; jj < jnum; jj++) {

      j = jlist[jj];

      j &= NEIGHMASK;

      // Check if atom j is in the group
      if (!(mask[j] & groupbit)) continue;

      // Molecule index (is shared between processors when using atom style "full")
      int j_index;
      j_index= atom->molecule[j] - 1;

      if (i_index == j_index) continue;
      factor_lj   = special_lj[sbmask(j)];
      factor_coul = special_coul[sbmask(j)];

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];

      rsq = delx * delx + dely * dely + delz * delz;
      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {
        eng = force->pair->single(i, j, itype, jtype, rsq, factor_coul, factor_lj, fpair);

        if (masstotal[i_index] > 0.0) {
          
          double com_j[3], delcom[3];

          // Get COM for molecule j
          for (int d = 0; d < 3; ++d) {
            com_j[d] = xcm_all[j_index][d] / masstotal[j_index];
          }

          // Get the intermolecular distance
          for (int d = 0; d < 3; ++d) {
            delcom[d] = com_i[d] - com_j[d];
          }

          //Rewrap the vector, so there are no discontinuities when molecules cross boundaries
          domain->minimum_image_big(delcom);

          //Reminder: delx * fpair is the interatomic force and dx is the intermolecular distance
          stress[i_index][0] += delx*fpair*delcom[0]*scale; //Txx
          stress[i_index][1] += delx*fpair*delcom[1]*scale; //Txy
          stress[i_index][2] += delx*fpair*delcom[2]*scale; //Txz
          stress[i_index][3] += dely*fpair*delcom[0]*scale; //Tyx
          stress[i_index][4] += dely*fpair*delcom[1]*scale; //Tyy
          stress[i_index][5] += dely*fpair*delcom[2]*scale; //Tyz
          stress[i_index][6] += delz*fpair*delcom[0]*scale; //Tzx
          stress[i_index][7] += delz*fpair*delcom[1]*scale; //Tzy
          stress[i_index][8] += delz*fpair*delcom[2]*scale; //Tzz

          if (newton_pair || j < nlocal) {
            stress[j_index][0] += delx*fpair*delcom[0]*scale; //Txx
            stress[j_index][1] += delx*fpair*delcom[1]*scale; //Txy
            stress[j_index][2] += delx*fpair*delcom[2]*scale; //Txz
            stress[j_index][3] += dely*fpair*delcom[0]*scale; //Tyx
            stress[j_index][4] += dely*fpair*delcom[1]*scale; //Tyy
            stress[j_index][5] += dely*fpair*delcom[2]*scale; //Tyz
            stress[j_index][6] += delz*fpair*delcom[0]*scale; //Tzx
            stress[j_index][7] += delz*fpair*delcom[1]*scale; //Tzy
            stress[j_index][8] += delz*fpair*delcom[2]*scale; //Tzz
          }
        } else {
          stress[i_index][0] = stress[i_index][1] = stress[i_index][2] = stress[i_index][3] = stress[i_index][4] = stress[i_index][5] = stress[i_index][6] = stress[i_index][7] = stress[i_index][8] = 0.0;
          if (newton_pair || j < nlocal) {
            stress[j_index][0] = stress[j_index][1] = stress[j_index][2] = stress[j_index][3] = stress[j_index][4] = stress[j_index][5] = stress[j_index][6] = stress[j_index][7] = stress[j_index][8] = 0.0;
          }
        }
      }
    }
  }

  //Gather all of the stress data
  MPI_Allreduce(&stress[0][0],
    &stress_all[0][0],
    nchunk * 9, MPI_DOUBLE, MPI_SUM, world);
}

double ComputeVirialcomstressChunk::memory_usage()
{
  double bytes = (double) (maxchunk*9*2+maxchunk*2+maxchunk*2*3)* sizeof(double);
  return bytes;
}
