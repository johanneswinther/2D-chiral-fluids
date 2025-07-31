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
#include "compute_virialcomstress_2dsum_chunk.h"
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
grained stress tensor for a 2D system of molecules.

---------------------------------------------------------------------- */

ComputeVirialcomstress2DSumChunk::ComputeVirialcomstress2DSumChunk(LAMMPS *lmp, int narg, char **arg) :
  ComputeChunk(lmp, narg, arg),
  id_temp(NULL), stress(nullptr),stress_all(nullptr), massproc(nullptr), masstotal(nullptr),xcm_all(nullptr),xcm(nullptr)
{
  if (narg < 4) error->all(FLERR,"Illegal compute virialcomstress/chunk command");

  vector_flag = 1;
  extvector = 0;
  size_vector = 4;


  perchunk_flag = 1;
  // presschunkflag = 1;
  // peratom_flag = 1;
  timeflag = 1;
  comm_reverse = 4;
  nmax = 0;
  keflag = 0; /*I changed something here*/
  pairflag = 0;
  bondflag = 0;

  ComputeVirialcomstress2DSumChunk::init();
  ComputeVirialcomstress2DSumChunk::allocate();
}

/* ---------------------------------------------------------------------- */

ComputeVirialcomstress2DSumChunk::~ComputeVirialcomstress2DSumChunk()
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

void ComputeVirialcomstress2DSumChunk::allocate()
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
  if (maxchunk == 0) error->all(FLERR, "No chunks found for virialcomstress2dsum/chunk");

  memory->create(stress,4,"virialcomstress2dsum/chunk:stress");
  memory->create(stress_all,4,"virialcomstress2dsum/chunk:stress_all");
  memory->create(massproc,maxchunk,"virialcomstress2dsum/chunk:massproc");
  memory->create(masstotal,maxchunk,"virialcomstress2dsum/chunk:masstotal");
  memory->create(xcm,maxchunk,2,"virialcomstress2dsum/chunk:xcm");
  memory->create(xcm_all,maxchunk,2,"virialcomstress2dsum/chunk:xcm_all");
  memory->create(stress_comp, 4, "virialcomstress2dsum/chunk:stress_comp");
  vector=stress_all;
}

void ComputeVirialcomstress2DSumChunk::init_list(int /*id*/, NeighList *ptr) {
  list = ptr;
}

void ComputeVirialcomstress2DSumChunk::init() {
  //neighbor_request_id = modify->request_neighbor(this);
  if (force->pair == nullptr)
    error->all(FLERR, "No pair style is defined for compute stress/spherical");
  if (force->pair->single_enable == 0)
    error->all(FLERR, "Pair style does not support compute stress/spherical");
  neighbor->add_request(this, NeighConst::REQ_OCCASIONAL);
}

/* ---------------------------------------------------------------------- */

void ComputeVirialcomstress2DSumChunk::setup()
{
  // one-time calculation of per-chunk mass
  // done in setup, so that ComputeChunkAtom::setup() is already called

  if (firstflag && cchunk->idsflag == ONCE) {
    compute_vector();
    firstflag = massneed = 0;
  }
}

/* ---------------------------------------------------------------------- */

void ComputeVirialcomstress2DSumChunk::compute_vector()
{
  ComputeChunk::compute_vector();
  size_array_rows = nchunk;

  invoked_vector = update->ntimestep;
    
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

  for (int c = 0; c < 4; c++) {
    stress[c]      = 0.0;
    stress_all[c]  = 0.0;
    stress_comp[c] = 0.0;   // zero the Kahan compensator
  }

  for (int i = 0; i < nchunk; i++) {
    for (int d = 0; d < 2; d++) {
      xcm[i][d]=0.0;
      xcm_all[i][d]=0.0;
    }
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
      if (massneed) massproc[index] += massone;
    }

  MPI_Allreduce(&xcm[0][0], &xcm_all[0][0], 2 * nchunk, MPI_DOUBLE, MPI_SUM, world);
  if (massneed) MPI_Allreduce(massproc, masstotal, nchunk, MPI_DOUBLE, MPI_SUM, world);

  for (int i = 0; i < nchunk; i++) { 
    if (masstotal[i] > 0.0) {
      xcm_all[i][0] /= masstotal[i];
      xcm_all[i][1] /= masstotal[i];
    } else
      xcm_all[i][0] = xcm_all[i][1] = 0.0;
  }

  double Lx = domain->boxhi[0] - domain->boxlo[0];
  double Ly = domain->boxhi[1] - domain->boxlo[1];

  nktv2p = force->nktv2p;
  inv_volume = 1.0 / (Lx*Ly);

  scale=nktv2p*inv_volume;

  // Get forces for each atom on the processor and bin them into their corresponding intermolecular interactions
  for (ii = 0; ii < inum; ii++) {

    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];

    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    // Check if atom i is in the group
    if (!(mask[i] & groupbit)) continue;

    // Molecule index (is shared between processors when using atom style "full")
    int i_index;
    i_index= atom->molecule[i] - 1;
    double com_i[2];

    //Get COM for molecule i
    for (int d = 0; d < 2; ++d) {
      com_i[d] = xcm_all[i_index][d];
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

      rsq = delx * delx + dely * dely;
      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {
        eng = force->pair->single(i, j, itype, jtype, rsq, factor_coul, factor_lj, fpair);

        if (masstotal[i_index] > 0.0) {
          
          double com_j[2], delcom[3];

          // Get COM for molecule j
          for (int d = 0; d < 2; ++d) {
            com_j[d] = xcm_all[j_index][d];
          }

          // Get the intermolecular distance
          for (int d = 0; d < 2; ++d) {
            delcom[d] = com_i[d] - com_j[d];
          }
          //zero out the third component of delcom
          delcom[2] = 0.0;

          //Get the minimum image distance between the two molecules
          domain->minimum_image_big(delcom);

          //print chunk
          // printf("time: %d chunk: %d to chunk: %d\n", update->ntimestep, i_index, j_index);
          // printf("  delcom: %f %f %f\n", delcom[0], delcom[1], delcom[2]);
          // printf("  com_i: %f %f %f\n", com_i[0], com_i[1], com_i[2]);
          // printf("  com_j: %f %f %f\n", com_j[0], com_j[1], com_j[2]);

          //Reminder: delx * fpair is the interatomic force and dx is the intermolecular distance
          //Normal summation
          // stress[0] += delx*fpair*delcom[0]; //Txx
          // stress[1] += delx*fpair*delcom[1]; //Txy
          // stress[2] += delx*fpair*delcom[2]; //Txz
          // stress[3] += dely*fpair*delcom[0]; //Tyx
          // stress[4] += dely*fpair*delcom[1]; //Tyy
          // stress[5] += dely*fpair*delcom[2]; //Tyz
          // stress[6] += delz*fpair*delcom[0]; //Tzx
          // stress[7] += delz*fpair*delcom[1]; //Tzy
          // stress[8] += delz*fpair*delcom[2]; //Tzz

          //Kahan summation
          double terms[4] = {
            delx*fpair*delcom[0],  // T_xx
            delx*fpair*delcom[1],  // T_xy
            dely*fpair*delcom[0],  // T_yx
            dely*fpair*delcom[1],  // T_yy
          };
          
          // Kahan‐compensated update for all 9 components
          for (int c = 0; c < 4; c++) {
            double y = terms[c] - stress_comp[c];
            double t = stress[c] + y;
            stress_comp[c] = (t - stress[c]) - y;
            stress[c]      = t;
          }
        }
      }
    }
  }

  //Gather all of the stress data
  MPI_Allreduce(stress, stress_all, 4, MPI_DOUBLE, MPI_SUM, world);

  for (int i = 0; i < 4; i++) {
    stress_all[i] = scale * stress_all[i];
  }

  // printf("stress_all: %f %f %f %f %f %f %f %f %f\n", stress_all[0], stress_all[1], stress_all[2], stress_all[3], stress_all[4], stress_all[5], stress_all[6], stress_all[7], stress_all[8]);
}

double ComputeVirialcomstress2DSumChunk::memory_usage()
{
  double bytes = (double) (4*2+maxchunk*2+maxchunk*2*2)* sizeof(double);
  return bytes;
}