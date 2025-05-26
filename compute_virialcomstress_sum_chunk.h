/* -*- c++ -*- ----------------------------------------------------------
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

#ifdef COMPUTE_CLASS

ComputeStyle(virialcomstresssum/chunk,ComputeVirialcomstressSumChunk)

#else

#ifndef LMP_COMPUTE_VIRIALCOMSTRESS_SUM_CHUNK_H
#define LMP_COMPUTE_VIRIALCOMSTRESS_SUM_CHUNK_H

#include "compute_chunk.h"

namespace LAMMPS_NS {

class ComputeVirialcomstressSumChunk : public ComputeChunk {
 public:
  ComputeVirialcomstressSumChunk(class LAMMPS *, int, char **);
  ~ComputeVirialcomstressSumChunk() override;

  void setup() override;
  void compute_vector() override;
  double memory_usage() override;
  void init_list(int, class NeighList *);
  void init();

 public:
  int keflag,pairflag,bondflag,angleflag,dihedralflag,improperflag;
  int kspaceflag,fixflag,biasflag,perchunk_flag,presschunkflag,timeflag;
  int size_perchunk_cols,invoked_perchunk,nchunk_local;
  Compute *temperature;
  char *id_temp;

  int nmax;
 private:
  class NeighList *list;  // Store neighbor list pointer
  double *massproc, *masstotal;
  double *stress,*stress_all;
  double **xcm, **xcm_all;
  double inv_volume,nktv2p, scale;
  int neighbor_request_id;

  void allocate() override;
};

}

#endif
#endif
