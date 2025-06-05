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

ComputeStyle(translationaltemp/chunk,ComputeTranslationaltempSumChunk)

#else

#ifndef LMP_COMPUTE_TRANSLATIONALTEMP_SUM_CHUNK_H
#define LMP_COMPUTE_TRANSLATIONALTEMP_SUM_CHUNK_H

#include "compute_chunk.h"

namespace LAMMPS_NS {

class ComputeTranslationaltempSumChunk : public ComputeChunk {
 public:
  ComputeTranslationaltempSumChunk(class LAMMPS *, int, char **);
  ~ComputeTranslationaltempSumChunk() override;


  void setup() override;
  void compute_vector() override;
  double memory_usage() override;

 public:
  int keflag,pairflag,bondflag,angleflag,dihedralflag,improperflag;
  int kspaceflag,fixflag,biasflag,perchunk_flag,presschunkflag,timeflag;
  int size_perchunk_cols,invoked_perchunk,nchunk_local;
  Compute *temperature;
  char *id_temp;

  int nmax;
 private:
  double *massproc, *masstotal;
  double **vcm, **vcmall, *stress, *stress_all;
  double inv_volume,scale;

  void allocate() override;
};

}

#endif
#endif
