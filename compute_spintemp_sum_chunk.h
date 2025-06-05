/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.

   This file is part of custom LAMMPS code for simulations of chiral matter and is based on compute_angmom_chunk.h.
    author: Johannes Winther
    e-mail: johannesiwk@gmail.com
    github: https://github.com/mandadapu-group/2D-chiral-fluids
------------------------------------------------------------------------- */

#ifdef COMPUTE_CLASS
// clang-format off
ComputeStyle(spintempsum/chunk,ComputeSpintempSumChunk);
// clang-format on
#else

#ifndef LMP_COMPUTE_SPINT_SUM_CHUNK_H
#define LMP_COMPUTE_SPINT_SUM_CHUNK_H

#include "compute_chunk.h"

namespace LAMMPS_NS {

class ComputeSpintempSumChunk : public ComputeChunk {
 public:
  ComputeSpintempSumChunk(class LAMMPS *, int, char **);
  ~ComputeSpintempSumChunk() override;
  void compute_vector() override;

  double memory_usage() override;

 private:
  double *massproc, *masstotal;
  double **com, **comall;
  double **angmom, **angmomall;
  double *spintempall;
  double inv_vol, scale, avg_mass, rho_bead;

  void allocate() override;
};
}    // namespace LAMMPS_NS
#endif
#endif
