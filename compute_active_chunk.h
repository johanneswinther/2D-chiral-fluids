/* -*- c++ -*- ----------------------------------------------------------
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

#ifdef COMPUTE_CLASS
// clang-format off
ComputeStyle(active/chunk,ComputeActiveChunk);
// clang-format on
#else

#ifndef LMP_COMPUTE_ACTIVE_CHUNK_H
#define LMP_COMPUTE_ACTIVE_CHUNK_H

#include "compute_chunk.h"

namespace LAMMPS_NS {

class ComputeActiveChunk : public ComputeChunk {
 public:
  ComputeActiveChunk(class LAMMPS *, int, char **);
  ~ComputeActiveChunk() override;
  void compute_vector() override;

  double memory_usage() override;

 private:
  double *massproc, *masstotal;
  double **com, **comall;
  double **stress, **stressall;
  double k_torque;
  double *summedstress;
  double inv_vol;
  double Lx, Ly;
  double fx, fy;
  double nktv2p;

  void allocate() override;
};

}    // namespace LAMMPS_NS

#endif
#endif
