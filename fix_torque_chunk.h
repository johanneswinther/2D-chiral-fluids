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


// ----------------------------------------------------------------------
// This fix is used to torque individual chunks in 2d, around their center of mass.
// ----------------------------------------------------------------------


#ifdef FIX_CLASS
// clang-format off
FixStyle(torque/chunk,FixTorqueChunk);
// clang-format on
#else

#ifndef LMP_FIX_TORQUE_CHUNK_H
#define LMP_FIX_TORQUE_CHUNK_H

#include "fix.h"

namespace LAMMPS_NS {

class FixTorqueChunk : public Fix {
 public:
  FixTorqueChunk(class LAMMPS *, int, char **);
  ~FixTorqueChunk() override;
  int setmask() override;
  void init() override;
  void setup(int) override;
  void min_setup(int) override;
  void post_force(int) override;
  void post_force_respa(int, int, int) override;
  void min_post_force(int) override;

 private:
  int ilevel_respa;
  double k_torque;
  char *idchunk, *idcom;

  int nchunk;

  class ComputeChunkAtom *cchunk;
  class ComputeCOMChunk *ccom;
};

}    // namespace LAMMPS_NS

#endif
#endif
