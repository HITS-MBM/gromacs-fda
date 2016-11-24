/*
 * Handling of arrays of pairwise forces.
 *
 * Copyright Bogdan Costescu 2010-2012
 */

#include "fda.h"
#include "pf_array.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

using namespace fda;

void pf_atom_add_bonded_nocheck(FDA *fda, int i, int j, int type, rvec force)
{
  fda->add_bonded_nocheck(i, j, type, force);
}

void pf_atom_add_nonbonded_single(FDA *fda, int i, int j, int type, real force, real dx, real dy, real dz)
{
  fda->add_nonbonded_single(i, j, type, force, dx, dy, dz);
}

void pf_atom_add_nonbonded(FDA *fda, int i, int j, real pf_coul, real pf_lj, real dx, real dy, real dz)
{
  fda->add_nonbonded(i, j, pf_coul, pf_lj, dx, dy, dz);
}

void pf_atom_virial_add(FDA *fda, int ai, tensor v, real s)
{
  fda->add_virial(ai, v, s);
}

void pf_atom_virial_bond(FDA *fda, int ai, int aj, real f, real dx, real dy, real dz)
{
  fda->add_virial_bond(ai, aj, f, dx, dy, dz);
}
