/*
 * $Id$
 * 
 *       This source code is part of
 * 
 *        G   R   O   M   A   C   S
 * 
 * GROningen MAchine for Chemical Simulations
 * 
 *               VERSION 2.0
 * 
 * Copyright (c) 1991-1999
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
 * Please refer to:
 * GROMACS: A message-passing parallel molecular dynamics implementation
 * H.J.C. Berendsen, D. van der Spoel and R. van Drunen
 * Comp. Phys. Comm. 91, 43-56 (1995)
 * 
 * Also check out our WWW page:
 * http://md.chem.rug.nl/~gmx
 * or e-mail to:
 * gromacs@chem.rug.nl
 * 
 * And Hey:
 * Green Red Orange Magenta Azure Cyan Skyblue
 */
static char *SRCID_names_c = "$Id$";

#include "typedefs.h"
#include "names.h"

/* note: these arrays should correspond to enums in include/types/enums.h */

char *eblock_names[ebNR+1]=
{
  "CGS","MOLS","SBLOCKS",NULL
};

char *epbc_names[epbcNR+1]=
{
  "xyz", "none", NULL
};

char *ens_names[enNR+1]=
{
  "Grid","Simple", NULL
};

char *ei_names[eiNR+1]=
{ 
  "md", "steep", "cg", "ld", "sd", NULL 
};

char *bool_names[BOOL_NR+1]=
{
  "FALSE","TRUE", NULL
};

char *yesno_names[BOOL_NR+1]=
{
  "no","yes", NULL
};

char *ptype_str[eptNR+1] = {
  "Atom", "Nucleus", "Shell", "Bond", "Dummy", NULL
};

char *eel_names[eelNR+1] = {
  "Cut-off", "Reaction-Field", "Generalized-Reaction-Field",
  "PME", "Ewald", "PPPM", "Poisson", "Switch", "Shift", "User", NULL
};

char *evdw_names[evdwNR+1] = {
  "Cut-off", "Switch", "Shift", "User", NULL
};

char *eshake_names[estNR+1] = {
  "Lincs", "Shake", NULL
};

char *egrp_nm[egNR+1] = { 
  "Coul-SR","LJ","Buck", "Coul-LR", "LJ-LR", "Coul-14", "LJ-14", NULL
};

char *epcoupl_names[epcNR+1] = {
  "No", "Isotropic", "Semiisotropic", "Anisotropic", "Surface-Tension", NULL
};

char *edisre_names[edrNR+1] = {
  "No", "Simple", "Ensemble", NULL
};

char *edisreweighting_names[edrwNR+1] = {
  "Equal", "Conservative", NULL
};

char *enbf_names[eNBF_NR+1] = {
  "", "LJ", "Buckingham", NULL
};

char *ecomb_names[eCOMB_NR+1] = {
  "", "Arithmetic", "Geometric", "ArithSigEps", NULL
};

char *gtypes[egcNR+1] = {
  "T-Coupling", "Energy Mon.", "Acceleration", "Freeze",
  "User1", "User2", "User3", "XTC", NULL
};

char *efep_names[efepNR+1] = {
  "no", "yes", NULL
};
