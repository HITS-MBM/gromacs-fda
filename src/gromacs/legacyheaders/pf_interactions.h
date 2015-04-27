/* 
 * Copyright Bogdan Costescu 2010-2012
 */

#ifndef pf_interactions_h
#define pf_interactions_h

/*
 * The following are #define-d rather than enum-ed, so that they can be used
 * for bitwise operations.
 */
#define PF_INTER_NONE		0
#define PF_INTER_BOND		0x01
#define PF_INTER_ANGLE		0x02
#define PF_INTER_DIHEDRAL	0x04
#define PF_INTER_POLAR		0x08
#define PF_INTER_COULOMB	0x10
#define PF_INTER_LJ			0x20
#define PF_INTER_NB14		0x40
/* defines for easy tests of bonded/nonbonded/all interactions */
#define PF_INTER_BONDED		PF_INTER_BOND + PF_INTER_ANGLE + PF_INTER_DIHEDRAL
#define PF_INTER_NONBONDED	PF_INTER_COULOMB + PF_INTER_LJ + PF_INTER_NB14
#define PF_INTER_ALL		PF_INTER_BONDED + PF_INTER_POLAR + PF_INTER_NONBONDED

/* used in pf_utils to convert the user provided string to a mask */
typedef struct {
  int val;
  char *str;
} t_pf_interactions_type;

#endif  /* pf_interactions_h */
