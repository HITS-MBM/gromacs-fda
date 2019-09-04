#!/bin/bash

declare -a arr=(
  "alagly_pairwise_forces_scalar"
  "alagly_pairwise_forces_scalar_atom_based"
  "alagly_pairwise_forces_scalar_no_residue_based"
  "alagly_pairwise_forces_scalar_detailed_no_residue_based"
  "alagly_pairwise_forces_vector"
  "alagly_punctual_stress"
  "alagly_pairwise_forces_scalar_detailed_nonbonded"
  "alagly_pairwise_forces_vector_detailed_nonbonded"
  "alagly_verlet_summed_scalar"
  "alagly_verlet_pbc_summed_scalar"
  "alagly_group_excl"
  "alagly_group_excl_uncomplete_cgs"
  "alagly_pairwise_forces_scalar_all"
  "glycine_trimer_group_excl1"
  "glycine_trimer_group_excl2"
  "glycine_trimer_group_excl3"
  "glycine_trimer_group_excl4"
  "glycine_trimer_group_excl5"
  "glycine_trimer_group_excl6"
  "glycine_trimer_group_bonded_excl1"
  "glycine_trimer_virial_stress"
  "glycine_trimer_virial_stress_von_mises"
  "alagly_deprecated_keywords"
  "alagly_unknown_option"
  "vwf_a2_domain_nframes1_pairwise_forces_scalar"
  "vwf_a2_domain_nframes1_punctual_stress"
  "vwf_a2_domain_nframes10_pairwise_forces_scalar"
  "vwf_a2_domain_nframes10_punctual_stress"
)

declare -a files=(
  "pfa"
  "pfr"
  "psa"
  "psr"
  "vsa"
  "vma"
)

for i in "${arr[@]}"
do
  for j in "${files[@]}"
  do
    echo "cp $1/*/$i/fda.$j $2/$i/fda.$j.ref"
    cp $1/*/$i/fda.$j $2/$i/fda.$j.ref
  done
done
