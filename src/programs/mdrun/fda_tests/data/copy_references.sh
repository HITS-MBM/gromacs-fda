#!/bin/bash

declare -a arr=(
  "alagly_pairwise_forces_scalar_summed"
  "alagly_pairwise_forces_scalar_summed_atom_based"
  "alagly_pairwise_forces_scalar_summed_no_residue_based"
  "alagly_pairwise_forces_scalar_detailed"
  "alagly_pairwise_forces_vector"
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
