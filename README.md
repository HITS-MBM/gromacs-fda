[![Build
Status](https://jenkins.h-its.org/buildStatus/icon?job=MBM/HITS-MBM/gromacs-fda/release-2019-fda)](https://jenkins.h-its.org/job/MBM/job/HITS-MBM/job/gromacs-fda/job/release-2019-fda/)
[![Build
Status](https://travis-ci.org/HITS-MBM/gromacs-fda.svg?branch=release-2019-fda)](https://travis-ci.org/HITS-MBM/gromacs-fda)


# Force Distribution Analysis (FDA)

## Background

Calculations of internal forces and stress distributions have long been used by
mechanical engineers in designing strong structures, from tall sky-scrapers and
kilometer-long bridges to Formula 1 cars and aircraft wings. The same kind of
analysis, but applied to atomic and molecular level, could shed light onto the
mechanical stability of (bio)molecules, allosteric mechanisms and signal
propagation in cells or mechanically induced processes such as blood
coagulation.

Molecular dynamics (MD) simulations apply Newton's equations of motion to
atomic systems. During each integration time step, atoms interact with each
other through pairwise forces. All pairwise forces acting on an atom are summed
up and the resulting force determines the displacement of the atom in the next
time step. These displacements are then typically used to analyze MD
simulations.

The analysis of individual atomic pairwise forces, which we call Force
Distribution Analysis, can offer a much more sensitive view of a (bio)molecular
system, be it in a stable state or during a transition between states. Such
analysis would be able to explain, for example, the functioning of an atomic
level Newton's cradle, where the application of force on an atom does not
translate directly into a displacement of all spheres but only of the outer
ones. Similarly, atoms can propagate forces in the absence of any significant
atomic displacements, for example through a protein's rigid core.

## Citations

If you use the FDA code, please cite: http://dx.doi.org/10.1186/2046-1682-6-5


# Installation

GROMACS-FDA base on the official GROMACS. The actual molecular dynamic
simulation will be performed with the official GROMACS, whereas the force
distribution analysis will be performed with GROMACS-FDA. To avoid
inconsistencies between different versions of GROMACS and GROMACS-FDA it is
crucial that the gromacs version numbers ar identical. Please find the
available FDA releases at https://github.com/HITS-MBM/gromacs-fda/releases

## Installation of GROMACS
The default installation of GROMACS can be build with following steps:

```
git clone https://github.com/gromacs/gromacs.git -b <gromacs-version>
cd gromacs
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=<installation path> ..
make -j <number of cores>
make check
make install
```

Please look at the official GROMACS documentation for details of the
installation procedure.

## Installation of GROMACS-FDA
GROMACS-FDA is also available on GitHub and can be build with following steps:

```bash
git clone https://github.com/HITS-MBM/gromacs-fda.git -b <gromacs-version>-fda<fda-version>
cd gromacs-fda
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=<installation path> -DGMX_BUILD_FDA=ON -DGMX_DEFAULT_SUFFIX=OFF -DGMX_BINARY_SUFFIX=_fda -DGMX_SIMD=NONE -DGMX_BUILD_UNITTESTS=ON -DGMX_GPU=OFF ..
make -j <number of cores>
make check
make install
```

This will generate an executable called gmx_fda which contains all the
functionality described in this document. To compile a double precision version
use following cmake parameters:

```
-DGMX_DOUBLE=ON -DGMX_BINARY_SUFFIX=_fda_d
```

which will generate an executable called gmx_fda_d. In the rest of this
document only gmx_fda (default name for a single precision version) is
mentioned, but everything applies to the double precision version as well.

## Installation of VMD visualization plugins
The installation is done by
unpacking the archive into the VMD plugins directory; normally this is
plugins/noarch/tcl. This operation will create 2 directories called
pf_loaduser1.0 and pf_draw1.0 which contain the corresponding Tcl files. To
test the installation, start VMD then at the console or in the Tk console type:

```
package require pf_loaduser package require pf_draw
```

# Principles
## Residue-based operation
For residue-based operation, the following steps are done in order:
1. building of an atom to residue correspondence table
1. if an atom pair is part of those selected for pairwise forces calculations,
   the residue numbers of the 2 atoms are obtained from the correspondence table
1. if the 2 residue numbers are equal (2 different atoms from the same residue),
   nothing is done because it makes no sense to calculate the interaction
   of a residue with itself and the code goes back to step 2 with a different atom pair
1. the atomic vector force is added to the vector force corresponding
   to the residue pair
   
So in the end, each residue pair will have a vector force calculated as:

where i is an atom of residue ri and j is an atom of residue rj, with ri and rj
being different.

Once calculated, the residue-based pairwise forces can be
written to output files just like the atom-based ones – as vector or scalar
values – or can be further used to calculate a per-residue
sum/average/minimum/maximum value.

## Residue (re)numbering
For GROMACS, residues are of secondary importance
during the actual MD run – they are useful only to write out a structure
file (f.e. PDB). The correspondence between numbering of atoms and residues is
made based on the information in the structure file where residues numbering
might not start from 1 and might not be contiguous – to deal with all these
cases, by default GROMACS doesn't make any attempt to renumber residues and
writes back whatever residue number was assigned to an atom in the input
structure file.

The PF2 code uses these residue numbers to identify which
atoms belong to which residues and to output residue-based data. Dealing with
residue numbers not starting from 1 or not being contiguous is easy – just
fill in the missing numbers starting from 1 and have the output contain zeros
for them. But residue numbers can also appear multiple times, f.e. if the
structure file contains several chains (say a protein, a ligand and solvent)
and in each of them residues start being numbered at 1 or when there are many
residues (more than 9999 which is the maximum number in a PDB file) and their
number wraps around. To get around this, GROMACS can do residue renumbering
which numbers residues with increasing numbers starting from 1. The PF2 code
can detect residue number collisions and, by default, automatically switches to
renumbering residues; this behavior can be controlled through a .pfi file
option.

Renumbering residues influences both the way the residue data is
stored in memory and the way it is written in the output files. The information
whether renumbering has taken place in the PF2 code (be it automatically or
enforced through the .pfi file option) has to be passed on to the VMD plugin;
failure to do so will result in wrong representations!

## Pairwise forces decomposition for 3- and 4-body potentials
For potentials
involving only two atoms (bonds, Coulomb, Lennard-Jones), the pairwise force is
the same as the atomic force, directly calculated from the potential. For
potentials involving more than two atoms (angles, dihedral angles – proper
and improper, as well as cross bond-bond and cross bond-angle), the atomic
forces need to be decomposed to reflect the interactions between each two
atoms. The vector sum of pairwise forces which result from this decomposition
should equal the atomic forces. The decomposition implemented in the PF2 code
is detailed in Appendix A.

# Running
## Additional mdrun options
FDA adds a few command line options to gmx_fda mdrun:

  * **-pfi input.pfi** – specifies an input file which controls how the
     code should run
  * **-pfn index.ndx** – is an index file containing one or
     more groups of atoms which should be used for pairwise forces calculations
  * **-pfa output.pfa** – output file, containing pairwise forces for atoms; the
     content and format depend on the choices in the .pfi file
  * **-pfr output.pfr** – output file, containing pairwise forces for residues; the content and
     format depend on the choices in the .pfi file
  * **-psa output.psa** – output file, containing punctual stress for atoms;
     the content and format depend on the choices in the .pfi file
  * **-psr output.psr** – output file, containing punctual stress for residues;
     the content and format depend on the choices in the .pfi file
  * **-vsa output.vsa** – output file, containing virial stress
     for atoms; the content and format depend on the choices in the .pfi file
  * **-vma output.vma** – output file, containing von Mises virial stress for atoms;
     the content and format depend on the choices in the .pfi file

`gmx_fda mdrun` should always be run with the -rerun option, requiring that an
initial trajectory is produced beforehand; apart from the PF2 specific options
described above, the options passed to `gmx_fda mdrun` should be the same as the
ones used when producing the original trajectory, with one exception: a single
CPU run should be enforced, as the PF2 code is not (yet) parallel. Please note
that a non-MPI version of GROMACS 4.5 runs by default as threaded which still
counts as parallel, so use `'-nt 1'` to enforce a non-parallel run:

```
gmx_fda mdrun -nt 1 -rerun traj.trr -pfi input.pfi -pfn index.ndx -pfa
output.pfa -pfr output.pfr <other options>
```

## The input file (.pfi)
The input file controls how the PF2 code should run;
it also specifies the content and format of the output files (.pfa and .pfr).
Like the .mdp file which GROMACS uses as input, the .pfi file is a simple text
file containing statements of form `'keyword = value'`, with comments
started by a semi-colon `';'`. The .pfi file is read
directly by gmx_fda mdrun, there is no need for a step equivalent to grompp.

The .pfi file keywords are:

**onepair** – describes how the pairwise forces between the same pair of atoms
i and j is handled. Typically, a forcefield can allow either several types of
non-bonded interactions or several types of bonded interactions between the
same pair of atoms. F.e. two distant atoms could interact through both Coulomb
and Lennard-Jones interactions, while two bonded atoms could interact through a
bond, several angle and several dihedral angle interactions. There are 2
possible values for this keyword:

*summed* – all pairwise forces between the
same pair of atoms are summed; if atoms i and j interact, the code stores
exactly one value of the pairwise force between them

*detailed* – all pairwise forces of the same type (f.e. angle interactions) between the same pair of atoms are stored separately; if atoms i and j interact, the code stores at least one, possibly many, pairwise forces between them, corresponding to
each type of interaction that exists between them.

Using *detailed* can lead to
a significant increase in the amount of memory needed compared with the summed
option. *detailed* should probably only be used in the rare cases where
simultaneous output of pairwise forces separated by potential type should be
obtained, like forcefield development/debugging. The default is *summed*. For
example, the following set of options will compute a sum of all interactions of
the same type (say dihedral angles) between the atoms i and j:

```
onepair = detailed
type = dihedral
```

**group1** and **group2** – name the groups which will be used for pairwise forces
calculations; if not present, they default to “Protein”. The 2 groups
are equivalent, there is no difference in how they are treated internally. The
groups names should exist in the index file specified with the -pfn command
line option; it's not necessary to have a separate index file for this purpose,
different from the one used for the **grompp** step of the original run – only
the groups named in group1 and group2 will be loaded, other groups will be
ignored. The 2 groups define which atoms are used for pairwise forces
calculations: only interactions for which an atom is part of group1 and the
other atom is part of group2 are being considered. It's fine to have one or
more atoms in both groups at the same time. It's also fine to have only one
group defined in the index file and the name of this group specified for both
group1 and group2 – this is equivalent to the original PF implementation
[1] where only one group was defined. Having 2 groups allows a greater
flexibility and efficiency: imagine a case where the interactions between a
ligand and a protein are of interest – group1 could be defined to contain
the protein and group2 could be defined to contain the ligand and the pairwise
interactions will only be computed for interactions between them; using only
one group as in the original PF implementation would mean including both the
protein and the ligand in this group and obtaining potentially useless
interactions protein-protein and ligand-ligand.

group1 and group2 always
contain atom indices, even when requesting only residue-based calculations.
Specifying all atoms of a residue in one or both groups remains the
responsibility of the user. If only some of the atoms of a residue are
specified and residue-based calculations are requested, the results
corresponding to that residue will only contain the interactions involving the
specified atoms. The code doesn't automatically fill in the missing atoms to
obtain full residues. This is done for flexibility: e.g. if only the
residue-based backbone interactions are of interest, the index file should
contain only backbone atoms and residue-based calculations should be requested.

**type** – specifies which types of interactions are of interest. It can have
one or more of the values:

```
bond angle dihedral polar coulomb lj nb14 bonded nonbonded all
```

*nb14* is defined as the combination of LJ14 and Coulomb14
interactions. LJ14 and Coulomb14 cannot be obtained separately.

*bonded* is defined as the combination of *bond*, *angle* and *dihedral*.

*nonbonded* is defined as the combination of coulomb, *lj* and *nb14*.

*all* is defined as the combination of *bonded*, *polar* and *nonbonded*.

If no type is specified, `'type = all'` is assumed.

GROMACS contains several functions which calculate potentials and forces for
bonds, angles and dihedral angles. Which ones are used depends on the
forcefield, choices in the .mdp file (f.e. “morse=yes”) and user
modifications to the topology files (f.e. to add extra harmonic potentials).
The PF2 code treats them based on the number of particles involved: all bonded
interactions involving 2 particles are considered bonds, all bonded
interactions involving 3 particles are considered angles and all bonded
interactions involving 4 particles are considered dihedral angles. This makes
impossible to differentiate between several types of interactions involving the
same number of atoms; f.e. in a protein topology using the OPLS-AA forcefield,
the dihedral angles are modeled using a Ryckaert-Bellemans form while the
improper dihedral angles are modeled using a classical dihedral angle form
– both these types of interactions will be selected by *dihedral*.

When requesting residue-based calculations, a pairwise residue-based interaction of
a certain type will exist if at least one atom from one residue and one atom
from the second residue are involved in an interaction of that type.

**atombased** and **residuebased** – specify how the pairwise forces are handled
and saved. The 2 keywords are independent: one can be set to “no” to
specify that there is no interest in storing, handling and saving that type of
data; the default is “no”. The following table shows the possible
values and their effect:

option | calculation | output data | output file 
-------|-------------|-------------|------------
no | Pairwise forces are not stored, handled or saved. | none |.pfa or .pfr
pairwise_forces_vector | Pairwise forces are saved as vectors. | vector | .pfa or .pfr
pairwise_forces_scalar | Pairwise forces are saved as signed scalars. | signed scalar | .pfa or .pfr
punctual_stress | For each atom/residue, the punctual stress is saved. Only works with `onepair=summed`. | scalar | .psa or .psr
virial_stress | For each atom, the virial stress is saved. Only works with atombased and `onepair=summed`. | tensor | .vsa
virial_stress_von_mises | For each atom, the von Mises virial stress is saved. Only works with `atombased` and `onepair=summed`. | scalar | .vma
compat_bin (deprecated option) | Pairwise forces are saved as signed scalars into a binary format compatible to the original PF implementation. Only works with `onepair=summed`. | signed scalar | .pfa or .pfr
compat_ascii (deprecated option) | Pairwise forces are saved as signed scalars into a text format compatible to the original PF implementation. Only works with `onepair=summed`. | signed scalar | .pfa or .pfr

**vector2scalar** – controls how the transformation between the vector pairwise
forces and their scalar representation is made. It can have the following
values:

*projection* – the scalar is the length of the projection of the
force vector on the vector determined by the position of the 2 atoms or the
position of the COM of the 2 residues.

*norm* – the scalar is the norm of the vector and the sign is determined based on the cosine of the angle between
the vector force and the vector determined by the positions of the 2 atoms or
the position of the COM of the 2 residues:

  * if it's positive, the vector force is oriented within -90 to 90 degrees with
    respect to the atoms/residues - considered to be same direction = attractive = negative
  * if it's negative, the vector force is oriented within 90 to 270 degrees with
    respect to the atoms/residues - considered to be in opposite direction =
    repulsive = positive
  * if it's zero, the vector force is oriented at exactly 90 or 270 degrees
    with respect to the atoms/residues - the force can be considered neither
    attractive nor repulsive, so it's set to zero
    
The default value is `'norm'`.

**residuesrenumber** – controls residue renumbering.
The possible values are:

*auto* – keep the residue numbering as in the structure file if there are not
collisions and do renumbering if collisions are present

*yes* – do residue renumbering independent of the presence of residue number collisions 

*no* – don't do residue renumbering. In case a collision is detected, a warning will
be written but renumbering will not take place. Using this option is dangerous
as residues cannot be uniquely identified! If residue number collisions are
present, pairwise forces of residues with the same residue number will be
summed up together – which is probably not desired.

This option is taken into account only if `residuebased` is set to something else
than `'no'`. The default value is `'auto'`.

**time_averages_period** – is the number of input trajectory frames between
averaging scalar data. This can only be used together with atombased and/or
residuebased set to `pairwise_forces_scalar`, `compat_bin` or `compat_ascii`. The
default value is 1 – no averaging is done. The output will contain one
frame for this number of input frames; if the input frames are not a multiple
of this number, a last frame is written after averaging over the remainder of
input frames. A value of 0 specifies that averaging should be done over all
frames present in the input trajectory; the output will then contain a single
frame.

**no_end_zeros** – controls the presence of zeros at the end of per
atom/residue data. The possible values are:

*yes* – suppress all zeros at the end of per atom/residue data

*no* – do not suppress zeros at the end of the per atom/residue data

The default value is no. This option is useful when
computing per atom/residue data for systems composed of the molecule(s) of
interest followed by lots of water; in this case, the per atom/stress only has
non-zero values for the molecule(s) of interest followed by lots of zeros,
increasing the storage space without bringing any useful information.

**energy_grp_exclusion** – controls the usage of energy group exclusions to
speedup the md rerun:

*yes* – unneeded forces between atoms not in the FDA groups will be not calculated

*no* – all forces will be calculated

The default is `'yes'`. The option is only for the case that the exclusion is not
working correctly for debugging.

**normalize_punctual_stress_per_residue** – Divide the per residue punctual
stress by the number of atoms in the residue. The possible values are:

*yes* – per residue punctual stress will be normalized

*no* – per residue punctual stress will not be normalized

The default is `'no'`.

**binary_result_file** – Store the FDA result file in a binary format. The
possible values are:

*yes* – FDA result file in binary format

*no* – FDA result file in text based format

The default is `'no'`.

**threshold** – is a real number for neglecting all pairwise forces which are
smaller than the threshold. The unit is kJ/mol/nm and the default value is
`'1e-7'`.

## Input file examples
An example .pfi file for only residue-based output in
text compatibility format, taking into consideration all interaction types and
with time averaging over all input frames:

```
; onepair could be detailed or summed
onepair = summed

; group1 and group2 as defined in the -pfn file
; if not defined, defaults to 'Protein'
group1 = Protein
group2 = Ligand

; no interest in atom-based information
atombased = no

; output residue-based information in compatibility format
residuebased = compat_ascii

; interactions type could be one of more of:
; bond angle dihedral polar coulomb lj nb14 bonded nonbonded all
type = all

; scalar summation over all steps time_averages_period = 0
```

Another example .pfi file, for only atom-based output containing a per-atom sum
of all bonded interactions:

```
; onepair could be detailed or summed
onepair = summed

; group1 and group2 as defined in the -pfn file
; if not defined, defaults to 'Protein'
group1 = Protein
group2 = Ligand

; sum of scalar pairwise forces per atom
atombased = per_atom_sum

; no interest in residue-based information
residuebased = no
 
; interactions type could be one of more of:
; bond angle dihedral polar coulomb lj nb14 bonded nonbonded all
type = bonded
```

Yet another example .pfi file, for atom-based output containing per-atom
maximum scalar pairwise interactions and residue-based output containing scalar
pairwise forces:

```
; onepair could be detailed or summed
onepair = summed 

; group1 and group2 as defined in the -pfn file
; if not defined, defaults to 'Protein'
group1 = Protein
group2 = Ligand

; maximum scalar pairwise force per atom
atombased = per_atom_max

; scalar pairwise forces between residues
residuebased = scalar

; interactions type could be one of more of:
; bond angle dihedral polar coulomb lj nb14 bonded nonbonded all
type = all
```

And finally an example .pfi file for only residue-based output in text
compatibility format, taking into consideration all interaction types and
performing residue renumbering, with time averaging over all input frames:

```
; onepair could be detailed or summed
onepair = summed

; group1 and group2 as defined in the -pfn file
; if not defined, defaults to 'Protein'
group1 = Protein
group2 = Ligand

; no interest in atom-based information
atombased = no

; output residue-based information in compatibility format
residuebased = compat_ascii

; ask for residues to be renumbered - starting from 1 and continuous
residuesrenumber = yes

; interactions type could be one of more of:
; bond angle dihedral polar coulomb lj nb14 bonded nonbonded all
type = all

; scalar summation over all steps
time_averages_period = 0
```

# Output file formats
The file formats are independent on whether they contain
atom-based or residue-based data (.pfa or .pfr file respectively). The
compatibility formats (atombased/residuebased set to compat_bin or
compat_ascii) are supposed to be as close as possible to the ones obtained from
the original PF implementation. The following mapping between PF2 and PF
interaction types is used:

PF2          | PF
------------ | -------------
bond         | iBond
angle        | iAngle
dihedral     | iDihedral
polar        | iPolar
lj           | iLJ
coulomb      | iCoul
nb14         | iDihedral

However, the compatibility formats contain further data: the original PF
implementation made the approximation that for an angle or dihedral angle only
the extreme atoms interact; PF2 correctly treats interactions between all atoms
in these cases, leading to more interactions present in the output file. These
formats are supposed to be used only to facilitate reading the data into the
current R FDA tools. The output formats containing pairwise forces
(atombased/residuebased set to scalar or vector) are composed of a set of
frames, each frame containing a variable set of interactions. A frame starts
with a line containing:

```
frame <frame_number>
```

where `<frame_number>` is an integer starting from 0, separated by a whitespace
from the word “frame”. This line is followed by zero or more lines of
form:

```
i j scalar_force interaction_type
```

for atombased/residuebased set to scalar or

```
i j force_x force_y force_z interaction_type
```

for atombased/residuebased set to vector. i and j are integer numbers which
represent the atom indeces (0-based); scalar_force is a floating point number
which represents the magnitude of the pairwise force vector; force_x, force_y
and force_z are floating point numbers which represent the X, Y and Z
components of the pairwise force vector; interaction_type is an integer number
which represents the type of interaction – for detailed interactions, this
is a constant defined in include/pf_interactions.h file and might change
depending on the PF2 version (the existing constants will normally be kept, new
ones will be added with higher values); for summed interactions, this is a sum
of the constants corresponding to all interactions present between atoms i and
j. All floating point numbers are output using the C printf format specifier
“%e”. All values on one line are separated by whitespaces.

The output
formats containing per-atom/per-residue data are composed of one line per frame
with values on the line separated by whitespaces. All values are floating point
numbers which represent the sum/average/minimum/maximum scalar pairwise
interaction; they are written using the C printf format specifier “%e”.
The first value on the line corresponds to the first atom/residue, the second
value on the line corresponds to the second atom/residue, etc.; on each line
there are as many per-atom values as the number of atoms and as many
per-residue values as the number of residues. If an atom/residue does not have
any pairwise interactions (f.e. because it is located further than the cutoff
from other atoms, because it doesn't have any interactions of the type
specified with type or because it was included neither in group1 nor in
group2), the corresponding value is 0.0.

For per-residue data, residue
numbering depends on whether residue renumbering has taken place. If the
residue renumbering has not taken place, the residue numbers start from 0 and
continue until the highest number found in the input structure file. F.e. if
the input structure contains 4 residues numbered 5, 6, 8 and 9, the output will
contain 10 values: the first 5 values and the 8th value will always be zero. If
the residue renumbering has taken place, the residue numbers start from 0 and
continue until the total number of residues minus 1.

## Binary file format
The FDA result files can also written in a binary format
to reduce and speed-up the disk usage. The binary format is used by setting the
option ‘binary_result_file’ to ‘yes’. All FDA result files
(pfa, pfr, psa, psr, vsa, vma) can be converted from text-based into binary
format and vice versa by using the conversion module:

```
gmx_fda fda_convert -i <input-file> -o <output-file>
```

If the input file is text-based to output file will be binary and the other way
round.

# Analysis modules
All analysis modules are integrated within GROMACS and can be executed by using:

```
gmx_fda <module> [options]
```

The documentation of every analysis module can be printed by using:

```
gmx_fda help <module>
```

## Module fda_graph
fda_graph converts a FDA force network into a PDB or
DIMACS graph. If the optional file -ipf-diff is used the differences of the
pairwise forces will be taken. The PDB graph allows an easy visualization with
a program of your choice. Only forces larger than the -t will be considered.
The default threshold is zero. Networks must contains at least the same number
of nodes as the the min-value (default: 2). If the option -big is used, only
the biggest network in term of number of nodes will be printed. Each network
will be determined and segment names will be assign to each of them, thus
coloring them by segment id will help the analysis (32 different colors). The
Bfactor column will be used for the value of the force and helps the coloring
as a function of the force magnitude. The CONNECT header will be used to create
bonds between nodes.

## Module fda_shortest_path
fda_shortest_path calculates
the k-shortest paths between a source node and a destination node of a FDA
force network. If the optional file -ipf-diff is used the differences of the
pairwise forces will be taken. The graph will printed in the PDB-format, which
allow an easy visualization with an program of your choice. Each path will be
determined and segment names will be assign to each of them, thus coloring them
by segment id will help the analysis (32 different colors). The Bfactor column
will be used for the value of the force and helps the coloring as a function of
the force magnitude. The CONNECT header will be used to create bonds between
nodes.

## Module fda_get_stress
fda_get_stress calculates the punctual stress
by pairwise forces. If the optional file -ipf-diff is used the differences of
the pairwise forces will be taken.

## Module fda_view_stress
fda_view_stress
illustrates the punctual or von Mises virial stress of FDA as xpm or pdb-file.
The x-axis of the xpm file represent the particle number and the y-axis the
frame number. For the xpm-file the number of different colors can be set with
the -nbColors option. For the pdb file the Bfactor column will be used for the
value of the stress and helps the coloring as a function of the stress
magnitude. The pdb-file can be visualized with an program of your choice.
