Force Distribution Analysis (FDA)
=================================
Jenkins: [![Build Status](https://jenkins.h-its.org/buildStatus/icon?job=MBM/HITS-MBM/gromacs-fda/master-fda)](https://jenkins.h-its.org/job/MBM/job/HITS-MBM/job/gromacs-fda/job/master-fda/)

Travis: [![Build Status](https://api.travis-ci.org/HITS-MBM/gromacs-fda.svg?branch=master-fda)](https://travis-ci.org/HITS-MBM/gromacs-fda)

Background
----------

Calculations of internal forces and stress distributions have long been used by mechanical engineers in designing strong structures, from tall sky-scrapers and kilometer-long bridges to Formula 1 cars and aircraft wings. The same kind of analysis, but applied to atomic and molecular level, could shed light onto the mechanical stability of (bio)molecules, allosteric mechanisms and signal propagation in cells or mechanically induced processes such as blood coagulation.

Molecular dynamics (MD) simulations apply Newton's equations of motion to atomic systems. During each integration time step, atoms interact with each other through pairwise forces. All pairwise forces acting on an atom are summed up and the resulting force determines the displacement of the atom in the next time step. These displacements are then typically used to analyze MD simulations.

The analysis of individual atomic pairwise forces, which we call Force Distribution Analysis, can offer a much more sensitive view of a (bio)molecular system, be it in a stable state or during a transition between states. Such analysis would be able to explain, for example, the functioning of an atomic level Newton's cradle, where the application of force on an atom does not translate directly into a displacement of all spheres but only of the outer ones. Similarly, atoms can propagate forces in the absence of any significant atomic displacements, for example through a protein's rigid core.

Installation and usage
----------------------

Please find the GROMACS-FDA manual at [fda-manual](fda-manual/fda-manual.pdf).

Citations
---------

If you use the FDA code, please cite: http://dx.doi.org/10.1186/2046-1682-6-5
