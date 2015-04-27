Force Distribution Analysis (FDA)
=================================

Background
----------

Calculations of internal forces and stress distributions have long been used by mechanical engineers in designing strong structures, from tall sky-scrapers and kilometer-long bridges to Formula 1 cars and aircraft wings. The same kind of analysis, but applied to atomic and molecular level, could shed light onto the mechanical stability of (bio)molecules, allosteric mechanisms and signal propagation in cells or mechanically induced processes such as blood coagulation.

Molecular dynamics (MD) simulations apply Newton's equations of motion to atomic systems. During each integration time step, atoms interact with each other through pairwise forces. All pairwise forces acting on an atom are summed up and the resulting force determines the displacement of the atom in the next time step. These displacements are then typically used to analyze MD simulations.

The analysis of individual atomic pairwise forces, which we call Force Distribution Analysis, can offer a much more sensitive view of a (bio)molecular system, be it in a stable state or during a transition between states. Such analysis would be able to explain, for example, the functioning of an atomic level Newton's cradle, where the application of force on an atom does not translate directly into a displacement of all spheres but only of the outer ones. Similarly, atoms can propagate forces in the absence of any significant atomic displacements, for example through a protein's rigid core.

Under the generic name of Force Distribution Analysis, you can find here both the initial implementation, simply called FDA, and the newest development, called Time-resolved FDA or TRFDA. In addition to tracking changes in pairwise forces, the implementation of TRFDA fixes several shortcomings of FDA and extends it to a unified treatment of atomic- and residue-level pairwise forces as well as stress calculations. Both implementations are based on the GROMACS MD package.

Stable states versus dynamic processes
--------------------------------------

FDA focuses on the analysis of averaged dynamical data. This makes sense only for simulations at quasi-equilibrium, i.e. for MD simulations sampling a single energy well. There should be no significant structural changes during the simulation and, therefore, pairwise forces are supposed to be relatively constant apart from thermal fluctuations. Such an MD simulation could be the equilibration of a protein, in the presence or absence of a ligand or a constant force. FDA can also be used to find differences in the distribution of internal forces between two or more stable states, like those between an apo and a ligand-bound state, which could be associated with an allosteric mechanism. Likewise, differences in force distribution between a molecule under high and low external force can give insight into structural motifs involved in bearing the external load.

In contrast, TRFDA focuses on the analysis of changes in pairwise forces over time. Such changes appear e.g. during MD simulations of protein unfolding or other conformational transitions, in atomic level signal propagation or upon molecular association. Starting from pairwise forces, TRFDA can also calculate an atomic level stress, which we call punctual stress; changes in the stress distribution offer a view of the mechanical stability of the simulated molecular system.

Citations
---------

If you use the FDA code, please cite: http://dx.doi.org/10.1186/2046-1682-6-5
