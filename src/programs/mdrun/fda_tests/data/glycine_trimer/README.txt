gmx grompp -f md.mdp -c glycine_trimer.gro -n index.ndx -p topol.top -o topol.tpr
gmx mdrun -s topol.tpr
