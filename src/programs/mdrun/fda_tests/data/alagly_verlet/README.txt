gmx grompp -f md.mdp -c conf.gro -p topol.top -o topol.tpr
gmx mdrun -s topol.tpr
