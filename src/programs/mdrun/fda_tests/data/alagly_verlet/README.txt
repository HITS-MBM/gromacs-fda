gmx grompp -f md.mdp -c conf.gro -p topol.top -o alagly.tpr
gmx mdrun -s alagly.tpr
