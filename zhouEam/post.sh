tail -n 1031 log-ZrNbAlTi-npt.lammps |head -n 1000 >dataAl03.txt
tail -n 1031 log-ZrNbCrTi-npt.lammps |head -n 1000 >dataCr03.txt
mv log-ZrNbAlTi-npt.lammps log-ZrNbAlTi-npt-03.lammps
mv log-ZrNbCrTi-npt.lammps log-ZrNbCrTi-npt-03.lammps 
gnumeric dataCr03.txt dataAl03.txt
