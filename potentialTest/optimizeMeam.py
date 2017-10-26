import os
os.environ["LAMMPS_COMMAND"] = "/opt/lmpizarro/GitHub/lammps/src/lmp_serial"

from ase.constraints import UnitCellFilter
from ase.eos import EquationOfState
from ase.optimize import FIRE
from ase.units import kJ, _e
from lammps import LAMMPS
from ase import Atoms
import numpy as np

files = [ "meamf", "meam.alsimgcufe" ]
mypc = [ "* * " + files[0] + " AlS SiS MgS CuS FeS " + files[1] ]
parameters = { "pair_style" : "meam", "pair_coeff" : mypc }
myext = "S"

def pick_elements(parameters, elems):
    for myel in elems:
        parameters["pair_coeff"][0]  += " " + myel + myext

species = ['Al', 'Mg', 'Fe']
pick_elements(parameters, species)

if __name__ == '__main__':

    calc = LAMMPS(parameters=parameters, files=files)
    a0 = 4.5
    b0 = a0 / 2.0
    bulk = Atoms(['Al']*4,
                 positions=[(0,0,0),(b0,b0,0),(b0,0,b0),(0,b0,b0)],
                 cell=[a0, a0, a0],
                 pbc=True)
    bulk.set_calculator(calc)

    bulk[0].set('symbol', 'Mg')
    bulk[1].set('symbol', 'Fe')

    # Ininit
    epa0 = bulk.get_potential_energy() / bulk.get_number_of_atoms()
    vpa0 = bulk.get_volume() / bulk.get_number_of_atoms()
    print "epa0:", epa0
    print "vpa0:", vpa0, "\n"
    print 'a0 %f\n'% (bulk.get_cell()[0][0])

   
    # Optimization
    FIRE(UnitCellFilter(bulk, mask=[1,1,1,0,0,0]),\
            logfile='FIRE.log').run(fmax=0.001)

    epaf = bulk.get_potential_energy() / bulk.get_number_of_atoms()
    vpaf = bulk.get_volume() / bulk.get_number_of_atoms()
    print "FIRE epaf:", epaf
    print "FIRE vpaf:", vpaf, "\n"

    opt_cell = bulk.get_cell()

    print 'a0 %f\n'% (opt_cell[0][0])

    # reoptimize/check volume
    #
    volumes = []
    energies = []
    for x in np.linspace(0.98, 1.02, 5):
        bulk.set_cell(opt_cell * x, scale_atoms=True)
        volumes.append(bulk.get_volume() / bulk.get_number_of_atoms())
        energies.append(bulk.get_potential_energy() / bulk.get_number_of_atoms())

    # fit EOS
    #
    eos = EquationOfState(volumes, energies)
    vpaf, epaf, B1 = eos.fit()
    print "EOS vpaf:", vpaf, "A^3"
    print "EOS epaf:", epaf, "eV"
    print "EOS B1:", B1 / kJ * 1.0e24, "GPa"

    strain = 0.001
    diag = ([1, 0, 0],
            [0, 1, 0],
            [0, 0, 1])

    # c44
    #
    defm1 = np.array([[0, strain, strain],
                  [strain, 0, strain],
                  [strain, strain, 0]])
    cell1 = np.dot(opt_cell, defm1 + diag)
    bulk.set_cell(cell1, scale_atoms=True)

    ene1 = bulk.get_potential_energy() / bulk.get_number_of_atoms()
    dele = (ene1 - epaf) / vpaf * _e / 1.0e-30  # eV/angstrom^3
    c44 = dele / 6 / strain / strain / 1e9  # GPa
    print "epaf:", epaf
    print "ene1:", ene1
    print "del:", ene1 - epaf
    print "c44:", c44, "GPa\n"
