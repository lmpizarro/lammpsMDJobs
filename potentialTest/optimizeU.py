from ase import Atoms
from lammps import LAMMPS
from ase.constraints import UnitCellFilter
from ase.optimize import FIRE

import os
os.environ["LAMMPS_COMMAND"] = "/opt/lmpizarro/GitHub/lammps/src/lmp_serial"


if __name__ == '__main__':

    pair_style = 'eam/alloy'
    Pd_eam_file = 'U2.eam'
    pair_coeff = [ '* * ' + Pd_eam_file + ' U' ]
    parameters = { 'pair_style' : pair_style, 'pair_coeff' : pair_coeff }
    files = [ Pd_eam_file ]

    calc = LAMMPS(parameters=parameters, files=files)
    a0 = 4.20
    b0 = a0 / 2.0
    bulk = Atoms(['U']*4,
                     positions=[(0,0,0),(b0,b0,0),(b0,0,b0),(0,b0,b0)],
                     cell=[a0]*3,
                     pbc=True)
    bulk.set_calculator(calc)

    FIRE(UnitCellFilter(bulk, mask=[1,1,1,0,0,0]), logfile='FIRE.log').run(fmax=0.0001)

    print bulk.get_cell()[0][0]
    print bulk.get_number_of_atoms()

    calc.clean()
