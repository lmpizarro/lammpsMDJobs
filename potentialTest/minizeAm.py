from ase.eos import EquationOfState
import atomman.lammps as lmp
from ase.units import kJ, _e
from ase import Atoms
import atomman as am            
import numpy as np
import parseLogs
import lmpsScripts as lmpS

def get_vars(sims, index):
    Lx = sims[0]['thermo']['Lx']
    Ly = sims[0]['thermo']['Ly']
    Lz = sims[0]['thermo']['Lz']
    PotEng = sims[0]['thermo']['PotEng']
    Atoms = sims[0]['thermo']['Atoms']
    index = len(Lx) -1
    Lx = Lx[index]
    Ly = Ly[index]
    Lz = Lz[index]
    PotEng = PotEng[index]
    Atoms = Atoms[index]
    return (Lx, Ly, Lz, PotEng, Atoms)

if __name__ == '__main__':
    a0 = 4.24
    b0 = a0 / 2.0
    size = 3

    potential = lmpS.pot_eam

    bulk = Atoms(['Al']*4,
                     positions=[(0,0,0),(b0,b0,0),(b0,0,b0),(0,b0,b0)],
                     cell=[a0]*3,
                     pbc=True)
    '''
    from ase.lattice.cubic import FaceCenteredCubic
    bulk = FaceCenteredCubic(directions=[[1,0,0], [0,1,0], [0,0,1]],\
            size=(size,size,size), symbol='U',\
            pbc=(1,1,1), latticeconstant=a0 )

    '''
    #bulk[0].set('symbol', 'Zr')
    #bulk[1].set('symbol', 'Fe')



    in_data = 'in.data'
    in_name = 'in.min'
    system, elements = am.convert.ase_Atoms.load(bulk)
    system_info = lmp.atom_data.dump(system, in_data)

    lammps_exe ='/opt/lmpizarro/GitHub/lammps/src/lmp_serial'
    #lammps_exe ='/usr/bin/lammps'

    lmpS.create_in(in_data, potential, lmpS.fix_min)
    output = lmp.run(lammps_exe, in_name, return_style='object')
    sims =  output.simulations
    print(len(sims), 'simulation run(s) in Log')
    status = parseLogs.get_vars2('log.lammps')
    Lx = status['Lx']
    PotEng = status['PotEng']
    Atoms = status['Atoms']

    # reoptimize/check volume
    #
    volumes = []
    energies = []
    cwll = bulk.get_cell()
    for x in np.linspace(0.75, 1.35, 40):
        cell = cwll * x
        bulk.set_cell(cell, scale_atoms=True)
        system, elements = am.convert.ase_Atoms.load(bulk)
        system_info = lmp.atom_data.dump(system, in_data)
        lmpS.create_in(in_data, potential, lmpS.fix_nve)
        output = lmp.run(lammps_exe, in_name, return_style='object')
        sims =  output.simulations
        status = parseLogs.get_vars2('log.lammps')
        Lx = status['Lx']
        Ly = status['Ly']
        Lz = status['Lz']
        PotEng = status['PotEng']
        Atoms = status['Atoms']

        volumes.append(Lx*Ly*Lz / Atoms)
        energies.append(PotEng / Atoms)

    import matplotlib.pyplot as plt

    plt.plot(volumes, energies)
    plt.show()
    '''
    eos = EquationOfState(volumes, energies)
    vpaf, epaf, B1 = eos.fit()
    print "EOS vpaf:", vpaf, "A^3"
    print "EOS epaf:", epaf, "eV"
    print "EOS B1:", B1 / kJ * 1.0e24, "GPa"
    '''

