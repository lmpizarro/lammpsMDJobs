from ase.eos import EquationOfState
import atomman.lammps as lmp
from ase.units import kJ, _e
from ase import Atoms
import atomman as am            
import numpy as np
import parseLogs

fix_min='''
fix 1 all box/relax iso 0.0 vmax 0.001
min_style cg 
minimize 1e-25 1e-25 5000 10000
'''

fix_nve ='''
fix 1 all nve
'''

pot_eam='''
### interactions 
pair_style eam/alloy 
pair_coeff * * U2.eam U 
'''

pot_meam='''
### interactions
pair_style meam 
pair_coeff * * meamf AlS SiS MgS CuS FeS meam.alsimgcufe AlS MgS FeS
'''

global in_frame
in_frame = '''
clear
units metal 
atom_style atomic
boundary p p p 
atom_modify sort 0 0.0 

### read atom configuration 
read_data %s 

### interactions 
%s
 
# Set Thermo 
thermo 10
thermo_style custom step temp press cpu pxx pyy pzz pxy pxz pyz ke pe etotal vol lx ly lz atoms 
thermo_modify flush yes

dump dump_all all custom 1 u.dump id type x y z vx vy vz fx fy fz
timestep 0.001

# run fix min o nve
%s

run 0
unfix 1
#fix 2 all nve
#run 1
#unfix 2
'''

def create_in(atoms_lmps, interaction, fix):
    in_lammps = in_frame % (atoms_lmps, interaction, fix)

    with open('in.min', 'w') as inscript:
         inscript.write( in_lammps )

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

    potential = pot_meam

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
    bulk[0].set('symbol', 'Mg')
    bulk[1].set('symbol', 'Fe')



    in_data = 'in.data'
    in_name = 'in.min'
    system, elements = am.convert.ase_Atoms.load(bulk)
    system_info = lmp.atom_data.dump(system, in_data)

    lammps_exe ='/opt/lmpizarro/GitHub/lammps/src/lmp_serial'
    #lammps_exe ='/usr/bin/lammps'

    create_in(in_data, potential, fix_min)
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
    for x in np.linspace(0.85, 1.15, 40):
        cell = cwll * x
        bulk.set_cell(cell, scale_atoms=True)
        system, elements = am.convert.ase_Atoms.load(bulk)
        system_info = lmp.atom_data.dump(system, in_data)
        create_in(in_data, potential, fix_nve)
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

    eos = EquationOfState(volumes, energies)
    vpaf, epaf, B1 = eos.fit()
    print "EOS vpaf:", vpaf, "A^3"
    print "EOS epaf:", epaf, "eV"
    print "EOS B1:", B1 / kJ * 1.0e24, "GPa"

    import matplotlib.pyplot as plt

    plt.plot(volumes, energies)
    plt.show()
