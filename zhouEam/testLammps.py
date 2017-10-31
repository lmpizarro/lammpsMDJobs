from ase.lattice.cubic import FaceCenteredCubic
import atomman.lammps as lmp
import atomman as am
import random

import genSetFlEamZhou as zhou
import sys
sys.path.append('../../pizza/src')

from log import log

keys_thermo = ['Step', 'Press','PotEng','TotEng', 'Lx', 'Ly','Lz','Atoms'] 

def get_vars2(LOG):
    lg = log(LOG)
    status = {}

    for k in keys_thermo:
        a =  lg.get(k)
        i = len(a) - 1 
        status[k] = a[i]

    Lx = status['Lx']
    PotEng = status['PotEng']
    Atoms = status['Atoms']

    return (Lx, PotEng / Atoms, Atoms)

fix_min='''
fix 1 all box/relax iso 0.0 vmax 0.001
min_style cg
minimize 1e-15 1e-15 5000 5000
'''

in_frame =''' 
clear
boundary    p p p
units       metal
atom_style  atomic

thermo 50

# 2 atoms positions
read_data %s	

# 1 interactions
%s

mass 1 1.0e-20

thermo_style custom step press pe etotal  lx ly lz atoms 
# 3 fix
%s

run 0
'''

def create_in(data_lmps, interaction, fix):
    in_lammps = in_frame % (data_lmps, interaction, fix)

    with open('in.min', 'w') as inscript:
         inscript.write( in_lammps )

import math

lammps_exe ='/opt/lmpizarro/GitHub/lammps/src/lmp_serial'


def test_01():
    ES = ['Al', 'Nb','Cr']

    potZhou = zhou.calcPotentials(ES)
    potZhou.createPot()

    atomsPerCell = 4
    a0 = 4.1

    for nL in range(len(ES)):
        nCells = math.pow(nL, 3)
        nAtoms = int(atomsPerCell * nCells)
        if nAtoms > len(ES):
            break

    n = int(nAtoms / len(ES) )
    npc = 100.0 * n / nAtoms

    mix = {}
    for k in ES:
        mix[k]={'pcw':0, 'pca':npc, 'nAtoms':n}

    print mix

    lNums = random.sample(range(nAtoms), nAtoms)

    for i,e in enumerate(ES):
        mix[e]['pos'] = lNums[i*n: (i + 1 )*n]

    formula = ''
    for e in mix:
        formula +=e+str(mix[e]['nAtoms'])
    lamp_data = formula + '.lmp'
    in_name = 'in.min'

    print formula, nL

    alloy = FaceCenteredCubic(directions=[[1,0,0], [0,1,0], [0,0,1]],\
            size=(nL,nL,nL), symbol=ES[0],\
            pbc=(1,1,1), latticeconstant=a0 )

    for e in mix:
        for p in mix[e]['pos']:
             alloy[p].symbol = e

    print alloy.get_chemical_formula()

    system, elements = am.convert.ase_Atoms.load(alloy)
    system_info = lmp.atom_data.dump(system, lamp_data , units='metal')
    
    create_in(lamp_data, potZhou.getEam(), fix_min)
    output = lmp.run(lammps_exe, in_name, return_style='object')

    (Lx, PotEng, Atoms) = get_vars2('log.lammps')
    

    print Lx / nCells, PotEng, Atoms, a0


if __name__ == '__main__':
    test_01()
