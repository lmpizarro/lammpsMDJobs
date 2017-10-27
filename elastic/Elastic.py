import eConstants
import sys
import atomman.lammps as lmp
import atomman as am

sys.path.append('../potentialTest')
import lmpsScripts as lmpS
import parseLogs

keys_pot = ['epsilon', 'sigma', 'a', 'lambda_', 'gamma',\
        'costheta0', 'A', 'B', 'p', 'q', 'tol']

keys_eC = ['C11', 'C12', 'C44', 'B', 'Sh1', 'Sh2', 'Pr']

sw_pot = {'name':'Si', 'parameters':{'epsilon':2.1683, 'sigma':2.0951,\
        'a':1.8, 'lambda_':10.0, 'gamma':1.2, 'costheta0':-0.333333,\
        'A':7.049556, 'B':0.602224, 'p':4.0, 'q':0.0, 'tol':0.0 }}

def set_epsilon(e):
    sw_pot['parameters']['epsilon'] = e

def eCToStr():
    space = ' ' 
    str_= '' 
    header = ''
    for k in keys_eC:
        str_ += str(eConstants.cxx[k]) + space
        header += k + space
    return str_, header


def getSwPot(sw_pot, filename):
    name = sw_pot['name']
    space = ' ' 
    str_= name + space + name + space + name + space
    for k in keys_pot:
        str_ += str(sw_pot['parameters'][k]) + space 

    with open(filename, 'w') as inscript:
         inscript.write( str_ )




def test01():
    getSwPot(sw_pot, "si.sw")

    values, header =  eCToStr()
    print header
    print values

    set_epsilon(2.0)
    getSwPot(sw_pot, "si.sw")



lammps_exe ='/opt/lmpizarro/GitHub/lammps/src/lmp_serial'

if __name__ == '__main__':
    potential = lmpS.pot_sw
    in_data = 'Al4.lmp'
    in_name = 'in.min'
    fix =  lmpS.fix_min

    lmpS.create_in(in_data, potential, fix)
    output = lmp.run(lammps_exe, in_name, return_style='object')
    status = parseLogs.get_vars2('log.lammps')
    Lx = status['Lx']
    Ly = status['Ly']
    Lz = status['Lz']
    Vol = Lx * Ly * Lz
    Atoms = status['Atoms']
    PotEng = status['PotEng'] / Atoms 

    print Lx, Ly, Lz, Vol, PotEng


    avect =  [ 4.050,  0.000,  0.000]
    bvect =  [ 0.000,  4.050,  0.000]
    cvect =  [ 0.0000001,  0.000,  4.050]
    origin = [ 0.000,  0.000,  0.000]

    prop = {'atype':1, 'pos':[[0.0,0.0,0.0], [0.5,0.5,0.0], [0.5,0.0,0.5],
        [0.0,0.5,0.5]]}
    atoms = am.Atoms(natoms=4, prop=prop)
    box = am.Box(avect=avect, bvect=bvect, cvect=cvect, origin=origin)

    fcc_cell = am.System(box=box, atoms=atoms, scale=True)
    print(fcc_cell)


    sys_info = lmp.atom_data.dump(fcc_cell, 'atom.dat')

