import sys
sys.path.append('../../pizza/src')

import numpy as np
from log import log
from dump import dump

LOG='log.lammps'
def get_vars(LOG):
    lg = log(LOG)
    PotEng = lg.get('PotEng')
    Atoms = lg.get('Atoms')
    Lx = lg.get('Lx')
    Ly = lg.get('Ly')
    Lz = lg.get('Lz')
    ix = len(Lx) -1
    iy = len(Ly) -1
    iz = len(Lz) -1
    iPE = len(PotEng) -1
    iAt = len(Atoms) -1
    return Lx[ix], Ly[iy], Lz[iz], PotEng[iPE], Atoms[iAt]


keys_thermo = ['Step', 'Press','PotEng','TotEng', 'Lx', 'Ly','Lz','Atoms'] 

def get_vars2(LOG):
    lg = log(LOG)
    status = {}

    for k in keys_thermo:
        a =  lg.get(k)
        i = len(a) - 1 
        status[k] = a[i]

    return status

def get_forces(dump_):
    d = dump(dump_)
    t = d.time()
    index = t[len(t) -1]
    d.aselect.all(index)
    return d.vecs(index, 'fx', 'fy', 'fz')

if __name__ == '__main__':
    print get_vars(LOG)
    fx, fy, fz = get_forces('u.dump')
    print fx
    print fy
    print fz

    d= dump('u.dump')

    fx,fy,fz = d.atom(1,"fx","fy",'fz')  

    print get_vars2(LOG) 
