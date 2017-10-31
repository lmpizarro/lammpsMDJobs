import sys
sys.path.append('../../pizza/src')

import numpy as np
from log import log
from dump import dump

keys_thermo = ['Step', 'Press','PotEng','TotEng', 'Lx', 'Ly','Lz','Atoms'] 

def get_vars2(LOG):
    lg = log(LOG)
    status = {}

    for k in keys_thermo:
        a =  lg.get(k)
        i = len(a) - 1 
        status[k] = a[i]

    return status

if __name__ == '__main__':
    LOG='log.lammps'
    print get_vars2(LOG) 
