import sys
import ase
from ase.units import _Nav

sys.path.append('../../pizza/src')

from log import log
import math

lammps_exe ='/opt/lmpizarro/GitHub/lammps/src/lmp_serial'


class RunLammps():
    keys_thermo = ['Step', 'Press','PotEng','TotEng', 'Lx', 'Ly','Lz','Atoms'] 
    def __init__(self, elements):

        self.atoms = []

        for e in elements:
            a = ase.Atom(e)
            mass = a.mass / _Nav
            self.atoms.append({'ase':a, 'mass': mass })

        self.in_frame =''' 
            clear
            boundary    p p p
            units       metal
            atom_style  atomic

            thermo 50

            # 1 atoms positions
            read_data %s	

            # 2 interactions
            %s

            # 3 mass 1 1.0e-20
            %s
            
            thermo_style custom step press pe etotal  lx ly lz atoms 
            # 4 fix
            %s

            run 0
        '''
        self.fix_min='''
            fix 1 all box/relax iso 0.0 vmax 0.001
            min_style cg
            minimize 1e-15 1e-15 5000 5000
        '''
        self.defInteraction = self.defInteraction()
        self.in_frame = self.formatMultiline(self.in_frame)
        self.fix = self.formatMultiline(self.fix_min)
        self.interaction = self.formatMultiline(self.defInteraction)


    def defInteraction (self):
        str_ = 'pair_style lj/cut 7\n'
     
        for i,e in enumerate(self.atoms):
            j = i + 1
            str_ += 'pair_coeff ' + str(j) + ' ' + str(j) + ' 1 2 3\n'

        return str_

    def getMasess(self):
        str_ =''
        i = 1
        for e in self.atoms:
            str_ += 'mass ' + str(i) + ' ' +str(e['mass']) + '\n'
            i+=1
        return  str_

    def get_vals(self, LOG):
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



    def formatMultiline(self, multiline):
        l = multiline.split('\n')
        str_=''
        for e in l:
           str_ +=e.lstrip() + '\n'
        return str_

    def create_in(self, fileName):
        mass = self.getMasess()

        in_lammps = self.in_frame % (self.data_lmps, self.interaction, mass, self.fix)

        with open(fileName, 'w') as inscript:
             inscript.write( in_lammps )

    def setDataLmp(self, name):
        self.data_lmps = name

    def setInteraction (self, interaction):
        self.interaction = interaction

    def setFix(self, fix):
        self.fix = fix

def test_01():
    elements = ['Zr', 'Fe', 'Al']
    rL = RunLammps(elements)
    print rL.getMasess()
    rL.setDataLmp('data.lmp')
    rL.create_in('in.min')
    pass


if __name__ == '__main__':
    test_01()
