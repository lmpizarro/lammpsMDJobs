import sys
import ase
from ase.units import _Nav, kB, kJ
import periodictable as pt
import mendeleev as mdlv

sys.path.append('../../pizza/src')

from log import log
import math

import random

lammps_exe ='/opt/lmpizarro/GitHub/lammps/src/lmp_serial'

JOULE_TO_EV = kJ / 1000
AVOGADRO = _Nav
EV_AT = JOULE_TO_EV * 1000 / AVOGADRO

def thermal_props(element, EV=True):

    element = mdlv.element(element)
    mp = element.melting_point
    bp = element.boiling_point

    if EV == True:
        ev_h = element.evaporation_heat * EV_AT
        h_o_f = element.heat_of_formation * EV_AT
        f_h = element.fusion_heat * EV_AT
    else:
        ev_h = element.evaporation_heat
        h_o_f = element.heat_of_formation
        f_h = element.fusion_heat

    data = {'evaporation_heat': ev_h,
            'heat_of_formation': h_o_f,
            'fusion_heat': f_h,
            'melting_point': mp,
            'boiling_point':bp
           }
    return data

'''
   REFS:
   Calculation of the parameters of the Lennard-Jones potential for pairs of
   identical atoms based on the properties of solid substances

   V. P. Filippova, S. A. Kunavin, M. S. Pugachev


'''
def calc_lj_01 (element):

    form = pt.formula(element)
    e_ = form.structure[0][1]

    data = thermal_props(element)

    if  e_.crystal_structure['symmetry'] == 'fcc' or \
        e_.crystal_structure['symmetry'] == 'hcp':

        sigma = e_.crystal_structure['a'] * 0.635
        epsilon = data['evaporation_heat'] * 0.172

    elif e_.crystal_structure['symmetry'] == 'BCC':
        sigma = e_.crystal_structure['a'] *0.788
        epsilon = data['evaporation_heat'] * 0.1819

    else:
        sigma = e_.crystal_structure['a'] * 0.7  # 0.648
        epsilon = data['evaporation_heat'] * 0.18
        #epsilon = data['fusion_heat']  * 0.843

    return ({'epsilon': epsilon, 'sigma':sigma})


def test_lj():
    elements = ['U', 'Fe', 'Mo','V', 'Cr','Cu', 'Ni', 'Al',\
                'Zr', 'Si', 'Mg', 'Co', 'Ag', 'Pt', 'W', \
                'Pb','Au', 'Pd', 'Ti', 'Nb', 'Ta' ]

    for e in elements:
        print (e, calc_lj_01(e))



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
            minimize 1e-15 1e-15 50000 50000
        '''
        self.defInteraction = self.defInteraction()
        self.in_frame = self.formatMultiline(self.in_frame)
        self.fix = self.formatMultiline(self.fix_min)
        self.interaction = self.formatMultiline(self.defInteraction)


    def defInteraction (self):
        str_ = 'pair_style lj/cut 7\n'
     
        for i,e in enumerate(self.atoms):
            el = e['ase'].symbol
            paramsLJ = calc_lj_01(el)
            # eps sigma cut
            coefPair = '  ' + str(paramsLJ['epsilon']) + ' ' +\
                    str(paramsLJ['sigma']) + '  ' + \
                    str(paramsLJ['sigma'] * 1.5) +  '\n'

            j = i + 1
            str_ += 'pair_coeff ' + str(j) + ' ' + str(j) + coefPair 

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

class DataLammps():
    def __init__(self, elements, mult, size=[10.0,10.0,10.0]):
        self.pos = None
        self.t1_ = None
        self.mult = mult
        self.nTypes = len(elements) 
        self.nAt = len(elements) * mult 
        self.box =[[0, size[0]],[0, size[1]],[0, size[2]]]

    def genFile(self, fileName):
        self.str_ = '\n'
        self.str_ += str(self.nAt) + ' atoms\n'
        self.str_ += str(self.nTypes) + ' atom types\n'
        xlo = self.box[0][0]
        xhi = self.box[0][1]
        ylo = self.box[1][0]
        yhi = self.box[1][1]
        zlo = self.box[2][0]
        zhi = self.box[2][1]
        self.str_ += str(float(xlo)) + ' ' +str(xhi)+ '  ' + ' xlo xhi\n'
        self.str_ += str(float(ylo)) + ' ' +str(yhi)+ '  ' + ' ylo yhi\n'
        self.str_ += str(float(zlo)) + ' ' +str(zhi)+ '  ' + ' zlo zhi\n'
        self.str_ +='\n'
        self.str_ +='Atoms\n'
        self.str_ +='\n'


        for i,e in enumerate(self.pos):
            self.str_ += str(i +1) + '  ' +   str(self.t1_[i]) + ' ' +\
                    str(e[0]) + ' ' + str(e[1])+ ' '  + str(e[2]) + '\n'
        with open(fileName, 'w') as inscript:
            inscript.write( self.str_)


    def setRandomPositions(self):
        Lx = self.box[0][1]
        Ly = self.box[1][1]
        Lz = self.box[2][1]

        self.t1_ = []
        [[self.t1_.append(e) for e in [i+1]*self.mult] \
                for i in range(self.nTypes)]

        self.pos =[]
        for i in range(self.nAt):
            x = random.random() * Lx
            y = random.random() * Ly
            z = random.random() * Lz
            self.pos.append([x,y,z])

        return self.pos, self.t1_, self.box

def test_01():

    data_lmp = 'data.lmp'
    in_lmp = 'in.min'
    elements = ['Zr', 'Fe', 'Al', 'Mo']
    rL = RunLammps(elements)
    rL.setDataLmp(data_lmp)
    rL.create_in(in_lmp)

    dL = DataLammps(elements, 50)
    dL.setRandomPositions()
    dL.genFile(data_lmp)

if __name__ == '__main__':
    test_01()
