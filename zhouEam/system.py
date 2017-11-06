import ase
from ase.units import _Nav#, kB, kJ
import periodictable as pt
import ljParameters
import EamZhouPotential as zhou
import sys
import random

class System():
    def __init__(self, setting):
        self.setting = setting
        self.atoms = []

        elements = self.setting['elements']

        self.calcAtoms2()

        per = self.setting['period']
        a = self.setting['a']
        self.box =[[0, per[0]*a],[0, per[1]*a],[0, per[2]*a]]


        for e in elements:
            a = ase.Atom(e)
            mass = a.mass / _Nav

            form = pt.formula(a.symbol)
            e_ = form.structure[0][1]
            crys = e_.crystal_structure['symmetry'] 
            a_ = e_.crystal_structure['a'] 
            self.atoms.append({'ase':a, 'mass': mass, 'structure': crys, 'a':a_ })

        if self.setting['structure'] == 'rnd':
            print 'rnd implemented'
            self.genRandomPositions()

            self.genStructure()
 
            self.setRandomStructure()
        if self.setting['structure'] == 'fcc':
            self.setCrystal('fcc')
            print 'fcc implemented'

            if self.setting['positions'] == 'rnd':
                self.setRandomStructure()
        if self.setting['structure'] == 'bcc':
            self.setCrystal('bcc')

            if self.setting['positions'] == 'rnd':
                self.setRandomStructure()
            print 'bcc no implemented'
            #sys.exit(0)
        if self.setting['structure'] == 'hcp':
            print 'hcp no implemented'
            sys.exit(0)

    def setCrystal(self, crys):
        px = self.setting['period'][0]
        py = self.setting['period'][1]
        pz = self.setting['period'][2]
        a =  self.setting['a']

        if crys == 'fcc':
            from ase.lattice.cubic import FaceCenteredCubic
            alloy = FaceCenteredCubic(directions=[[1,0,0], [0,1,0], [0,0,1]],
                                        size=(px,py,pz), symbol='Cu',
                    pbc=(1,1,1), latticeconstant=a)

        if crys == 'bcc':
            from ase.lattice.cubic import BodyCenteredCubic
            alloy = BodyCenteredCubic(directions=[[1,0,0], [0,1,0], [0,0,1]],
                                        size=(px,py,pz), symbol='Cu',
                    pbc=(1,1,1), latticeconstant=a)

        self.setting['nAtoms'] =  alloy.get_number_of_atoms()

        self.calcAtoms2()

        self.genStructure()

        self.pos = alloy.get_positions()
        self.nAt = self.setting['nAtoms']

    def genStructure(self):
        self.t1_ = []
        for i,e in enumerate(self.setting['nAt']):
            [self.t1_.append(i+1) for j in range(e)]
 
    def setRandomStructure(self):
        x = [int(random.random()*len(self.t1_)) for i in range(len(self.t1_))]

        for i in range(len(x) - 1):
            t1 = self.t1_[x[i]]
            t2 = self.t1_[x[i+1]]

            self.t1_[x[i]] = t2
            self.t1_[x[i+1]] = t1
            
        return self.pos, self.t1_, self.box

    def genRandomPositions(self):
        Lx = self.sys.box[0][1]
        Ly = self.sys.box[1][1]
        Lz = self.sys.box[2][1]

        self.pos =[]
        for i in range(self.nAt):
            x = random.random() * Lx
            y = random.random() * Ly
            z = random.random() * Lz
            self.pos.append([x,y,z])


    def getAtoms(self):
        return self.atoms

    def Interaction (self):
        if self.setting['pot'] == 'lj':
            ljp = ljParameters.LjParameters()
            str_ = ljp.lammpsInteraction(self.atoms)
        if self.setting['pot'] == 'zhou':
            gz = zhou.calcPotentials(self.setting['elements'])
            print "TODO: not yet implemented"
            gz.createPot()
            str_ = gz.getEam()

        return str_

    def calcAtoms2(self):
        nAt = []

        sum_ = 0
        sumpc = 0.0000001
        len_elements = len(self.setting['elements'])
        len_pca = len(self.setting['pca'])

        print len_pca, len_elements

        if len_elements != len_pca + 1:
            print 'error len'
            sys.exit(0)

        for p in self.setting['pca']:
            sumpc +=p

        if sumpc >= 100 or sumpc <= 0:
            print 'error pc'
            sys.exit(0)

        sumpc = 0

        for i in range(len(self.setting['elements']) - 1):
            p = self.setting['pca'][i]
            sumpc +=p
            t = int(p * self.setting['nAtoms'] / 100.0)
            sum_ +=t
            nAt.append(t)

        nAt.append(self.setting['nAtoms'] - sum_)
        self.setting['nAt'] = nAt

        #print 'pca', self.setting['pca']
        #print 'nAt', self.setting['nAt']


    def getMasess(self):
        str_ =''
        i = 1
        for e in self.atoms:
            str_ += 'mass ' + str(i) + ' ' +str(e['mass']) + '\n'
            i+=1
        return  str_

def test_01():

    lammps_setting = {'data_lmp':'data.lmp', 
                      'in_lmp':'in.min',
                       'lammps_exe' :'/opt/lmpizarro/GitHub/lammps/src/lmp_serial'}

    setting ={'elements':['Al'], 'pot':'zhou', \
              'pca':[], 'nAtoms':250,\
              #'structure':'bcc',\
              #'positions':'rnd','a':3.0, 'period':[5,5,5]}

              'structure':'fcc',\
              'positions':'rnd','a':4.2, 'period':[5,5,5],\
              'lammps_setting':lammps_setting }


    
    sys = System(setting)
    print sys.setting['elements']
    print sys.setting['nAt']
    print sys.setting['pca']
    print sys.atoms

    print sys.Interaction()

def test_02():

    lammps_setting = {'data_lmp':'data.lmp', 
                      'in_lmp':'in.min',
                       'lammps_exe' :'/opt/lmpizarro/GitHub/lammps/src/lmp_serial'}

    setting ={'elements':['Al'], 'pot':'zhou', \
              'pca':[], 'nAtoms':250,\
              #'structure':'bcc',\
              #'positions':'rnd','a':3.0, 'period':[5,5,5]}

              'structure':'fcc',\
              'positions':'rnd','a':4.2, 'period':[5,5,5],\
              'lammps_setting':lammps_setting }


    sys = System(setting)
    print sys.setting['elements']
    print sys.setting['nAt']
    print sys.setting['pca']
    print sys.atoms

    print sys.Interaction()

if __name__ == '__main__':
    test_01()
