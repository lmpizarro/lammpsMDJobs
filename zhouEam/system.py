import ase
from ase.units import _Nav#, kB, kJ
import periodictable as pt
import ljParameters
import EamZhouPotential as zhou
import sys

class System():
    def __init__(self, setting):
        self.setting = setting
        self.atoms = []

        elements = self.setting['elements']

        self.calcAtoms2()

        for e in elements:
            a = ase.Atom(e)
            mass = a.mass / _Nav

            form = pt.formula(a.symbol)
            e_ = form.structure[0][1]
            crys = e_.crystal_structure['symmetry'] 
            a_ = e_.crystal_structure['a'] 
            self.atoms.append({'ase':a, 'mass': mass, 'structure': crys, 'a':a_ })

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
    setting ={'elements':['Al', 'Fe', 'Cr'], 'pot':'zhou',\
              'pca':[20, 70 ], 'nAtoms':250,\
              #'structure':'bcc',\
              #'positions':'rnd','a':3.0, 'period':[5,5,5]}

              'structure':'fcc',\
              'positions':'rnd','a':4.0, 'period':[4,4,4]}


    sys = System(setting)
    print sys.setting['elements']
    print sys.setting['nAt']
    print sys.setting['pca']
    print sys.atoms

    print sys.Interaction()

def test_02():
    setting ={'elements':['Zr', 'Fe', 'Al', 'Mo'], 'pot':'zhou',\
              'pca':[10, 10, 10], 'nAtoms':250,\
              #'structure':'bcc',\
              #'positions':'rnd','a':3.0, 'period':[5,5,5]}

              'structure':'fcc',\
              'positions':'rnd','a':4.0, 'period':[4,4,4]}


    sys = System(setting)
    print sys.setting['elements']
    print sys.setting['nAt']
    print sys.setting['pca']
    print sys.atoms

    print sys.Interaction()

if __name__ == '__main__':
    test_01()
