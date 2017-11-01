import ase
from ase.units import _Nav#, kB, kJ
import periodictable as pt
import ljParameters

class System():
    def __init__(self, setting):
        self.setting = setting
        self.atoms = []

        elements = self.setting['elements']

        self.calcAtoms()

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

    def basicInteraction (self):
        ljp = ljParameters.LjParameters()

        str_ = ljp.lammpsInteraction(self.atoms)

        return str_

    def calcAtoms(self):
        nAt = []

        sum_ = 0
        sumpc = 0
        for p in self.setting['pca']:
            sumpc +=p
            t = int(p * self.setting['nAtoms'] / 100.0)
            sum_ +=t
            nAt.append(t)
        nAt.append(self.setting['nAtoms'] - sum_)
        self.setting['pca'].append(100-sumpc)
        self.setting['nAt'] = nAt

def test_01():
    setting ={'elements':['Zr', 'Fe', 'Al', 'Mo'],\
              'pca':[10, 10, 10], 'nAtoms':250,\
              'structure':'rnd',\
              'positions':'rnd','a':3.0, 'period':[4,4,4]}

    sys = System(setting)
    print sys.setting['elements']
    print sys.setting['nAt']
    print sys.setting['pca']
    print sys.atoms


    pass

if __name__ == '__main__':
    test_01()
