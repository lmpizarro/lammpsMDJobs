from ase.calculators.calculator import PropertyNotImplementedError
import scriptsLammps as sLmps

import sys

sys.path.append('../zhouEam')
sys.path.append('../../pizza/src')

import system as atSys
import wDataLmp as dLmps

from log import log
import numpy as np
from ase.units import GPa

class Lammps_():
    def __init__(self,settings,  label='lammps'):
        self.settings = settings
        self.label = label
        self.calls = 0

        if 'in_lmp' in settings:
            self.in_lmp = settings['in_lmp']
        else:
            self.in_lmp = 'in.min'

        if 'log' in settings:
            self.log = settings['log']
        else:
            self.log = 'log.lammps'

        self.lammps_exe = settings['lammps_exe']

        self.scriptLammps = sLmps.SLammps(settings)

        self.system = settings['sys']


    def get_forces(self, atoms):
        self.update(atoms)

        # TODO
        raise PropertyNotImplementedError

    def get_stress(self, atoms):
        self.update(atoms)
        # TODO
        #raise PropertyNotImplementedError
        pss = ['Pxx', 'Pyy', 'Pzz','Pyz',  'Pxz', 'Pxy' ]
        # TODO
        return np.array([self.run_results[i] for i in pss])*(-1e-4*GPa)



    def get_potential_energy(self, atoms):
        self.update(atoms)
        # TODO
        #raise PropertyNotImplementedError
        self.calculate(atoms)
        return self.run_results['PotEng']

    def update(self, atoms):
        if not hasattr(self,'atoms') or self.atoms != atoms:
            #or 'minimize' in self.parameters:
            self.calculate(atoms)

    def calculate(self, atoms):
       self.atoms = atoms.copy()
       self.run()

    def run(self):
        self.calls += 1

        dL = dLmps.DataLammps(self.settings)
        dL.genFile() # create data.lamp
        self.scriptLammps.create_in() # create in.min


        import atomman.lammps as lmp
        output = lmp.run(self.lammps_exe, self.in_lmp, return_style='object')

        self.parse_log()


    def parse_log(self):
        self.keys_thermo = ['Step', 'Press','CPU', 'Pxx', 'Pyy', 'Pzz', 'Pxy',
        'Pxz', 'Pyz', 'KinEng', 'PotEng','TotEng', 'Lx', 'Ly','Lz','Atoms'] 


        lg = log(self.log)
        status = {}

        for k in self.keys_thermo:
            a =  lg.get(k)
            i = len(a) - 1 
            status[k] = a[i]

        Step = status['Step']
        Press = status['Press']
        CPU = status['CPU']

        Pxx = status['Pxx']
        Pyy = status['Pyy']
        Pzz = status['Pzz']
        Pxy = status['Pxy']
        Pyz = status['Pyz']
        Pxz = status['Pxz']

        KinEng = status['KinEng']
        PotEng = status['PotEng']
        TotEng = status['TotEng']
        Atoms = status['Atoms']


        Lx = status['Lx']
        Ly = status['Ly']
        Lz = status['Lz']

        self.run_results = {'Step':Step, 'Press':Press,'CPU':CPU, 
                'Pxx':Pxx, 'Pyy':Pyy, 'Pzz':Pzz, 'Pxy':Pxy, 'Pxz':Pxz, 'Pyz':Pyz, 
                'KinEng':KinEng, 'PotEng':PotEng,'TotEng':TotEng, 
                'Lx':Lx, 'Ly':Ly,'Lz':Lz,'Atoms':Atoms}

        #print self.run_results


def fitEos(volumes, energies):
    # fit EOS
    #
    from ase.utils.eos import EquationOfState
    from ase.units import kJ, _e

    eos = EquationOfState(volumes, energies)
    vpaf, epaf, B1 = eos.fit()

    print ("vpaf: %f, epaf: %f, B1: %f") % (vpaf, epaf, 1.0e24 * B1 / kJ )



def test_01():
    sys_setting ={'elements':['Al'], 'pot':'zhou', \
                  'pca':[], 'nAtoms':250,\
              #'structure':'bcc',\
              #'positions':'rnd','a':3.0, 'period':[5,5,5]}

              'structure':'fcc',\
              'positions':'rnd','a':4.5, 'period':[5,5,5]}

    sys1 = atSys.System(sys_setting)

    lammps_setting = {'data_lmp':'data.lmp', 
                      'in_lmp':'in.min',
                      'lammps_exe':\
                              '/opt/lmpizarro/GitHub/lammps/src/lmp_serial',
                              'log': 'log.lammps', 'sys':sys1, 'minimize':True}

    calc = Lammps_(lammps_setting)

    atoms = sys1.bulk 
    atoms.set_calculator(calc)

    print (atoms.get_stress())
    epa0 = atoms.get_potential_energy() / atoms.get_number_of_atoms()

    init_cell = np.eye(3) * calc.run_results['Lx']

    print epa0, init_cell

    volumes = []
    energies = []
    for x in np.linspace(0.99, 1.16, 25):

        atoms = sys1.bulk 
        atoms.set_cell(init_cell * x)
        lammps_setting['minimize'] = False
        calc = Lammps_(lammps_setting)
        atoms.set_calculator(calc)

        energies.append(atoms.get_potential_energy() /
                atoms.get_number_of_atoms())

        volumes.append(atoms.get_volume() /
                atoms.get_number_of_atoms())

    print volumes, energies

    fitEos(volumes, energies)
    import matplotlib.pyplot as plt

    plt.plot(volumes, energies)
    plt.show()

def test_elastic():
    sys_setting ={'elements':['Fe'], 'pot':'zhou', \
                  'pca':[], 'nAtoms':250,\
              #'structure':'bcc',\
              #'positions':'rnd','a':3.0, 'period':[5,5,5]}

              'structure':'fcc',\
              'positions':'rnd','a':4.5, 'period':[2,2,2]}

    sys1 = atSys.System(sys_setting)

    lammps_setting = {'data_lmp':'data.lmp', 
                      'in_lmp':'in.min',
                      'lammps_exe': '/opt/lmpizarro/GitHub/lammps/src/lmp_serial',
                      'log': 'log.lammps', 'sys':sys1, 'minimize':'elastic'}

    calc = Lammps_(lammps_setting)
    atoms = sys1.bulk 
    atoms.set_calculator(calc)

    for i in range(2):
        pe = atoms.get_potential_energy()

        print pe
        import eConstants as eco

        print ("c44: %f c11 %f c12 %f B: %f Cp: %f")%\
                (eco.cxx['C44'], eco.cxx['C11'],\
                 eco.cxx['C12'], eco.cxx['B'],\
                 eco.cxx['Pr'])



if __name__ == '__main__':
    test_elastic()
