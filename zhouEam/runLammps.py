import sys
sys.path.append('../../pizza/src')
from log import log

from ase.calculators.calculator import PropertyNotImplementedError

import system

import numpy as np
from ase.units import GPa

import dataLammps as dLmps

class RLammps():
    def __init__(self, system, label='lammps'):
        self.label = label

        self.system = system
        self.keys_thermo = ['Step', 'Press','CPU', 'Pxx', 'Pyy', 'Pzz', 'Pxy',
        'Pxz', 'Pyz', 'KinEng', 'PotEng','TotEng', 'Lx', 'Ly','Lz','Atoms'] 

        self.atoms = system.getAtoms()

        self.data_lmps = self.system.setting['lammps_setting']['data_lmp']
        self.in_lmp =  self.system.setting['lammps_setting']['in_lmp']
        self.lammps_exe = self.system.setting['lammps_setting']['lammps_exe']


        self.log = 'log.lammps'

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
            
            thermo_style custom step press cpu pxx pyy pzz pxy pxz pyz ke pe etotal  vol lx ly lz atoms 
            # 4 fix
            %s

            run 0
        '''
        self.fix_min='''
            fix 1 all box/relax iso 0.0 vmax 0.001
            min_style cg
            minimize 1e-15 1e-15 50000 50000
        '''
        self.defInteraction = self.system.Interaction()
        self.in_frame = self.formatMultiline(self.in_frame)
        self.fix = self.formatMultiline(self.fix_min)
        self.interaction = self.formatMultiline(self.defInteraction)

    def update(self):
        self.atoms = self.system.getAtoms()

    def get_vals(self):

        status = self.run_results 

        nx = self.system.setting['period'][0]
        ny = self.system.setting['period'][1]
        nz = self.system.setting['period'][2]

        Lx = status['Lx']/ nx
        Ly = status['Ly']/ ny
        Lz = status['Lz']/ nz
        PotEng = status['PotEng']
        Atoms = status['Atoms']


        Pxx = status['Pxx']
        Pyy = status['Pyy']
        Pzz = status['Pzz']
        Pxy = status['Pxy']
        Pyz = status['Pyz']
        Pxz = status['Pxz']


        return (Lx, Ly, Lz, PotEng / Atoms, Atoms)


    def parse_log(self):
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


    def get_forces(self):
        # TODO
        raise PropertyNotImplementedError

    def get_stress(self):
        pss = ['Pxx', 'Pyy', 'Pzz','Pyz',  'Pxz', 'Pxy' ]
        # TODO
        return np.array([self.run_results[i] for i in pss])*(-1e-4*GPa)

    def get_potential_energy(self):
        return self.run_results['PotEng']

    def formatMultiline(self, multiline):
        l = multiline.split('\n')
        str_=''
        for e in l:
           str_ +=e.lstrip() + '\n'
        return str_

    def create_in(self):
        mass = self.system.getMasess()
        in_lammps = self.in_frame % (self.data_lmps, self.interaction, mass, self.fix)

        with open(self.in_lmp, 'w') as inscript:
             inscript.write( in_lammps )

    def setInteraction (self, interaction):
        self.interaction = interaction

    def setFix(self, fix):
        self.fix = fix

    def run(self):
        self.update()
        self.system.update()
        dL = dLmps.DataLammps(self.system)
        dL.genFile()
        self.create_in()


        import atomman.lammps as lmp
        output = lmp.run(self.lammps_exe, self.in_lmp, return_style='object')

        self.parse_log()


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

    sys = system.System(setting)
    rL = RLammps(sys)
    rL.run()
    (Lx, Ly, Lz, PotEng, atoms) = rL.get_vals()

    stress =  rL.get_stress()
    print rL.get_potential_energy()
    print Lx, Ly, Lz, PotEng, atoms

    ##############################################

    cell  = sys.bulk.get_cell()
    vol0 = sys.bulk.get_volume()

    sxx0, syy0, szz0, syz0, szx0, sxy0  = rL.get_stress()

    ## C11
    eps = 0.00001
    T = np.diag( [ eps, 0.0, 0.0 ] )
    sys.bulk.set_cell( np.dot(np.eye(3)+T, cell.T).T, scale_atoms=True )

    rL.run()
    (Lx, Ly, Lz, PotEng, atoms) = rL.get_vals()
    print Lx, Ly, Lz, PotEng, atoms


    sxx11, syy11, szz11, syz11, szx11, sxy11  = rL.get_stress()
    from ase.units import GPa


    C11  = (sxx11-sxx0)/eps

    print C11 / GPa

    vol2 =  sys.bulk.get_volume()

    print vol2 -vol0

    ##############################################

def test_02():

    lammps_setting = {'data_lmp':'data.lmp', 'in_lmp':'in.min' ,
                       'lammps_exe' :'/opt/lmpizarro/GitHub/lammps/src/lmp_serial'}

    setting ={'elements':['Zr', 'Fe', 'Al'], 'pot':'zhou', \
              'pca':[10, 10], 'nAtoms':250,\
              'structure':'bcc',\
              'positions':'rnd','a':3.0, 'period':[5,5,5], \
              'lammps_setting':lammps_setting}

              #'structure':'fcc',\
              #'positions':'rnd','a':4.2, 'period':[4,4,4]}

    sys = system.System(setting)
    rL = RLammps(sys)
    rL.run()

    (Lx, Ly, Lz, PotEng, Atoms) = rL.get_vals()
    print Lx, PotEng, Atoms

def test_03():

    lammps_setting = {'data_lmp':'data.lmp', 'in_lmp':'in.min' ,
                       'lammps_exe' :'/opt/lmpizarro/GitHub/lammps/src/lmp_serial'}


    import dbEamZhou as dbz
    parameters = dbz.parameters

    fccs = []
    bccs = []
    hcps = []
    for p in parameters:
        if parameters[p]['struct'] == 'fcc':
            fccs.append(p)
        if parameters[p]['struct'] == 'bcc':
            bccs.append(p)
        if parameters[p]['struct'] == 'hcp':
            hcps.append(p)

    for e in bccs:
        fccs.append(e)

    import periodictable as pt
    str_ = ''
    for e in fccs:

        form = pt.formula(e)
        a = form.structure[0][1].crystal_structure['a']


        setting ={'elements':[e], 'pot':'zhou', \
              'pca':[], 'nAtoms':250,\
              'structure':parameters[e]['struct'],\
              'positions':'rnd','a':4.2, 'period':[5,5,5], \
              'lammps_setting':lammps_setting}

        sys = system.System(setting)
        rL = RLammps(sys)
        rL.run()

        (Lx, Ly, Lz, PotEng, Atoms) = rL.get_vals()
        pc = 100  * (Lx - a) / a
        str_ += ('%s %f %f %f %f %f %f %f  \n')%(e, a, pc, Lx, Ly, Lz, PotEng, Atoms)

    print str_

def test_04():

    lammps_setting = {'data_lmp':'data.lmp', 
                      'in_lmp':'in.min',
                       'lammps_exe' :'/opt/lmpizarro/GitHub/lammps/src/lmp_serial'}

    setting ={'elements':['Al'], 'pot':'zhou', \
              'pca':[], 'nAtoms':250,\
              #'structure':'bcc',\
              #'positions':'rnd','a':3.0, 'period':[5,5,5]}

              'structure':'rnd',\
              'positions':'rnd','a':4.2, 'period':[5,5,5],\
              'lammps_setting':lammps_setting }

    sys = system.System(setting)
    calc = RLammps(sys)
    calc.run()

    (Lx, Ly, Lz, PotEng, Atoms) = calc.get_vals()


    print Lx, Ly, Lz, PotEng, Atoms



if __name__ == '__main__':
    #test_04()
    test_03()
    test_02()
    #test_01()
