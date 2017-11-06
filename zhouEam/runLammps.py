import sys
sys.path.append('../../pizza/src')
from log import log

from ase.calculators.calculator import PropertyNotImplementedError

import system

import numpy as np
from ase.units import GPa



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



    def get_vals(self):
        lg = log(self.log)
        status = {}

        for k in self.keys_thermo:
            a =  lg.get(k)
            i = len(a) - 1 
            status[k] = a[i]

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

    def setDataLmp(self, name):
        self.data_lmps = name

    def setInteraction (self, interaction):
        self.interaction = interaction

    def setFix(self, fix):
        self.fix = fix

    def run(self):
        import atomman.lammps as lmp
        output = lmp.run(self.lammps_exe, self.in_lmp, return_style='object')


class DataLammps():
    def __init__(self, sys):

        self.sys = sys
        self.settings = sys.setting

        self.data_lmp = self.settings['lammps_setting']['data_lmp']

        elements = self.settings['elements']
        self.nTypes = len(elements) 
      
        self.nAt = self.settings['nAtoms']

        self.pos = None
        self.t1_ = None


    def genFile(self):
        self.str_ = '\n'
        self.str_ += str(self.nAt) + ' atoms\n'
        self.str_ += str(self.nTypes) + ' atom types\n'
        xlo = self.sys.box[0][0]
        xhi = self.sys.box[0][1]
        ylo = self.sys.box[1][0]
        yhi = self.sys.box[1][1]
        zlo = self.sys.box[2][0]
        zhi = self.sys.box[2][1]
        self.str_ += str(float(xlo)) + ' ' +str(xhi)+ '  ' + ' xlo xhi\n'
        self.str_ += str(float(ylo)) + ' ' +str(yhi)+ '  ' + ' ylo yhi\n'
        self.str_ += str(float(zlo)) + ' ' +str(zhi)+ '  ' + ' zlo zhi\n'
        self.str_ +='\n'
        self.str_ +='Atoms\n'
        self.str_ +='\n'

        for i,e in enumerate(self.sys.pos):
            self.str_ += str(i +1) + '  ' +   str(self.sys.t1_[i]) + ' ' +\
                    str(e[0]) + ' ' + str(e[1])+ ' '  + str(e[2]) + '\n'
        with open(self.data_lmp, 'w') as inscript:
            inscript.write( self.str_)

        print 'Generated: ', self.data_lmp

def test_04():
    data_lmp = 'data.lmp'
    in_lmp = 'in.min'

    setting ={'elements':['Al'], 'pot':'zhou', \
              'pca':[], 'nAtoms':250,\
              #'structure':'bcc',\
              #'positions':'rnd','a':3.0, 'period':[5,5,5]}

              'structure':'fcc',\
              'positions':'rnd','a':4.2, 'period':[5,5,5]}

    sys = system.System(setting)

    calc = RLammps(sys)


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
    rL.create_in()

    dL = DataLammps(sys)
    dL.genFile()

    rL.run()


    (Lx, Ly, Lz, PotEng, Atoms) = rL.get_vals()

    rL.parse_log()
    print rL.get_stress()
    print rL.get_potential_energy()
    print Lx, Ly, Lz, PotEng, Atoms


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

    rL.create_in()

    dL = DataLammps(sys)
    dL.genFile()

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
        rL.create_in()

        dL = DataLammps(sys)
        dL.genFile()

        rL.run()

        (Lx, Ly, Lz, PotEng, Atoms) = rL.get_vals()
        pc = 100  * (Lx - a) / a
        str_ += ('%s %f %f %f %f %f %f %f  \n')%(e, a, pc, Lx, Ly, Lz, PotEng, Atoms)

    print str_


if __name__ == '__main__':
    test_03()
