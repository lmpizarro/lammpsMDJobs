import sys
sys.path.append('../../pizza/src')
from log import log

import system
import random

lammps_exe ='/opt/lmpizarro/GitHub/lammps/src/lmp_serial'

class RunLammps():
    def __init__(self, system):
        self.system = system
        self.keys_thermo = ['Step', 'Press','PotEng','TotEng', 'Lx', 'Ly','Lz','Atoms'] 

        self.atoms = system.getAtoms()

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

        return (Lx, Ly, Lz, PotEng / Atoms, Atoms)

    def formatMultiline(self, multiline):
        l = multiline.split('\n')
        str_=''
        for e in l:
           str_ +=e.lstrip() + '\n'
        return str_

    def create_in(self, fileName):
        mass = self.system.getMasess()

        self.in_lmp = fileName

        in_lammps = self.in_frame % (self.data_lmps, self.interaction, mass, self.fix)

        with open(fileName, 'w') as inscript:
             inscript.write( in_lammps )

    def setDataLmp(self, name):
        self.data_lmps = name

    def setInteraction (self, interaction):
        self.interaction = interaction

    def setFix(self, fix):
        self.fix = fix

    def run(self):
        import atomman.lammps as lmp
        output = lmp.run(lammps_exe, self.in_lmp, return_style='object')


class DataLammps():
    def __init__(self, sys):

        self.sys = sys
        self.settings = sys.setting

        elements = self.settings['elements']
        self.nTypes = len(elements) 
      
        self.nAt = self.settings['nAtoms']

        self.pos = None
        self.t1_ = None

        per = self.settings['period']
        a = self.settings['a']
        self.box =[[0, per[0]*a],[0, per[1]*a],[0, per[2]*a]]

        if self.settings['structure'] == 'rnd':
            print 'rnd implemented'
            self.genRandomPositions()

            self.genStructure()
 
            self.setRandomStructure()
        if self.settings['structure'] == 'fcc':
            self.setCrystal('fcc')
            print 'fcc implemented'

            if self.settings['positions'] == 'rnd':
                self.setRandomStructure()
        if self.settings['structure'] == 'bcc':
            self.setCrystal('bcc')

            if self.settings['positions'] == 'rnd':
                self.setRandomStructure()
            print 'bcc no implemented'
            #sys.exit(0)
        if self.settings['structure'] == 'hcp':
            print 'hcp no implemented'
            sys.exit(0)

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

        print 'Generated: ', fileName


    def genStructure(self):
        self.t1_ = []
        for i,e in enumerate(self.settings['nAt']):
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
        Lx = self.box[0][1]
        Ly = self.box[1][1]
        Lz = self.box[2][1]

        self.pos =[]
        for i in range(self.nAt):
            x = random.random() * Lx
            y = random.random() * Ly
            z = random.random() * Lz
            self.pos.append([x,y,z])

    def setCrystal(self, crys):
        px = self.settings['period'][0]
        py = self.settings['period'][1]
        pz = self.settings['period'][2]
        a =  self.settings['a']

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

        self.settings['nAtoms'] =  alloy.get_number_of_atoms()

        self.sys.calcAtoms2()

        self.genStructure()

        self.pos = alloy.get_positions()
        self.nAt = self.settings['nAtoms']


def test_01():

    data_lmp = 'data.lmp'
    in_lmp = 'in.min'

    setting ={'elements':['Al'], 'pot':'zhou', \
              'pca':[], 'nAtoms':250,\
              #'structure':'bcc',\
              #'positions':'rnd','a':3.0, 'period':[5,5,5]}

              'structure':'fcc',\
              'positions':'rnd','a':4.2, 'period':[5,5,5]}

    sys = system.System(setting)

    rL = RunLammps(sys)
    rL.setDataLmp(data_lmp)
    rL.create_in(in_lmp)

    dL = DataLammps(sys)
    dL.genFile(data_lmp)

    rL.run()

    (Lx, Ly, Lz, PotEng, Atoms) = rL.get_vals()
    print Lx, Ly, Lz, PotEng, Atoms


def test_02():

    data_lmp = 'data.lmp'
    in_lmp = 'in.min'

    setting ={'elements':['Zr', 'Fe', 'Al'], 'pot':'zhou', \
              'pca':[10, 10], 'nAtoms':250,\
              'structure':'bcc',\
              'positions':'rnd','a':3.0, 'period':[5,5,5]}

              #'structure':'fcc',\
              #'positions':'rnd','a':4.2, 'period':[4,4,4]}

    sys = system.System(setting)

    rL = RunLammps(sys)

    rL.setDataLmp(data_lmp)
    rL.create_in(in_lmp)

    dL = DataLammps(sys)
    dL.genFile(data_lmp)

    rL.run()

    (Lx, PotEng, Atoms) = rL.get_vals()
    print Lx, PotEng, Atoms

def test_03():
    import dbEamZhou as dbz

    data_lmp = 'data.lmp'
    in_lmp = 'in.min'


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
              'positions':'rnd','a':4.2, 'period':[5,5,5]}

        sys = system.System(setting)

        rL = RunLammps(sys)
        rL.setDataLmp(data_lmp)
        rL.create_in(in_lmp)

        dL = DataLammps(sys)
        dL.genFile(data_lmp)

        rL.run()

        (Lx, Ly, Lz, PotEng, Atoms) = rL.get_vals()
        pc = 100  * (Lx - a) / a
        str_ += ('%s %f %f %f %f %f %f %f  \n')%(e, a, pc, Lx, Ly, Lz, PotEng, Atoms)

    print str_


if __name__ == '__main__':
    test_03()
