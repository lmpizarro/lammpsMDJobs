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
        self.defInteraction = self.system.basicInteraction()
        self.in_frame = self.formatMultiline(self.in_frame)
        self.fix = self.formatMultiline(self.fix_min)
        self.interaction = self.formatMultiline(self.defInteraction)


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

        for k in self.keys_thermo:
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
    def __init__(self, settings):

        self.settings = settings
        elements = settings['elements']
        self.nTypes = len(elements) 

        self.nAt = settings['nAtoms']

        self.pos = None
        self.t1_ = None

        per = settings['period']
        a = settings['a']
        self.box =[[0, per[0]*a],[0, per[1]*a],[0, per[2]*a]]
        print self.settings['structure']
        if self.settings['structure'] == 'rnd':
            print 'rnd implemented'
            self.setRandomStructure()
        if self.settings['structure'] == 'fcc':
            self.setFcc()
            print 'fcc no implemented'
            sys.exit(0)
        if self.settings['structure'] == 'bcc':
            print 'bcc no implemented'
            sys.exit(0)
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


    def setRandomStructure(self):

        self.genRandomPositions()

        self.t1_ = []

        for i,e in enumerate(self.settings['nAt']):
            [self.t1_.append(i+1) for j in range(e)]
                

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

    def setFcc(self):
        px = self.settings['period'][0]
        py = self.settings['period'][1]
        pz = self.settings['period'][2]
        a =  self.settings['a']
        from ase.lattice.cubic import FaceCenteredCubic
        alloy = FaceCenteredCubic(directions=[[1,0,0], [0,1,0], [0,0,1]],
                                        size=(px,py,pz), symbol='Cu',
                pbc=(1,1,1), latticeconstant=a)

        print alloy
        print self.settings['pca']



def test_01():

    data_lmp = 'data.lmp'
    in_lmp = 'in.min'

    setting ={'elements':['Zr', 'Fe', 'Al', 'Mo'],\
              'pca':[10, 10, 10], 'nAtoms':250,\
              'structure':'fcc',\
              'positions':'rnd','a':3.0, 'period':[4,4,4]}

    sys = system.System(setting)

    rL = RunLammps(sys)
    rL.setDataLmp(data_lmp)
    rL.create_in(in_lmp)

    dL = DataLammps(setting)
    dL.genFile(data_lmp)

    import atomman.lammps as lmp
    output = lmp.run(lammps_exe, in_lmp, return_style='object')

    (Lx, PotEng, Atoms) = rL.get_vals('log.lammps')
    print Lx, PotEng, Atoms


if __name__ == '__main__':
    test_01()
