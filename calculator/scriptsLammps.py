import sys

sys.path.append('../zhouEam')
import ljParameters
import EamZhouPotential as zhou
import system as sys


class SLammps():
    def __init__(self, settings):

        self.setting = settings

        self.data_lmps = self.setting['data_lmp']
        self.in_lmp =  self.setting['in_lmp']
        self.system = self.setting['sys']


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

        self.fix_nve='''
            fix 1 all nve
        '''



        self.defInteraction = self.system.Interaction()
        self.in_frame = self.formatMultiline(self.in_frame)

        if settings['minimize'] == True:
            self.fix = self.formatMultiline(self.fix_min)
        else:
            self.fix = self.formatMultiline(self.fix_nve)


        self.interaction = self.formatMultiline(self.defInteraction)

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

def test_01():

    sys_setting ={'elements':['Al'], 'pot':'zhou', \
              'pca':[], 'nAtoms':250,\
              #'structure':'bcc',\
              #'positions':'rnd','a':3.0, 'period':[5,5,5]}

              'structure':'rnd',\
              'positions':'rnd','a':4.2, 'period':[5,5,5]}

    sys1 = sys.System(sys_setting)

    lammps_setting = {'data_lmp':'data.lmp', 
                      'in_lmp':'in.min',
                      'lammps_exe':\
                              '/opt/lmpizarro/GitHub/lammps/src/lmp_serial',
                              'log': 'log.lammps',
                              'sys':sys1}


    sLmp = SLammps(lammps_setting)
    sLmp.create_in()
    pass

if __name__ == '__main__':
    test_01()