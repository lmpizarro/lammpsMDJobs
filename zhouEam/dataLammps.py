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
        print self.settings['nAtoms']
        self.str_ += str(self.settings['nAtoms']) + ' atoms\n'
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

        print '----', len(self.sys.t1_)
        for i,e in enumerate(self.sys.pos):
            self.str_ += str(i +1) + '  ' +   str(self.sys.t1_[i]) + ' ' +\
                    str(e[0]) + ' ' + str(e[1])+ ' '  + str(e[2]) + '\n'
        with open(self.data_lmp, 'w') as inscript:
            inscript.write( self.str_)

        print 'Generated: ', self.data_lmp


