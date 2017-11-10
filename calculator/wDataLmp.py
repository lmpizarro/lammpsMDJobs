class DataLammps():
    def __init__(self, settings):

        self.data_lmp = settings['data_lmp']
        self.system = settings['sys']
        elements = self.system.setting['elements']

        self.nTypes = len(elements)
        self.nAtoms = self.system.bulk.get_number_of_atoms()

        self.system.update()

        self.pos = None
        self.t1_ = None


    def genFile(self):
        self.str_ = '\n\n'
        print self.nAtoms
        self.str_ = '\nxy xz yz\n'
        self.str_ += str(self.nAtoms) + ' atoms\n'
        self.str_ += str(self.nTypes) + ' atom types\n'
        xlo = self.system.box[0][0]
        xhi = self.system.box[0][1]
        ylo = self.system.box[1][0]
        yhi = self.system.box[1][1]
        zlo = self.system.box[2][0]
        zhi = self.system.box[2][1]
        self.str_ += str(float(xlo)) + ' ' +str(xhi)+ '  ' + ' xlo xhi\n'
        self.str_ += str(float(ylo)) + ' ' +str(yhi)+ '  ' + ' ylo yhi\n'
        self.str_ += str(float(zlo)) + ' ' +str(zhi)+ '  ' + ' zlo zhi\n'
        self.str_ +='\n'
        self.str_ +='Atoms\n'
        self.str_ +='\n'

        for i,e in enumerate(self.system.pos):
            self.str_ += str(i +1) + '  ' +   str(self.system.t1_[i]) + ' ' +\
                    str(e[0]) + ' ' + str(e[1])+ ' '  + str(e[2]) + '\n'
        with open(self.data_lmp, 'w') as inscript:
            inscript.write( self.str_)


        print 'Generated: ', self.data_lmp


