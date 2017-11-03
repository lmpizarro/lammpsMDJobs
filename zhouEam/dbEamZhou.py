import periodictable as pt
'''
re: r_e,   fe: f_e,    rhoe: rho_e_,  alfa: gamma_****, beta:omega_****,
A: A_,     B:B_,       kappa: kappa_, lambda: lambda_,   Fn0: F_n0_, 
Fn1:F_n1_, Fn2: F_n2_, Fn3:F_n3_,     F0:F_0_,          F1:F_1_, 
F2:F_2_,   F3:F_3_,    eta:eta_,      Fo:F_e

P(r): A, B, gamma_ (alfa), omega_ (beta), r_e, lambda, kappa_
f(r): f_e, omega_, r_e, lambda
F(f): rho_e, eta_, F_n0_,  F_n1_, F_n2_, F_n3_, F_0_, F_1_, F_2_, F_3_, Fo 

Atomic scale structure of sputtered metal multilayers
Acta Mater. 49 (2001) 4005-4015
X.W. Zhou, H.N.G.Wadley ...

Cu Ag Au Ni Pd Pt Al Pb
Fe Mo Ta W Mg Co Ti Zr

Computational study of the generation of crystal
defects in a bcc metal target irradiated 
by short laser pulses
Zhibin Lin, Robert A. Johnson ..

Physical Review B 77, 214108 (2008)

Cr

An n-body potential for Zr-Nb system based on the
embedded-atom method
De-ye Lin, S S Wang ...
J.Phys.: Condens. Matter 25 (2013) 105404 (14pp)
Zr-Nb
_
'''
'''
A_: afects crystal parameter in fcc
'''

'''
Interatomic potentials of the binary transition metal systems and
some applications in materials physics
J.H. Li, X.D. Dai,...
PhysicsReports455(2008)1-134
Pag. 50

http://www.knowledgedoor.com/2/elements_handbook/cohesive_energy.html

E:young modulus B: bulk modulus 
Sh: rigidity modulus (shear modulus)
Pr: Poisson ratio
'''
exp_props ={'Fe': {'a': 2.87, 'Ec': 4.28, 'Evf':1.79, 
                   'cijexp':{'C11':226, 'C12': 140, 'C44':116}},
            'V':{'a':3.03, 'Ec':5.31 , 
                'cijexp':{'C11':229, 'C12':140, 'C44':116}, 
                'Evf':2.2 },
            'Mo':{'a':3.1472, 'Ec':6.82, 'Evf':3.10, 
                'cijexp':{'C11':463.7, 'C12':157.8 , 'C44':109.2}},
            'Nb':{'a':3.3, 'Ec':7.57, 'Evf':2.75 , 
                'cijexp':{'C11':247, 'C12':135 , 'C44':28.7}},
            'Ta':{'a':3.3, 'Ec':8.1, 'Evf':2.18, 
                  'cijexp':{'C11':266.3, 'C12':158.2 , 'C44':87.4}},
            'W':{'a':3.16, 'Ec':8.9, 'Evf':3.95, 
                 'cijexp':{'C11':532, 'C12':204.9 , 'C44':163.1}},
            'Cu':{'a':3.61, 'Ec':3.49, 'Evf':1.28, 
                'cijexp':{'C11':168.4, 'C12':121.4 , 'C44':75.4}},
            'Ag':{'a':4.09 , 'Ec':2.95, 'Evf':1.1, 
                'cijexp':{'C11':124, 'C12':93.7 , 'C44':46.1}},
            'Au':{'a':4.08 , 'Ec':3.81, 'Evf':0.9, 
                'cijexp':{'C11':192.3, 'C12':163.1 , 'C44':42.0}},
            'Ni':{'a':3.52 , 'Ec':4.44, 'Evf':1.6, 
                  'cijexp':{'C11':245, 'C12':140 , 'C44':125}},
            'Pd':{'a':3.89, 'Ec':3.89, 'Evf':1.4, 
                  'cijexp':{'C11':227.1, 'C12':176.1 , 'C44':71.7}},
            'Pt':{'a':3.92, 'Ec':5.84, 'Evf':1.5, 
                  'cijexp':{'C11':347.0, 'C12':251 , 'C44':76.5}},
            'Al':{'a':4.05, 'Ec':3.34, 'Sh':26, 'B':76, 'Pr':0.35, 'E':70,
                'cijcalc':{'fcc':{'C11':101.0,'C12':61,'C44':25.4}},
                'cijexp':{'C11':114.0,'C12':61.9,'C44':31.6}}, 
            'Mg':{'Ec':1.53,
                'cijexp':{'C11':59.3,'C33':61.5, 'C44':16.4, 'C12':25.7,
                    'C13':21.4}
                },
            'Ti':{'Ec':1.87,
                'cijexp':{'C11':160,'C33':181, 'C44':46.5, 'C12':90, 'C13':66}
                },
            'Zr':{'Ec':6.136,
                'cijexp':{'C11':144,'C33':166, 'C44':33.4, 'C12':74, 'C13':67}
                },
            'Co':{'Ec':4.387,
                'cijexp':{'C11':295,'C33':335, 'C44':71, 'C12':159, 'C13':111}
                },
            'Pb':{'Ec':2.04,
                'cijexp':{'C11':55.5,'C12':45.4,'C44':19.4}
                },
            'Cr':{'Ec':4.1,
                'cijexp':{'C11':391,'C12':90.0,'C44':103.2}
                    },
            'Y':{'Ec':4.387,
                'cijexp':{'C11':77.9,'C33':76.9, 'C44':24.3, 'C12':29.2,
                    'C13':20}
                    },
            'Si':{'a':5.430099,'Ec':4.64, 'Sh':'no-data', 'B':100, 
                'Pr':'no-data', 'E':47, 
                'cijcalc':{'fcc':{'C11':64.4,'C12':87.2,'C44':4.7}},
                'cijexp':{'C11':166.0,'C12':63.9,'C44':79.6}
                    },
            'U':{'Ec':5.405,

                    'cijexp':{'C11':13.6, 'C44':165.6, 'C12':20.2},
                    'cijcalc':{'fcc':{'C11':13.6, 'C44':165.6, 'C12':20.2}}
                }
            }


parameters={
        'Ni': {'eta_': 0.469, 'number': 28, 'omega_': 4.471175, 'F_1_': 0.0, 'B_': 0.633531, 
            'F_0_': -2.7, 'index': 3, 'F_3_': -0.152856, 'kappa_': 0.443599, 'r_e': 2.488746, 
            'F_n2_': 0.241442, 'symbol': 'Ni', 'rho_e_': 27.562015, 'F_2_': 0.26539, 
            'F_n0_': -2.693513, 'F_e_': -2.699486, 'A_': 0.429046, 'rho_s_': 27.93041, 
            'f_e': 2.007018, 'F_n3_': -2.375626, 'F_n1_': -0.076445, 'name': u'Nickel', 
            'gamma_': 8.383453, 'lambda_': 0.820658, 'mass': 58.6934,
            'struct':'fcc'},
        'Mg': {'eta_': 0.431425, 'number': 12, 'omega_': 5.455311, 'F_1_': 0.0, 'B_': 0.22593, 
            'F_0_': -0.9, 'index': 12, 'F_3_': -0.22601, 'kappa_': 0.5, 'r_e': 3.196291, 
            'F_n2_': 0.162232, 'symbol': 'Mg', 'rho_e_': 7.1326, 'F_2_': 0.122838, 
            'F_n0_': -0.896473, 'F_e_': -0.899702, 'A_': 0.137518, 'rho_s_': 7.1326, 
            'f_e': 0.544323, 'F_n3_': -0.68995, 'F_n1_': -0.044291, 'name': u'Magnesium', 
            'gamma_': 10.228708, 'lambda_': 1.0, 'mass': 24.305, 'struct':'hcp'},
        'Co': {'eta_': 0.694608, 'number': 27, 'omega_': 4.629134, 'F_1_': 0.0, 'B_': 0.640107, 
            'F_0_': -2.56, 'index': 13, 'F_3_': -0.68714, 'kappa_': 0.5, 'r_e': 2.505979, 
            'F_n2_': 0.733381, 'symbol': 'Co', 'rho_e_': 27.206789, 'F_2_': 0.705845, 
            'F_n0_': -2.541799, 'F_e_': -2.559307, 'A_': 0.421378, 'rho_s_': 27.206789, 
            'f_e': 1.975299, 'F_n3_': -1.589003, 'F_n1_': -0.219415, 'name': u'Cobalt', 
            'gamma_': 8.679625, 'lambda_': 1.0, 'mass': 58.933195, 'struct':'hcp'},
        'Ag': {'eta_': 0.783924, 'number': 47, 'omega_': 4.870405, 'F_1_': 0.0, 'B_': 0.419611, 
            'F_0_': -1.75, 'index': 1, 'F_3_': -1.15065, 'kappa_': 0.33971, 'r_e': 2.891814, 
            'F_n2_': 0.91205, 'symbol': 'Ag', 'rho_e_': 14.6041, 'F_2_': 0.744561, 
            'F_n0_': -1.729364, 'F_e_': -1.748423, 'A_': 0.277758, 'rho_s_': 14.604144, 
            'f_e': 1.106232, 'F_n3_': -0.561432, 'F_n1_': -0.255882, 'name': u'Silver', 
            'gamma_': 9.13201, 'lambda_': 0.750758, 'mass': 107.8682,
            'struct':'fcc'},
        'Pt': {'eta_': 1.45, 'number': 78, 'omega_': 3.78975, 'F_1_': 0.0, 'B_': 0.696037, 
            'F_0_': -4.17, 'index': 5, 'F_3_': -2.420128, 'kappa_': 0.385255, 'r_e': 2.771916, 
            'F_n2_': 0.528491, 'symbol': 'Pt', 'rho_e_': 33.367564, 'F_2_': 3.010561, 
            'F_n0_': -1.455568, 'F_e_': -4.145597, 'A_': 0.556398, 'rho_s_': 35.205357, 
            'f_e': 2.336509, 'F_n3_': 1.222875, 'F_n1_': -2.149952, 'name': u'Platinum', 
            'gamma_': 7.105782, 'lambda_': 0.77051, 'mass': 192.084, 'struct':'fcc'},
        'W': {'eta_': -0.582714, 'number': 74, 'omega_': 4.746728, 'F_1_': 0.0, 'B_': 1.394592, 
            'F_0_': -4.96, 'index': 11, 'F_3_': 0.348147, 'kappa_': 0.139209, 'r_e': 2.74084, 
            'F_n2_': 0.365057, 'symbol': 'W', 'rho_e_': 37.234847, 'F_2_': 0.661935, 
            'F_n0_': -4.946281, 'F_e_': -4.961306, 'A_': 0.882435, 'rho_s_': 37.234847, 
            'f_e': 3.48734, 'F_n3_': -4.432406, 'F_n1_': -0.148818, 'name': u'Tungsten', 
            'gamma_': 8.900114, 'lambda_': 0.278417, 'mass': 183.84,
            'struct':'bcc'},
        'Mo': {'eta_': 0.790879, 'number': 42, 'omega_': 4.47655, 'F_1_': 0.0, 'B_': 1.120373, 
            'F_0_': -3.71, 'index': 9, 'F_3_': 0.776222, 'kappa_': 0.13764, 'r_e': 2.7281, 
            'F_n2_': 0.38045, 'symbol': 'Mo', 'rho_e_': 29.354065, 'F_2_': 0.875874, 
            'F_n0_': -3.692913, 'F_e_': -3.712093, 'A_': 0.708787, 'rho_s_': 29.354065, 
            'f_e': 2.72371, 'F_n3_': -3.13365, 'F_n1_': -0.178812, 'name': u'Molybdaenum', 
            'gamma_': 8.393531, 'lambda_': 0.27528, 'mass': 95.96, 'struct':'bcc'},
        'Al': {'eta_': 0.785902, 'number': 13, 'omega_': 3.527021, 'F_1_': 0.0, 'B_': 0.365551, 
            'F_0_': -2.83, 'index': 6, 'F_3_': -2.488244, 'kappa_': 0.379846, 'r_e': 2.863924, 
            'F_n2_': 1.258562, 'symbol': 'Al', 'rho_e_': 20.418205, 'F_2_': 0.622245, 
            'F_n0_': -2.807602, 'F_e_': -2.824528, 'A_': 0.314873, 'rho_s_': 23.19574, 
            'f_e': 1.403115, 'F_n3_': -1.247604, 'F_n1_': -0.301435, 'name': u'Aluminium', 
            'gamma_': 6.613165, 'lambda_': 0.759692, 'mass': 26.9815386,
            'struct':'fcc'},
        'Pb': {'eta_': 0.93538, 'number': 82, 'omega_': 5.212457, 'F_1_': 0.0, 'B_': 0.236884, 
            'F_0_': -1.44, 'index': 7, 'F_3_': -0.538766, 'kappa_': 0.250805, 'r_e': 3.499723, 
            'F_n2_': 0.682886, 'symbol': 'Pb', 'rho_e_': 8.450154, 'F_2_': 0.702726, 
            'F_n0_': -1.42237, 'F_e_': -1.439436, 'A_': 0.161219, 'rho_s_': 8.450063, 
            'f_e': 0.647872, 'F_n3_': -0.529378, 'F_n1_': -0.210107, 'name': u'Lead', 
            'gamma_': 9.121799, 'lambda_': 0.764955, 'mass': 207.2,
            'struct':'fcc'},
        'Zr': {'eta_': 0.597133, 'number': 40, 'omega_': 4.564902, 'F_1_': 0.0, 'B_': 0.640054, 
             'F_0_': -4.51, 'index': 15, 'F_3_': -0.98187, 'kappa_': 0.5, 'r_e': 3.199978, 
             'F_n2_': 0.990148, 'symbol': 'Zr', 'rho_e_': 30.879991, 'F_2_': 0.928602, 
             'F_n0_': -4.485793, 'F_e_': -4.509025, 'A_': 0.424667, 'rho_s_': 30.879991, 
             'f_e': 2.230909, 'F_n3_': -3.202516, 'F_n1_': -0.293129, 'name': u'Zirkonium', 
             'gamma_': 8.55919, 'lambda_': 1.0, 'mass': 91.224, 'struct':'hcp'},
        'Au': {'eta_': 1.021095, 'number': 79, 'omega_': 5.075228, 'F_1_': 0.0, 'B_': 0.356666, 
             'F_0_': -2.98, 'index': 2, 'F_3_': -1.134778, 'kappa_': 0.35657, 'r_e': 2.885034, 
             'F_n2_': 1.601954, 'symbol': 'Au', 'rho_e_': 19.991632, 'F_2_': 1.706587, 
             'F_n0_': -2.937772, 'F_e_': -2.978815, 'A_': 0.229762, 'rho_s_': 19.991509, 
             'f_e': 1.529021, 'F_n3_': -0.83553, 'F_n1_': -0.500288, 'name': u'Gold', 
             'gamma_': 9.516052, 'lambda_': 0.748798, 'mass': 196.966569,
             'struct':'fcc'},
        'Fe': {'eta_': 0.39175, 'number': 26, 'omega_': 5.236411, 'F_1_': 0.0, 'B_': 0.646243, 
             'F_0_': -2.54, 'index': 8, 'F_3_': -0.14877, 'kappa_': 0.170306, 'r_e': 2.481987, 
             'F_n2_': 0.193065, 'symbol': 'Fe', 'rho_e_': 20.041463, 'F_2_': 0.200269, 
             'F_n0_': -2.534992, 'F_e_': -2.539945, 'A_': 0.392811, 'rho_s_': 20.041463, 
             'f_e': 1.885957, 'F_n3_': -2.282322, 'F_n1_': -0.059605, 'name': u'Iron', 
             'gamma_': 9.81827, 'lambda_': 0.340613, 'mass': 55.845,
             'struct':'bcc'},
        'Pd': {'eta_': 1.13, 'number': 46, 'omega_': 4.638612, 'F_1_': 0.0, 'B_': 0.59888, 
             'F_0_': -2.36, 'index': 4, 'F_3_': -1.675615, 'kappa_': 0.397263, 'r_e': 2.750897, 
             'F_n2_': 1.615343, 'symbol': 'Pd', 'rho_e_': 21.335246, 'F_2_': 1.481742, 
             'F_n0_': -2.321006, 'F_e_': -2.352753, 'A_': 0.406763, 'rho_s_': 21.940073, 
             'f_e': 1.595417, 'F_n3_': -0.231681, 'F_n1_': -0.473983, 'name': u'Palladium', 
             'gamma_': 8.697397, 'lambda_': 0.754799, 'mass': 106.42,
             'struct':'fcc'},
        'Ti': {'eta_': 0.558572, 'number': 22, 'omega_': 4.68023, 'F_1_': 0.0, 'B_': 0.570968, 
             'F_0_': -3.22, 'index': 14, 'F_3_': -0.75071, 'kappa_': 0.5, 'r_e': 2.933872, 
             'F_n2_': 0.683779, 'symbol': 'Ti', 'rho_e_': 25.565138, 'F_2_': 0.608587, 
             'F_n0_': -3.203773, 'F_e_': -3.219176, 'A_': 0.373601, 'rho_s_': 25.565138, 
             'f_e': 1.8632, 'F_n3_': -2.321732, 'F_n1_': -0.198262, 'name': u'Titanium', 
             'gamma_': 8.775431, 'lambda_': 1.0, 'mass': 47.867, 'struct':'hcp'},
        'Cr': {'eta_': 0.456, 'number': 24, 'omega_': 7.170494, 'F_1_': 0.0, 'B_': 1.827556, 
             'F_0_': -2.02, 'F_n2_': -0.183611, 'F_3_': 0.439144, 'kappa_': 0.18533, 
             'r_e': 2.493879, 'symbol': 'Cr', 'rho_e_': 17.641302, 'F_2_': -0.056517, 
             'F_n0_': -2.022754, 'F_e_': -2.020038, 'A_': 1.551848, 'rho_s_': 19.60545, 
             'f_e': 1.793835, 'F_n3_': -2.245972, 'F_n1_': 0.039608, 'name': u'Chromium', 
             'gamma_': 8.604593, 'lambda_': 0.277995, 'mass': 51.9961,
             'struct':'bcc', 'index':20},
        'Nb': {'eta_': 0.890133, 'number': 41, 'omega_': 4.5, 'F_1_': 0.0, 'B_': 0.905874, 
             'F_0_': -4.975703, 'F_n2_': 1.680064, 'F_3_': -0.765504, 'kappa_': 0.136112, 
             'r_e': 2.85823, 'symbol': 'Nb', 'rho_e_': 31.7449, 'F_n1_': -0.549044, 
             'F_n0_': -4.92855, 'F_e_': -4.975568, 'A_': 0.56965, 'rho_s_': 31.7449, 
             'f_e': 2.889832, 'F_n3_': -2.699442, 'F_2_': 1.980875, 'name': u'Niobium', 
             'gamma_': 7.514546, 'lambda_': 0.388893, 'mass': 92.90638,
             'index':21, 'struct': 'bcc'},
        'Cu': {'eta_': 0.31049, 'number': 29, 'omega_': 4.334731, 'F_1_': 0.0, 'B_': 0.548085, 
             'F_0_': -2.19, 'index': 0, 'F_3_': -2.100595, 'kappa_': 0.308782, 'r_e': 2.556162, 
             'F_n2_': 1.088878, 'symbol': 'Cu', 'rho_e_': 21.175871, 'F_2_': 0.56183, 
             'F_n0_': -2.170269, 'F_e_': -2.186568, 'A_': 0.39662, 
             'rho_s_': 21.175395, 'f_e': 1.554485, 'F_n3_': -0.817603, 'F_n1_': -0.263788, 
             'name': u'Copper', 'gamma_': 8.12762, 'lambda_': 0.756515, 'mass':
             63.546, 'struct':'fcc'},
        'Ta': {'eta_': 0.848843, 'number': 73, 'omega_': 4.527748, 'F_1_': 0.0, 'B_': 1.032101, 
             'F_0_': -5.14, 'index': 10, 'F_3_': 0.221375, 'kappa_': 0.176977, 'r_e': 2.860082, 
             'F_n2_': 1.112997, 'symbol': 'Ta', 'rho_e_': 33.787168, 'F_2_': 1.640098, 
             'F_n0_': -5.103845, 'F_e_': -5.141526, 'A_': 0.611679, 'rho_s_': 33.787168, 
             'f_e': 3.086341, 'F_n3_': -3.585325, 'F_n1_': -0.405524, 'name': u'Tantalum', 
             'gamma_': 8.489528, 'lambda_': 0.353954, 'mass': 180.94788,
             'struct':'bcc'}}

class EamZhou():
    def __init__(self):
        self.keys_pot =['eta_', 'number', 'omega_', 'F_1_', 'B_',\
                'F_0_', 'index', 'F_3_', 'kappa_', 'r_e',\
                'F_n2_', 'symbol', 'rho_e_', 'F_n1_', 'F_n0_', 'F_e_',\
                'A_', 'rho_s_', 'f_e', 'F_n3_', 'F_2_',\
                'name','gamma_','lambda_', 'mass',]


    def __str__(self):
        str_ = ' '*18
        for s in self.symbols:
            str_ += '{:14}'.format(s)
        str_ +='\n' 

        for c in self.keys_pot:
            if c != 'name' and c != 'symbol':
                str_ +='{:8}'.format(c)
                for s in self.symbols:
                    str_ +='{:14f}'.format(float(parameters[s][c]))
                str_ +='\n'

        return str_

    def printAll(self):
        symbols =[k for k in parameters]

        t = []
        c = 0
        for s in symbols:
            t.append(s)
            c +=1
            if c > 7:
                self.symbols = t
                print self
                self.symbols = None
                c=0
                t=[]
        self.symbols = t
        print self



def test2():
    ez = EamZhou()

    ez.printAll()

def latticeA(el):
    form = pt.formula(el)
    e_ = form.structure[0][1]
    crys = e_.crystal_structure['symmetry'] 
    a_ = e_.crystal_structure['a'] 

    return  a_

def miningEpsSig(element):
    cijsi = exp_props[element]

    a0 = latticeA(element)
    Ec0 = cijsi['Ec']

    E=''
    S=1e6
    for e in parameters:
        sum_ = 0
        Ac = ((a0 - latticeA(e)) / a0) **2
        Ec = ((exp_props[e]['Ec'] - Ec0) / Ec0) ** 2
        sum_ += Ec + Ac
        if sum_ < S:
            S = sum_
            E = e
    return E


def miningCij(element, cijs=[ 'C11','C12','C44']):
    cijsi = exp_props[element]

    keys = cijs 

    E=''
    S=1e6
    for e in parameters:
        sum_ = 0
        for k in keys:
            sum_ += ((exp_props[e]['cijexp'][k] - cijsi['cijexp'][k]) / \
                    cijsi['cijexp'][k])**2
        if sum_ < S:
            S = sum_
            E = e
    return   E

def test3():
    element = 'Si'
    ec = miningCij(element)
    ees = miningEpsSig(element)

    print element, ' elastic :', ec, 'sigma epsilon: ', ees


def test4():
    from periodictable import formula
    elsInOrder =['Mg', 'Al', 'Ti', 'Cr', 'Fe','Ni', 'Co', 
                 'Cu', 'Zr', 'Nb', 'Mo', 'Pd', 'Ag', 'Ta', 
                 'W', 'Pt', 'Au', 'Pb']
    print len(elsInOrder)
    for e in elsInOrder:
        print e, formula(e).mass, parameters[e]['A_'], parameters[e]['r_e']

if __name__ == '__main__':
    test3()
    test4()
