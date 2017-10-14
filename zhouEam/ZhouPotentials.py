from periodic.table import element
from atsim.potentials import writeSetFL
from atsim.potentials import Potential
from atsim.potentials import EAMPotential

import math


def parametersDB():
  niobiumParams = {'symbol': 'Nb', \
                 'r_e':2.858230,'f_e':2.889832, 'rho_e_':31.744900,'rho_s_':31.744900, 'kappa_':0.136112,\
                 'lambda_':0.388893, 'A_':0.569650, 'B_':0.905874, 'gamma_':7.514546,  'omega_':4.5 ,\
                 'F_0_': -4.975703,'F_1_':0.0, 'F_2_':1.980875,'F_3_':-0.765504, 'eta_':0.890133,\
                 'F_n0_': -4.928550,'F_n1_':-0.549044, 'F_n2_':1.680064,'F_n3_':-2.699442,'F_e_': -4.975568
           }

  chromiunParams = {'symbol': 'Cr', \
                 'r_e':2.493879,'f_e':1.793835, 'rho_e_':17.641302,'rho_s_':19.60545, 'kappa_':0.18533,\
                 'lambda_':0.277995, 'A_':1.551848, 'B_':1.827556, 'gamma_':8.604593,  'omega_':7.170494 ,\
                 'F_n0_': -2.022754,'F_n1_':0.039608, 'F_n2_':-0.183611,'F_n3_':-2.245972, 'eta_':0.456,\
                 'F_0_': -2.02,'F_1_':0.0, 'F_2_':-0.056517,'F_3_':0.439144,'F_e_': -2.020038
           }

  matriz = ''' 
element Cu Ag Au Ni Pd Pt Al Pb  Fe Mo Ta W Mg Co Ti Zr
r_e     2.556162  2.891814 2.885034 2.488746 2.750897 2.771916 2.863924 3.499723 2.481987 2.728100 2.860082 2.740840 3.196291 2.505979 2.933872 3.199978
f_e     1.554485  1.106232 1.529021 2.007018 1.595417 2.336509 1.403115 0.647872 1.885957 2.723710 3.086341 3.487340 0.544323 1.975299 1.863200 2.230909
rho_e_  21.175871 14.604100 19.991632 27.562015 21.335246 33.367564 20.418205 8.450154 20.041463 29.354065 33.787168 37.234847 7.132600 27.206789 25.565138 30.879991
rho_s_  21.175395 14.604144 19.991509 27.930410 21.940073 35.205357 23.195740 8.450063 20.041463 29.354065 33.787168 37.234847 7.132600 27.206789 25.565138 30.879991
gamma_   8.127620  9.132010  9.516052  8.383453  8.697397  7.105782 6.613165  9.121799 9.818270 8.393531 8.489528 8.900114 10.228708 8.679625 8.775431 8.559190
omega_   4.334731  4.870405  5.075228  4.471175  4.638612  3.789750 3.527021  5.212457 5.236411 4.476550 4.527748 4.746728  5.455311 4.629134 4.680230 4.564902 
A_      0.396620  0.277758  0.229762  0.429046  0.406763  0.556398 0.134873  0.161219 0.392811 0.708787 0.611679 0.882435 0.137518 0.421378 0.373601 0.424667 
B_      0.548085  0.419611  0.356666  0.633531  0.598880  0.696037 0.365551  0.236884 0.646243 1.120373 1.032101 1.394592 0.225930 0.640107 0.570968 0.640054 
kappa_   0.308782  0.339710  0.356570  0.443599  0.397263  0.385255 0.379846  0.250805 0.170306 0.137640 0.176977 0.139209 0.5 0.5 0.5 0.5
lambda_  0.756515  0.750758  0.748798  0.820658  0.754799  0.770510 0.759692  0.764955 0.340613 0.275280 0.353954 0.278417 1.0 1.0 1.0 1.0
F_n0_   -2.170269 -1.729364 -2.937772 -2.693513 -2.321006 -1.455568 -2.807602 -1.422370 -2.534992 -3.692913 -5.103845 -4.946281 -0.896473 -2.541799 -3.203773 -4.485793
F_n1_   -0.263788 -0.255882 -0.500288 -0.076445 -0.473983 -2.149952 -0.301435 -0.210107 -0.059605 -0.178812 -0.405524 -0.148818 -0.044291 -0.219415 -0.198262 -0.293129 
F_n2_    1.088878  0.912050 1.601954 0.241442 1.615343 0.528491 1.258562 0.682886 0.193065 0.380450 1.112997 0.365057 0.162232 0.733381 0.683779 0.990148
F_n3_   -0.817603 -0.561432 -0.835530 -2.375626 -0.231681 1.222875 -1.247604 -0.529378 -2.282322 -3.133650 -3.585325 -4.432406 -0.689950 -1.589003 -2.321732 -3.202516
F_0_    -2.19     -1.75 -2.98 -2.70 -2.36 -4.17 -2.83 -1.44  -2.54 -3.71 -5.14 -4.96 -0.90 -2.56 -3.22 -4.51
F_1_     0         0 0 0 0 0 0 0  0 0 0 0 0 0 0 0 
F_2_     0.561830  0.744561 1.706587 0.265390 1.481742 3.010561 0.622245 0.702726 0.200269 0.875874 1.640098 0.661935 0.122838 0.705845 0.608587 0.928602 
F_3_    -2.100595 -1.150650 -1.134778 -0.152856 -1.675615 -2.420128 -2.488244 -0.538766 -0.148770 0.776222 0.221375 0.348147 -0.226010 -0.687140 -0.750710 -0.981870
eta_    0.310490  0.783924 1.021095 0.469000 1.130000 1.450000 0.785902 0.935380 0.391750 0.790879 0.848843 -0.582714 0.431425 0.694608 0.558572 0.597133
F_e_    -2.186568 -1.748423 -2.978815 -2.699486 -2.352753 -4.145597 -2.824528 -1.439436 -2.539945 -3.712093 -5.141526 -4.961306 -0.899702 -2.559307 -3.219176 -4.509025
'''

  data = matriz.split('\n')

  lines = []
  for d in data:
    line = (d.split())
    if len(line) > 0:
        lines.append(line)
        
  parameters = {}

  for i,el in enumerate(lines[0][1:]):
    parameters[el] = {'index':i}


  for k in parameters.keys():
    parameters[k]['symbol'] = k
    index = parameters[k]['index']
    for l in lines[1:]:
        parameters[k][l[0]] = float(l[index  + 1])
        
  parameters['Nb'] = niobiumParams
  parameters['Cr'] = chromiunParams

  for k in parameters:
    a_1 = element(k) 
    parameters[k]['mass'] = a_1.mass
    parameters[k]['number'] = a_1.atomic
    parameters[k]['name'] = a_1.name
    
  return parameters

#parameters = parametersDB()

from dbEamZhou import parameters

def makeFunc(a, b, r_e, c):
  # Creates functions of the form used for density function.
  # Functional form also forms components of pair potential.
  def func(r):
    return (a * math.exp(-b*(r/r_e -1)))/(1+(r/r_e - c)**20.0)
  return func


def makePairPotAA(A, gamma, r_e, kappa,
                  B, omega, lamda):
  # Function factory that returns functions parameterised for homogeneous pair interactions
  f1 = makeFunc(A, gamma, r_e, kappa)
  f2 = makeFunc(B, omega, r_e, lamda)
  def func(r):
    return f1(r) - f2(r)
  return func


def makePairPotAB(dens_a, phi_aa, dens_b, phi_bb):
  # Function factory that returns functions parameterised for heterogeneous pair interactions
  def func(r):
    return 0.5 * ( (dens_b(r)/dens_a(r) * phi_aa(r)) + (dens_a(r)/dens_b(r) * phi_bb(r)) )
  return func

def makeEmbed(rho_e, rho_s, F_ni, F_i, F_e, eta):
  # Function factory returning parameterised embedding function.
  rho_n = 0.85*rho_e
  rho_0 = 1.15*rho_e

  def e1(rho):
    return sum([F_ni[i] * (rho/rho_n - 1)**float(i) for i in range(4)])

  def e2(rho):
    return sum([F_i[i] * (rho/rho_e - 1)**float(i) for i in range(4)])

  def e3(rho):
    return F_e * (1.0 - eta*math.log(rho/rho_s)) * (rho/rho_s)**eta

  def func(rho):
    if rho < rho_n:
      return e1(rho)
    elif rho_n <= rho < rho_0:
      return e2(rho)
    return e3(rho)
  return func



def makePotentialSimple(E1):
    dens_E1   = makeFunc(E1['f_e'], E1['omega_'], E1['r_e'], E1['lambda_'])

    F_ni_E1 = [E1['F_n0_'], E1['F_n1_'], E1['F_n2_'], E1['F_n3_']]
    F_i_E1 = [E1['F_0_'], E1['F_1_'], E1['F_2_'], E1['F_3_']]
    embed_E1 = makeEmbed(E1['rho_e_'], E1['rho_s_'], F_ni_E1, F_i_E1, E1['F_e_'], E1['eta_'])
    
    # Wrap them in EAMPotential objects
    eamPotentials = [
        EAMPotential(E1['symbol'], E1['number'], E1['mass'], embed_E1, dens_E1)
    ]
    
    # Define pair functions
    pair_E1E1 = makePairPotAA(E1['A_'], E1['gamma_'], E1['r_e'], E1['kappa_'], E1['B_'], E1['omega_'], E1['lambda_'])
    
       # Wrap them in Potential objects
    pairPotentials = [
      Potential(E1['symbol'], E1['symbol'], pair_E1E1)
    ]
    
    return eamPotentials, pairPotentials


def makePotentialObjects(E1, E2):
    print 'Potenciales para: ', E1['symbol'], E2['symbol']  
    # Define the density functions
    dens_E1   = makeFunc(E1['f_e'], E1['omega_'], E1['r_e'], E1['lambda_'])
    dens_E2   = makeFunc(E2['f_e'], E2['omega_'], E2['r_e'], E2['lambda_'])
    
    # Finally, define embedding functions for each species
    F_ni_E1 = [E1['F_n0_'], E1['F_n1_'], E1['F_n2_'], E1['F_n3_']]
    F_i_E1 = [E1['F_0_'], E1['F_1_'], E1['F_2_'], E1['F_3_']]
    embed_E1 = makeEmbed(E1['rho_e_'], E1['rho_s_'], F_ni_E1, F_i_E1, E1['F_e_'], E1['eta_'])
    
    F_ni_E2 = [E2['F_n0_'], E2['F_n1_'], E2['F_n2_'], E2['F_n3_']]
    F_i_E2 = [E2['F_0_'], E2['F_1_'], E2['F_2_'], E2['F_3_']]
    embed_E2 = makeEmbed(E2['rho_e_'], E2['rho_s_'], F_ni_E2, F_i_E2, E2['F_e_'], E2['eta_'])
    
    # Wrap them in EAMPotential objects
    eamPotentials = [
        EAMPotential(E1['symbol'], E1['number'], E1['mass'], embed_E1, dens_E1),
        EAMPotential(E2['symbol'], E2['number'], E2['mass'], embed_E2, dens_E2)
    ]
    
    
    # Define pair functions
    pair_E1E1 = makePairPotAA(E1['A_'], E1['gamma_'], E1['r_e'], E1['kappa_'], E1['B_'], E1['omega_'], E1['lambda_'])
    pair_E2E2 = makePairPotAA(E2['A_'], E2['gamma_'], E2['r_e'], E2['kappa_'], E2['B_'], E2['omega_'], E2['lambda_'])
    pair_E1E2 = makePairPotAB(dens_E2, pair_E2E2, dens_E1, pair_E1E1)
    
    
    # Wrap them in Potential objects
    pairPotentials = [
      Potential(E1['symbol'], E1['symbol'], pair_E1E1),
      Potential(E2['symbol'], E2['symbol'], pair_E2E2),
      Potential(E1['symbol'], E2['symbol'], pair_E1E2)
    ]

    return eamPotentials, pairPotentials
    
   
def test(): 
  E1 = 'Al'
  E2 = 'Nb'

  eamPotentials, pairPotentials = makePotentialObjects (parameters[E1], parameters[E2])

  E3 = 'Zr'
  eamPotentials3, pairPotentials3 = makePotentialSimple(parameters[E3])


  # Perform tabulation
  # Make tabulation
  nrho = 5000
  drho = 0.02

  nr = 5000
  dr = 0.0015

  with open("Zhou_"+E1+E2+".setfl", 'wb') as outfile:
    writeSetFL(
      nrho, drho,
      nr, dr,
      eamPotentials,
      pairPotentials,
      out= outfile,
      comments = ['Zhou '+ E1 + ' ' + E2, "", ""]) # <-- Note: title lines given as list of three strings
    

  with open("Zhou_"+E3+".setfl", 'wb') as outfile:
    writeSetFL(
      nrho, drho,
      nr, dr,
      eamPotentials3,
      pairPotentials3,
      out= outfile,
      comments = ['Zhou '+ E3, "", ""]) # <-- Note: title lines given as list of three strings
    
  pars = ['Mg','Al','Ti','Cr','Fe','Co','Ni','Cu','Zr','Nb','Mo', 'Pd', 'Ag','Ta','W', 'Pt','Au', 'Pb']
  for p in pars:
    print parameters[p]['r_e'], parameters[p]['gamma_'], parameters[p]['omega_'], parameters[p]['number'], p

def makePotentials(ES):
    eamPotentials = []
    densities = []
    pairsE1E1 = []
    pairsE1E2 = []
    Nelements =len(ES)

    for e in ES:
        E1 = parameters[e]
        dens_E1   = makeFunc(E1['f_e'], E1['omega_'], E1['r_e'], E1['lambda_'])
        densities.append(dens_E1)

        F_ni_E1 = [E1['F_n0_'], E1['F_n1_'], E1['F_n2_'], E1['F_n3_']]
        F_i_E1 = [E1['F_0_'], E1['F_1_'], E1['F_2_'], E1['F_3_']]
        embed_E1 = makeEmbed(E1['rho_e_'], E1['rho_s_'], F_ni_E1, F_i_E1, \
                E1['F_e_'], E1['eta_'])
 
        eamPotentials.append(EAMPotential(E1['symbol'], E1['number'], \
                E1['mass'], embed_E1, dens_E1))

        pair_E1E1 = makePairPotAA(E1['A_'], E1['gamma_'], E1['r_e'], E1['kappa_'], E1['B_'], E1['omega_'], E1['lambda_'])

        pairsE1E1.append(pair_E1E1)
  
    for i in range(Nelements):
        for j in range (i+1, Nelements):
           pair_E1E2 = makePairPotAB(densities[i], pairsE1E1[i], densities[j], pairsE1E1[j])
           pairsE1E2.append(pair_E1E2)


    # Wrap them in Potential objects
    #pairPotentials = [
    #  Potential(E1['symbol'], E1['symbol'], pair_E1E1),
    #  Potential(E2['symbol'], E2['symbol'], pair_E2E2),
    #  Potential(E1['symbol'], E2['symbol'], pair_E1E2)
    #]

    return eamPotentials
 

def test2():
    ES = ['Zr', 'Nb','Al', 'Ti']
    makePotentials(ES)

if __name__ == '__main__':
    test2()

    print parameters['Zr']

