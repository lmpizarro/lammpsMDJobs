from periodic.table import element
from atsim.potentials import writeSetFL
from atsim.potentials import Potential
from atsim.potentials import EAMPotential
import  dbEamZhou as dbZhou 
#import parameters

import math

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

class calcPotentials():
    def __init__(self, ES):
        self.ES = ES
        self.eamPotentials = []
        self.pairsE1E2 = []
        self.Nelements =len(ES)

        self.myParameters={}

        for e in ES:
            self.myParameters[e] = dbZhou.parameters[e]

        self.symbols = self.myParameters.keys()
        self.k_p = self.myParameters[self.symbols[0]].keys()

        for e in ES:
            E1 =  self.myParameters[e]
            dens_E1 = makeFunc(E1['f_e'], E1['omega_'], E1['r_e'], E1['lambda_'])
            self.myParameters[e]['dens'] = dens_E1

            F_ni_E1 = [E1['F_n0_'], E1['F_n1_'], E1['F_n2_'], E1['F_n3_']]
            F_i_E1 = [E1['F_0_'], E1['F_1_'], E1['F_2_'], E1['F_3_']]
            embed_E1 = makeEmbed(E1['rho_e_'], E1['rho_s_'], F_ni_E1, F_i_E1, \
                E1['F_e_'], E1['eta_'])
            self.myParameters[e]['embed'] = embed_E1 
 
            self.eamPotentials.append(EAMPotential(E1['symbol'], E1['number'], \
                E1['mass'], embed_E1, dens_E1))

            pair_E1E1 = makePairPotAA(E1['A_'], E1['gamma_'], E1['r_e'], E1['kappa_'], E1['B_'], E1['omega_'], E1['lambda_'])

            self.myParameters[e]['pair']=pair_E1E1

        #print self.myParameters

    # pyformat.info
    def __str__(self):
        str_ = ' '*18
        for s in self.symbols:
            str_ += '{:14}'.format(s)
        str_ +='\n' 

        for c in self.k_p:
            if c != 'name' and c != 'symbol':
                str_ +='{:8}'.format(c)
                for s in self.symbols:
                    str_ +='{:14f}'.format(float(self.myParameters[s][c]))
                str_ +='\n'

        return str_ 

 
if __name__ == '__main__':
    ES = ['Zr', 'Nb','Al', 'Ti']

    c = calcPotentials(ES)

    print c


