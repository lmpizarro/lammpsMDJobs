#! /usr/bin/env python
import math
from atsim.potentials import EAMPotential
from atsim.potentials import Potential
from atsim.potentials import potentialforms
from atsim.potentials import writeSetFL

def embed(rho):
  #return -math.sqrt(rho)
  return 0

def density(rij):
  #if rij == 0:
  #  return 0.0
  #return (2.928323832 / rij) ** 6.0
  return 0


def pair_AgAg(rij):
    if rij == 0:
      return 0.0
    return (2.485883762/rij) ** 12


def buck(A, rho, C):
  """Returns a Buckingham potential function for a given set of A, rho and C parameters

  .. math ::

    U(r_{ij}) = A \exp \left( \\frac{-r_{ij}}{\\rho} \\right) - \\frac{C}{r_{ij}^6}

  :param A: Buckingham A parameter
  :param rho: Buckingham rho parameter :math:`\\rho`
  :param C: Buckingham C parameter
  :return: Function that will evaulate Buckingham potential for given A, rho and C"""

  def potential(r):
    return A * math.exp(-r/rho) #- C/r**6
  return potential

def main():
  # Create EAMPotential
  f_UU = buck(294.640000, 0.327022, 0.0)
  eamPotentials = [ EAMPotential("U", 92, 238.03, embed, density) ]
  pairPotentials = [ Potential('U', 'U', f_UU) ]

  nrho = 50000
  drho = 0.001

  nr = 12000
  dr = 0.001

  from atsim.potentials import writeFuncFL

  with open("U.eam", 'wb') as outfile:
    writeSetFL(
            nrho, drho,
            nr, dr,
            eamPotentials,
            pairPotentials,
            out= outfile,
            comments=['Buckingham Uranium as U.eam','',''])

if __name__ == "__main__":
  main()
