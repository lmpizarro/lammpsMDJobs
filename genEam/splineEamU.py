import eamU as eu

from atsim.potentials import EAMPotential
from atsim.potentials import Potential
from atsim.potentials import potentialforms
from atsim.potentials import writeSetFL


plot = False
if plot:
    import matplotlib.pyplot as plt
    xinit = 0
    xnew = np.linspace(xinit, 8, num=41, endpoint=True)
    plt.grid(True)
    plt.plot(rij, phi_rij, 'X')
    plt.plot(xnew, phi(xnew), '--')
    xnew2 = np.linspace(xinit, 8, num=41, endpoint=True)
    plt.legend(['data', 'linear', 'cubic'], loc='best')
    plt.show()

    plt.plot(rij, rho_rij, 'X')
    plt.plot(xnew, rho(xnew), '--')
    plt.grid(True)
    plt.show()

    xnew = np.linspace(0, 1.353, num=41, endpoint=True)
    plt.plot(rhoj, F_rhoj, 'X')
    plt.plot(xnew, F(xnew), '--')
    plt.grid(True)
    plt.show()


def embed(rho):
  if rho > 1.353:
      return 0.0
  return eu.F(rho)

def density(rij):
  if rij >8.0:
    return 0.0
  return eu.rho(rij) 


def pair_UU(rij):
  if rij >8.0:
    return 0.0
  return eu.phi(rij)

def main():
  # Create EAMPotential
  eamPotentials = [ EAMPotential("U", 92, 238.03, embed, density) ]
  pairPotentials = [ Potential('U', 'U', pair_UU) ]

  cutoff = 8.0
  #nrho = 962 # n points density
  nrho = 5000 # n points density
  #drho = 0.00103950 # distance points density
  drho = 0.02 # distance points density

  #nr = 1462 # n points embed & Vpair
  nr = 5000 # n points embed & Vpair
  dr = cutoff / nr # distance points embed & Vpair

  from atsim.potentials import writeFuncFL

  with open("U2.eam", 'wb') as outfile:
    writeSetFL(
            nrho, drho,
            nr, dr,
            eamPotentials,
            pairPotentials,
            out= outfile,
            comments=['Spline Uranium as U.eam','',''])

if __name__ == "__main__":
  main()

