# http://periodictable.readthedocs.io/en/latest/
from periodictable import elements

'''
https://sites.google.com/site/eampotentials/Home/MgY#TOC-Latest:-PdSi.lammps.eam-10-10-2011-
'''

fileName = 'MgY.lammps.eam'
potentials = {'properties':{}}



def readFile():
    with open(fileName) as f:
        content = f.readlines()
    content = [x.strip() for x in content]

    return content

def getRhoF(content, init, nrho, nr):

    e = content[init].split()
    aNumber = int(e[0])
    symbol = str(elements[aNumber])
    potentials[symbol]['aMass'] = float(e[1])
    potentials[symbol]['a0'] = float(e[2])
    potentials[symbol]['lat'] = e[3]

    F1 = []
    NRHO = 0
    RHO1 = []
    NR = 0
    cLines = 0
    for  c in content[init+1:]:
         cLines +=1
         dataL = c.split()
         for d in dataL:
             F1.append(float(d))
             NRHO +=1
         if NRHO == nrho:
             print "F1 ready"
             break

    print len(F1), cLines, F1[len(F1)-1]

    print content[cLines+init + 2]

    r = cLines
    for i, c in enumerate(content[r+init + 1:]):
         cLines +=1
         dataL = c.split()
         for d in dataL:
             RHO1.append(float(d))
             NR +=1
         if NR == nr:
             print "RHO1 ready"
             break

    print len(RHO1), cLines, RHO1[len(RHO1)-1]

    return F1, RHO1, cLines

def main():
    content = readFile()

    e = content[3].split()
    potentials['properties']['nEls'] = int(e[0])
    for i in range(potentials['properties']['nEls']):
        potentials[e[1+i]] ={}

    e = content[4].split()
    potentials['properties']['nrho'] =  int(e[0])
    potentials['properties']['drho'] =  float(e[1])
    potentials['properties']['nr'] =  int(e[2])
    potentials['properties']['dr'] =  float(e[3])
    potentials['properties']['cutoff'] = float(e[4])

    nrho = potentials['properties']['nrho']
    nr = potentials['properties']['nr']

    F1, RHO1, init = getRhoF(content, 5, nrho, nr)

    print potentials
    print init 

    import matplotlib.pyplot as plt

    #plt.plot(F1)
    #plt.show()
    #plt.plot(RHO1)
    #plt.show()

    F2, RHO2, init = getRhoF(content, init + 6, nrho, nr)
    plt.plot(F2)
    plt.show()
    plt.plot(RHO2)
    plt.show()


if __name__ == "__main__":
  main()
