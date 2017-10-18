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


def getElement(content, cLines):
    e = content[cLines].split()
    aNumber = int(e[0])
    symbol = str(elements[aNumber])
    potentials[symbol]['aMass'] = float(e[1])
    potentials[symbol]['a0'] = float(e[2])
    potentials[symbol]['lat'] = e[3]

    return symbol

def getRhoV(symbol, content, cLines, nr):
    V = []
    NR = 0

    for  c in content[cLines+1:]:
         cLines +=1
         dataL = c.split()
         for d in dataL:
             V.append(float(d))
             NR +=1
         if NR == nr:
             break

    return V, cLines


def getRhoF(symbol, content, cLines, nrho, nr):
    F1 = []
    NRHO = 0
    RHO1 = []
    NR = 0

    for  c in content[cLines+1:]:
         cLines +=1
         dataL = c.split()
         for d in dataL:
             F1.append(float(d))
             NRHO +=1
         if NRHO == nrho:
             break

    potentials[symbol]['F1'] = F1


    r = cLines
    for i, c in enumerate(content[r + 1:]):
         cLines +=1
         dataL = c.split()
         for d in dataL:
             RHO1.append(float(d))
             NR +=1
         if NR == nr:
             break

    potentials[symbol]['F1'] = RHO1

    return cLines

def main():
    content = readFile()

    cLine = 3
    e = content[cLine].split()
    potentials['properties']['nEls'] = int(e[0])
    for i in range(potentials['properties']['nEls']):
        potentials[e[1+i]] ={}

    cLine +=1
    e = content[cLine].split()
    potentials['properties']['nrho'] =  int(e[0])
    potentials['properties']['drho'] =  float(e[1])
    potentials['properties']['nr'] =  int(e[2])
    potentials['properties']['dr'] =  float(e[3])
    potentials['properties']['cutoff'] = float(e[4])

    nrho = potentials['properties']['nrho']
    nr = potentials['properties']['nr']

    cLine +=1
    symbol = getElement(content, cLine)
    cLine = getRhoF(symbol, content, cLine, nrho, nr)
    print symbol

    cLine +=1
    symbol = getElement(content, cLine)
    cLine = getRhoF(symbol, content, cLine, nrho, nr)
    print symbol

    print cLine
    #cLine +=1
    #def getRhoV(symbol, content, cLines, nr):
    V1, cLine = getRhoV(symbol, content, cLine, nr)
    V2, cLine = getRhoV(symbol, content, cLine, nr)
    V3, cLine = getRhoV(symbol, content, cLine, nr)

    import matplotlib.pyplot as plt

    #plt.plot(V1)
    #plt.show()
    #plt.plot(V2)
    #plt.show()
    #plt.plot(V3)
    #plt.show()


if __name__ == "__main__":
  main()
