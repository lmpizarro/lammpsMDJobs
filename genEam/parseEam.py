# http://periodictable.readthedocs.io/en/latest/
from periodictable import elements

'''
https://sites.google.com/site/eampotentials/Home/MgY#TOC-Latest:-PdSi.lammps.eam-10-10-2011-
'''
fileName = 'MgY.lammps.eam'
potentials = {'properties':{}}
symbols =[]

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

    potentials[symbol]['V'] = V

    return cLines

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

    potentials[symbol]['F'] = F1

    r = cLines
    for i, c in enumerate(content[r + 1:]):
         cLines +=1
         dataL = c.split()
         for d in dataL:
             RHO1.append(float(d))
             NR +=1
         if NR == nr:
             break

    potentials[symbol]['RHO'] = RHO1

    return cLines

def main():
    content = readFile()

    cLine = 3
    e = content[cLine].split()
    potentials['properties']['nEls'] = int(e[0])
    s12 = ""
    for i in range(potentials['properties']['nEls']):
        s12 +=e[1+i]
        potentials[e[1+i]] ={}
    potentials[s12] = {}



    cLine +=1
    e = content[cLine].split()
    potentials['properties']['nrho'] =  int(e[0])# n points density
    potentials['properties']['drho'] =  float(e[1])# distance points density
    potentials['properties']['nr'] =  int(e[2])# n points embed & Vpair
    potentials['properties']['dr'] =  float(e[3])# distance points embed & Vpair
    potentials['properties']['cutoff'] = float(e[4])

    nrho = potentials['properties']['nrho']
    nr = potentials['properties']['nr']

    cLine +=1
    symbol = getElement(content, cLine)
    cLine = getRhoF(symbol, content, cLine, nrho, nr)
    symbols.append(symbol)

    cLine +=1
    symbol = getElement(content, cLine)
    cLine = getRhoF(symbol, content, cLine, nrho, nr)
    symbols.append(symbol)

    #cLine +=1
    #def getRhoV(symbol, content, cLines, nr):
    cLine = getRhoV(symbols[0], content, cLine, nr)
    symbol = symbols[0]+symbols[1]

    symbols.append(symbol)
    cLine = getRhoV(symbol, content, cLine, nr)
    cLine = getRhoV(symbols[1], content, cLine, nr)

    for e in symbols:
        for k in  potentials[e]:
            print e, k
        print ""

    rhoY =  potentials['Y']['RHO']
    fY =  potentials['Y']['F']
    vY =  potentials['Y']['V']

    nrho = potentials['properties']['nrho']# n points density
    drho = potentials['properties']['drho']# distance points density
    nr = potentials['properties']['nr']# n points embed & Vpair
    dr = potentials['properties']['dr']# distance points embed & Vpair
    cutoff = potentials['properties']['cutoff']

    print 'nrho*drho: ', nrho * drho
    print 'nr*dr:     ', nr * dr
    print 'cutoff:    ', cutoff

    import numpy as np
    xrho = np.linspace(0, nrho*drho, num=nrho, endpoint=True)
    xr = np.linspace(0, nr*dr, num=nr, endpoint=True)

    tr = 250
    import matplotlib.pyplot as plt
    plt.plot(xr[tr:], vY[tr:])
    plt.show()

    plt.plot(xr[tr:], rhoY[tr:])
    plt.show()

    plt.plot(xrho, fY)
    plt.show()

if __name__ == "__main__":
  main()
