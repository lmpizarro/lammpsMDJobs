import math

dbMishin={'Fe':{'rc':0.5055, 'hc':0.6202, 'V0':-35674, \
        'r1': 0.1769, 'b1':5.4824e-2, 'b2':7.7124, \
        'delta': 3.6665e4, 'm': 7.0735e2, 'y': 20.201, 'gamma': 12.981,\
        'C0': 0.1391, 'B0':1.1171e5,'r0': -0.4502, 'beta':0.0, 'd1':1.9315, \
        'd2': -10.796, 'd3': -0.8928, 'q1':-5.8954, 'q2': -13.872, 'q3':247.9},

        'Ni':{'rc':0.5168, 'hc':0.3323, 'V0':-3516, \
        'r1': 3.8673e-5, 'b1':4.7067e-3, 'b2':0.1511, \
        'delta': 3.6046e3, 'm': 0.0, 'y': 19.251, 'gamma': 16.802,\
        'C0': 0.2033, 'B0':1.1914e5,'r0': -0.3138, \
        'beta':0.4890e-2, 'd1':4.4657e-2, \
        'd2': -1.3702e1, 'd3': -0.9611, 'q1':-6.4502e2, 'q2': 0.628, 'q3':602.08}
        }

def makeDens(el):
    def func (r):
        r /=10
        z = r - el['r0']
        x = (r - el['rc']) / el['hc']

        X_X = 0
        if x < 0:
            X_X = math.pow(x, 4) / (1 + math.pow(x,4))


        f = X_X * (math.pow(z, el['y'])*math.exp(-el['gamma']*z)*\
                (1+ el['B0']*math.exp(-el['gamma']*z)) + el['C0'])
        return f
    return func

def makePot(el):
    def func(r):
        r /=10
        z = r / el['r1']
        x = (r - el['rc']) / el['hc']
        B1 = el['b1'] / math.pow(z, el['b2'])
        B2 = el['b2'] / math.pow(z, el['b1'])
        BB = B2 - B1 

        X_X = 0
        if x < 0:
            X_X = math.pow(x, 4) / (1 + math.pow(x,4))

        t = el['V0'] / (el['b1'] + el['b2'])

        f = X_X * (t * BB + el['delta']) 
        return f
    return func

import numpy as np
import matplotlib.pyplot as plt

def test_01():


    r = np.linspace(0.0,5.6, 100)

    el = dbMishin['Fe']

    for e in el:
        print e, el[e]
    dens = makeDens(el)

    densvals = []
    for e in r:
        densvals.append( dens(e))

    
    plt.plot(r, densvals)
    plt.show()


def test_02():
    r = np.linspace(.01,5.6, 100)

    el = dbMishin['Fe']
    dens =  makeDens(el)
    pot = makePot(el) 

    densvals = []
    for e in r:
        #print pot(e)
        #densvals.append( el['m'] * dens(e) + pot(e))
        densvals.append( pot(e))

    
    plt.plot(r, densvals)
    plt.show()

import dbEamZhou as dbA

def makeRose(el, beta):
    Ec = el['Ec']
    a0 = el['a']
    B= el['B']
    V0 = math.pow(a0,3)

    alfa = math.sqrt(V0 * B / Ec / 160.217662)

    def func(a):
        x = -1 + a / a0
        f = -Ec * (1 + alfa * x) * math.exp(-alfa * x)
        return f
    return func

def test_rose_01():
    Al = dbA.exp_props['Fe']


    f = makeRose(Al, 1e-3)

    a0 = Al['a']
    a = a0 + np.linspace(-1, 1, 100) * a0 / 10
    
    valsF = []
    for e in a:
        valsF.append( f(e))


    plt.plot(a, valsF)
    plt.show()

def test_rose_02():
    elements = ['Al', 'Fe', 'Cr']

    for el in elements:
        print el
        Al = dbA.exp_props[el]


        f = makeRose(Al, 1e-6)

        a0 = Al['a']
        a = a0 + np.linspace(-1, 1, 100) * a0 / 10
    
        valsF = []
        for e in a:
            valsF.append( f(e))
        plt.plot(a, valsF)


    plt.show()




if __name__ == '__main__':
    test_01()
    test_02()

