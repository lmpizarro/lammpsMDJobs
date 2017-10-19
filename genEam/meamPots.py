# meam data from vax files fcc,bcc,dia
potKeys = [['elt', 'lat', 'z', 'ielement', 'atwt'],
 ['alpha', 'b0', 'b1', 'b2', 'b3', 'alat', 'esub', 'asub'],
['t0', 't1', 't2', 't3', 'rozero', 'ibar']]


elements={'U':{'elt':'U','lat':'fcc', 'z':12, 'ielement':92, 'atwt':238.03,\
'alpha':5.1, 'b0':4.80, 'b1':6.0, 'b2':6, 'b3':6, 'alat':4.280, 'esub':5.27,\
'asub':0.98, 't0':1,  't1':2.50, 't2':4, 't3':1.0, 'rozero':1, 'ibar':0},
'Zr':{'elt':'Zr', 'lat':'bcc', 'z':8, 'ielement':40, 'atwt':91.224,\
'alpha':4.10, 'b0':2.80, 'b1':2.0, 'b2':7.00, 'b3':1.0, 'alat':3.535,\
'esub':6.20, 'asub':0.48, 't0':1.0, 't1':3.00, 't2':2.0, 't3':-7.0,\
'rozero':1, 'ibar':0},
'Fe':{'elt':'Fe', 'lat':'bcc', 'z':8, 'ielement':26, 'atwt':55.845,\
'alpha':5.07, 'b0':2.94, 'b1':1.0, 'b2':1.00, 'b3':1.0, 'alat':2.8665,\
'esub':4.29, 'asub':0.89, 't0':1.0, 't1':3.94, 't2':4.12, 't3':-1.5,\
'rozero':1, 'ibar':0},
'Y':{'elt':'Y', 'lat':'hcp', 'z':12, 'ielement':39, 'atwt':88.90585,\
'alpha':5.07, 'b0':4.74, 'b1':1.0, 'b2':2.50, 'b3':1.0, 'alat':2.8665,\
'esub':5.30, 'asub':0.9, 't0':1.0, 't1':3.30, 't2':3.2, 't3':-2.0,\
'rozero':1, 'ibar':0},


}
'''
for e in elements:
    line1 =''
    line2 =''
    for ls in potKeys:
        for p in ls:
            line2 += '  ' + str(elements[e][p])
            line1 += '  ' + str(p)
        line2 +='\n'
        line1 +='\n'
    print line1[:-1]
    print line2[:-1]
'''
def headerMeam():
    a =  '{:8} {:8} {:8} {:8} {:8}\n'.format(potKeys[0][0], potKeys[0][1],\
                potKeys[0][2], potKeys[0][3], potKeys[0][4])
    a +=  '{:8} {:8} {:8} {:8} {:8} {:8} {:8} {:8}\n'.format(potKeys[1][0],\
                potKeys[1][1], potKeys[1][2], potKeys[1][3], potKeys[1][4],\
                potKeys[1][5], potKeys[1][6], potKeys[1][7])
    a +=  '{:8} {:8}          {:8} {:8}          {:8} {:8}\n'.format(potKeys[2][0], potKeys[2][1],\
                potKeys[2][2], potKeys[2][3], potKeys[2][4],\
                potKeys[2][5])
    return a

def paramV(element):
    d = elements[element]
    keys = potKeys
    a = '{:8} {:} {:8} {:8} {:8}\n'.format(d[keys[0][0]], d[keys[0][1]],\
                d[keys[0][2]], d[keys[0][3]], d[keys[0][4]])
    a += '{:} {:8} {:8} {:8} {:8} {:8} {:8} {:8}\n'.format(d[keys[1][0]],\
                d[keys[1][1]], d[keys[1][2]], d[keys[1][3]], d[keys[1][4]],\
                d[keys[1][5]], d[keys[1][6]], d[keys[1][7]])
    a += '{:} {:8}          {:8} {:8}          {:8} {:8}\n'.\
                format(d[keys[2][0]], d[keys[2][1]],\
                d[keys[2][2]], d[keys[2][3]], d[keys[2][4]],\
                d[keys[2][5]])
    return a


def genMeamf(els):
    comment = ''
    for e in els:
        comment += e + ' '
    comment +='\n'

    pzr = '# meam data for: ' + comment
    for e in els:
        if e in elements.keys():
            pzr +=  headerMeam()
            pzr +=  paramV(e)
    return pzr


def main ():
    meamf = genMeamf(['U', 'Zr', 'Fe'])
    print meamf

if __name__ == "__main__":
  main()
