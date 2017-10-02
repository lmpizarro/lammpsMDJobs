import sys

sys.path.append('/home/lmpizarro/ENVS/pizza/src')

import numpy as np
from log import log


folder ="./"
Time  = []
Press = []
Temp = []
Toteng = []
Enthalpy= []
#print "tempCalc press temp toteng, volume, poteng"
EquilTime =  110

def analisis(LOG):
    lg = log(LOG)
    time,temp,press, toteng, volume, poteng, density, enthalpy= lg.get("Step","Temp","Press", \
            "TotEng", "Volume", "PotEng", "Density", 'Enthalpy')
    time = time[EquilTime:]
    [Time.append(e) for e in time[EquilTime:]]
    [Press.append(e) for e in press[EquilTime:]]
    [Temp.append(e) for e in temp[EquilTime:]] 
    [Toteng.append(e) for e in toteng[EquilTime:]]
    [Enthalpy.append(e) for e in enthalpy[EquilTime:]]
    Mpress = np.array(press[EquilTime:]).mean()
    Vpress = np.array(press[EquilTime:]).var()
    Mtemp = np.array(temp[EquilTime:]).mean()
    Vtemp = np.array(temp[EquilTime:]).var()
    Mtoteng = np.array(toteng[EquilTime:]).mean()
    Vtoteng = np.array(toteng[EquilTime:]).var()
    Mvolume = np.array(volume[EquilTime:]).mean()
    Vvolume = np.array(volume[EquilTime:]).var()
    Mpoteng = np.array(poteng[EquilTime:]).mean()
    Vpoteng = np.array(poteng[EquilTime:]).var()
    Mdensity = np.array(density[EquilTime:]).mean()
    Vdensity = np.array(density[EquilTime:]).var()

    Menthalpy = np.array(enthalpy[EquilTime:]).mean()
    Venthalpy = np.array(enthalpy[EquilTime:]).var()
    return  time[len(time)-1], Mpress, Vpress, Mtemp, Vtemp, Mtoteng, Vtoteng, Mvolume, Vvolume, Mpoteng, Vpoteng, Menthalpy, Venthalpy  

print ("%s %.4e %.4e %.4e %.4e %.4e %.4f %.4f %.4e %.4e %.4e %.4f %.4f")%analisis("log.lammps")
