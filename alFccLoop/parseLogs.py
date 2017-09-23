import sys

sys.path.append('/pathto/pizza/src')

import numpy as np
from log import log


folder ="./logs"
Time  = []
Press = []
Temp = []
Toteng = []
print "tempCalc press temp toteng"
EquilTime =  53
for t in range(50,700,25):
    LOG = folder + "/log_" + str(t) + ".lammps"
    #print LOG
    lg = log(LOG)
    time,temp,press, toteng, volume, poteng, density= lg.get("Step","Temp","Press", \
            "TotEng", "Volume", "PotEng", "Density")
    time = time[EquilTime:]
    [Time.append(e) for e in time[EquilTime:]]
    [Press.append(e) for e in press[EquilTime:]]
    [Temp.append(e) for e in temp[EquilTime:]] 
    [Toteng.append(e) for e in toteng[EquilTime:]]
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



    print  t, Mpress, Vpress, Mtemp, Vtemp, Mtoteng, Vtoteng, Mvolume, Vvolume, Mpoteng, Vpoteng
