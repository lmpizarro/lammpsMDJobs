import sys

sys.path.append('/home/lmpizarro/ENVS/pizza/src')

import numpy as np
from log import log


folder ="./logs"
ensembles = ["NPT", "NVE", "NVT"]

Time  = []
Press = []
Temp = []
Toteng = []
Volume = []
Poteng = []
Density = []


print "ensem tref time press temp toteng, volume, poteng"
for tref in range(250,350,25):
    for ens in ensembles:
        LOG = folder + "/log_" + ens + "_"+ str(tref) + ".lammps"
        lg = log(LOG)
        time,temp,press, toteng, volume, poteng, density= lg.get("Step","Temp","Press", \
                "TotEng", "Volume", "PotEng", "Density")

        [Time.append(e) for e in time]
        [Press.append(e) for e in press]
        [Temp.append(e) for e in temp] 
        [Toteng.append(e) for e in toteng]
        [Volume.append(e) for e in volume]
        [Poteng.append(e) for e in poteng]
        [Density.append(e) for e in density]


        for i,t in enumerate(Time):
            print ens, tref, t, Press[i], Temp[i], Toteng[i], Volume[i], Poteng[i], Density[i]
