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


format_ = ("%s %f %d %f %f %f %f %f %f\n")
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

        fileName = "./post/data_" + str(tref) + "_" + ens + ".dat" 
        f = open(fileName, 'w')
        f.write( "ensem tref time press temp toteng, volume, poteng\n")
        for i,t in enumerate(Time):
            str_ = format_ %(ens, tref, t, Press[i], Temp[i], Toteng[i], Volume[i], Poteng[i], Density[i])
            f.write (str_) 
        f.close()
        Time  = []
        Press = []
        Temp = []
        Toteng = []
        Volume = []
        Poteng = []
        Density = []


