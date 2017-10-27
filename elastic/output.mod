
print "cxx={'C11':${C11all},'C12':${C12all},\" file eConstants.py
print "'C44': ${C44all},'B':${bulkmodulus},\" append eConstants.py
print "'Sh1':${shearmodulus1},'Sh2':${shearmodulus2},\" append eConstants.py
print "'Pr':${poissonratio}}" append eConstants.py

print "C11 C12 C44 B Sh1 Sh2 Pr" file info.dat

print "${C11all} ${C12all} ${C44all} ${bulkmodulus} ${shearmodulus1} ${shearmodulus2} ${poissonratio}}" append info.dat

print "=========================================" append info.dat
print "Components of the Elastic Constant Tensor" append info.dat
print "========================================="
print "Elastic Constant C11all = ${C11all} ${cunits}" append info.dat
print "Elastic Constant C22all = ${C22all} ${cunits}" append info.dat
print "Elastic Constant C33all = ${C33all} ${cunits}" append info.dat

print "Elastic Constant C12all = ${C12all} ${cunits}" append info.dat
print "Elastic Constant C13all = ${C13all} ${cunits}" append info.dat
print "Elastic Constant C23all = ${C23all} ${cunits}" append info.dat

print "Elastic Constant C44all = ${C44all} ${cunits}" append info.dat
print "Elastic Constant C55all = ${C55all} ${cunits}" append info.dat
print "Elastic Constant C66all = ${C66all} ${cunits}" append info.dat

print "Elastic Constant C14all = ${C14all} ${cunits}" append info.dat
print "Elastic Constant C15all = ${C15all} ${cunits}" append info.dat
print "Elastic Constant C16all = ${C16all} ${cunits}" append info.dat

print "Elastic Constant C24all = ${C24all} ${cunits}" append info.dat
print "Elastic Constant C25all = ${C25all} ${cunits}" append info.dat
print "Elastic Constant C26all = ${C26all} ${cunits}" append info.dat

print "Elastic Constant C34all = ${C34all} ${cunits}" append info.dat
print "Elastic Constant C35all = ${C35all} ${cunits}" append info.dat
print "Elastic Constant C36all = ${C36all} ${cunits}" append info.dat

print "Elastic Constant C45all = ${C45all} ${cunits}" append info.dat
print "Elastic Constant C46all = ${C46all} ${cunits}" append info.dat
print "Elastic Constant C56all = ${C56all} ${cunits}" append info.dat

print "=========================================" append info.dat
print "Average properties for a cubic crystal" append info.dat
print "=========================================" append info.dat

print "Bulk Modulus = ${bulkmodulus} ${cunits}" append info.dat
print "Shear Modulus 1 = ${shearmodulus1} ${cunits}" append info.dat
print "Shear Modulus 2 = ${shearmodulus2} ${cunits}" append info.dat
print "Poisson Ratio = ${poissonratio}" append info.dat

shell rm restart.equil
shell rm log.lammps
