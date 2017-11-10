import sys

sys.path.append('../zhouEam')

import system as sysAt
global potential_mod

class Elastic():
    def __init__(self, lammps_settings):
        self.lammps_settings = lammps_settings
        self.in_lmp = self.lammps_settings['in_lmp']

        if self.lammps_settings['minimize'] == 'elastic':
            print 'calc elastic'

        self.potential_mod ='''
# NOTE: This script can be modified for different pair styles 
# See in.elastic for more info.

# Choose potential


%s

# Setup neighbor style
neighbor 1.0 nsq
neigh_modify once no every 1 delay 0 check yes

# Setup minimization style
min_style	     cg
min_modify	     dmax ${dmax} line quadratic

# Setup output
thermo		1
thermo_style custom step temp pe press pxx pyy pzz pxy pxz pyz lx ly lz vol
thermo_modify norm no
'''

        self.displace_mod='''
# NOTE: This script should not need to be
# modified. See in.elastic for more info.
#
# Find which reference length to use

if "${dir} == 1" then &
   "variable len0 equal ${lx0}" 
if "${dir} == 2" then &
   "variable len0 equal ${ly0}" 
if "${dir} == 3" then &
   "variable len0 equal ${lz0}" 
if "${dir} == 4" then &
   "variable len0 equal ${lz0}" 
if "${dir} == 5" then &
   "variable len0 equal ${lz0}" 
if "${dir} == 6" then &
   "variable len0 equal ${ly0}" 

# Reset box and simulation parameters

clear
box tilt large
read_restart restart.equil
include potential.mod

# Negative deformation

variable delta equal -${up}*${len0}
variable deltaxy equal -${up}*xy
variable deltaxz equal -${up}*xz
variable deltayz equal -${up}*yz
if "${dir} == 1" then &
   "change_box all x delta 0 ${delta} xy delta ${deltaxy} xz delta ${deltaxz} remap units box"
if "${dir} == 2" then &
   "change_box all y delta 0 ${delta} yz delta ${deltayz} remap units box"
if "${dir} == 3" then &
   "change_box all z delta 0 ${delta} remap units box"
if "${dir} == 4" then &
   "change_box all yz delta ${delta} remap units box"
if "${dir} == 5" then &
   "change_box all xz delta ${delta} remap units box"
if "${dir} == 6" then &
   "change_box all xy delta ${delta} remap units box"

# Relax atoms positions

minimize ${etol} ${ftol} ${maxiter} ${maxeval}

# Obtain new stress tensor
 
variable tmp equal pxx
variable pxx1 equal ${tmp}
variable tmp equal pyy
variable pyy1 equal ${tmp}
variable tmp equal pzz
variable pzz1 equal ${tmp}
variable tmp equal pxy
variable pxy1 equal ${tmp}
variable tmp equal pxz
variable pxz1 equal ${tmp}
variable tmp equal pyz
variable pyz1 equal ${tmp}

# Compute elastic constant from pressure tensor

variable C1neg equal ${d1}
variable C2neg equal ${d2}
variable C3neg equal ${d3}
variable C4neg equal ${d4}
variable C5neg equal ${d5}
variable C6neg equal ${d6}

# Reset box and simulation parameters

clear
box tilt large
read_restart restart.equil
include potential.mod

# Positive deformation

variable delta equal ${up}*${len0}
variable deltaxy equal ${up}*xy
variable deltaxz equal ${up}*xz
variable deltayz equal ${up}*yz
if "${dir} == 1" then &
   "change_box all x delta 0 ${delta} xy delta ${deltaxy} xz delta ${deltaxz} remap units box"
if "${dir} == 2" then &
   "change_box all y delta 0 ${delta} yz delta ${deltayz} remap units box"
if "${dir} == 3" then &
   "change_box all z delta 0 ${delta} remap units box"
if "${dir} == 4" then &
   "change_box all yz delta ${delta} remap units box"
if "${dir} == 5" then &
   "change_box all xz delta ${delta} remap units box"
if "${dir} == 6" then &
   "change_box all xy delta ${delta} remap units box"

# Relax atoms positions

minimize ${etol} ${ftol} ${maxiter} ${maxeval}

# Obtain new stress tensor
 
variable tmp equal pe
variable e1 equal ${tmp}
variable tmp equal press
variable p1 equal ${tmp}
variable tmp equal pxx
variable pxx1 equal ${tmp}
variable tmp equal pyy
variable pyy1 equal ${tmp}
variable tmp equal pzz
variable pzz1 equal ${tmp}
variable tmp equal pxy
variable pxy1 equal ${tmp}
variable tmp equal pxz
variable pxz1 equal ${tmp}
variable tmp equal pyz
variable pyz1 equal ${tmp}

# Compute elastic constant from pressure tensor

variable C1pos equal ${d1}
variable C2pos equal ${d2}
variable C3pos equal ${d3}
variable C4pos equal ${d4}
variable C5pos equal ${d5}
variable C6pos equal ${d6}

# Combine positive and negative 

variable C1${dir} equal 0.5*(${C1neg}+${C1pos})
variable C2${dir} equal 0.5*(${C2neg}+${C2pos})
variable C3${dir} equal 0.5*(${C3neg}+${C3pos})
variable C4${dir} equal 0.5*(${C4neg}+${C4pos})
variable C5${dir} equal 0.5*(${C5neg}+${C5pos})
variable C6${dir} equal 0.5*(${C6neg}+${C6pos})

# Delete dir to make sure it is not reused

variable dir delete
'''


        self.init_mod='''
# NOTE: This script can be modified for different atomic structures, 
# units, etc. See in.elastic for more info.
#

# Define the finite deformation size. Try several values of this
# variable to verify that results do not depend on it.
variable up equal 0.5e-6
 
# Define the amount of random jiggle for atoms
# This prevents atoms from staying on saddle points
variable atomjiggle equal 1.0e-5

# Uncomment one of these blocks, depending on what units
# you are using in LAMMPS and for output

# metal units, elastic constants in eV/A^3
#units		metal
#variable cfac equal 6.2414e-7
#variable cunits string eV/A^3

# metal units, elastic constants in GPa
units		metal
variable cfac equal 1.0e-4
variable cunits string GPa

# real units, elastic constants in GPa
#units		real
#variable cfac equal 1.01325e-4
#variable cunits string GPa

# Define minimization parameters
variable etol equal 1.0e-15 
variable ftol equal 1.0e-15
variable maxiter equal 10000
variable maxeval equal 10000
variable dmax equal 1.0e-2

# generate the box and atom positions using a diamond lattice
#variable a equal 4.240
variable a equal 5.43

boundary	p p p

#lattice         diamond $a
#region		box prism 0 2.0 0 3.0 0 4.0 0.0 0.0 0.0
#create_box	1 box
#create_atoms	1 box

read_data data.lmp
#read_data atom.dat

# Need to set mass to something, just to satisfy LAMMPS
mass 1 1.0e-20
'''



        self.output_mod='''
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
'''

        self.in_elastic='''
echo both
# Compute elastic constant tensor for a crystal
#
# Written by Aidan Thompson (Sandia, athomps@sandia.gov)
#
#  This script uses the following three include files.
#
#   init.mod      (must be modified for different crystal structures)
# 	       	  Define units, deformation parameters and initial
#		  configuration of the atoms and simulation cell.  
#
#
#   potential.mod    (must be modified for different pair styles)
# 		     Define pair style and other attributes 
#		     not stored in restart file
#
#
#   displace.mod    (displace.mod should not need to be modified)
# 		    Perform positive and negative box displacements 
# 		    in direction ${dir} and size ${up}. 
# 		    It uses the resultant changes 
#		    in stress to compute one
# 		    row of the elastic stiffness tensor
#		    
#		    Inputs variables:
#		    	   dir = the Voigt deformation component 
#		    		    (1,2,3,4,5,6)  
#		    Global constants:
#       	    	   up = the deformation magnitude (strain units)
#       		   cfac = conversion from LAMMPS pressure units to 
#               	   output units for elastic constants 
#
#
#  To run this on a different system, it should only be necessary to 
#  modify the files init.mod and potential.mod. In order to calculate
#  the elastic constants correctly, care must be taken to specify
#  the correct units in init.mod (units, cfac and cunits). It is also
#  important to verify that the minimization of energy w.r.t atom
#  positions in the deformed cell is fully converged.
#  One indication of this is that the elastic constants are insensitive
#  to the choice of the variable ${up} in init.mod. Another is to check
#  the final max and two-norm forces reported in the log file. If you know
#  that minimization is not required, you can set maxiter = 0.0 in 
#  init.mod. 
#

include init.mod
include potential.mod

# Compute initial state
fix 3 all box/relax  aniso 0.0
minimize ${etol} ${ftol} ${maxiter} ${maxeval}

variable tmp equal pxx
variable pxx0 equal ${tmp}
variable tmp equal pyy
variable pyy0 equal ${tmp}
variable tmp equal pzz
variable pzz0 equal ${tmp}
variable tmp equal pyz
variable pyz0 equal ${tmp}
variable tmp equal pxz
variable pxz0 equal ${tmp}
variable tmp equal pxy
variable pxy0 equal ${tmp}

variable tmp equal lx
variable lx0 equal ${tmp}
variable tmp equal ly
variable ly0 equal ${tmp}
variable tmp equal lz
variable lz0 equal ${tmp}

# These formulas define the derivatives w.r.t. strain components
# Constants uses $, variables use v_ 
variable d1 equal -(v_pxx1-${pxx0})/(v_delta/v_len0)*${cfac}
variable d2 equal -(v_pyy1-${pyy0})/(v_delta/v_len0)*${cfac}
variable d3 equal -(v_pzz1-${pzz0})/(v_delta/v_len0)*${cfac}
variable d4 equal -(v_pyz1-${pyz0})/(v_delta/v_len0)*${cfac}
variable d5 equal -(v_pxz1-${pxz0})/(v_delta/v_len0)*${cfac}
variable d6 equal -(v_pxy1-${pxy0})/(v_delta/v_len0)*${cfac}

displace_atoms all random ${atomjiggle} ${atomjiggle} ${atomjiggle} 87287 units box

# Write restart
unfix 3
write_restart restart.equil

# uxx Perturbation

variable dir equal 1
include displace.mod

# uyy Perturbation

variable dir equal 2
include displace.mod

# uzz Perturbation

variable dir equal 3
include displace.mod

# uyz Perturbation

variable dir equal 4
include displace.mod

# uxz Perturbation

variable dir equal 5
include displace.mod

# uxy Perturbation

variable dir equal 6
include displace.mod

# Output final values

variable C11all equal ${C11}
variable C22all equal ${C22}
variable C33all equal ${C33}

variable C12all equal 0.5*(${C12}+${C21})
variable C13all equal 0.5*(${C13}+${C31})
variable C23all equal 0.5*(${C23}+${C32})

variable C44all equal ${C44}
variable C55all equal ${C55}
variable C66all equal ${C66}

variable C14all equal 0.5*(${C14}+${C41})
variable C15all equal 0.5*(${C15}+${C51})
variable C16all equal 0.5*(${C16}+${C61})

variable C24all equal 0.5*(${C24}+${C42})
variable C25all equal 0.5*(${C25}+${C52})
variable C26all equal 0.5*(${C26}+${C62})

variable C34all equal 0.5*(${C34}+${C43})
variable C35all equal 0.5*(${C35}+${C53})
variable C36all equal 0.5*(${C36}+${C63})

variable C45all equal 0.5*(${C45}+${C54})
variable C46all equal 0.5*(${C46}+${C64})
variable C56all equal 0.5*(${C56}+${C65})

# Average moduli for cubic crystals

variable C11cubic equal (${C11all}+${C22all}+${C33all})/3.0
variable C12cubic equal (${C12all}+${C13all}+${C23all})/3.0
variable C44cubic equal (${C44all}+${C55all}+${C66all})/3.0

variable bulkmodulus equal (${C11cubic}+2*${C12cubic})/3.0
variable shearmodulus1 equal ${C44cubic}
variable shearmodulus2 equal (${C11cubic}-${C12cubic})/2.0
variable poissonratio equal 1.0/(1.0+${C11cubic}/${C12cubic})
  
# For Stillinger-Weber silicon, the analytical results
# are known to be (E. R. Cowley, 1988):
#               C11 = 151.4 GPa
#               C12 = 76.4 GPa
#               C44 = 56.4 GPa


include output.mod
'''


    def write_in(self):
        print self.in_lmp 
        with open(self.in_lmp, 'w') as inscript:
             inscript.write( self.in_elastic )


    def write_output_mod(self):
        with open('output.mod', 'w') as inscript:
             inscript.write( self.output_mod )

    def write_init_mod(self):
        with open('init.mod', 'w') as inscript:
             inscript.write( self.init_mod )

    def write_displace_mod(self):
        with open('displace.mod', 'w') as inscript:
             inscript.write( self.displace_mod )


    def write_potential_mod(self, potential):
        self.potential_mod = self.potential_mod % potential
        with open('potential.mod', 'w') as inscript:
             inscript.write( self.potential_mod )


    def write_scripts(self, potential):
        self.write_in()
        self.write_output_mod()
        self.write_potential_mod(potential)
        self.write_init_mod()
        self.write_displace_mod()



def test_01():
    sys_setting ={'elements':['Al'], 'pot':'zhou', \
                  'pca':[], 'nAtoms':250,\
              #'structure':'bcc',\
              #'positions':'rnd','a':3.0, 'period':[5,5,5]}

              'structure':'fcc',\
              'positions':'rnd','a':4.5, 'period':[5,5,5]}

    sys1 = sysAt.System(sys_setting)

    lammps_setting = {'data_lmp':'data.lmp', 
                      'in_lmp':'in.min',
                      'lammps_exe':\
                              '/opt/lmpizarro/GitHub/lammps/src/lmp_serial',
                              'log': 'log.lammps', 'sys':sys1, 'minimize':True,
                              'elastic': True}

    elasTic = Elastic(lammps_setting)
    mi_potential ='sooy un potencial'
    elasTic.write_scripts(mi_potential)

if __name__ == '__main__':
    test_01()
