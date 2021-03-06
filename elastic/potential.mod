# NOTE: This script can be modified for different pair styles 
# See in.elastic for more info.

# Choose potential
#pair_style	sw
#pair_coeff * * Si.sw Si

#pair_style meam 
#pair_coeff * * library.meam U U.meam U

#pair_style meam 
#pair_coeff * * library.meam  SiS AlSiMgCuFe.meam SiS 

pair_style eam/alloy
pair_coeff * * Zhou_Al.setfl Al 



#pair_style  	meam/spline
#pair_coeff   	* * Ti.meam.spline Ti 

# https://sites.google.com/site/eampotentials/Au
#pair_style  	eam/alloy
#pair_coeff   	* *  Au.lammps.eam Au

# https://sites.google.com/site/eampotentials/Al
#pair_style  	eam/alloy
#pair_coeff   	* *  Al.lammps.eam Al

# https://sites.google.com/site/eampotentials/Cu
#pair_style  	eam/alloy
#pair_coeff   	* *  Cu.lammps.eam Cu

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
