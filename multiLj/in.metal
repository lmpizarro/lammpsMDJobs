# 3d Lennard-Jones melt
echo both
units		metal
atom_style	atomic

read_data Al.lmp

variable tinit equal 300
variable tend equal 300
variable nruns1 equal 1000
variable nruns2 equal ${nruns1}*2
variable nruns3 equal ${nruns2}*10
variable tdump equal .1
variable pdump equal 1
variable pdrag equal 1


variable Cu equal 1
variable Al equal 2
variable Ni equal 3

variable aCu equal 3.6147
variable mCu  equal 63.546
variable rcCu equal ${aCu}*1.3
variable sCu  equal 2.29726
variable eCu  equal 0.52031

variable aAl equal 4.0496
variable mAl  equal 26.98
variable rcAl equal ${aAl}*0.85
variable sAl  equal 2.57366
variable eAl  equal 0.50722

variable aNi equal 3.5238
variable mNi  equal 58.6934
variable rcNi equal ${aNi}*0.85
variable sNi  equal 2.223949
variable eNi  equal 0.66092

#  Lorentz-Berthelot rules
variable rcAlCu equal (${rcAl}+${rcCu})/2.0
variable sAlCu equal (${sAl}+${sCu})/2.0
variable eAlCu equal sqrt(${eAl}*${eCu})/2.0

variable rcNiCu equal (${rcNi}+${rcCu})/2.0
variable sNiCu equal (${sNi}+${sCu})/2.0
variable eNiCu equal sqrt(${eNi}*${eCu})

variable rcNiAl equal (${rcAl}+${rcNi})/2.0
variable sNiAl equal (${sAl}+${sNi})/2.0
variable eNiAl equal sqrt(${eAl}*${eNi})

mass		${Cu} ${mCu} 
mass		${Al} ${mAl} 
mass		${Ni} ${mNi} 

pair_style	lj/cut 7
pair_coeff	${Cu} ${Cu} ${eCu} ${sCu} ${rcCu}
pair_coeff	${Al} ${Al} ${eAl} ${sAl} ${rcAl} 
pair_coeff	${Cu} ${Al} ${eAlCu} ${sAlCu} ${rcAlCu} 
pair_coeff	${Ni} ${Ni} ${eNi} ${sNi} ${rcNi} 
pair_coeff	${Al} ${Ni} ${eNiAl} ${sNiAl} ${rcNiAl} 
pair_coeff	${Cu} ${Ni} ${eNiCu} ${sNiCu} ${rcNiCu} 



# molecule           ε[×10-22J]       σ[×10-10m]
# ------------------------------------------------ 
# hydrogen           3.37          3.290                                
# oxygen             16.6          3.369 
# A Modeling of Thermal Properties of Hydrogen/Oxygen 
# System Using Molecular Simulations
# Shin-ichi Tsuda and Nobuhiro Yamanishi

                         

neighbor	4 bin
neigh_modify    delay 5 every 1

compute eng all pe/atom                                 ## compute potential energy for each atom
compute eatoms all reduce sum c_eng                     ## compute total energy for whole system
timestep 0.001
# Set thermo output
thermo 100
thermo_style custom step  temp pe ke  press vol enthalpy  c_eatoms 


dump    1  all xyz 100 output.xyz
#dump_modify 1 element Cu  2 element Al

min_style cg
minimize 1e-15 1e-15 5000 5000
run 0

velocity	all create ${tinit} 87287 dist gaussian
fix		1 all nve
run		${nruns1}
unfix 1

compute myRDF all rdf 50
fix rdf all ave/time 100 1 100 c_myRDF[*] file tmp.rdf mode vector

fix     2 all nvt temp ${tinit} ${tend} 0.001
run     ${nruns2}
unfix 2

fix 3 all npt temp ${tinit} ${tend} ${tdump} iso 0 0 ${pdump} drag ${pdrag}
run ${nruns3}

unfix 3
