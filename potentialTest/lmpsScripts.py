fix_min='''
fix 1 all box/relax iso 0.0 vmax 0.001
min_style cg 
minimize 1e-25 1e-25 5000 10000
'''

fix_nve ='''
fix 1 all nve
'''

pot_eam='''
### interactions 
pair_style eam/alloy 
pair_coeff * * U2.eam U 
'''
pot_eam_zhou='''
### interactions 
pair_style eam/alloy 
pair_coeff * * Zhou_ZrAl.setfl Zr 
'''

pot_sw='''
### interactions 
pair_style	sw
pair_coeff * * si.sw Si
'''

pot_meam='''
### interactions
pair_style meam 
pair_coeff * * meamf AlS SiS MgS CuS FeS meam.alsimgcufe AlS MgS FeS
'''

global in_frame
in_frame = '''
clear
units metal 
atom_style atomic
boundary p p p 
atom_modify sort 0 0.0 

### read atom configuration 
read_data %s 

mass 1 1.0e-20
### interactions 
%s
 
# Set Thermo 
thermo 10
thermo_style custom step temp press cpu pxx pyy pzz pxy pxz pyz ke pe etotal vol lx ly lz atoms 
thermo_modify flush yes

dump dump_all all custom 1 u.dump id type x y z vx vy vz fx fy fz
timestep 0.001

# run fix min o nve
%s

run 0
unfix 1
#fix 2 all nve
#run 1
#unfix 2
'''

def create_in(atoms_lmps, interaction, fix):
    in_lammps = in_frame % (atoms_lmps, interaction, fix)

    with open('in.min', 'w') as inscript:
         inscript.write( in_lammps )
