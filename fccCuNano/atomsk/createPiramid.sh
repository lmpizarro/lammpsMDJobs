atomsk --create fcc 5.431 Al -duplicate 10 10 10 Al_supercell.cfg
atomsk Al_supercell.cfg -cut above 40 [111] Al_cut.cfg

atomsk --create fcc 5.431 U -duplicate 10 10 10 U_supercell.cfg
atomsk U_supercell.cfg -cut below 40 [111] U_cut.cfg

atomsk --merge 2 Al_cut.cfg U_cut.cfg piramid.cfg

atomsk piramid.cfg -select out box 0.25*BOX 0.25*BOX  0.25*BOX 0.75*BOX  0.75*BOX 0.75*BOX -rmatom select double.cfg


#atomsk --create fcc 5.431 Al -duplicate 10 10 10 Al_supercell.cfg
atomsk Al_supercell.cfg -select in box 0.25*BOX 0.75*BOX  0.25*BOX 0.75*BOX  0.25*BOX 0.75*BOX -rmatom select hollow.cfg

atomsk --merge 2 double.cfg hollow.cfg sys.cfg
Atomeye sys.cfg
