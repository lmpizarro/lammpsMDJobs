#!/bin/bash

# This example shows how to build an interface
# between fcc copper and bcc iron using atomsk.

rm -f interface.* *.atsk

# First, create the bottom crystal of fcc Cu
atomsk --create fcc 4.24 Cu orient [100] [010] [001] -duplicate 30 30 30 Cu_supercell.xsf


atomsk Cu_supercell.xsf -select out box 0.25*BOX 0.75*BOX  0.25*BOX 0.75*BOX  0.25*BOX 0.75*BOX -rmatom select lammps