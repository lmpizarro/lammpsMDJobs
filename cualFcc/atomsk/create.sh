#!/bin/bash

# This example shows how to build an interface
# between fcc copper and bcc iron using atomsk.

rm -f interface.* *.atsk

# First, create the bottom crystal of fcc Cu
atomsk --create fcc 3.61 Cu orient [100] [010] [001] -duplicate 10 10 5 bottom_cu.atsk

# Second, create the top crystal of bcc Fe
atomsk --create fcc 3.61 Al orient [100] [010] [001] -duplicate 10 10 5 top_al.atsk

# Finally, merge the two systems and output to XYZ, XSF, CFG and LAMMPS data format
atomsk --merge z 2 bottom_cu.atsk top_al.atsk interface.xyz xsf cfg lmp

# Remove temporary files
rm -f *.atsk

