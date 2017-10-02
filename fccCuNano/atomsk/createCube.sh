#!/bin/bash

# This example shows how to build an interface
# between fcc copper and bcc iron using atomsk.

rm -f interface.* *.atsk

# First, create the bottom crystal of fcc Cu
atomsk --create fcc 4.24 U orient [100] [010] [001] -duplicate 10 10 10 cfg
atomsk U.cfg -select out box 0.25*BOX 0.75*BOX  0.25*BOX 0.75*BOX  0.25*BOX 0.75*BOX -rmatom select U_r.cfg

atomsk --create fcc 4.24 Al orient [100] [010] [001] -duplicate 10 10 10 cfg
atomsk Al.cfg -select in box 0.25*BOX 0.75*BOX  0.25*BOX 0.75*BOX  0.25*BOX 0.75*BOX -rmatom select Al_r.cfg

atomsk --merge 2 Al_r.cfg U_r.cfg Al_U_r.cfg

Atomeye Al_U_r.cfg
