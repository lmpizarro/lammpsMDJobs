#!/bin/bash

# This example shows how to build an interface
# between fcc copper and bcc iron using atomsk.

rm -f interface.* *.atsk

# First, create the bottom crystal of fcc Cu
atomsk --create fcc 4.24 Al orient [100] [010] [001] -duplicate 18 18 8 Al_.xsf
atomsk --create fcc 4.24 U orient [100] [010] [001] -duplicate 18 18 3 U_.xsf

atomsk --merge Z 3 Al_.xsf U_.xsf Al_.xsf U_Al.cfg

atomsk U_Al.cfg -select out box 0.3*BOX 0.3*BOX  0.3*BOX 0.70*BOX  0.70*BOX 0.7*BOX -rmatom select U_Al_nano.cfg

Atomeye U_Al_nano.cfg
