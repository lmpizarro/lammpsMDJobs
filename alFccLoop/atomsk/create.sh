#!/bin/bash

# This example shows how to build an interface
# between fcc copper and bcc iron using atomsk.

rm -f interface.* *.atsk

# First, create the bottom crystal of fcc Cu
atomsk --create fcc 4.24 Al orient [100] [010] [001] -duplicate 10 10 10 lammps


