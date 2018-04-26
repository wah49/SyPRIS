#!/bin/usr/env python

import math, sys, os
from transform import *
from eutilities import *

def external_rotator(pdb_file, rotation_degrees):
    degree_rotate = float(rotation_degrees)
    pdb = read_file(pdb_file)
    t_pdb = Transform(pdb)
    t_pdb.rotate_object([0.0,0.0,0.0],[0.0,0.0,10.0], math.radians(degree_rotate))
    write_file(t_pdb.return_to_pdb(), "klab2internally_rotated" + str(degree_rotate) + ".pdb")
    t_pdb.translate([100.0,0.0,0.0])
    write_this = t_pdb.return_to_pdb()
    write_file(write_this,"klab2translated_pdb.pdb")
    t_pdb.rotate_object([0.0,0.0,0.0],[0.0,0.0,10.0], math.radians(degree_rotate))
    write_that = t_pdb.return_to_pdb()
    write_file(write_that, "klab2externally_rotated_"+str(rotation_degrees)+".pdb")

def main():
    args = sys.argv
    pdb_file = args[1]
    rotation_amount = args[2]
    external_rotator(pdb_file, rotation_amount)

if __name__ == "__main__":
    main()
