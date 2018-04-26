#!/bin/usr/env python

import math, sys, os
from AssemB_pkgs import *
from eutilities import *

def assemb_tester(pdb_file_1, pdb_file_2):
    two_component = TwoComponentAssem(pdb_file_1, pdb_file_2)
    print two_component.get_axis_coord_by_name('A')
    print two_component.get_vector_difference()
    two_comp_assem = Cn_Cn_AssemMaker(pdb_file_1, pdb_file_2)
    two_comp_assem.generate_alignments()
    #two_comp_assem.generate_conformation(0.,0.,0.,100.,False)
    #two_comp_assem.generate_conformation(100.,0.,0.,0.,False)
    #two_comp_assem.generate_conformation(100.,0.,0.,100.,False)
    #two_comp_assem.generate_conformation(100.,-60.,0.,100.,False)
    two_comp_assem.generate_conformation(20.,-60.,90.,20.,True)
    ten_grid = two_comp_assem.mobile_obj.get_atom_grid(10.0)
    create_arbPDB_from_xyz(convert_grid_string_to_list(ten_grid), 'ten.pdb')
    create_arbPDB_from_xyz(convert_grid_string_to_list(two_comp_assem.mobile_obj.get_atom_grid(5.0)), 'five.pdb')
    create_arbPDB_from_xyz(convert_grid_string_to_list(two_comp_assem.mobile_obj.get_atom_grid(3.0)), 'three.pdb')
    sys.exit()
    two_comp_assem.static_obj.write_to_file('STATIC.pdb')
    print two_comp_assem.mobile_obj.get_symmetric_axes()
    create_arbPDB_from_xyz([x for i,x in two_comp_assem.mobile_obj.get_symmetric_axes().iteritems()], 'HIHI.pdb')
    two_comp_assem.mobile_obj.write_to_file('LATER.pdb')
    for name, thingy in two_comp_assem.ensemble.iteritems():
        thingy.write_to_file('%s.pdb' % name)
        create_arbPDB_from_xyz([x for i,x in thingy.get_symmetric_axes().iteritems()], 'YOOO.pdb')
    sys.exit()


def main():
    args = sys.argv
    pdb1 = args[1]
    pdb2 = args[2]
    assemb_tester(pdb1, pdb2)

if __name__ == "__main__":
    main()
