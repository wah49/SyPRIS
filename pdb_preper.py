#!/usr/bin/env python

import math, csv, sys, optparse, os
from klab4 import *

def main(argv):

    parser = optparse.OptionParser(usage="\n\nTake SyPRIS output and mutate matched residues.")

    parser.add_option('--pdb-path', dest = 'pdb_path',
        help = 'The path to the pdbs listed in the outcsv \n \
                Example: /home/user/my_pdbs/')

    parser.add_option('--list-path', dest = 'list_path',
        help = 'The path to the pdbs listed in the outcsv \n \
                Example: /home/user/my_pdbs/')
    
    parser.add_option('-o', dest = 'final_out_path',
        default = "./",
        help = 'The path to where the output INPUT files should be stored. \n \
                Example: /home/user/my_pdbs/')

    parser.add_option('--chain-present',
        help = "If passed the program will look driectly for the chain of the pdb indicated",
        default = False,
        action = "store_true", dest='chain')

    parser.add_option('--keep_center',
        help = "If passed the program will not center the pdb indicated",
        default = False,
        action = "store_true", dest='keep_center')
    
    parser.add_option('--keep-resnum',
        help = "If passed the program will keep the wild type numbering of pdb indicated",
        default = False,
        action = "store_true", dest='wt_res')

    parser.add_option('--keep-remarks',
        help = "If passed the program will keep the remark lines of pdb indicated",
        default = False,
        action = "store_true", dest='remarks')
    
    parser.add_option('--keep-hetatm',
        help = "If passed the program will keep the hetatm lines of pdb indicated",
        default = False,
        action = "store_true", dest='hetatm')

    parser.add_option('--keep-water',
        help = "If passed the program will keep the HOH atoms of pdb indicated",
        default = False,
        action = "store_true", dest='HOH')

    (options,args) = parser.parse_args()
    
    list_of_pdbs = read_file(options.list_path)
    
    pdb_line_starters = ["ATOM  "]
    add_ons = []
    if options.remarks == True:
        pdb_line_starters.append("REMARK")
        add_ons.append("R")
    if options.hetatm == True:
        pdb_line_starters.append("HETATM")
        add_ons.append("H")
    
    for pdb_line in list_of_pdbs:
        pdb_line = pdb_line.strip('\r\n').split()
        name = pdb_line[0]
        if pdb_line[0][:-4] != '.pdb':
            name += '.pdb'
        try:
            pdb = read_file('%s%s' % (options.pdb_path, name))
            print name
        except IOError:
            print name, "Does not exist, passing and moving on."
            continue

        if options.HOH == True:
            if "W" not in add_ons:
                add_ons.append("W")
            pdb = [x for x in pdb if x[0:6] in pdb_line_starters]
        else:
            pdb = [x for x in pdb if x[0:6] in pdb_line_starters and x[17:20] != "HOH"]
            
        if "REMARK" in pdb_line_starters:
            remarks = [x for x in pdb if x[0:6] == "REMARK"]
        if options.chain == True:
            #pdb_line.split()
            #pdb = read_file('%s%s' % (options.pdb_path, pdb_line[0]))
            chain = pdb_line[1]
            pdb = [x for x in pdb if x[21] == chain]
        '''elsedd:
            pdb = read_file('%s%s' % (options.pdb_path, pdb_line.strip('\r\n')))'''

        #pdb = [x for x in pdb if x[0:6] in pdb_line_starters]
        #print pdb[0:100]
        #sys.exit()
        if options.wt_res == False:
            if "P" not in add_ons:
                add_ons.append("P")
            pdb = convert_to_pose_num(pdb)
        
        if options.keep_center == False: 
            if "C" not in add_ons:
                add_ons.append("C")

            pdb_obj = Transform(pdb)
            pdb_obj.translate([0.0,0.0,0.0])
            pdb = pdb_obj.return_to_pdb()
        
        if "REMARK" in pdb_line_starters:
            pdb = remarks + pdb
        name = name[:-4]
        if not add_ons:
            continue
        else:
            name += '_prep'
            for x in add_ons:
                name += x

        with open('%s%s.pdb' % ( options.final_out_path, name ), 'w') as myoutfile:
            myoutfile.writelines(pdb)


if __name__ == '__main__':
    main(sys.argv[1:])
