#!/bin/env python

from transform import *
from rSyPRIS_pkgs import *
import sys, math, os, time, ast

def write_dictionary_file(my_dict):
    with open('rotlib_cont.txt', 'a') as mf:
        mf.write(str(my_dict))

def average_rot_dict_items2(rot_dict):
    avg_rot_list = {}
    for key, vals in rot_dict.iteritems():
        length = len(vals)
        avg_rot_list[key] = [round(sum(x)/length,4) for x in zip(*vals)]
    return avg_rot_list
                    
def dunbrack_rotLib_reader2(my_file):
    print "in db_rotlib_reader"
    if type(my_file) != list:
        my_file = read_file(my_file)
    rot_dict = {}
    type_rot = {}
    for line in my_file:
        column_sep = line.split()
        if line == my_file[0]:
            last_resType = column_sep[0]
        if column_sep[0] == last_resType:
            rot_id = int(column_sep[4]+\
                         column_sep[5]+\
                         column_sep[6]+\
                         column_sep[7]) 
            if rot_id not in type_rot.keys():
                #type_rot[rot_id] = [[x for i,x in enumerate(column_sep) if i > 8]]
                type_rot[rot_id] = [[float(x) for x in column_sep[9:]]]
            else:
                current_data = type_rot[rot_id][:]
                #current_data.append([x for i,x in enumerate(column_sep) if i > 8])
                current_data.append([float(x) for x in column_sep[9:]])
                type_rot[rot_id] = current_data
        
        else:
            final_old_rot_avgs = average_rot_dict_items(type_rot)           
            
            rot_dict[last_resType] = final_old_rot_avgs
            type_rot = {}
            last_resType = column_sep[0]

    print rot_dict['PHE']
    write_dictionary_file(rot_dict)

def main(mf):
    rot_dict = dunbrack_rotLib_reader(mf)
    #print rot_dict['PHE']
    write_dictionary_file(rot_dict)
    #print read_dictionary_file(mf)

if __name__ == "__main__":
    s = time.clock()
    my_file = read_file(sys.argv[1])
    main(my_file)
    print time.clock() -s 
