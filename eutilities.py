#!/usr/bin/env python

from transform import *
import math, sys, os, numpy, time

#writes a pdb file to a pdb file. Enter the new name you want, and it will be written to the current directory.

def vector_difference(a,b):
    return [a[0]-b[0],a[1]-b[1],a[2]-b[2]]


#Dihedral angle calculator that accepts the 4 sets of xyz coordinates
def dihedral_point_calculator(a,b,c,d):
    v1 = vector_difference(a,b)
    v2 = vector_difference(b,c)
    v3 = vector_difference(c,d)

    v12 = cross_product(v1,v2)
    v23 = cross_product(v2,v3)

    dot_over_mag = dot_product(v12,v23)/(magnitude_calc(v12)*magnitude_calc(v23))

    dihedral_angle = math.acos(dot_over_mag)
    
    return dihedral_angle


#input res string
#input initial string type (3 or 1)
def res3_to_res1(res_matrix):
    new_res_string = []
    string_length = 1
    if len(res_matrix[0]) == 3:
        string_length = 3
    
    res1_dictionary = {"A" : "ALA",\
                     "R" : "ARG",\
                     "N" : "ASN",\
                     "D" : "ASP",\
                     "C" : "CYS",\
                     "Q" : "GLN",\
                     "E" : "GLU",\
                     "G" : "GLY",\
                     "H" : "HIS",\
                     "I" : "ILE",\
                     "L" : "LEU",\
                     "K" : "LYS",\
                     "M" : "MET",\
                     "F" : "PHE",\
                     "P" : "PRO",\
                     "S" : "SER",\
                     "T" : "THR",\
                     "W" : "TRP",\
                     "Y" : "TYR",\
                     "V" : "VAL"}
    res3_dictionary = {v: k for k, v in res1_dictionary.iteritems()}
    if string_length == 3:
        for residue in res_matrix:
            new_res_string.append(res3_dictionary[residue])
    else:
        for residue in res_matrix:
            new_res_string.append(res1_dictionary[residue])
    print new_res_string
    print len(new_res_string)

    return new_res_string


#renames chains in pdbs and residue numbers.
def chain_renamer(pdb, chain_name):
    renamed = []
    for line in pdb:
        if line.startswith(("ATOM","HETATM")):
            line_out = line[:21] + chain_name + line[22:]
            renamed.append(line_out)

    return renamed

#renumbers the residues in a pdb file. Ensures rosetta corresponds correctly with the pdb.
#plug in a read pdb file into pdb and declare if you want residues to continue to increment through the different chains or reset for each chain
def residue_renumber(pdb,res_increment = False):
    
    res = 0
    atom_count = 0
    renumpdb = []
    for line in pdb:
        if line.startswith(("ATOM","HETATM")):
            atom_count += 1
            if res == 0:
                cur_chain = line[21]
            if line[12:16].strip() == 'N':
                if line[21] == cur_chain:
                    res += 1
                else:
                    if res_increment == False:
                        res = 1
            cur_chain = line[21]
            renumpdb.append(line[:6]+str(atom_count).rjust(5,' ')+line[11:22]+str(res).rjust(4,' ')+line[30:].rjust(55,' '))
        else:
            renumpdb.append(line)
    return renumpdb

###TAKE THIS####
#finds the names of chains in the pdb
def chaingrabber(pdb):
    chains = []
    for line in pdb:
        if line.startswith(("ATOM", "HETATM")):
            curchain = line[21]
            if curchain not in chains:
                chains.append(curchain)
    return chains


###TAKE THIS####
def gen_alphabet():
    return "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"

def gen_bb_atomtype():
    return ['N','CA','C','O','OXT']

###TAKE THIS####
def atom_renumber(pdb):
    new_pdb = []
    atom_count = 0
    for lines in pdb:
        if lines.startswith(("ATOM","HETATM")):
            atom_count += 1
            new_pdb.append(lines[:6]+str(atom_count).rjust(5,' ') + lines[11:])
    return new_pdb

###TAKE THIS####
def res_sequence_extractor(pdb):
    res_sequence = []
    for lines in pdb:
        if lines.startswith(("ATOM","HETATM")):
            if lines[11:17].strip() == "N":
                res_sequence.append(lines[17:20].strip())
                #print lines[17:20].strip()
    return res_sequence


###TAKE THIS####  modify for all termini (breaks included) write in class
def first_last_pdb_residues(pdb,chain):
    first_res = []
    last_res = []
    for lines in pdb:
        if lines.startswith(("ATOM","HETATM")):
            if lines[21] == chain:
                cur_res = int(lines[23:27].strip())
                if not first_res or cur_res < first_res:
                    first_res = cur_res
                if not last_res or cur_res > last_res:
                    last_res = cur_res
    return first_res, last_res

def standard_deviation(x):
    mean = 0.0
    length = float(len(x))
    print x
    for y in x:
        mean += y/length
    
    deviation_from_mean = 0.0
    for y in x:
        deviation_from_mean += abs(y-mean)**2.0

    stdev = (deviation_from_mean/length) ** 0.5
    return stdev

#This code takes a container and a domain.
#It finds the indices in the container that contains the domain
#if the entire domain is in the continer, it will find the entire domain
#it will find the largest chunk of domain in the container.
#eg. Container: ASLKDENGLIKDJLDLKNE Domain: EFDENGLIKDFEASF
#Similar sequence:efDENGLIKDfeasf  caps in the container.
#the indices will be for DENGLIKD in the container. So: 5 and 12 

####TAKE THIS#####
def domain_index(container, domain):
    if len(container) < len(domain):
        print "No domain in this container. And You broke the code"
        return
    
    container_res = res_sequence_extractor(container)
    domain_res = res_sequence_extractor(domain)
    res_container = ''.join(res3_to_res1(''.join(container_res),3))
    res_domain =    ''.join(res3_to_res1(''.join(domain_res), 3))

    print "Single Residue Container:" + str(len(res_container))
    print res_container
    print "Single Residue Domain:" + str(len(res_domain))
    print res_domain

    index_check = False
    test_domain = res_domain
    cutter = checker = trim_count = 0

    while index_check == False:
        checker += 1
        if test_domain in res_container:
            start_index = res_container.index(test_domain)
            end_index = start_index + len(test_domain)
            index_check = True

        else:
            if float(len(test_domain)) < float(len(res_domain))/ 5.:
                print "Domain not found with high enough accuracy. \nThe testing domain is less than 1/5 the original size. \nThe smallest it should be is 1/3 or \nwe are finding small chunks of the protein that gives bad results.\n Especially if the chunks are repeats."
                start_index = 0
                end_index = 0
                index_check = True
                break

            if len(test_domain) + trim_count == len(res_domain):
                cutter += 1
                trim_count = 0
                trim = len(res_domain) - cutter
                test_domain = res_domain[:trim]
                print str(cutter)
            
            else:
                print str(trim_count)
                trim_count += 1
                test_domain = res_domain[trim_count:trim+trim_count]

    print "Total times through the while loop: " + str(checker)

    print "Start Index: " + str(start_index)
    print "End Index: " + str(end_index)

    return start_index, end_index


###TAKE THIS####  modify
def mutant_indexer(pdb1, pdb2, res_type = []):
    pdb1_sequence = res_sequence_extractor(pdb1)
    pdb2_sequence = res_sequence_extractor(pdb2)
    mutant_indices = []
    for index, res in enumerate(pdb1_sequence):
        #print "index: " + str(index) + "  pdb1 : " + str(pdb1_sequence[index]) + "  pdb2: " + str(pdb2_sequence[index])
        if pdb1_sequence[index] != pdb2_sequence[index]:
            if not res_type:
                mutant_indices.append(index)
            elif res_type == pdb2_sequence[index]:
                mutant_indices.append(index)
    #print "mutant count: " + str(mutant_count)
    return mutant_indices
            


### TAKE THIS #####
def vector_extension(vector, distance):
    initial_mag = magnitude_calc(vector)
    final_mag = initial_mag + distance
    square_initial = initial_mag * initial_mag
    square_final = final_mag * final_mag
    vector_squares = [x*x for x in vector]
    center_dir = [x/abs(x) for x in vector]
    #print center_dir
    vector_proportions = [x/square_initial*square_final for x in vector_squares]
    vector_final = [x ** 0.5 for x in vector_proportions]
    finale = []
    for x in xrange(len(vector_final)):
        finale.append(vector_final[x] * center_dir[x])
    print finale
    return finale
    

