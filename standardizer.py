#!/usr/bin/env python
"""
Multigen.py is property of the Khare Lab at Rutegers University.
Created by: William Hansen
Please refrence the lab and creator if you wish to use this program.

Multigen.py is a program that, when given a sub unit amount, will 
search through a pdb file for all oligomeric biomolecules of the 
chosen subtype. A file will be created for each biomolecule using 
the BIOMT transformation matrices provided.
"""

import os, sys, math

#----------------------------------------------------------
#Program will run to main.
#multimerize is the main definition and takes 3 arguements:
#pdb = pdb file
#mer_num = sub unit number input from main
#mer_name = name equivlent to sub unit number
#----------------------------------------------------------      
def initial_matrix(seq1, seq2):
    '''This function take two sequences, find their length m and n, return two (m+1)*(n+1) matrices (nested list) filled with value of 0 (score_matrix) and "none" (pointer_matrix).'''

    m = len(seq1)+1
    n = len(seq2)+1
    score_matrix = []
    pointer_matrix = []
    for i in range(m):
        nested_list = []
        nested_list2 = []
        for i2 in range(n):
            nested_list.append(0)
            if i == 0 and i2 > 0:
                nested_list2.append('left')
            elif i2 == 0 and i > 0:
                nested_list2.append('up')
            else:
                nested_list2.append('none')
        score_matrix.append(nested_list)
        pointer_matrix.append(nested_list2)
    return score_matrix, pointer_matrix #score_matrix stores the computed scores, pointer_matrix stores the information ("diagonal", "up", "left" or "none") for tracking back
def compute_matrix(seq1, seq2):
    '''This function creates score_matrix and pointer_matrix'''
    #call take_parameter() for scores
    match = +1
    mismatch = -1
    gap = 0

    #call initial_matrix() for initializaiton of score_matrix and pointer_matrix
    (score_matrix,pointer_matrix) = initial_matrix(seq1, seq2)

    #initialization of the diagonal, up and left scores.
    diagonal_score = 0 # Score for the "diagonal" step when tracking back, either a match or a mismatch
    up_score = 0 # Score for the "up" step when tracking back, insert a gap in the (first or second?) sequence
    left_score = 0 # Score for the "left" step when tracking back, insert a gap in the (first or second?) sequence

    #initialization of the max score and its position in matrix
    max_score = 0
    iMax = 0
    jMax = 0

    #fill in the score_matrix and keep track ("diagonal", "up", "left" or "none") in the pointer_matrix. Also keep track of the max_score and its position.
    #use nested for loops to accomplish your task --> return what you are asked for

    for index, letter in enumerate(seq1):
        for index2, letter2 in enumerate(seq2):
            if letter == letter2:
                score_matrix[index+1][index2+1] += int(match) + score_matrix[index][index2]
                pointer_matrix[index+1][index2+1] = 'diagonal'
                if score_matrix[index+1][index2+1] > max_score:
                    max_score = score_matrix[index+1][index2+1]
                    iMax = index2+1
                    jMax = index+1
            if letter != letter2:
                if score_matrix[index][index2] >= score_matrix[index+1][index2] >0 and \
                   score_matrix[index][index2] >= score_matrix[index][index2+1] >0:
                    score_matrix[index+1][index2+1] += int(mismatch) + score_matrix[index][index2]
                    pointer_matrix[index+1][index2+1] = 'diagonal'
                else:
                    if score_matrix[index+1][index2] > score_matrix[index][index2+1] >= 0:
                        score_matrix[index+1][index2+1] += int(gap) + score_matrix[index+1][index2]
                        pointer_matrix[index+1][index2+1] = 'left'
                    elif score_matrix[index][index2+1] > score_matrix[index+1][index2] >= 0:
                        score_matrix[index+1][index2+1] += int(gap) + score_matrix[index][index2+1]
                        pointer_matrix[index+1][index2+1] = 'up'
                    else:
                        score_matrix[index+1][index2+1] = 0
                        pointer_matrix[index+1][index2+1] = 'None'

    return score_matrix, pointer_matrix, iMax, jMax


def print_score_matrix(score_matrix):
    '''This function prints score matrix'''

    for item in score_matrix:
        print ' '.join(str(x) for x in item)
    return

def track_back(pointer_matrix, seq1, seq2, iMax, jMax):
    '''Tracks back to create the aligned sequence pair'''
    #initialization of aligned sequences
    aligned_seq1 = ''
    aligned_seq2 = ''

    #start from the position where the max score is.
    i = iMax
    j = jMax
    point = ''
    while pointer_matrix[j][i] != 'none':
        if pointer_matrix[j][i] == 'diagonal':
            i -= 1
            j -= 1
            aligned_seq1 += seq1[j]
            aligned_seq2 += seq2[i]
            point = pointer_matrix[j][i]
        elif pointer_matrix[j][i] == 'up':
            i -= 0
            j -= 1
            aligned_seq1 += '_'
            #aligned_seq2 += seq2[i]
            point = pointer_matrix[j][i]
        else:
            i -= 1
            j -= 0
            #aligned_seq1 += seq1[j]
            aligned_seq2 += '_'
            point = pointer_matrix[j][i]

    #track backwards in the pointer_matrix, based on the info ("diagonal", "up", "left" or "none") decide for each position whether it is a match/mismatch, a gap in the first/second sequence, or to stop.
    aligned_seq1 = aligned_seq1[::-1]
    aligned_seq2 = aligned_seq2[::-1]

    return (aligned_seq1, aligned_seq2)

def three_to_one(seq):
    d = {'ALA': 'A',\
         'ARG': 'R',\
         'ASN': 'N',\
         'ASP': 'D',\
         'CYS': 'C',\
         'GLN': 'Q',\
         'GLU': 'E',\
         'GLY': 'G',\
         'HIS': 'H',\
         'ILE': 'I',\
         'LEU': 'L',\
         'LYS': 'K',\
         'MET': 'M',\
         'PHE': 'F',\
         'PRO': 'P',\
         'SER': 'S',\
         'THR': 'T',\
         'TRP': 'W',\
         'TYR': 'Y',\
         'VAL': 'V',\
         'MLY': 'B',\
         'PTR': 'J',\
         'TPO': 'O',\
         'SEP': 'U'}
    
    new_seq = ''
    for code in seq:
        new_seq += d[code]
    
    return new_seq

def one_to_three(seq):
    d = {'A': 'ALA',\
         'R': 'ARG',\
         'N': 'ASN',\
         'D': 'ASP',\
         'C': 'CYS',\
         'Q': 'GLN',\
         'E': 'GLU',\
         'G': 'GLY',\
         'H': 'HIS',\
         'I': 'ILE',\
         'L': 'LEU',\
         'K': 'LYS',\
         'M': 'MET',\
         'F': 'PHE',\
         'P': 'PRO',\
         'S': 'SER',\
         'T': 'THR',\
         'W': 'TRP',\
         'Y': 'TYR',\
         'V': 'VAL',\
         'B': 'MLY',\
         'J': 'PTR',\
         'O': 'TPO',\
         'U': 'SEP'}
    
    new_seq = []
    for code in seq:
        if code == '_':
            new_seq.append('___')
        else:
            new_seq.append(d[code])
    
    return new_seq


def mag_checker(xyz1, xyz2):
    difx_sq = (xyz1[0] - xyz2[0])**2
    dify_sq = (xyz1[1] - xyz2[1])**2
    difz_sq = (xyz1[2] - xyz2[2])**2
    mag = math.sqrt(difx_sq + dify_sq + difz_sq)
    return mag

def polar_coord(xyz, xyz2):

    radius_xyz = math.sqrt((xyz[0]**2)+(xyz[1]**2)+(xyz[2]**2))
    radius_xyz2 = math.sqrt((xyz2[0]**2)+(xyz2[1]**2)+(xyz2[2]**2))

    angle = math.acos(((xyz[0]*xyz2[0]) + (xyz[1]*xyz2[1]) + (xyz[2]*xyz2[2])) / \
            (radius_xyz * radius_xyz2))

    return angle

def vector_finder(xyz, xyz2):

    vector = [xyz2[0]-xyz[0], xyz2[1]-xyz[1], xyz2[2]-xyz[2]]
    
    return vector


def chain_builder(chain, final_ca):
    final_chain = []
    for line in chain:
        if int(line[22:30]) in final_ca:
            final_chain.append(line)
    return final_chain        

def block_create(chain, res_list):
    resi_blocks = []
    res_index = []
    for index, line in enumerate(chain):
        #wordlist = line.split()
        prev_line = chain[index -1]
        #prev_wordlist = prev_line.split()
        prev_resid = prev_line[22:27]
        if prev_resid != line[22:27]:
            
            res_index.append(index)
    
    res_index.append(len(chain))
    print res_index
    print ""
    print len(res_index)
    for ind in range(len(res_index) -1):
        if ind != len(res_index)-1:
            begin = res_index[ind]
            end = res_index[ind+1]
            residue = []
            for j in range(begin, end):
                line = chain[j]
                residue.append(line)
        resi_blocks.append(residue)
    
    num = 0
    res_list2 = [x[0][17:20] for x in resi_blocks]
    #print resi_blocks
    print len(resi_blocks)
    print res_list2
    print ""
    print res_list
    print len(res_list)
    new_resi_blocks = []
    for res_ind, res in enumerate(res_list):
        if res == '___':
            num += 1
            continue
        else:
            new_resi_blocks.append(resi_blocks[res_ind])
    #print len(resi_blocks)
    master = []
    for residues in new_resi_blocks:
        for line in residues:
            master.append(line)
    #print master
    return master, num

def chain_alignment(chains,names):

    num_chains = len(chains)
    chaincounts = []

#counting the chain length in the chains
    for x in range(0,num_chains):
        chaincounts.append(len(chains[x]))
    
    sequences = residues = final_sequences = [[] for n in range(0,num_chains)]
    original_lengths = []
    delta_lengths = []
    for x in range(0, num_chains):
        for line in chains[x]:
            if line[13:15] == 'CA':
                residues[x].append(line[17:20])
        
        sequences[x] = three_to_one(residues[x])
        original_lengths.append(len(sequences[x]))

        print "The " + names[x] + " sequence is: " + sequences[x]
        print "    the length is " + str(len(sequences[x]))


    #Call compute_matrix() and calculate the matrix.
    for x in range(0,num_chains):
        for y in range(x+1, num_chains):
            (score_matrix, pointer_matrix, iMax, jMax) = compute_matrix(sequences[x], sequences[y])
            (sequences[x], sequences[y]) = track_back(pointer_matrix, sequences[x], sequences[y], iMax, jMax)
    
    #check these lengths to see if the sequences are good.
    length_check = []

    for x in range(0, num_chains):
        delta_lengths.append(original_lengths[x]-len(sequences[x]))
        length_check.append(len(sequences[x]))
        if delta_lengths[x] > 0:
            for y in range(delta_lengths[x]):
                sequences[x] +='_'
    
    final_lengths = []
    #convert to 3 letter amino acids
#Final conversions to final sequences, getting final lengths of sequences
    for y in range(0, num_chains):
        print "Delta Length of chain " + names[y] + " is " + str(delta_lengths[y])
        sequenced = list(sequences[y])
        residues[y] = one_to_three(sequenced)
        #print residues[y]
        final_sequences[y], number_of_cutaway_res = block_create(chains[y], residues[y])
        #print final_sequences[y]
        final_lengths.append(len(final_sequences[y]))   
    
    print final_lengths
    print '\n\n'
    print number_of_cutaway_res, "Residues found to be different"
    print '\n\n'

#Comparing the sequence lengths equal eachother

    if len(set(length_check)) <= 1:# == D_fin_len:
        last_change = True
    else:
        last_change = False

    return final_sequences, last_change, number_of_cutaway_res

def chain_sep(pdb):
    standard = []
    with open(pdb, 'r') as pdb_file:
        pdblist = list(pdb_file)
        pdb_file.close()
        atom_list = []
        chains = []
        chainindex = []
        indexcount = -1
        for line in pdblist:
            wordlist = line.split()
            if wordlist[0] == "ATOM":
                indexcount += 1
                atom_list.append(line)
                if line[21] not in chains:
                    chains.append(line[21])
                    chainindex.append(indexcount)
        
        chainindex.append(len(pdblist)+10)
        #print chains
        new_atomlist = []
        
        for index, line in enumerate(atom_list):
            for k in range(0,len(chainindex)):
                if index >= chainindex[k] and index < chainindex[k+1]:
                    newline = line[0:21] + chains[k] + line[22:]
                    new_atomlist.append(newline)

        num_chains = len(chains)

        if not new_atomlist:
            print("will use original atom_list")
        else:
            atom_list = new_atomlist
            chainletters = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z']
        
        chainnames = chainletters[0:num_chains]

        #print chains
        ca_num = []
        chain = [[] for k in range(0, num_chains)]
        num = 0
        for index, line in enumerate(atom_list):
            for k in range(0,num_chains):
                if index >= chainindex[k] and index < chainindex[k+1]:
                    newline = line[0:21] + chainnames[k] + line[22:]
                    chain[k].append(newline)
                    if line[13:15] == "CA":
                        num += 1
                    else:
                        continue
                    ca_num.append(num)
        
        #These sets below are coords for the nitrogen atoms of the first residue.
        #Used to compare the distnace of the individual chains.
        
        last_change = False
        while last_change == False:
            (chain, last_change, number_of_cutaway_res) = chain_alignment(chain, chainnames) 
        
        if last_change == True:
            final_chains = chain
        
        for x in range(0,num_chains):
            standard += final_chains[x] 
    
    #print standard
    with open('%s_standard.pdb' % (pdb[:-4]), 'w') as myfile:
        for line in standard:
            myfile.write(line)

def main(argv):
    args = sys.argv
    pdb = args[1]
    chain_sep(pdb)
                
if __name__ == "__main__":
    main(sys.argv[1:])
