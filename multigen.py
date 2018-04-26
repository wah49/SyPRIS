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

import os, sys, glob, re

#----------------------------------------------------------
#Program will run to main.
#multimerize is the main definition and takes 3 arguements:
#pdb = pdb file
#mer_num = sub unit number input from main
#mer_name = name equivlent to sub unit number
#----------------------------------------------------------
def multimerize(pdb, mer_num, mer_name):
    with open(pdb, 'r') as pdb_file:
        pdblist = list(pdb_file)
        pdb_file.close()
        remark_list = []
        atom_list = []
        
        rem_nums = []
        ecoli = True
        for line in pdblist:
            wordlist = line.split()
            if wordlist[0] == "SOURCE":
                if line[11:46] == "EXPRESSION_SYSTEM: ESCHERICHIA COLI":
                    ecoli = True
            if wordlist[0] == "REMARK":
                rem_nums.append(wordlist[1])
        if "350" not in rem_nums:
            print("This pdb does not contain a REMARK 350 section.")
            return

        #------------------------------------------------
        #Will iterate through all lines in pdb file and 
        #create lists of remark and atom lines.
        #-----------------------------------------------
        for line in pdblist:
            wordlist = line.split()
            if wordlist[0] == "REMARK" and\
                wordlist[1] == "300"   and\
                len(wordlist) > 2      and\
                wordlist[2] == "BIOMOLECULE:":
                    num_of_biomolecules = len(wordlist) - 3
                    
            elif wordlist[0] == "REMARK":
                remark_list.append(line)

            elif wordlist[0] == "ATOM"    or\
                 wordlist[0] == "MODEL"   or\
                 wordlist[0] == "ENDMDL":
                atom_list.append(line)

            elif wordlist[0] == "HETATM":
                if str(line[17:20]) != "HOH":
                    atom_list.append(line)
                else:
                    continue
        
        atom_models = []
        model_indices = []
        has_model = False
        for line_num, line in enumerate(atom_list):
            if line.startswith("MODEL"):
                model_indices.append(line_num+1)
                has_model = True
        if not model_indices:
            atom_models.append(atom_list)
        else:
            for model_num in range(len(model_indices)):
                if model_num != len(model_indices)-1:
                    begin = model_indices[model_num]
                    end = model_indices[model_num+1]
                    model = []
                    for line_num in range(begin,end):
                        line = atom_list[line_num]
                        model.append(line)
                atom_models.append(model[:-2])
        #--------------------------------------------------------
        #bio_indices list will contain line numbers of individual 
        #biomoecules within a pdb file. Every two list items are 
        #the beginning and end for each biomolecule. Biomolecules
        #are found in the REMARK 350 section of every pdb file.
        #--------------------------------------------------------
        bio_indices=[]
        last = False
        for line in remark_list:
            words = line.split()
            if int(words[1]) > 350:
                remark_end = words[1]
                break
        
        for line_number in range(len(remark_list)):
            for biomolecule_number in range(num_of_biomolecules):
                if remark_list[line_number].rstrip() == 'REMARK 350 BIOMOLECULE: %s' % str(biomolecule_number+1):
                    bio_indices.append(line_number)
                elif remark_list[line_number].rstrip() == 'REMARK 350 BIOMOLECULE:  %s' % str(biomolecule_number+1):
                    bio_indices.append(line_number)
                if remark_list[line_number].rstrip() == ('REMARK %s' % remark_end) and not last:
                    last=True
                    bio_indices.append(line_number)

        #--------------------------------------------------------
        #chain_indices list will contain the chains in which each
        #transformation must be made. Every list item are the 
        #chain letters for each biomolecule.
        #--------------------------------------------------------
        chain_indices=[]
        for index, line in enumerate(remark_list):    
            remove_list = ['REMARK', '350', 'APPLY', 'THE', 'FOLLOWING', 'TO', 'CHAINS:']
            if remark_list[index].startswith('REMARK 350 APPLY THE FOLLOWING TO CHAINS: ') == True:
                chain_indices.append(index)
                line.rstrip()
                chain_list = line.split()
                x = ' '.join([j for j in chain_list if j not in remove_list])
                x = x.split(', ')
                print 'x is: ', x
        #--------------------------------------------------------
        #mybiomolecules list items are all lines (as strings). The 
        #bio_indices items
        saved_mers = []
        mybiomolecules =[]
        for i in range(len(bio_indices)-1):
            if i != len(bio_indices)-1:
                begin = bio_indices[i]
                end = bio_indices[i+1]
                mymol=[]
                for j in range(begin,end):
                    line = remark_list[j]
                    mymol.append(line)
            mybiomolecules.append(mymol)
        
        #print "Mybiomolecules are: ", mybiomolecules
        for mol in mybiomolecules: 
            if mer_name in mol[1] and num_of_biomolecules > 0:
                saved_mers.append(mol)
                
        #print("saved mers are: ", saved_mers)
        for mol in saved_mers:
            apply_indices = []
            for apply_i, apply_l in enumerate(mol):
                if apply_l.startswith('REMARK 350 APPLY THE FOLLOWING TO CHAINS: ') == True:
                    apply_indices.append(apply_i)
                if apply_i == len(mol)-1:
                    apply_indices.append(apply_i+1)
            
            if len(apply_indices) == 1:
                apply_cases = []
                myapply = []
                for line in mol:
                    myapply.append(line)
                apply_cases.append(myapply)
            
            elif len(apply_indices) > 1:
                apply_cases = []
                for obj_i in range(len(apply_indices)-1):
                    if obj_i != len(apply_indices)-1:
                        begin2 = apply_indices[obj_i]
                        end2 = apply_indices[obj_i+1]
                        myapply=[]
                        for j in range(begin2,end2):
                            line = mol[j]
                            myapply.append(line)
                    apply_cases.append(myapply)
            
            #print("apply cases are: ", apply_cases)
            biomol_num = mol[0].split(': ')[-1].rstrip()
            if len(biomol_num) > 1:
                biomol_num = biomol_num.split(' ')[-1].rstrip()
            if num_of_biomolecules > 1:
                biomol_name = ("_biomol_" + biomol_num)
            else:
                biomol_name = ""
            for model_num, atom_list_item in enumerate(atom_models):
                if has_model == True:
                    model_name = ("_model_" + str(model_num +1))
                else:
                    model_name = ""
                if ecoli == True:
                    with open('%s_oligo%s%s.pdb' % (pdb[:-4], model_name, \
                             biomol_name), 'w') as myfile:
                        for line in remark_list:
                            words = line.split()
                            if int(words[1]) < 350:
                                myfile.write(line)

                        for index, line in enumerate(mol):
                            myfile.write(line)
                            #need to grab chain ids here to apply biomt to atoms in last for loop
                        
                        for line in remark_list:
                            words = line.split()
                            if int(words[1]) > 350:
                                myfile.write(line)

                        for apply_block_index, apply_block in enumerate(apply_cases):    
                            chain_id_list = []
                            print "apply_block is: ", apply_block
                            for line in apply_block:
                                if line.startswith("REMARK 350 APPLY THE FOLLOWING TO CHAINS: "):
                                    new_line_ls = line.split(': ')
                                    chain_id = new_line_ls[-1].rstrip()
                                    chain_id_list_item = chain_id.split(', ')
                                    chain_id_list.extend(chain_id_list_item)

                            for x, chain in enumerate(chain_id_list):
                                print "chain is: ", chain
                                biomt_set = []
                                for biomts in xrange(1, 15):
                                    biomt_block = []
                                    for mol_line in apply_block:
                                        #print "mol_line is: ", mol_line
                                        biomt_words = mol_line.split()
                                        if len(biomt_words) > 3:
                                            pull = biomt_words[2]
                                            if re.match(pull[0:5], "BIOMT") and int(biomt_words[3]) == biomts:
                                                #print biomt_words[3]
                                                #print biomts
                                                biomtnums = []
                                                if biomt_words[2] == "BIOMT1":
                                                    biomtnums.append(float(biomt_words[4]))
                                                    biomtnums.append(float(biomt_words[5]))
                                                    biomtnums.append(float(biomt_words[6]))
                                                    biomtnums.append(float(biomt_words[7]))
                                                if biomt_words[2] == "BIOMT2":
                                                    biomtnums.append(float(biomt_words[4]))
                                                    biomtnums.append(float(biomt_words[5]))
                                                    biomtnums.append(float(biomt_words[6]))
                                                    biomtnums.append(float(biomt_words[7]))
                                                if biomt_words[2] == "BIOMT3":   
                                                    biomtnums.append(float(biomt_words[4]))
                                                    biomtnums.append(float(biomt_words[5]))
                                                    biomtnums.append(float(biomt_words[6]))
                                                    biomtnums.append(float(biomt_words[7]))
                                                biomt_block.append(biomtnums)
                                    if not biomt_block:
                                        continue
                                    else:
                                        biomt_set.append(biomt_block) 
                            
                                print "biomt set: ", biomt_set
                                for bio_ind, bio_set in enumerate(biomt_set): 
                                    xset = bio_set[0]
                                    yset = bio_set[1]
                                    zset = bio_set[2]
                                    for index, line in enumerate(atom_list_item):
                                        atom_word_list = line.split()
                                        string = atom_list_item[index]
                                        if line[21] == chain:
                                            
                                            x_comp = "%.3f" % ((float(line[30:38]) * xset[0]) + (float(line[38:46]) * xset[1]) + (float(line[46:54]) * xset[2]) + xset[3]) 
                                            y_comp = "%.3f" % ((float(line[30:38]) * yset[0]) + (float(line[38:46]) * yset[1]) + (float(line[46:54]) * yset[2]) + yset[3])
                                            z_comp = "%.3f" % ((float(line[30:38]) * zset[0]) + (float(line[38:46]) * zset[1]) + (float(line[46:54]) * zset[2]) + zset[3])
                                            
                                            x_comp = ('{0:>8}').format(x_comp)
                                            y_comp = ('{0:>8}').format(y_comp)
                                            z_comp = ('{0:>8}').format(z_comp)
                                            

                                            alphabet = 'A B C D E F G H I J K L M N O P Q R S T U V W X Y Z'
                                            alphabet_list = alphabet.split(' ')
                                            
                                            if len(apply_cases) > 1:
                                                new_x = x*2
                                            else:
                                                new_x = x
                                            if len(chain_id_list) > 1:
                                                new_bio_ind = bio_ind*2
                                            else:
                                                new_bio_ind = bio_ind
                                            new_chain = alphabet_list[new_bio_ind + new_x + apply_block_index]
                                            newline = line[0:21] + new_chain + line[22:30] + str(x_comp) + str(y_comp) + str(z_comp) + line[54:]
                                            myfile.write(newline)
                                       
                                        else:
                                            continue

def main(argv):
    args = sys.argv
    mer_num = int(args[2]) 
    if mer_num == 1:
        mer_name = "MONOMERIC"
    elif mer_num == 2:
        mer_name = "DIMERIC"
    elif mer_num == 3:
        mer_name = "TRIMERIC"
    elif mer_num == 4:
        mer_name = "TETRAMERIC"
    elif mer_num == 5:
        mer_name = "PENTAMERIC"
    elif mer_num == 6:
        mer_name = "HEXAMERIC"
    elif mer_num == 7:
        mer_name = "HEPTAMERIC"
    elif mer_num == 8:
        mer_name = "OCTAMERIC"
    else:
        print("Sorry Will has only created multigen to work\nfrom monomer to octamer at this time.\nPlease reexecute program.")
        sys.exit()
    
    pdb_file = args[1] 
    print(pdb_file)
    multimerize(pdb_file, mer_num, mer_name)
                
if __name__ == "__main__":
    main(sys.argv[1:])
