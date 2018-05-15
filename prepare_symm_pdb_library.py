import os, sys, requests, string

mer_type = { 1:'MONOMERIC', 2:'DIMERIC', 3:'TRIMERIC', 4:'TETRAMERIC', 5:'PENTAMERIC',\
             6:'HEXAMERIC', 7:'HEPTAMERIC', 8:'OCTAMERIC', 9:'NONAMERIC', 10:'DECAMERIC' }

def convert_res_type(code):
    three_to_one = { 'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q',\
                     'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',\
                     'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W',\
                     'TYR': 'Y', 'VAL': 'V', 'MLY': 'B', 'PTR': 'J', 'TPO': 'O', 'SEP': 'U' }
    one_to_three = {v:k for k,v in three_to_one.iteritems()}
    try:
        return three_to_one[code]
    except KeyError:
        return one_to_three[code]

def parse_pdb_for_biomt(pdb_lines, symmetry_num):
    '''Locates the biomt lines in the REMARK 350 section and outputs a dictionary
    of the chains and transformations needed to be made to make the biological unit'''
    #Each index contains a dictionary of symmetry to apply to a given chain
    biological_units = [] 
    for l_index, line in enumerate(pdb_lines):
        if ["REMARK", "350", "BIOMOLECULE:"] == line.split()[:3]:
            if pdb_lines[l_index+1].split()[-1] != mer_type[symmetry_num]:
                continue
            biomts = {}
            for line2 in pdb_lines[l_index+1:]:
                sline2 = line2.split()
                if ["REMARK", "350", "BIOMOLECULE:"] == sline2[:3] or int(sline2[1]) > 350:
                    biological_units.append({chain_letters: biomts})
                    break
                elif ["APPLY", "THE", "FOLLOWING", "TO", "CHAINS:"] == sline2[2:7]:
                    chain_letters = ''.join(line2.rstrip().split(': ')[-1].split(', '))
                if sline2[2].startswith("BIOMT"):
                    try:
                        biomts[sline2[3]].append([float(x) for x in sline2[4:8]])
                    except KeyError:
                        biomts[sline2[3]] = [[float(x) for x in sline2[4:8]]]
    return biological_units

def split_by_model(pdb_lines):
    '''Split a pdb file into sub-files if more than one model solved'''
    split_models = []
    pdb = []
    for line_index, line in enumerate(pdb_lines):
        if line[0:6] == 'ENDMDL' or line_index == len(pdb_lines)-1:
            split_models.append(pdb[:])
            pdb = []
        else:
            pdb.append(line)
    return split_models
    

def apply_biomt_transform(pdb_line, biomt_block):
    '''Applies the biomt matrix transformation to a pdb_line, formats, and returns new pdb line'''
    x = '{0:>8.3f}'.format(float(pdb_line[30:38]) * biomt_block[0][0] +\
                           float(pdb_line[38:46]) * biomt_block[0][1] +\
                           float(pdb_line[46:54]) * biomt_block[0][2] + biomt_block[0][3] )

    y = '{0:>8.3f}'.format(float(pdb_line[30:38]) * biomt_block[1][0] +\
                           float(pdb_line[38:46]) * biomt_block[1][1] +\
                           float(pdb_line[46:54]) * biomt_block[1][2] + biomt_block[1][3] )

    z = '{0:>8.3f}'.format(float(pdb_line[30:38]) * biomt_block[2][0] +\
                           float(pdb_line[38:46]) * biomt_block[2][1] +\
                           float(pdb_line[46:54]) * biomt_block[2][2] + biomt_block[2][3] )
    return pdb_line[:30] + x + y + z + pdb_line[54:]

def apply_biomts(pdb_lines, relevant_biomols):
    '''Applies the relevant biological molecule transformation BIOMT lines to the the pdb lines
    with the corresponding chain identification and ouputs all generated symmetric proteins'''
    all_symm_pdbs = []
    for bio_index, biomol in enumerate(relevant_biomols):
        new_pdb = {}
        for chains, biomt_blocks in biomol.iteritems():
            chain_count = 0
            for block_number, biomt in biomt_blocks.iteritems():
                for chain in chains:
                    for line in pdb_lines:
                        if line[0:6] in ["ATOM  ", "HETATM"] and line[21] == chain and str(line[17:20]) != 'HOH':
                            try:
                                new_pdb[chain_count].append(apply_biomt_transform(line, biomt))
                            except KeyError:
                                new_pdb[chain_count] = [apply_biomt_transform(line, biomt)]
                    chain_count +=1
        final_pdb = []
        for num in xrange(len(new_pdb)):
            final_pdb += [x[:21] + string.ascii_uppercase[num] + x[22:] for x in new_pdb[num]]
        all_symm_pdbs.append(final_pdb)
    return all_symm_pdbs

def sequence_alignment(seq1, seq2, match=-1, missmatch=20, gap=1):
    '''Global alignment algorithm will match two symmetric chains well despite
    missing residues (electron density) at the termini'''
    direction_map = {}
    #Create starting matrix where columns are seq1 and rows are seq2
    dynamic_matrix = [range(0,gap*(len(seq1)+1),gap)] + [[x] for x in range(1,gap*(len(seq2)+1),gap)]
    #print dynamic_matrix
    for row_i in xrange(1, len(seq2)+1):
        for column_j in xrange(1, len(seq1)+1):
            if seq2[row_i-1] == seq1[column_j-1]:
                mod = match
            else:
                mod = missmatch
            diagonal = dynamic_matrix[row_i-1][column_j-1] + mod
            up = dynamic_matrix[row_i-1][column_j] + gap
            left = dynamic_matrix[row_i][column_j-1] + gap
            min_val = min(diagonal, up, left)
            dynamic_matrix[row_i].append(min_val)
            if diagonal == min_val:
                direction_map[(row_i,column_j)] = (row_i-1,column_j-1)
            elif up == min_val:
                direction_map[(row_i,column_j)] = (row_i-1,column_j)
            elif left == min_val:
                direction_map[(row_i,column_j)] = (row_i,column_j-1)
    
    current_pos = (row_i,column_j)
    aligned_seq1 = ''
    aligned_seq2 = ''
    while current_pos != (0,0):
        try:
            next_pos = direction_map[current_pos]
        except KeyError:
            if current_pos[0] == 0:
                next_pos = (0, current_pos[1]-1)
            elif current_pos[1] == 0:
                next_pos = (current_pos[0]-1,0)
            else:
                raise KeyError, "current_pos is not a key in direction_map and is not an epsilon"
        if next_pos[0] < current_pos[0] and next_pos[1] < current_pos[1]:
            aligned_seq1 += seq1[current_pos[1]-1]
            aligned_seq2 += seq2[current_pos[0]-1]
        elif next_pos[0] < current_pos[0]:
            aligned_seq1 += '_'
            aligned_seq2 += seq2[current_pos[0]-1]
        elif next_pos[1] < current_pos[1]:
            aligned_seq1 += seq1[current_pos[1]-1]
            aligned_seq2 += '_'

        current_pos = next_pos
    print aligned_seq1[::-1]
    print aligned_seq2[::-1]
    #start with max location and work back with dictionary
    return aligned_seq1[::-1], aligned_seq2[::-1]

def compartmentalize_pdb(pdb_lines):
    #{chain_key: {resi: list_of_atom_lines}}
    atomlines_blocked_by_chain = {}
    sequence_by_chain = {}
    xyz_by_chain = {}
    for line in pdb_lines:
        if line[0:6] == "ATOM  ":
            try:
                atomlines_blocked_by_chain[line[21]][int(line[22:26])].append(line)
                xyz_by_chain[line[21]][int(line[22:26])].append([float(line[30:38]),float(line[38:46]),float(line[46:54])])
            except KeyError:
                try:
                    atomlines_blocked_by_chain[line[21]][int(line[22:26])] = [line]
                    sequence_by_chain[line[21]] += convert_res_type(line[17:20])
                    xyz_by_chain[line[21]][int(line[22:26])] = [[float(line[30:38]),float(line[38:46]),float(line[46:54])]]
                except KeyError:
                    atomlines_blocked_by_chain[line[21]] = {int(line[22:26]): [line]}
                    sequence_by_chain[line[21]] = convert_res_type(line[17:20])
                    xyz_by_chain[line[21]] = {int(line[22:26]): [float(line[30:38]),float(line[38:46]),float(line[46:54])]}

    print 'length of chains', len(atomlines_blocked_by_chain)
    for chain, seq in sequence_by_chain.iteritems():
        print 'sequence:', chain, '\n', seq

    print atomlines_blocked_by_chain['A'][6]
    print xyz_by_chain['A'][6]

    return atomlines_blocked_by_chain, sequence_by_chain, xyz_by_chain

sequence_alignment('ABCDEFGHHIJKLM','BCDHIJL')
sequence_alignment('BCDHIJL','ABCDEFGHHIJKLM')


#Using these lines for testing purposes but will be turned into main at some point
#my_url = "http://www.rcsb.org/pdb/files/4HWG.pdb"
my_url = "http://www.rcsb.org/pdb/files/1KNV.pdb"
r = requests.get(my_url)
myfile = r.content
mylines = [x+'\n' for x in myfile.split('\n')]
#mylines = myfile.split('\n')
#with open('4HWG.pdb', 'r') as myFH:
#mylines = myFH.readlines()

relevant_biomols = parse_pdb_for_biomt(mylines, 4)
print len(relevant_biomols)
for v in relevant_biomols:
    #print 'k:\n', k, '\nv:\n', v
    print 'v:\n', v
pdb_by_model = split_by_model(mylines)
print len(pdb_by_model)
count = 0
for model in pdb_by_model:
    all_symm_pdbs = apply_biomts(model, relevant_biomols)
    for symm_pdb in all_symm_pdbs:
        #I want to maximize the occupancy of the residues in this step
        compartmentalize_pdb(symm_pdb)
        #exit(0)
        with open('../blah_%i.pdb' % count, 'w') as myNewFH:
            myNewFH.writelines(symm_pdb)
        count += 1
        


                
