import os, sys, requests, string

mer_type = { 1:'MONOMERIC', 2:'DIMERIC', 3:'TRIMERIC', 4:'TETRAMERIC', 5:'PENTAMERIC',\
             6:'HEXAMERIC', 7:'HEPTAMERIC', 8:'OCTAMERIC', 9:'NONAMERIC', 10:'DECAMERIC' }

def parse_pdb_for_biomt(pdb_lines, symmetry_num):
    '''Locates the biomt lines in the REMARK 350 section and output a dictionary
    of the chains and transformations needed to be made to make a biological unit'''
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
            for line in mylines:
                if line[0:6] in ["ATOM  ", "HETATM"] and line[21] in chains and str(line[17:20]) != 'HOH':
                    for chain_num, biomt in biomt_blocks.iteritems():
                        try:
                            new_pdb[chain_num].append(apply_biomt_transform(line, biomt))
                        except KeyError:
                            new_pdb[chain_num] = [apply_biomt_transform(line, biomt)]
        final_pdb = []
        for num in xrange(len(new_pdb)):
            final_pdb += [x[:21] + string.ascii_uppercase[num] + x[22:] for x in new_pdb[str(num+1)]]
        all_symm_pdbs.append(final_pdb)
    return all_symm_pdbs


#Using these lines for testing purposes but will be turned into main at some point
with open('4HWG.pdb', 'r') as myFH:
    mylines = myFH.readlines()

relevant_biomols = parse_pdb_for_biomt(mylines, 6)
print len(relevant_biomols)
pdb_by_model = split_by_model(mylines)
print len(pdb_by_model)
count = 0
for model in pdb_by_model:
    all_symm_pdbs = apply_biomts(model, relevant_biomols)
    for symm_pdb in all_symm_pdbs:
        with open('blah_%i.pdb' % count, 'w') as myNewFH:
            myNewFH.writelines(symm_pdb)
        count += 1
        



                
