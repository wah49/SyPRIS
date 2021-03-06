#!/usr/bin/env/ python
from transform import *
import sys, os, math, time, ast
from collections import OrderedDict as OD

torsion_map = { #conical set
                'ARG' : ['CZ', 'NE', 'CD', 'CG', 'CB', 'CA', 'N'],\
                'ASN' : ['OD1', 'CG', 'CB', 'CA', 'N'],\
                'ASP' : ['OD1', 'CG', 'CB', 'CA', 'N'],\
                'CYS' : ['SG' , 'CB', 'CA', 'N'],\
                'GLN' : ['OE1', 'CD', 'CG', 'CB', 'CA', 'N'],\
                'GLU' : ['OE1', 'CD', 'CG', 'CB', 'CA', 'N'],\
                'HIS' : ['ND1', 'CG', 'CB', 'CA', 'N'],\
                'ILE' : ['CD1', 'CG1', 'CB', 'CA', 'N'],\
                'LEU' : ['CD1', 'CG', 'CB', 'CA', 'N'],\
                'LYS' : ['NZ' , 'CE', 'CD', 'CG', 'CB', 'CA', 'N'],\
                'MET' : ['CE', 'SD', 'CG', 'CB', 'CA', 'N'],\
                'PHE' : ['CD1', 'CG', 'CB', 'CA', 'N'],\
                'PRO' : ['CD', 'CG', 'CB', 'CA', 'N'],\
                'SER' : ['OG' , 'CB', 'CA', 'N'],\
                'THR' : ['OG1', 'CB', 'CA', 'N'],\
                'TRP' : ['CD1', 'CG', 'CB', 'CA', 'N'],\
                'TYR' : ['CD1', 'CG', 'CB', 'CA', 'N'],\
                'VAL' : ['CG1', 'CB', 'CA', 'N'],\
                #ncaa set
                'BPY' : ['CD1', 'CG', 'CB', 'CA', 'N'],\
                'HQA' : ['CD1', 'CG', 'CB', 'CA', 'N'] }

metals = ['LI','NA','K','MG','CA','MN','FE','CO','NI','CU','ZN','MO','RU','PD','PT']

def read_dictionary_file(my_file):
    with open(my_file, 'r') as f:
        s = f.read()
        my_dict = ast.literal_eval(s)
    return my_dict

# rcp = read, classify, and prepare
# ligands will be labeled as hetatms
# potential branches labeled as atoms
def rcp_cofactor(cofactor_path):
    cofactor = read_file(cofactor_path)
    cofactor = clean_pdb(cofactor)
    ordered_chains = []
    chain_lines = []
    
    for line in cofactor:
        if line == cofactor[0]:
            last_chain = line[21].upper()
        elif line == cofactor[-1]:    
            chain_lines.append(line)
            ordered_chains.append(chain_lines[:])
            break
        if line[21].upper() == last_chain:
            chain_lines.append(line)
        elif line[21].upper() != last_chain:    
            ordered_chains.append(chain_lines[:])
            chain_lines = []
            chain_lines.append(line)
            last_chain = line[21].upper()
    
    symmetric_chains = {}
    ligand_residues = []
    ligand_atom_names = []
    branch_names_by_chain = {}
    all_branch_resn = set()
    metal_ions = []
    for chain in ordered_chains:
        blocked_chain = block_pdb_by_res(chain)
        residues_with_bb = []
        branches = []
        for res_block in blocked_chain:
            res_block_atom_names = []
            N = False
            CA = False
            CB = False
            for line in res_block:
                res_block_atom_names.append(str(int(line[22:30])).upper() + '_'+\
                                            line[21].upper() + '_' +\
                                            line[17:20].strip().upper() + '_' +\
                                            line[11:17].strip().upper())
                if line[11:17].strip().upper() == 'N':
                    N = True
                if line[11:17].strip().upper() == 'CA':
                    CA = True
                    resname_of_ca = str(int(line[22:30])).upper() + line[17:20].upper()
                if line[11:17].strip().upper() == 'CB':
                    CB = True
            if N == CA == CB == True:
                residues_with_bb.append(res_block[:])
                branches.append(res_block_atom_names[:])
                all_branch_resn.add(resname_of_ca)
            else:
                ligand_residues.append(res_block[:])
                ligand_atom_names += res_block_atom_names[:]
                for atom_full_name in res_block_atom_names:
                    atom_name = atom_full_name.split('_')[-1]
                    if atom_name in metals:
                        metal_ions.append(atom_name)

        #print residues_with_bb
        if residues_with_bb:
            if chain == ordered_chains[0]:
                first_chain_number_of_residues = len(residues_with_bb)
            if len(residues_with_bb) != first_chain_number_of_residues:
                print "The input cofactor available backbones were not considered symmetric"
                sys.exit()
            unblocked_bb =  unblock_pdb_by_res(residues_with_bb[:])
            symmetric_chains[chain[0][21].upper()] = unblocked_bb[:]
        if branches:
            branch_names_by_chain[chain[0][21].upper()] = branches[:]
    
    subunits = len(symmetric_chains.keys())
    unblocked_ligand = unblock_pdb_by_res(ligand_residues)
    ordered_by_chain = OD(sorted(symmetric_chains.items()))
    
    final_cofactor = []
    for k,v in ordered_by_chain.iteritems():
        final_cofactor += ["ATOM  " + x[6:] for x in v]
        #this will essentially put copies of all ligands in each chain to insure symmetry
        final_cofactor += ["HETATM" + x[6:21] + k + x[22:] for x in unblocked_ligand]
    
    final_cofactor_obj = Transform(final_cofactor)
    #add and make ligand and branch all same resn/resi

    branch_names_plus_cofactor_names = {}
    for k,v in branch_names_by_chain.iteritems():
        updated_branches = []
        for branch in v:
            branch_names = branch[:]
            for name in ligand_atom_names:
                name_split = name.split('_')
                name_split[1] = k
                new_name =  '_'.join(name_split)
                branch_names.append(new_name)
        
            updated_branches.append(branch_names[:])
        branch_names_plus_cofactor_names[k] = updated_branches[:]

    return final_cofactor_obj, branch_names_plus_cofactor_names, all_branch_resn, subunits, metal_ions

def get_first_three_chis_of_resnum(obj, resnum):
    obj.get_xyz_by_resnum(resnum, update=True)
    obj.get_xyz_by_chain('A', update=True)
    xyz_by_atom = {}
    atom_map = []
    for i, xyz_name in enumerate(obj.xyz_names):
        split_name = xyz_name.split('_')
        if i == 0:
            if split_name[-2] not in torsion_map:
                return ['X','X','X']
            atom_map = torsion_map[split_name[-2]][:]
        if split_name[-1] in torsion_map[split_name[-2]]:
            xyz_by_atom[split_name[-1]] = obj.xyz_total[i][:]
    chi_set = []
    #print 'atom map', atom_map
    #print 'xyz by atom', xyz_by_atom
    for i in xrange(len(atom_map)-3):
        dihedral = get_dihedral( xyz_by_atom[atom_map[-(i+1)]][:],\
                                 xyz_by_atom[atom_map[-(i+2)]][:],\
                                 xyz_by_atom[atom_map[-(i+3)]][:],\
                                 xyz_by_atom[atom_map[-(i+4)]][:] )
        chi_set.append(math.degrees(dihedral))
    while len(chi_set) < 3:
        chi_set.append('X')
    
    return chi_set


def read_SyPRIS_table(file_path):
    with open(file_path, 'r') as f:
        my_file = f.readlines()
        my_file = [x.strip('\r\n').split(',') for x in my_file]
    pdbs = [x[0] for x in my_file]
    residues = [x[1] for x in my_file]
    rigid_rots = [x[2] for x in my_file]
    rmsds = [x[3] for x in my_file]
    anglelogs = [x[4] for x in my_file]
    chis = [x[5:] for x in my_file]
    
    return pdbs, residues, rigid_rots, rmsds, anglelogs, chis

def average_dict_lists(my_dict):
    avg_rot_list = {}
    for key, vals in my_dict.iteritems():
        length = len(vals)
        avg_rot_list[key] = [round(sum(x)/length,4) for x in zip(*vals)]
    return avg_rot_list
                   
def dunbrack_rotLib_reader(my_file):
    if type(my_file) != list:
        my_file = read_file(my_file)
    rot_dict = {}
    type_rot = {}
    append_items = False
    cont_vals = [-121.,-61.,-1.0,59.,119.,179.]
    for line in my_file:
        column_sep = line.split()
        if line == my_file[0]:
            last_resType = column_sep[0]
        elif line == my_file[-1]:
            append_items = True
            
        if column_sep[0] == last_resType:
            if last_resType in ['ASP', 'ASN']:
                rot_ids = []
                rep_index = 1
                for i, cont_val in enumerate(cont_vals):
                    rot_ids.append(int(column_sep[4]+\
                                       str(i)+\
                                       column_sep[6]+\
                                       column_sep[7]))
            elif last_resType in ['GLU', 'GLN']:
                rot_ids = []
                rep_index = 2
                for i, cont_val in enumerate(cont_vals):
                    rot_ids.append(int(column_sep[4]+\
                                       column_sep[5]+\
                                       str(i)+\
                                       column_sep[7]))
            else:    
                rot_ids = [int(column_sep[4]+\
                             column_sep[5]+\
                             column_sep[6]+\
                             column_sep[7])]
            for rot_ind, rot_id in enumerate(rot_ids):
                if rot_id not in type_rot.keys():
                    corrected_vals = []
                    for v_ind, value in enumerate(column_sep[9:]):
                        value = float(value)
                        if value < 0.0:
                            value = 360. + value
                        if len(rot_ids) > 1:
                            if v_ind == rep_index:
                                corrected_vals.append(cont_vals[rot_ind])
                            elif v_ind == rep_index+4:
                                corrected_vals.append(15.0)
                            else:
                                corrected_vals.append(value)
                        else:
                            corrected_vals.append(value)
                    type_rot[rot_id] = [corrected_vals[:]]

                else:
                    current_data = type_rot[rot_id][:]
                    corrected_vals = []
                    for v_ind, value in enumerate(column_sep[9:]):
                        value = float(value)
                        if value < 0.0:
                            value = 360. + value
                        if len(rot_ids) > 1:
                            if v_ind == rep_index:
                                corrected_vals.append(cont_vals[rot_ind])
                            elif v_ind == rep_index+4:
                                corrected_vals.append(15.0)
                            else:
                                corrected_vals.append(value)
                        else:
                            corrected_vals.append(value)
                    current_data.append(corrected_vals)
                    type_rot[rot_id] = current_data[:]
        
        else:
            append_items = True
        
        if append_items == True:
            final_old_rot_avgs = average_dict_lists(type_rot)
            final_rot_avgs = {}
            for rot_key, rot_vals in final_old_rot_avgs.iteritems():
                recorrected_rot_vals = []
                for val in rot_vals:
                    if val > 180.0:
                        val = round((val - 360.),4)
                    recorrected_rot_vals.append(val)
                final_rot_avgs[rot_key] = recorrected_rot_vals
            
            rot_dict[last_resType] = final_rot_avgs
            type_rot = {}
            last_resType = column_sep[0]
            append_items = False
    
    return rot_dict

def rotate_grid_about_axis(grid, axis, theta, grid_size):
    rotate_vector = axisangle_to_q(axis, theta)
    grid_list = convert_grid_string_to_list(grid)

    rotated_grid = []
    for grid_point in grid_list:
        mag = get_mag(grid_point)

        Rxyz = q_xyz_mult(rotate_vector, grid_point)
        Rxyz = ((Rxyz[0] * mag),\
                (Rxyz[1] * mag),\
                (Rxyz[2] * mag))
        
        sRxyz = str([int(round(x/math.ceil(grid_size), 0)*math.ceil(grid_size)) for x in Rxyz])
        rotated_grid.append(sRxyz)
    return rotated_grid
    

def angle_avg(test, scaffold, ratio=(math.pi/6.0)):
    def zero_vector_pair(v_pair):
        trans = [0. - v_pair[0][0], 0. - v_pair[0][1], 0. - v_pair[0][2]]
        new_pair = []
        for point in v_pair:
            new_point = add_vector_difference(point, trans)
            new_pair.append(new_point)
        return new_pair

    compare_vectors = zip(test, scaffold)
    vector_ang_sum = 0.0
    for vector_pairs in compare_vectors:
        vector1 = zero_vector_pair(vector_pairs[0])
        vector2 = zero_vector_pair(vector_pairs[1])
        vector_ang_sum += vector_angle(vector1[1], vector2[1])
    ang_avg = ((vector_ang_sum/ratio)**2 / len(compare_vectors))**0.5
    return ang_avg
#for angle_log to work it requires a list of vectors
#example: [ [[1,2,3],[4,5,6]], [[],[]] ]
#where paired atom coords would produce a vector
#from atom1 to atom2 threshold default is 20 degrees, 
#and must be supplied as radians
def angleLog(test, scaffold, threshold=(math.pi/6.0)):
    def zero_vector_pair(v_pair):
        trans = [0. - v_pair[0][0], 0. - v_pair[0][1], 0. - v_pair[0][2]]
        new_pair = []
        for point in v_pair:
            new_point = add_vector_difference(point, trans)
            new_pair.append(new_point)
        return new_pair

    compare_vectors = zip(test, scaffold)
    vector_ang_sum = 0.00001
    for vector_pairs in compare_vectors:
        vector1 = zero_vector_pair(vector_pairs[0])
        vector2 = zero_vector_pair(vector_pairs[1])
        vector_ang_sum += vector_angle(vector1[1], vector2[1])
    angle_log = math.log10((vector_ang_sum) / (len(compare_vectors) * threshold))
    return angle_log

class SyPRISCofactor(object):
    
    def __init__(self, cofactor_path):
        if type(cofactor_path) == str:
            #self.cofactor_file = read_file(cofactor_path)
            self.cofactor, \
            self.branch_names_by_chain, \
            self.branch_residue_names, \
            self.subunits, \
            self.metal_ions = rcp_cofactor(cofactor_path)
        else:
            print "SyPRISCofactor arguement is not of type str\n \
                   #and thus not a valid path name."
            sys.exit()
        
        self.zero_tor_center = []
        self.branch_connectivity_map = {}
        self.downstream_atom_connections = {}
        self.branch_objs = {}
        self.bb_rot_ensemble = {}
        self.cof_grid_size = 1.0
        self.bb_grid_ensemble = {}
        self.grid_center = []
        self.lasa = ''
        self.lasa_geo_centers = {}
        self.lasa_grid_centers = {}
        self.lasa_vectors = {}
        self.lasa_radii = {}
        self.frozen_atom_names = []
        self.dunbrack_rotamers = {}
        self.use_input_rotamer = False
        self.original_rotamers = {}
        self.lock_rotamer_addition = False
        self.std_mult = 1.0
        self.match_cofactors = {}
        self.match_data_set = {}

    #Need to give this a rot lib made by dunbrack_rotLib_reader
    def read_in_rotlib(self, my_file):
        with open(my_file, 'r') as f:
            s = f.read()
            my_dict = ast.literal_eval(s)
            if type(my_dict) == dict:
                self.dunbrack_rotamers = my_dict
            else:
                self.dunbrack_rotamers = dunbrack_rotLib_reader(my_file)
        if self.use_input_rotamer:
            if self.original_rotamers:
                for k,v in self.original_rotamers.iteritems():
                    restype = k.split('_')[-1]
                    for i,tors in enumerate(v):
                        self.dunbrack_rotamers[restype].update({'ex%s' % i : tors })
                #print "The read_in_rotlib did not evaluate a dictionary"
                #sys.exit()

    #generate atom map
    #not sure if you would ever want to do this several times but it is distance based
    #if something was moved and atoms are in the incorrect place it will be incorrect
    def generate_branch_connectivity_map(self):
        self.branch_connectivity_map = {}
        for chain_name, branch_set in self.branch_names_by_chain.iteritems():
            #branch_map = {}
            for branch in branch_set:
                #atom_connection_dict = {}
                for atom_name1 in branch:
                    #if atom_name1 == branch[0]:
                    #branch_name = '_'.join([x for i,x in enumerate(atom_name1.split('_')) if i != 3])
                    atom1 = self.cofactor.get_xyz_from_name_list(atom_name1)[0]
                    connected_atoms = []
                    for atom_name2 in branch:
                        if atom_name1 != atom_name2:
                            atom2 = self.cofactor.get_xyz_from_name_list(atom_name2)[0]
                            if get_mag_difference(atom1, atom2) < 2.0:
                                connected_atoms.append(atom_name2)
                    if connected_atoms:
                        #atom_connection_dict[atom_name1] = connected_atoms
                        self.branch_connectivity_map[atom_name1] = connected_atoms[:]
                #branch_map[branch_name] = atom_connection_dict.copy()
            #self.branch_connectivity_map[chain_name] = branch_map.copy()
                                
    def get_downstream_connected_atoms(self, starting_atoms):
        if type(starting_atoms) != list:
            if type(starting_atoms) == str:
                name_precursor = starting_atoms
                starting_atoms = [starting_atoms]
            else:
                starting_atoms = list(starting_atoms)

        included_atom_names = []
        old_length = -1
        while len(included_atom_names) > old_length:
            if len(included_atom_names) == 0:
                atoms_to_check = starting_atoms[:]
            else:
                atoms_to_check = included_atom_names[:]
                old_length = len(atoms_to_check)
            for atom_precursor in atoms_to_check:
                for connected_atom in self.branch_connectivity_map[atom_precursor]:
                    if connected_atom not in included_atom_names and \
                       connected_atom not in self.frozen_atom_names:
                        included_atom_names.append(connected_atom)
        
        self.downstream_atom_connections[name_precursor] = included_atom_names[:]
        #sys.exit()
        #return included_atom_names
                       #connected_atom not in static_tor_atoms[:-1] and \
        

    #if looking to work on a specific chain, input chain letter
    def generate_branch_objects(self, chain=[]):
        self.branch_objs = {}
        for k,v in self.branch_names_by_chain.iteritems():
            if k in chain:
                for branch in v:
                    branch_obj = copy_transform_object(self.cofactor)
                    branch_obj.get_xyz_from_name_list(branch[:], True)
                    branch_name = '_'.join([x for i,x in enumerate(branch[0].split('_')) if i != 3])
                    self.branch_objs[branch_name] = copy_transform_object(branch_obj)

    #freeze_atoms
    def freeze_atoms(self, ligand=True, extra=[]):
        if ligand == True:
            for i, line in enumerate(self.cofactor.pdb[:]):
                if line.startswith('HETATM'):
                    self.frozen_atom_names.append(self.cofactor.xyz_names[i])
        if extra:
            for line in extra:
                if line in self.cofactor.xyz_names:
                    self.frozen_atom_names.append(line)

    def apply_tor_to_branch(self, branch_obj='', branch_name='', torsions=0):
        if torsions != 0:
            #we multiply by -1 because we are rotating in reverse
            rotamer_set = [math.radians(-1.0*float(x)) for x in torsions[0:4] if float(x) != 0.0]
        if not self.branch_connectivity_map:
            self.generate_branch_connectivity_map()

        if type(branch_obj) == str:
            if not branch_name:
                branch_name = branch_obj
            branch_restype = branch_name.split('_')[2]
            branch_obj = self.branch_objs[branch_name]
        else:
            branch_restype = branch_name.split('_')[2]

        if branch_restype in torsion_map.keys():
            tor_set = torsion_map[branch_restype][:]
        else:
            print 'tor_set unkown from given cofactor'
            sys.exit()
        
        number_of_chis = len(tor_set)-3
        if torsions != 0:
            print 'rotamer set of torsions: ', rotamer_set
            if number_of_chis != len(rotamer_set):
                print '#ofchis', number_of_chis
                print 'branch_restype:', branch_restype
                print 'tor_set', tor_set
                print 'rotamer set',rotamer_set
                tor_set.insert(0, self.metal_ions[0])
                number_of_chis += 1
                print "rotamer set provided is not the same size as number of chis to sample"
        saved_dihedrals = []
        for i in xrange(number_of_chis):
            #the first 3 atoms do not move, we will later define the movement
            #as all atoms connected to the last of these 3 (LASA)
            pre_static_tor_atoms = [tor_set[i], tor_set[i+1], tor_set[i+2], tor_set[i+3]]
            #static_tor_atoms = [ ''.join([branch_name, '_', x]) for x in static_tor_atoms]
            static_tor_atoms = []
            for tor in pre_static_tor_atoms:
                if tor not in self.metal_ions:
                    static_tor_atoms.append(''.join([branch_name, '_', tor]))
                else:
                    chain = branch_name.split('_')[1]
                    metal_copy = copy_transform_object(branch_obj)
                    metal_copy.get_xyz_by_atomtype(self.metal_ions, update=True)
                    metal_copy.get_xyz_by_chain(chain, update=True)
                    static_tor_atoms.append(metal_copy.xyz_names[0])
                    
            lasa = branch_name + '_' + tor_set[i+2]
            static_tor_xyz = branch_obj.get_xyz_from_name_list(static_tor_atoms)
            self.freeze_atoms(True, static_tor_atoms[:-1])
            if lasa not in self.downstream_atom_connections.keys():
                self.get_downstream_connected_atoms(lasa)
            included_atom_names = self.downstream_atom_connections[lasa][:]
            self.frozen_atom_names = []
            #self.unfreeze_atoms( static_tor_atoms[:-1] )
            
            #sys.exit()
            truncated_branch = copy_transform_object(branch_obj)
            truncated_branch.get_xyz_from_name_list(included_atom_names, True)
            #translate it be tor_set[i+1]
            #truncated_branch.write_to_file('first')
            truncated_branch.translate([0.,0.,0.], static_tor_xyz[1])
            if torsions == 0:
                print 'here'
                print static_tor_atoms
                print static_tor_xyz
                dihedral_angle = get_dihedral(static_tor_xyz[0],\
                                              static_tor_xyz[1],\
                                              static_tor_xyz[2],\
                                              static_tor_xyz[3])
                if self.use_input_rotamer:
                    if not self.lock_rotamer_addition:
                        saved_dihedrals.insert(0, math.degrees(dihedral_angle))
                truncated_branch.rotate_about_axis(get_vector_difference(static_tor_xyz[2],static_tor_xyz[1]),\
                                                   dihedral_angle, True)
            else:
                #need to apply dunbrack in reverse
                dihedral_angle = rotamer_set[-(i+1)]
                print 'dihedral angle is:', dihedral_angle
                truncated_branch.rotate_about_axis(get_vector_difference(static_tor_xyz[2],static_tor_xyz[1]),\
                                                   dihedral_angle, True)
            
            truncated_branch.translate(static_tor_xyz[1], [0.,0.,0.])
            branch_obj.update_coordinates_from_truncated_obj(copy_transform_object(truncated_branch))
            if type(branch_obj) == list:
                self.branch_objs[branch_name] = copy_transform_object(branch_obj)
        
        if self.use_input_rotamer:
            if not self.lock_rotamer_addition:
                orig_tor = [0.,0.,0.,0.]
                orig_std = [0.,0.,0.,0.]
                for i in xrange(len(saved_dihedrals)):
                    orig_tor[i] += saved_dihedrals[i]
                    orig_std[i] += 10.0
                orig_rotamer = orig_tor + orig_std
                if branch_name in self.original_rotamers.keys():
                    current_rots = self.original_rotamers[branch_name][:]
                    current_rots.append(orig_rotamer)
                    self.original_rotamers[branch_name] = current_rots[:]
                else:
                    self.original_rotamers[branch_name] = [orig_rotamer]
                    
        return copy_transform_object(branch_obj)
    

    def generate_cof_center(self):
        branch_restype_in_cof = set()
        for branch_name in self.branch_objs.keys():
            branch_restype_in_cof.add(branch_name.split('_')[-1])
        list(branch_restype_in_cof)
        copy_of_cof = copy_transform_object(self.cofactor)
        copy_of_cof.get_xyz_by_restype(branch_restype_in_cof, update=True)
        self.zero_tor_center = list(copy_of_cof.get_geo_center())

    #this definition should go through and zero all torsion angles
    def zero_branch_torsion_angles(self):
        for chain_name, branch_set in self.branch_names_by_chain.iteritems():
            self.generate_branch_objects(chain_name)
            for branch_name in self.branch_objs.keys():
                print 'branch name', branch_name
                self.apply_tor_to_branch(branch_obj=branch_name)
                branch_obj = self.branch_objs[branch_name]
                self.cofactor.update_coordinates_from_truncated_obj(branch_obj)
        self.generate_cof_center()

    def generate_inverse_rot_ensemble(self, grid_size=1.0, std_periodicity=0, lasa='CB',\
                                      my_rotlib='~/Rosetta/main/database/rotamer/bbdep02.May.sortlib-correct.12.2010'):
        #for every branch we want to apply the normal rotamers from a read-in library
        if not self.lasa:
            self.lasa = lasa
        def apply_std(std, torsion_angles, start_ind):
            additional_torsions = {}
            print 'torsion_angles:', torsion_angles
            for ind, value in enumerate(std):
                if ind < start_ind:
                    continue
                for term in xrange(std_periodicity):
                    print 'term', term
                    plus_torsion = torsion_angles[:]
                    plus_torsion[ind] += (value * (term+1)/(std_periodicity))
                    additional_torsions['a_%s_%s' % (ind, term)] = plus_torsion[:]
                    
                    minus_torsion = torsion_angles[:]
                    minus_torsion[ind] -= (value * (term+1)/(std_periodicity))
                    additional_torsions['b_%s_%s' % (ind, term)] = minus_torsion[:]
                
            return additional_torsions

        self.cof_grid_size = grid_size
        additional_rotamers = {}
        if not self.dunbrack_rotamers:
            self.read_in_rotlib(my_rotlib)
        for residue_name in self.branch_residue_names:
            print residue_name
            restype = residue_name[-3:]
            if restype in self.dunbrack_rotamers.keys():
                rotamers = self.dunbrack_rotamers[restype]
            else:
                print "Rotamer library supplied did not have:", restype
                sys.exit()
            
            for rot_type, torsion_angles in rotamers.iteritems():
                if std_periodicity != 0:
                    std_values = [x*self.std_mult for x in torsion_angles[4:] if float(x) != 0.0]
                    first_name = residue_name + '_' + str(rot_type)
                    print residue_name
                    start_ind = len(torsion_map[residue_name[-3:]]) - (torsion_map[residue_name[-3:]].index(self.lasa) +1)
                    print start_ind
                    dict_of_torsion_angles = {first_name : torsion_angles[:]}
                    add_tors = apply_std(std_values[:], torsion_angles[:], start_ind)
                    
                    for type_name, tor in add_tors.iteritems():
                        split_tname = type_name.split('_')
                        my_type = split_tname[0] + split_tname[2]
                        index_in_type_name = int(split_tname[1])
                        list_of_rot_type = list(str(rot_type))
                        list_of_rot_type[index_in_type_name] += my_type 
                        name = residue_name + '_' + ''.join(list_of_rot_type)
                        dict_of_torsion_angles[name] = tor[:]
                        if residue_name not in additional_rotamers.keys():
                            additional_rotamers[residue_name] = {''.join(list_of_rot_type) : tor[:] }
                        else:
                            additional_rotamers[residue_name][''.join(list_of_rot_type)] = tor[:]
                else:
                    name = residue_name + '_' + str(rot_type)
                    dict_of_torsion_angles = {name : torsion_angles[:]}
                
                for name, torsion_angle_set in dict_of_torsion_angles.iteritems():
                    for chain_name, branch_set in self.branch_names_by_chain.iteritems():
                        self.generate_branch_objects(chain_name)
                        for branchtype in self.branch_objs.keys():
                            split_name = branchtype.split('_')
                            split_name.pop(1)
                            if ''.join(split_name) != residue_name:
                                continue
                            self.apply_tor_to_branch(branchtype, '', torsion_angle_set[:])
                            branch_obj = self.branch_objs[branchtype]
                            self.cofactor.update_coordinates_from_truncated_obj(branch_obj)
                            
                        
                    self.bb_grid_ensemble[name] = { 'N' : self.cofactor.get_atom_grid(grid_size, 'N'),\
                                                    'CA' : self.cofactor.get_atom_grid(grid_size, 'CA'),\
                                                    'CB' : self.cofactor.get_atom_grid(grid_size, 'CB') }
                    self.bb_rot_ensemble[name] = copy_transform_object(self.cofactor)
                    
                    lasa_xyz_total = self.cofactor.get_xyz_by_atomtype(lasa)
                    atom_sets_by_name = {}
                    for lasa_xyz in lasa_xyz_total:
                        atom_name_list = self.cofactor.get_name_from_xyz(lasa_xyz).split('_')
                        atom_name_list.pop(1)
                        if name.split('_')[0][-3:] == atom_name_list[1]:
                            atom_name = '_'.join(atom_name_list)
                            if atom_name not in atom_sets_by_name.keys():
                                atom_sets_by_name[atom_name] = [lasa_xyz]
                            else:
                                current_points = atom_sets_by_name[atom_name]
                                current_points.append(lasa_xyz)
                                atom_sets_by_name[atom_name] = current_points[:]

                    for a_name, val_set in atom_sets_by_name.iteritems(): 
                        lasa_center = get_geo_center(val_set)
                        lasa_grid_center = gridify_xyz(lasa_center, grid_size)
                        self.lasa_geo_centers[name] = lasa_center[:]
                        vector = get_vector_difference(val_set[0], lasa_center)
                        radius = get_mag(vector)
                        self.lasa_vectors[name] = vector
                        self.lasa_radii[name] = radius

                    self.grid_center = gridify_xyz(self.zero_tor_center[:], grid_size)
                    self.zero_branch_torsion_angles()
        if additional_rotamers:
            for res, rot_set in additional_rotamers.iteritems():
                for rot_type, torsion_angles in rot_set.iteritems():
                    self.dunbrack_rotamers[res[-3:]][rot_type] = torsion_angles[:]

    
    def separate_by_branch(self, rotamer_res_name):
        copy_of_cof = copy_transform_object(self.cofactor)
        print copy_of_cof.xyz_names
        atom_namelist = []
        for cof_name in copy_of_cof.xyz_names:
            if cof_name.split('_')[0] == rotamer_res_name[:-3]:
                atom_namelist.append(cof_name)
            if cof_name.split('_')[-1] in self.metal_ions\
            and cof_name.split('_')[-2] in self.metal_ions[-3:]:
                atom_namelist.append(cof_name)
        #copy_of_cof.get_xyz_by_resnum(rotamer_res_name[:-3], update=True)
        copy_of_cof.get_xyz_from_name_list(atom_namelist[:], update=True)
        branch_namelists = []
        same_name = []
        print 'rot_res_name', rotamer_res_name
        for name_type in copy_of_cof.xyz_names:
            print name_type
            name_list = name_type.split('_')
            new_name_type = '_'.join(name_list[:-1])
            if name_type == copy_of_cof.xyz_names[0]:
                previous_name = new_name_type
            elif name_type == copy_of_cof.xyz_names[-1]:
                same_name.append(name_type)
                branch_namelists.append(same_name)
                break 
            if new_name_type == previous_name:
                same_name.append(name_type)
            else:
                print 'split_name: ', name_type.split('_')
                print '\n'
                print self.metal_ions
                if new_name_type.split('_')[-1] in self.metal_ions:
                    same_name.append(name_type)
                else:
                    branch_namelists.append(same_name)
                    same_name = []
                    previous_name = new_name_type
                    same_name.append(name_type)

        print 'branch_namelists', branch_namelists
        branches = {}
        for branch_namelist in branch_namelists:
            branch_name = '_'.join(branch_namelist[0].split('_')[:-1])
            branch_obj = copy_transform_object(copy_of_cof)
            branch_obj.get_xyz_from_name_list(branch_namelist[:], update=True)
            branches[branch_name] = copy_transform_object(branch_obj)
        
        return branches

    def sample_std_deviations(self, scaffold, rotamer_res_name, scaff_atom_set, rotamer_torsion):
        rmsd_sum = 0.0
        vec_sqr_sum = 0.0
        sum_of_chis = {}
        print 'rot_name', rotamer_res_name
        cof_restype = rotamer_res_name[-3:]
        print cof_restype
        print 'rot_tor', rotamer_torsion
        std_values = [x * self.std_mult for x in rotamer_torsion[4:]]
        print std_values
        #cofactor is already setup
        operating_branches = self.separate_by_branch(rotamer_res_name)
        
        
        for branch_name, branch_obj in operating_branches.iteritems():
        
            atom_xyz_by_name = branch_obj.pair_name_and_xyz()
            print 'atom_xyz_by_name', atom_xyz_by_name
            metal_xyz_by_name = { x.split('_')[-1] : atom_xyz_by_name[x]\
                                      for x in atom_xyz_by_name.keys() if x.split('_')[-2] != cof_restype }
            atom_xyz_by_name = { x.split('_')[-1] : atom_xyz_by_name[x]\
                                      for x in atom_xyz_by_name.keys() if x.split('_')[-2] == cof_restype }
            print 'atom_xyz_by_name', atom_xyz_by_name
            print '\n\n'
            print 'metal_xyz_by_name', metal_xyz_by_name
            for atom_name, atom_xyz in atom_xyz_by_name.iteritems():
                split_atom_name = atom_name.split('_')
                print split_atom_name
                if split_atom_name[-1] == self.lasa:
                    for scaff_atom_xyz in scaff_atom_set:
                        lasa_distance = get_mag_difference(scaff_atom_xyz, atom_xyz)
                        if scaff_atom_xyz == scaff_atom_set[0]:
                            previous_length = lasa_distance
                            match_scaff_xyz = scaff_atom_xyz
                        else:
                            if lasa_distance < previous_length:
                                previous_length = lasa_distance
                                match_scaff_xyz = scaff_atom_xyz
            matched_chain = scaffold.get_name_from_xyz(match_scaff_xyz).split('_')[1]
            scaff_bb_obj = copy_transform_object(scaffold)
            scaff_bb_obj.get_xyz_by_chain(matched_chain, update=True)
            #scaff_bb_obj.get_xyz_by_resnum(  , update=True)
            scaff_atom_xyz_by_name = scaff_bb_obj.pair_name_and_xyz()
            scaff_atom_xyz_by_name = { x.split('_')[-1] : scaff_atom_xyz_by_name[x][:]\
                                      for x in scaff_atom_xyz_by_name.keys() }
            
            #do a map check first
            tor_atomtypes = torsion_map[rotamer_res_name[-3:]]
            for index, atomtype in enumerate(tor_atomtypes):
                if atomtype == self.lasa:
                    tor_lasa_ind = index
                if atomtype == 'N':
                    #this value is one less than number_of_chis
                    number_of_chis = index - tor_lasa_ind
            
            upstream_lasa = tor_atomtypes[:tor_lasa_ind+1]
            print 'upstream lasa', upstream_lasa
            downstream_lasa = tor_atomtypes[tor_lasa_ind:]
            print '\n', downstream_lasa, '\n'
            if tor_lasa_ind < 2:
                upstream_lasa.insert(0, self.metal_ions[0])
                if tor_lasa_ind == 0:
                    upstream_lasa.insert(0, 'zero')

            for chi_num in xrange(number_of_chis):
                if upstream_lasa[chi_num] == 'zero':
                    prefinal_atom = metal_xyz_by_name[self.metal_ions[0]][:]
                    final_anchor_atom = [0.0,0.0,0.0,]
                elif upstream_lasa[chi_num] in self.metal_ions:
                    prefinal_atom = atom_xyz_by_name[upstream_lasa[chi_num+1]][:]
                    final_anchor_atom = metal_xyz_by_name[self.metal_ions[0]][:]
                    rotamer_torsion[1] += 0.000000000000001
                else:
                    prefinal_atom = atom_xyz_by_name[upstream_lasa[(chi_num+1)]][:]
                    final_anchor_atom = atom_xyz_by_name[upstream_lasa[chi_num]][:]

                print final_anchor_atom
                print prefinal_atom
                print atom_xyz_by_name[downstream_lasa[chi_num]]
                print atom_xyz_by_name[downstream_lasa[chi_num+1]]
                print scaff_atom_xyz_by_name[downstream_lasa[chi_num]]
                print scaff_atom_xyz_by_name[downstream_lasa[chi_num+1]]
                print 'there'
                first_chi = get_dihedral( final_anchor_atom[:],\
                                          prefinal_atom[:],\
                                          atom_xyz_by_name[downstream_lasa[chi_num]][:],\
                                          atom_xyz_by_name[downstream_lasa[chi_num+1]][:] )
                print 'everywhere'
                new_chi = get_dihedral( final_anchor_atom[:],\
                                        prefinal_atom[:],\
                                        scaff_atom_xyz_by_name[downstream_lasa[chi_num]][:],\
                                        scaff_atom_xyz_by_name[downstream_lasa[chi_num+1]][:] )
                #with this chi apply at approp. position
                first_chi = math.degrees(first_chi)
                new_chi = math.degrees(new_chi)
                if first_chi < 0.0:
                    first_chi += 360.0
                if new_chi < 0.0:
                    new_chi += 360.0
                
                print 'first_chi degrees', first_chi
                print 'new_chi degrees', new_chi
                tor_values = rotamer_torsion[:4]
                print tor_values
                delta_chi = (new_chi - first_chi)
                if abs(delta_chi) > 180.0:
                    if first_chi > 180.:
                        new_chi += 360.
                    else:
                        first_chi += 360.
                    delta_chi = (new_chi - first_chi)

                #check here to see std
                std_for_this_chi_num = std_values[number_of_chis - (chi_num+1)]
                print delta_chi
                print std_for_this_chi_num
                if std_for_this_chi_num == 0.0:
                    additional_rotation = delta_chi
                else:
                    if abs(delta_chi) > abs(std_for_this_chi_num):
                        additional_rotation = std_for_this_chi_num * (delta_chi/abs(delta_chi))
                    else:
                        additional_rotation = delta_chi
                print additional_rotation
                if chi_num in sum_of_chis.keys():
                    sum_of_chis[chi_num] += (additional_rotation + tor_values[number_of_chis - (chi_num+1)])
                else:
                    sum_of_chis[chi_num] = (additional_rotation + tor_values[number_of_chis - (chi_num+1)])

                new_torsion_to_apply = [0.,0.,0.,0.,0.,0.,0.,0.]
                non_zero_tor = [x for x in tor_values if x != 0]
                for i in xrange(len(non_zero_tor)):
                    new_torsion_to_apply[i] += 0.000000000000001

                new_torsion_to_apply[number_of_chis - (chi_num+1)] += additional_rotation
                print 'new tor to apply', new_torsion_to_apply
                print 'number of chis', number_of_chis
                print 'additional rotation', additional_rotation
                branch_obj = self.apply_tor_to_branch(branch_obj, branch_name, new_torsion_to_apply[:])


            self.cofactor.update_coordinates_from_truncated_obj(copy_transform_object(branch_obj))

            print 'atom_xyz before:', atom_xyz_by_name
            updated_atom_xyz_by_name = branch_obj.pair_name_and_xyz()
            print '\nfirst_updated_atom_xyz_by_name:', updated_atom_xyz_by_name
            updated_atom_xyz_by_name = { x.split('_')[-1] : updated_atom_xyz_by_name[x][:]\
                                      for x in updated_atom_xyz_by_name.keys() if x.split('_')[-2] == rotamer_res_name[-3:]}
            print '\nupdated_atom_xyz_by_name:', updated_atom_xyz_by_name
            print '\nscaff_atom_xyz_by_name:', scaff_atom_xyz_by_name
            #determine rmsd and anglelog
            #sqrt of the sum of all distances for like atoms
            mag_sum = 0.0
            branch_two_point_vectors = []
            scaff_bb_two_point_vectors = []
            for atom_index, atomtype in enumerate(downstream_lasa):
                print 'atoms that will be for anglelog and rmsd', atomtype
                mag_sum += get_mag_difference(updated_atom_xyz_by_name[atomtype][:],
                                              scaff_atom_xyz_by_name[atomtype][:])**2
                if atomtype != 'N':
                    branch_two_point_vectors.append([updated_atom_xyz_by_name[atomtype][:],\
                                                     updated_atom_xyz_by_name[downstream_lasa[atom_index+1]][:]])
                    scaff_bb_two_point_vectors.append([scaff_atom_xyz_by_name[atomtype][:],\
                                                       scaff_atom_xyz_by_name[downstream_lasa[atom_index+1]][:]])
            rmsd = (mag_sum/len(downstream_lasa))**0.5
            print 'indv rmsd:', rmsd
            rmsd_sum += rmsd
            
            #anglelog just requires the xyzs of the branch starting from lasa
            vec_ang = angle_avg(branch_two_point_vectors[:], scaff_bb_two_point_vectors[:])
            vec_sqr_sum += vec_ang
            #anglelog = angleLog(branch_two_point_vectors[:], scaff_bb_two_point_vectors[:])
            #print 'indv anglelog:', anglelog

            #anglelog_sum += anglelog
        
        print 'rmsd:', rmsd_sum
        print 'vec_sqr_sum:', vec_sqr_sum
        #avg the chis then swap old chis with new chis in reverse order
        avg_chis = {k : v/self.subunits for k,v in sum_of_chis.iteritems()}
        chi_position = 0
        for i in xrange(len(avg_chis)-1, -1, -1):
            if avg_chis[i] > 180.0:
                avg_chis[i] -= 360.0
            elif avg_chis[i] < -180.0:
                avg_chis[i] += 360
            output_torsion = rotamer_torsion[:]
            output_torsion[chi_position] = avg_chis[i]
            chi_position += 1

        rmsd_out = rmsd_sum/self.subunits
        vec_sqr_out = vec_sqr_sum/self.subunits
        score_out = (rmsd_out + vec_sqr_out)/2.0

        return { 'rmsd' : rmsd_out,\
                 'rmsa' : vec_sqr_out,\
                 'score' : score_out,\
                 'torsions' : output_torsion[:4] }


class CsymmCofactor(SyPRISCofactor):
    
    def __init__(self, cofactor_path, scaffold_path):
        if type(cofactor_path) == str:
            self.cofactor_file_path = cofactor_path
            self.scaffold_file_path = scaffold_path
            self.scaffold = Transform(clean_pdb(read_file(scaffold_path)))
            #self.branch_connectivity_map = {}
            #self.branch_objs = {}
            #self.bb_grid_ensemble = {}
            #self.frozen_atom_names = []
            super(CsymmCofactor, self).__init__(self.cofactor_file_path)

        else:
            print "SyPRISCofactor arguement is not of type str\n \
                   #and thus not a valid path name."
            sys.exit()
        self.reduced_scaffold = []
        self.scaff_axes = {}
        self.scaff_lasa_geo_centers = {}
        self.scaff_lasa_grid_centers = {}
        self.scaff_lasa_vectors = {}
        self.scaff_lasa_radii = {}
        self.scaff_lasa_atom_sets = {}
        self.scaff_CA_names = set()
        self.scaff_bb_grid = {}
        self.axis_centers = set()
        #self.rigid_body_xyz_ensemble = {}
        self.coarse_grid_matches = {}
    
    def get_scaff_lasa_info(self, lasa='', radius_max=0):
        if not lasa:
            if not self.lasa:
                lasa = 'CB'
            else:
                lasa = self.lasa
        lasa_xyz_total = self.scaffold.get_xyz_by_atomtype(lasa)
        broken_into_chains = zip(*(iter(lasa_xyz_total),)*(len(lasa_xyz_total)/self.subunits))
        atom_sets_by_chain = zip(*broken_into_chains)
        residues_in_range = []
        for atom_set in atom_sets_by_chain:
            lasa_center = get_geo_center(atom_set)
            lasa_xyz = atom_set[0] 
            pre_name = self.scaffold.get_name_from_xyz(lasa_xyz)
            name = '_'.join([pre_name.split('_')[0], pre_name.split('_')[2], lasa])
            vector = get_vector_difference(lasa_xyz, lasa_center)
            radius = get_mag(vector)
            if radius_max == 0 or radius <= radius_max:
                self.scaff_lasa_geo_centers[name] = lasa_center
                self.scaff_lasa_grid_centers[name] = gridify_xyz(lasa_center, self.cof_grid_size)
                self.scaff_lasa_vectors[name] = vector
                self.scaff_lasa_radii[name] = radius
                self.scaff_lasa_atom_sets[name] = atom_set
                residues_in_range.append(pre_name.split('_')[0])
        copy_scaff_obj = copy_transform_object(self.scaffold)
        copy_scaff_obj.get_xyz_by_resnum(residues_in_range, update=True)
        self.reduced_scaffold = copy_scaff_obj
    
    def generate_scaff_grid(self, grid_size=1.0):
        if self.scaff_CA_names:
            reduced_scaff = copy_transform_object(self.scaffold)
            reduced_scaff.get_xyz_from_name_list(self.scaff_CA_names, True)
        else:
            reduced_scaff = self.scaffold
        self.scaff_bb_grid = { 'N' : reduced_scaff.get_atom_grid(grid_size, 'N'),\
                               'CA' : reduced_scaff.get_atom_grid(grid_size, 'CA'),\
                               'CB' : reduced_scaff.get_atom_grid(grid_size, 'CB') }


    def create_grid_axis(self, grid_size=1.0):
        scaffold_axes = self.scaffold.get_symmetric_axes(bb=True)
        scaff_center = self.scaffold.get_geo_center()
        pointer = TransformVector([0.0,0.0,0.0])
        pointer.translate(scaff_center)
        max_distance = int(max(self.scaffold.geo_center_mags[:]) + 2)
        axis_center_set = set()
        for i in xrange(max_distance):
            pointer.translate_along_axis(scaffold_axes['major'], i, update=True)
            axis_center_set.add(str(gridify_xyz(pointer.xyz[:], grid_size)))
            pointer.translate(scaff_center)
            pointer.translate_along_axis(scaffold_axes['minor'], i, update=True)
            axis_center_set.add(str(gridify_xyz(pointer.xyz[:], grid_size)))
            pointer.translate(scaff_center)
        self.axis_centers = convert_grid_string_to_list(axis_center_set)

    def origin_align(self, flip=False, align_major_to=[]):
        self.zero_branch_torsion_angles()
        self.cofactor.write_to_file('post_zero.pdb')
        if any(i > 0.001 for i in list(self.scaffold.get_geo_center())):
            self.scaffold.translate([0.,0.,0.])
        self.cofactor.translate([0.,0.,0.])#, self.zero_tor_center)
        scaff_axes = self.scaffold.get_symmetric_axes(bb=True)
        self.scaff_axes = scaff_axes
        cof_axes = self.cofactor.get_symmetric_axes()
        #check to see if aligned already
        self.cofactor.write_to_file('post_zero_translated.pdb')
        cof_axes = self.cofactor.get_symmetric_axes()
        meh = []
        for k,v in cof_axes.iteritems():
            meh.append(v)
        create_arbPDB_from_xyz(meh, 'fail_cof_axes')
        #self.cofactor.rotate_about_axis(self.scaff_axes['A'], math.pi/2.0)
        self.cofactor.major_symmetric_axis_align(self.scaffold, flip)
        self.cofactor.write_to_file('post_zero_translated_aligned.pdb')
        cof_axes = self.cofactor.get_symmetric_axes()
        meh = []
        for k,v in cof_axes.iteritems():
            meh.append(v)
        create_arbPDB_from_xyz(meh, 'fail_cof2_axes')
        meh2 = []
        for k2,v2 in scaff_axes.iteritems():
            meh2.append(v2)
        create_arbPDB_from_xyz(meh2, 'fail_scaff_axes')

    

    def coarse_grid_enum_scan(self, grid_size=1.0, flip=False, std=False, rot_interval=1, my_rotlib='', lasa='CB'):
        #self.origin_align(flip)
        #self.cofactor.write_to_file('cof3.pdb')
        #self.scaffold.write_to_file('scaff3.pdb')
        #sys.exit()
        if not self.scaff_bb_grid:
            self.generate_scaff_grid(grid_size)
        if not self.axis_centers:
            self.create_grid_axis(grid_size)

        planar_sector = int(360/self.subunits)
        scaff_axes = self.scaffold.get_symmetric_axes(bb=True)
        if flip == False:
            axis_of_rotation = scaff_axes['major']
        else:
            axis_of_rotation = scaff_axes['minor']
        matches = {}
        if not my_rotlib:
            self.generate_inverse_rot_ensemble(grid_size, std, lasa)
        else:    
            self.generate_inverse_rot_ensemble(grid_size, std, lasa, my_rotlib)
        self.scaffold.get_atom_grid(grid_size, 'CA')
        for rotation in xrange(0, planar_sector, int(rot_interval)):
            truncated_grid = {}
            for rot_name, grids_by_atom in self.bb_grid_ensemble.iteritems():
                CB_grids = grids_by_atom['CB']
                rotated_CB_grids = rotate_grid_about_axis(CB_grids, axis_of_rotation, math.radians(rotation), grid_size)
                replacement = set(rotated_CB_grids)
                
                CA_grids = grids_by_atom['CA']
                rotated_CA_grids = rotate_grid_about_axis(CA_grids, axis_of_rotation, math.radians(rotation), grid_size)
                replacement_CA = set(rotated_CA_grids)
                truncated_grid[rot_name] = [replacement,replacement_CA]
            
            for grid_center in self.axis_centers:
                addition_factor = get_vector_difference(grid_center, self.grid_center)
                match_boolean = []
                res_rotamer_table = {}
                for residue in self.branch_residue_names:
                    branch_found_match = False
                    for rotamer, my_set in truncated_grid.iteritems():
                        if rotamer.split('_')[0] == residue:
                            CB_set = my_set[0]
                            CB_lists = convert_grid_string_to_list(CB_set)
                            new_CB_set= set([str(add_vector_difference(x, addition_factor)) for x in CB_lists])
                            overlap_CB = self.scaff_bb_grid['CB'] & new_CB_set
                            
                            CA_set = my_set[1]
                            CA_lists = convert_grid_string_to_list(CA_set)
                            new_CA_set= set([str(add_vector_difference(x, addition_factor)) for x in CA_lists])
                            overlap_CA = self.scaff_bb_grid['CA'] & new_CA_set
                            
                            if not overlap_CB:
                                continue
                            else:
                                if not overlap_CA:
                                    continue
                                else:
                                    #if len(self.bb_grid_ensemble[rotamer]['N'] & overlap_N) > 0:
                                    #print 'yay'
                                    branch_found_match = True
                                    for grid_point in overlap_CA:
                                        full_CA_name = self.scaffold.grid_name_dict[grid_point]
                                        truncated_resn_resi = full_CA_name.split('_')
                                        scaffold_residue = '_'.join([truncated_resn_resi[0],\
                                                                     truncated_resn_resi[2]])
                                        if scaffold_residue in res_rotamer_table.keys():
                                            previous_rotamers = res_rotamer_table[scaffold_residue]
                                            previous_rotamers.append(rotamer)
                                        else:
                                            res_rotamer_table[scaffold_residue] = [rotamer]
                    match_boolean.append(branch_found_match)
                if all(residues_matched == True for residues_matched in match_boolean) == True:
                    matches['%s_%s' % (rotation, grid_center)] = res_rotamer_table
        if matches:
            if flip==False:
                flip_name = 'major'
            else:
                flip_name = 'minor'
            self.coarse_grid_matches[flip_name] = matches


    def inverse_rotamer_sampler(self, primary_score=1.0, secondary_score=1.0):
        
        if not self.bb_rot_ensemble:
            print 'No bb_rot_ensemble found, recursive_half_angle_sampler failing'
            sys.exit()

        if not self.scaff_lasa_geo_centers:
            radius = self.lasa_radii[max(self.lasa_radii, key=lambda i: self.lasa_radii[i])]
            radius += (radius**0.666667)   #add scaling factor to include more residues
            self.get_scaff_lasa_info(self.lasa, radius)
            self.scaffold.write_to_file('scaff_centered.pdb')
        
        reference_cof = copy_transform_object(self.cofactor)
        for res_match_name, scaff_lasa_vector in self.scaff_lasa_vectors.iteritems():
            scaff_res = res_match_name.split('_')[0]
            copy_of_reduced_scaff = copy_transform_object(self.reduced_scaffold)
            copy_of_reduced_scaff.get_xyz_by_resnum(scaff_res, update=True)
            for res_rot_name, rot_cof_obj in self.bb_rot_ensemble.iteritems():
                copy_rot_cof_obj = copy_transform_object(rot_cof_obj)
                
                cof_center = self.lasa_geo_centers[res_rot_name]
                cof_lasa_vector = self.lasa_vectors[res_rot_name]

                scaff_center = self.scaff_lasa_geo_centers[res_match_name]
                
                scaff_atom_set = self.scaff_lasa_atom_sets[res_match_name]
                vector_angle_diff = vector_angle(cof_lasa_vector, scaff_lasa_vector)
                rotate_axis = cross(cof_lasa_vector, scaff_lasa_vector)
                #need to compare rotate_axis to major and minor
                copy_rot_cof_obj.rotate_about_axis(rotate_axis, vector_angle_diff)
                copy_rot_cof_obj.write_to_file('rotated_cofactor.pdb')
                copy_rot_cof_obj.translate(scaff_center, cof_center)
                copy_rot_cof_obj.write_to_file('translated_cofactor.pdb')
                create_arbPDB_from_xyz([scaff_lasa_vector, [0.,0.,0.]], 'scaff_vec.pdb')
                create_arbPDB_from_xyz([cof_lasa_vector], 'cof_vec.pdb')
                rotamer_res_name = res_rot_name.split('_')[0]
                rotamer_type = res_rot_name.split('_')[1]
                            
                #for main residue: we need to evaluate
                try:
                    rot_lookup = int(rotamer_type)
                except ValueError:
                    rot_lookup = rotamer_type
                primary_rotamer_torsion = self.dunbrack_rotamers[rotamer_res_name[-3:]][rot_lookup][:]
                
                self.cofactor = copy_transform_object(copy_rot_cof_obj)
                print 'restype in map', rotamer_res_name[-3:]
                print 'rotamer about to enter std_devs', rot_lookup
                print 'res_match_name', res_match_name
                print 'scaff_atom_set', scaff_atom_set
                print 'straight from primary', primary_rotamer_torsion
                self.cofactor.write_to_file('rot_trans_pre_adjusted.pdb')
                primary_data_dict = self.sample_std_deviations(copy_of_reduced_scaff, rotamer_res_name,\
                                                               scaff_atom_set[:], primary_rotamer_torsion[:])
                
                self.cofactor.write_to_file('rot_trans_adjusted.pdb')
                print primary_data_dict
                print rotamer_res_name
                print res_match_name
                print res_rot_name
                '''
                if res_match_name == '79_HIS_CB':
                    if res_rot_name.split('_')[0] == '79HIS':
                        if primary_data_dict['score'] < 0.0:
                            sys.exit()
                '''
                if primary_data_dict['score'] > primary_score:
                    self.cofactor = copy_transform_object(reference_cof)
                    continue
                elif math.isnan(float(primary_data_dict['rmsa'])) == True:
                    self.cofactor = copy_transform_object(reference_cof)
                    continue
                print 'found one'
                
                bb_scaff_resname_primary = ''.join(res_match_name.split('_')[:-1])
                comb_res_rot_name = ''.join(res_rot_name.split('_'))
                prime_name = bb_scaff_resname_primary + '_' + comb_res_rot_name
                prime_data_dict = { prime_name : primary_data_dict.copy() }
                self.match_cofactors[prime_name] = copy_transform_object(self.cofactor)
                self.match_data_set[prime_name] = {prime_name : primary_data_dict.copy()}
                # { 133HIS_133HIS2100 : { 'torsions' : [10.0, 3.0, 0.0, 0.0], 'alog': -0.5, 'rmsd': 0.3 }, }
                all_sec_data_dict = {}
                #self.cofactor.write_to_file('final_match.pdb')
                pre_secondary_ref_obj = copy_transform_object(self.cofactor)
                for secondary_res_rot_name in self.bb_rot_ensemble.keys():
                    secondary_res = secondary_res_rot_name.split('_')[0]
                    secondary_rotamer = secondary_res_rot_name.split('_')[1]
                    if secondary_res != rotamer_res_name:
                        print secondary_res_rot_name
                        
                        try:
                            sec_rot_lookup = int(secondary_rotamer)
                        except ValueError:
                            sec_rot_lookup = secondary_rotamer
                        
                        secondary_rotamer_torsion = self.dunbrack_rotamers[secondary_res[-3:]][sec_rot_lookup][:]
                        print 'secondary_rotamer_torsion', secondary_rotamer_torsion
                         
                        operating_branches = self.separate_by_branch(secondary_res)
                        sec_scaff_match_atoms = set()
                        for sec_branch_name, sec_branch_obj in operating_branches.iteritems():
                            #locate each branch with res name and num and apply tor to it
                            sec_branch_obj = self.apply_tor_to_branch(copy_transform_object(sec_branch_obj), sec_branch_name,\
                                                     secondary_rotamer_torsion[:])
                            self.cofactor.update_coordinates_from_truncated_obj(copy_transform_object(sec_branch_obj))
                            #get grid here and match it to self.reduced_scaffold
                            #if pass, update a copy of primary_cofactor then run std sampler
                            cof_lasa_grid = sec_branch_obj.get_atom_grid(grid_size=5.0,\
                                                                         atom_name=self.lasa)
                            scaff_lasa_grid = self.reduced_scaffold.get_atom_grid(grid_size=5.0,\
                                                                                  atom_name=self.lasa)
                            matching_atoms = cof_lasa_grid & scaff_lasa_grid
                            if matching_atoms:
                                print matching_atoms
                                print sec_branch_name
                                print secondary_res_rot_name
                                #problem with grabing grid xyz is that many other things share
                                #the same value
                                name_xyz_from_grid = self.reduced_scaffold.get_xyz_by_grid(matching_atoms)
                                for key in name_xyz_from_grid.keys():
                                    split_key = key.split('_')
                                    split_key.pop(1)
                                    sec_scaff_match_atoms.add('_'.join(split_key))
                        self.cofactor.write_to_file('rot_trans_first_secondary.pdb')
                        pre_fine_adjustment_ref = copy_transform_object(self.cofactor)
                        for sec_bb_match in sec_scaff_match_atoms:
                            sec_scaff_atom_set = self.scaff_lasa_atom_sets[sec_bb_match][:]
                            
                            sec_scaff_res = sec_bb_match.split('_')[0]
                            sec_copy_of_reduced_scaff = copy_transform_object(self.reduced_scaffold)
                            sec_copy_of_reduced_scaff.get_xyz_by_resnum(sec_scaff_res, update=True)
                            
                            secondary_data_dict = self.sample_std_deviations(\
                                                 sec_copy_of_reduced_scaff, secondary_res,\
                                                 sec_scaff_atom_set[:], secondary_rotamer_torsion[:])
                            
                            self.cofactor.write_to_file('rot_trans_post_secondary.pdb')
                            '''if res_match_name == '86_ASP_CB':
                                if res_rot_name == '86ASP_3200':
                                    if sec_bb_match == '174_HIS_CB':
                                        if secondary_res_rot_name == '174HIS_2100':
                                            print secondary_data_dict['alog']
                                            print secondary_data_dict['rmsd']
                                            #sys.exit()'''
                            print 'sec_data', secondary_data_dict
                            if secondary_data_dict['score'] > secondary_score:
                                self.cofactor = copy_transform_object(pre_fine_adjustment_ref)
                                continue
                            
                            elif math.isnan(float(secondary_data_dict['rmsa'])) == True:
                                self.cofactor = copy_transform_object(pre_fine_adjustment_ref)
                                continue
                           
                            comb_sec_res_rot_name = ''.join(secondary_res_rot_name.split('_'))
                            bb_scaff_resname_secondary = ''.join(sec_bb_match.split('_')[:-1])
                            
                            secondary_name = bb_scaff_resname_secondary + '_' + comb_sec_res_rot_name
                            sec_data_dict = {secondary_name : secondary_data_dict.copy()}
                            prime_data_dict.update(sec_data_dict.copy())
                            self.match_cofactors[prime_name + '_' + secondary_name] = copy_transform_object(self.cofactor)
                            self.match_data_set[prime_name + '_' + secondary_name] = {prime_name : primary_data_dict.copy(),\
                                                                                      secondary_name : secondary_data_dict.copy()}
                            
                        #for other residues in the cofactor beside the rotamer we are on
                        #we need to generate an ensemble of positions while cofactor is in place
                        #select rotamers that seem promising and refine
                    self.cofactor = copy_transform_object(pre_secondary_ref_obj)
                
                self.cofactor = copy_transform_object(reference_cof)


class DysmmCofactor(SyPRISCofactor):
    
    def __init__(self, cofactor_path, scaffold_path):
        if type(cofactor_path) == str:
            self.cofactor_file_path = cofactor_path
            self.scaffold_file_path = scaffold_path
            self.scaffold_obj = Transform(read_file(scaffold_path))
            self.branch_connectivity_map = {}
            self.branch_objs = {}
            self.bb_grid_ensemble = {}
            self.frozen_atom_names = []
            super(DysmmCofactor, self).__init__(self.cofactor_file_path)

        else:
            print "SyPRISCofactor arguement is not of type str\n \
                   #and thus not a valid path name."
            sys.exit()
        self.rigid_body_rotational_ensemble = {}

    def generate_RB_rot_ensemble(self):
        return
        
