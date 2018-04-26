#!/usr/bin/env python
"""
SyPRIS v1.3

This program is intended to generate backbone matches for a given cofactor
within an oligomeric interface of any scaffold provided. The tolerance for
matches can be provided by the user.

AUTHOR: WILLIAM A. HANSEN
PROPERTY OF: KAHRE LAB

For general questions about SyPRIS email:
hawi03@gmail.com

>>>last changed 12/17/2015

changes:
SyPRIS v1.0 (9/18/15) 
1) Cofactors now only match chain A of the scaffold (makes symmetric design in Roestta simple)
2) Commented in Dsymm and avg_atom_center
3) Removed some temporary lines

SyPRIS v1.1 (12/17/2015)
1) Added user set sample_num
2) Added user set degree_factor
3) Commented in Csymm lines
4) Issue with Csymm vector to LASA fixed (Csymm fully operational) -- fixed in avg_atom_center

SyPRIS v1.2 (1/21/2016)
1) Added user option gen_scaff
2) SyPRIS will generate a logfile of pdbs that make matches under the given tolerance

SyPRIS v1.3(2/2/2016)
1) Fixed bug with added atom line
2) Included a check for if first chi has only 1 input
3) Added ForkManager to SyPRIS take -j # as input and runs on that many processors
4) Removed the os.remove(empty SyPRIS files) was not working

"""
from transform import *
import os, sys, math, numpy, time, csv, optparse, fork_manager, blargs

class FileForkManager :
    def __init__( self ) :
        self.file_for_job = {}
        self.all_files_ran = True
        self.files_that_failed = []
    def handle_successful_run( self, fm, pid ) :
        if pid not in self.file_for_job :
            print "Critical error.  Could not find file assigned to process ", pid
            for pid in self.file_for_job :
                print "Process ", pid, "responsible for", self.file_for_job[ pid ]
                sys.exit(1)
        else :
            del self.file_for_job[pid]
    def handle_failed_run( self, fm, pid ) :
        if pid not in self.file_for_job :
            print "Critical error.  Could not find file assigned to process ", pid
            for pid in self.file_for_job :
                print "Process ", pid, "responsible for", self.file_for_job[ pid ]
                sys.exit(1)
        else :
            self.files_that_failed.append( self.file_for_job[ pid ] )
            self.all_files_beautified = False
            del self.file_for_job[ pid ]                                     


def run_SyPRIS(run_list, num_cpu, options):

    print "Preparing to run SyPRIS on", len(run_list), "files"
    
    ffm = FileForkManager()
    fm = fork_manager.ForkManager(num_cpu)
    fm.error_callback = ffm.handle_failed_run
    fm.success_callback = ffm.handle_successful_run

    for file_path in run_list:
        pid = fm.mfork()
        pdb = read_file(file_path.strip('\r\n'))
        if pid == 0:
            ion_sampler( options[0], pdb, file_path.strip('\r\n'), \
                         options[1], \
                         options[2], \
                         options[3], \
                         options[4], \
                         options[5], \
                         options[6], \
                         options[7], \
                         options[8], \
                         options[9], \
                         options[10], \
                         options[11], \
                         options[12], \
                         options[13], \
                         options[14] )
            sys.exit(0)
        else:
            ffm.file_for_job[pid] = file_path
    
    fm.wait_for_remaining_jobs()
    return ffm

def exit_SyPRIS( ffm ):
    if ffm.all_files_ran:
        sys.exit(0)
    else:
        for file_name in ffm.files_that_failed:
            print "File", file_name, "did not finish"
        sys.exit(1)
             

#this is assumed that your placed atom is closer to the backbone than CB
def score_placement(prospect, axis_atom, placed_atom, scaffold_pdb, residue, scaff_chain, cof_res_ID):
    prosp_orig = atom_coordinate(prospect[0], axis_atom, int(cof_res_ID))
    prosp_atom = Transform([prosp_orig] + [atom_coordinate(prospect[0], placed_atom, int(cof_res_ID))])
    scaff_orig = atom_coordinate(scaffold_pdb, axis_atom, residue, scaff_chain)
    scaff_atom = Transform([scaff_orig] + [atom_coordinate(scaffold_pdb, placed_atom, residue, scaff_chain)])

    angle_log = angleLog([prosp_atom.xyz_total], \
                         [scaff_atom.xyz_total])
    placed_atom_mag = (scaff_atom.xyz_total[1][0] - prosp_atom.xyz_total[1][0])**2 + \
                      (scaff_atom.xyz_total[1][1] - prosp_atom.xyz_total[1][1])**2 + \
                      (scaff_atom.xyz_total[1][2] - prosp_atom.xyz_total[1][2])**2
    
    return angle_log, placed_atom_mag


def locate_half_angle(score_list, list_min, chi):
    for num, item in enumerate(score_list):
        
        if item[0] == list_min:
            if num == 0:
                half_angle = (float(item[1]) + float(score_list[1][1])) / 2.0
                second = num+1
            
            elif num == len(score_list)-1 and len(score_list) > 2:
                half_angle = (float(item[1]) + float(score_list[-2][1])) / 2.0
                second = num-1
            
            else:
                compare = [score_list[num-1], score_list[num+1]]
                compare_min = min(x[0] for x in compare)
                
                if compare[0][0] == compare_min:
                    half_angle = (float(item[1]) + float(compare[0][1])) / 2.0
                    second = num-1
                else:
                    half_angle = (float(item[1]) + float(compare[1][1])) / 2.0
                    second = num+1
            
            return half_angle, num, second 

##################################################################################
#   This definition is used to pair up all atoms that are similar in oligomers   #
#   with EXACTLY the same number of atoms at the same indices. It produces the   #
#   atoms geometric centers and their respective radii from that center. This    #
#   must be given an atom type to search.                                        #
##################################################################################
def avg_atom_center(pdb, atom_type, subunits, tolerance, base_distance = 0.0):
    ## A quick hardcoded list of possible static atoms in the scaffold
    static_atom_types = ['N','CB','CA'] 
    
    ## Atom_type_list contains the cart. coordinates of the LASA (last static atom)
    atom_type_list = []
    
    ## This is for the cases where the LASA is not apart of a traditional static atom and is in pdb still
    ## Currently not in use 
    ## **IN PROGRESS**
    added_atom_distance = 0.0
    
    ## Iterate over the pdb lines and pull out the LASA based on atom_type given
    for line in pdb:
        
        ## Check to see if the LASA is a traditional static atom stor all occurances of atom_type
        if atom_type in static_atom_types:
            
            ## Check pdb atom name against atom_type, if same store in atom_type_list
            if line[12:16].strip(' ') == str(atom_type):
                '''print 'atom equal to atom_type'
                #print '\n\n\n%s\n\n\n' % line[12:16]'''
                
                ## XYZ coordinates as first three arguements in atom_type_list, 4th is the residue number
                atom_type_list.append([float(line[30:38]), float(line[38:46]), float(line[46:54]), str(line[22:26].strip(' '))])
        
        ## If the LASA is not a traditional static atom, just pull all occurances of CB (works for cofactor)
        else:
            
            ## Check pdb atom name against 'CB', if same store in atom_type_list
            if line[12:16].strip(' ') == str('CB'):
                
                ## XYZ coordinates as first three arguements in atom_type_list, 4th is the residue number
                atom_type_list.append([float(line[30:38]), float(line[38:46]), float(line[46:54]), str(line[22:26].strip(' '))])

            ## Useful for cofactor LASA's that do not have a CB 
            ## **IN PROGRESS**
            '''elif line[12:16].strip(' ') == str(atom_type):
                q_atom = [float(line[30:38]), float(line[38:46]), float(line[46:54]), str(line[22:26].strip(' '))]
                CB_to_q = math.sqrt((CB[0] - q_atom[0])**2 +\
                                    (CB[1] - q_atom[1])**2 +\
                                    (CB[2] - q_atom[2])**2)
                if added_atom_distance == 0.0:
                    added_atom_distance += CB_to_q
                atom_type_list.append(q_atom)'''

    ## Sets is an integer that represents the number of residues with the corresponding atom_type devided by the number of subunits
    sets = len(atom_type_list)/int(subunits)
    
    ## Pairs all symmetrically corresponding residue and places them into a list further placed into set_list.
    set_list = [ atom_type_list[i::sets] for i in xrange(sets) ]
    
    ## Each list is filled with information corresponding to each unique symmetric residue set
    ## Coord is filled with non-match information while match is filled only with information based on matches found
    coord_centers = []
    coord_radii = []
    coord_vector = [] #vector from geometric center
    match_centers = []
    match_radii = []
    match_vector = [] #vector from geometric center
    
    ## For each unique residue set in set_list
    for my_set in set_list:
        
        ## Pull pdb residue from first residue in the set (last arguement)
        residue = my_set[0][-1]
        
        ## Remove the Residue from each list item in my_set
        ## This leaves just three members [x, y, z] coordinates
        my_set = [ x[:-1] for x in my_set ]

        ## Center is determined from the sum of each componenet devided by number of residues in my_set 
        center = [x/len(my_set) for x in map(sum, zip(*my_set))]

        ## Radius determined from the first residue atom_type
        ## This radius reflects the length of the atom_type coordinate from the symmetric axis
        radius = math.sqrt((my_set[0][0] - center[0])**2 + \
                           (my_set[0][1] - center[1])**2 + \
                           (my_set[0][2] - center[2])**2)
        
        ## This is from the non_static atom_type information that has yet to be integrated
        ## **IN PROGRESS**
        '''if added_atom_distance > 0.0:
            radius += added_atom_distance'''
        
        ## Coordinate vector of the atom_type with respect to the atom_type center
        vector = [my_set[0][0], \
                  my_set[0][1], \
                  my_set[0][2]]
        
        ## base_distance is the distance of some acceptable threshold (tunable by 
        ## If base_distance is greater than 0.0, add tolerance
        ## Base normally 0 unless provided
        if base_distance > 0.0:
            if radius <= (base_distance + tolerance) and radius >= (base_distance - tolerance): #tolerance was 1.0
                match_centers.append(center + [str(residue)])
                match_radii.append(radius)
                match_vector.append(vector)
        else:
            coord_centers.append(center)
            coord_radii.append(radius)
            coord_vector.append(vector)

    ## If match_centers is empty return coordinate outputs
    if not match_centers:
        
        ## Will return [], [], [] if base_distance not matched
        ## If base is 0.0 this will return a list of coordinate_centers, a list of radii, 
        ## and a list of vectors per unique symetric residue set
        return coord_centers, coord_radii, coord_vector
    
    ## If match_centers has members, return matched outputs
    else:

        ## Will only return a list of centers, radii, and vectors per residue set if base > 0.0
        return match_centers, match_radii, match_vector

def half_angle_sampler(score_set, score_min, chi, i, branch, base, scaff_chain, origin_atom_point, axis_atom_point, \
                       axis_atom, placed_atom, scaff_atom_type, scaffold_pdb, residue, score_type, static_atom, cof_res_ID):
    
    print '\n'
    print branch
    print '\n'
    print base
    print '\n'
    half_angle, first_ind, second_ind = locate_half_angle(score_set, score_min, chi)
    
    branch_obj = Transform(branch)
    print branch_obj.xyz_total
    branch_obj.rotate_object(origin_atom_point, axis_atom_point, math.radians(half_angle))
    print branch_obj.xyz_total
    print branch_obj.return_to_pdb()
    half_cofactor = [branch_obj.return_to_pdb() + base, half_angle] 
    print '\nhalf_cofactor'
    print half_cofactor
    if static_atom == True:
        new_anglelog, new_atom_mag = score_placement(half_cofactor, axis_atom, \
                                                     placed_atom, scaffold_pdb, residue, scaff_chain, int(cof_res_ID))
    else:
        prosp_coords = atom_coordinate(half_cofactor[0], placed_atom, int(cof_res_ID))
        scaff_coords = atom_coordinate(scaffold_pdb, scaff_atom_type, residue, scaff_chain)
        new_atom_mag = (scaff_coords[0] - prosp_coords[0])**2 + \
                       (scaff_coords[1] - prosp_coords[1])**2 + \
                       (scaff_coords[2] - prosp_coords[2])**2

    if score_type == 'anglelog':
        three_scores = [score_set[first_ind], score_set[second_ind], [new_anglelog, half_angle]]
    else:
        three_scores = [score_set[first_ind], score_set[second_ind], [new_atom_mag, half_angle]]
    
    score_max = max(x[0] for x in three_scores)
    score_min = min(x[0] for x in three_scores)
    score_set = [x for x in three_scores if x[0] != score_max]
    while len(score_set) < 2:
        if score_type == 'anglelog':
            score_set.append([new_anglelog, half_angle])
        else:
            score_set.append([new_atom_mag, half_angle])
    
    best_rotation = [x[1] for x in three_scores if x[0] == score_min][0]
    score_set = sorted([x[0],x[1]] for x in score_set)
    
    return score_set, score_min, half_cofactor
    
def sample_chi(tolerance, cofactor_pdb, scaffold_pdb, scaff_atom_type, subunits, chi_file, residue, cof_res_ID, sample_num, degree_factor, scaff_chain = 'A'):
    master_out = []
    master_cof_out = []
    #the last static atom is within the scaffold
    lasa_in_pdb = True
    for master, chi_set in enumerate(chi_file):
        group_of_accepted_chis = []
        group_of_accepted_params = []
        #print len(chi_set)
        for chi in xrange(len(chi_set)):
            #print 'new chi', chi
            origin_atom = chi_set[chi][1]
            axis_atom = chi_set[chi][2]
            placed_atom = chi_set[chi][3]
            if not group_of_accepted_chis:
                if chi > 0:
                    #print 'higher chi and not accepted_chi'
                    #print residue
                    #sys.exit()
                    break
                else:
                    cofactor_set = [cofactor_pdb]
                    old_degrees = []
                    old_atom_mag = []
                    old_anglelog = []
            else:
                print 'current chi ', str(chi)
                usable_group = group_of_accepted_chis[chi-1]
                usable_params = group_of_accepted_params[chi-1]
                print '\n\n\n\n\n'
                print usable_group
                print '\n\n\n\n\n'
                cofactor_set = [x[0] for x in usable_group]
                old_degrees = [x[1:] for x in usable_group]
                old_atom_mag = [x[0] for x in usable_params]
                old_anglelog = [x[1] for x in usable_params]
                #print old_degrees
            chis_for_each_independant = []
            atom_params_for_each_independant = []
            #print cofactor_set
            for indice, independant_cofactor in enumerate(cofactor_set):
                print independant_cofactor
                print origin_atom
                print 'cof_res_id', cof_res_ID
                #sys.exit()
                bushy_branch = []
                base = []
                for ind, line in enumerate(independant_cofactor):
                    if int(line[22:26]) == int(cof_res_ID):
                        bushy_branch.append(line)
                        if line[12:16].strip(' ') == axis_atom:
                            print int(cof_res_ID)
                            axis_atom_point = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
                        if line[12:16].strip(' ') == 'N':
                            start_of_block = ind + 1
                        if line[12:16].strip(' ') == origin_atom:
                            '''if chi == 0:
                                branch_start = (ind +1)/cof_res_ID
                            else:
                                branch_start = ind + 1'''
                            branch_start = ind + 1 - start_of_block
                            print branch_start
                            print int(cof_res_ID)
                            origin_atom_point = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
                    else:
                        base.append(line)
                
                print len(bushy_branch)
                branch = bushy_branch[:branch_start]
                base = bushy_branch[branch_start:] + base
                print '\n\n\n\n'
                print int(cof_res_ID)
                print branch
                print '\n'
                print base
                print '\ndone\n'
                #sys.exit()
                #print base
                first_chi_block = []
                for degree in xrange(int(chi_set[chi][-3]), int(chi_set[chi][-2]) + int(chi_set[chi][-1]), int(chi_set[chi][-1])):
                    #print degree
                    branch_obj = Transform(branch)
                    print axis_atom_point
                    print origin_atom_point
                    branch_obj.rotate_object(origin_atom_point, axis_atom_point, math.radians(degree))
                    first_chi_block.append([branch_obj.return_to_pdb() + base, degree])
                
                #print len(first_chi_block)
                '''right before appending them we want to only append the best groups'''
                #This is if CB is moving generally
                if placed_atom == scaff_atom_type:
                    lasa_in_pdb = False
                    atom_mags = []
                    for prospect in first_chi_block:
                        prosp_coords = atom_coordinate(prospect[0], placed_atom, int(cof_res_ID))
                        scaff_coords = atom_coordinate(scaffold_pdb, scaff_atom_type, residue, scaff_chain)
                        atom_mags.append([((scaff_coords[0] - prosp_coords[0])**2 + \
                                         (scaff_coords[1] - prosp_coords[1])**2 + \
                                         (scaff_coords[2] - prosp_coords[2])**2), prospect[-1]])
                    atom_mags_min = min(x[0] for x in atom_mags)
                    #print atom_mags_min
                    #didn't take a sqrt so 0.64 is equiv to 0.8 mag diff
                    if math.sqrt(atom_mags_min) < float(tolerance):

                        for i in xrange(10):
                            atom_mags, atom_mags_min, best_atom_mag = half_angle_sampler(atom_mags, atom_mags_min, chi, i, \
                                                                         branch, base, scaff_chain, origin_atom_point, \
                                                                         axis_atom_point, axis_atom, placed_atom, scaff_atom_type, \
                                                                         scaffold_pdb, residue, 'atom_mags', False, int(cof_res_ID))
                            if i == 9:
                                print best_atom_mag
                                print 'made it here'
                                #sys.exit()
                                rotation_atom_mags = [best_atom_mag[1]]
                                if sample_num > 0:
                                    for rotation_index in xrange(sample_number): #(0,10,2): and change 1.0 to 2.0
                                        additional_rotation = degree_factor * (float(rotation_index) + 1.0)
                                        if best_atom_mag[1]+additional_rotation <= float(chi_set[chi][-2]):
                                            rotation_atom_mags.append(best_atom_mag[1]+additional_rotation)
                                        if best_atom_mag[1]-additional_rotation >= float(chi_set[chi][-3]):
                                            rotation_atom_mags.append(best_atom_mag[1]-additional_rotation)
                                saved_atom_mag_cofactors = []
                                saved_atom_mag_params = []
                                for rotation in rotation_atom_mags:
                                    branch_obj = Transform(branch)
                                    branch_obj.rotate_object(origin_atom_point, axis_atom_point, math.radians(rotation))
                                    
                                    #half_cofactor = [branch_obj.return_to_pdb() + base, rotation] 
                                    half_cofactor = [branch_obj.return_to_pdb() + base, rotation]
                                    
                                    #first_chi_block.append([branch_obj.return_to_pdb() + base, degree])
                                    saved_atom_mag_cofactors.append(half_cofactor)
                                    prosp_coords = atom_coordinate(half_cofactor[0], placed_atom, int(cof_res_ID))
                                    scaff_coords = atom_coordinate(scaffold_pdb, scaff_atom_type, residue, scaff_chain)
                                    new_atom_mag = (scaff_coords[0] - prosp_coords[0])**2 + \
                                                   (scaff_coords[1] - prosp_coords[1])**2 + \
                                                   (scaff_coords[2] - prosp_coords[2])**2
                                    saved_atom_mag_params.append([new_atom_mag, 0.0])
                    
                    else:
                        #print '\n\n\n\nchi before break.\n\n\n\n', chi

                        break
                    #print anglelogs
                    #print atom_mags
                    chis_for_each_independant += saved_atom_mag_cofactors
                    atom_params_for_each_independant += saved_atom_mag_params
               
                elif len(first_chi_block) == 1:
                    chis_for_each_independant += first_chi_block
                    atom_params_for_each_independant += [[0.0,0.0]]

                else:
                    anglelogs = []
                    atom_mags = []
                    for prospect in first_chi_block:
                        angle_log, placed_atom_mag = score_placement(prospect, axis_atom, placed_atom, scaffold_pdb, residue, scaff_chain, int(cof_res_ID))
                        anglelogs.append([angle_log, prospect[-1]])
                        atom_mags.append([placed_atom_mag, prospect[-1]])
                    
                    anglelog_min = min(x[0] for x in anglelogs)
                    atom_mags_min = min(x[0] for x in atom_mags)
                    #print '\n'
                    #print anglelog_min
                    #print math.sqrt(atom_mags_min)
                    if anglelog_min < (math.log10(float(tolerance))/1.5)  and \
                       math.sqrt(atom_mags_min) < (float(tolerance)):
                        for i in xrange(10):
                            anglelogs, anglelog_min, best_anglelog = half_angle_sampler(anglelogs, anglelog_min, chi, i, \
                                                                         branch, base, scaff_chain, origin_atom_point, \
                                                                         axis_atom_point, axis_atom, placed_atom, scaff_atom_type, \
                                                                         scaffold_pdb, residue, 'anglelog', True, int(cof_res_ID))
                            atom_mags, atom_mags_min, best_atom_mag = half_angle_sampler(atom_mags, atom_mags_min, chi, i, \
                                                                         branch, base, scaff_chain, origin_atom_point, \
                                                                         axis_atom_point, axis_atom, placed_atom, scaff_atom_type, \
                                                                         scaffold_pdb, residue, 'atom_mags', True, int(cof_res_ID))

                            print anglelogs
                            print atom_mags
                            #sys.exit()
                            if i == 9:

                                rotation_anglelogs = [best_anglelog[1]]
                                if sample_num > 0:
                                    for rotation_index in xrange(sample_number): #(0,10,2): and change 1.0 to 2.0
                                        additional_rotation = degree_factor * (float(rotation_index) + 1.0)
                                        if best_anglelog[1]+additional_rotation <= float(chi_set[chi][-2]):
                                            rotation_anglelogs.append(best_anglelog[1]+additional_rotation)
                                        if best_anglelog[1]-additional_rotation >= float(chi_set[chi][-3]):
                                            rotation_anglelogs.append(best_anglelog[1]-additional_rotation)
                                saved_anglelog_cofactors = []
                                saved_anglelog_params = []
                                for rotation in rotation_anglelogs:
                                    branch_obj = Transform(branch)
                                    branch_obj.rotate_object(origin_atom_point, axis_atom_point, math.radians(rotation))
                                    if not old_degrees:
                                        half_cofactor = [branch_obj.return_to_pdb() + base, rotation] 
                                    else:
                                        half_cofactor = [branch_obj.return_to_pdb() + base]
                                        for deg in old_degrees[indice]:
                                            half_cofactor.append(deg)
                                        half_cofactor.append(rotation)
                                    
                                    saved_anglelog_cofactors.append(half_cofactor)
                                    new_anglelog, new_atom_mag = score_placement(half_cofactor, axis_atom, \
                                                                                 placed_atom, scaffold_pdb, \
                                                                                 residue, scaff_chain, int(cof_res_ID))

                                    if not old_anglelog:
                                        saved_anglelog_params.append([new_atom_mag, new_anglelog])
                                    else:    
                                        saved_anglelog_params.append([new_atom_mag + old_atom_mag[indice], new_anglelog + old_anglelog[indice]])
                                    
                                    if chi == len(chi_set)-1:
                                        #place in here the RMSD and ANGLELOG average
                                        if lasa_in_pdb == False:
                                            bb_rmsd = math.sqrt((new_atom_mag + old_atom_mag[indice])/len(chi_set))
                                            avg_anglog = (new_anglelog + old_anglelog[indice])/(len(chi_set)-1)
                                            master_out.append([bb_rmsd] + [avg_anglog] + half_cofactor[1:])
                                            master_cof_out.append(half_cofactor[0])
                                        else:
                                            bb_rmsd = math.sqrt((new_atom_mag + old_atom_mag[indice])/len(chi_set))
                                            avg_anglog = (new_anglelog + old_anglelog[indice])/len(chi_set)
                                            master_out.append([bb_rmsd] + [avg_anglog] + half_cofactor[1:])
                                            master_cof_out.append(half_cofactor[0])
                                            
                                #print best_anglelog
                                rotation_atom_mags = [best_atom_mag[1]]
                                if sample_num > 0:
                                    for rotation_index in xrange(sample_number): #(0,10,2): and change 1.0 to 2.0
                                        additional_rotation = degree_factor * (float(rotation_index) + 1.0)
                                        if best_atom_mag[1]+additional_rotation <= float(chi_set[chi][-2]):
                                            rotation_atom_mags.append(best_atom_mag[1]+additional_rotation)
                                        if best_atom_mag[1]-additional_rotation >= float(chi_set[chi][-3]):
                                            rotation_atom_mags.append(best_atom_mag[1]-additional_rotation)
                                saved_atom_mag_cofactors = []
                                saved_atom_mag_params = []
                                for rotation in rotation_atom_mags:
                                    branch_obj = Transform(branch)
                                    branch_obj.rotate_object(origin_atom_point, axis_atom_point, math.radians(rotation))
                                    #half_cofactor = [branch_obj.return_to_pdb() + base, rotation] 
                                    
                                    if not old_degrees:
                                        half_cofactor = [branch_obj.return_to_pdb() + base, rotation] 
                                    else:
                                        half_cofactor = [branch_obj.return_to_pdb() + base]
                                        for deg in old_degrees[indice]:
                                            half_cofactor.append(deg)
                                        half_cofactor.append(rotation)
                                    
                                    saved_atom_mag_cofactors.append(half_cofactor)
                                    new_anglelog, new_atom_mag = score_placement(half_cofactor, axis_atom, \
                                                                                 placed_atom, scaffold_pdb, \
                                                                                 residue, scaff_chain, int(cof_res_ID))
                                    if not old_atom_mag:
                                        saved_atom_mag_params.append([new_atom_mag, new_anglelog])
                                    else:
                                        saved_atom_mag_params.append([new_atom_mag + old_atom_mag[indice], new_anglelog + old_anglelog[indice]])

                                    if chi == len(chi_set)-1:
                                        if lasa_in_pdb == False:
                                            bb_rmsd = math.sqrt((new_atom_mag + old_atom_mag[indice])/len(chi_set))
                                            avg_anglog = (new_anglelog + old_anglelog[indice])/(len(chi_set)-1)
                                            master_out.append([bb_rmsd] + [avg_anglog] + half_cofactor[1:])
                                            master_cof_out.append(half_cofactor[0])
                                        else:
                                            bb_rmsd = math.sqrt((new_atom_mag + old_atom_mag[indice])/len(chi_set))
                                            avg_anglog = (new_anglelog + old_anglelog[indice])/len(chi_set)
                                            master_out.append([bb_rmsd] + [avg_anglog] + half_cofactor[1:])
                                            master_cof_out.append(half_cofactor[0])
                        
                    else:
                        continue
                    chis_for_each_independant += saved_anglelog_cofactors + saved_atom_mag_cofactors
                    atom_params_for_each_independant += saved_anglelog_params + saved_atom_mag_params
            
            if not chis_for_each_independant:
                break
            
            else:
                group_of_accepted_chis.append(chis_for_each_independant)
                group_of_accepted_params.append(atom_params_for_each_independant)
        group_of_accepted_chis = []
        group_of_accepted_params = []

    if not master_out:
        print 'no master out in sample chi'
        return [], []
    else:
        return master_out, master_cof_out 

def Csymm_fitter(cofactor_obj, pdb_obj, cof_rot, pdb_center, pdb_atomtype_radius, pdb_vector, cof_atom_type, subunits, tolerance, sample_num, clash_factor=2.8):
    
    ## generate the pdb axes of symmetry
    pdb_axes = pdb_obj.generate_axes()
    #pdb_geo_center = pdb_obj.get_geo_center()
    center = [0.0,0.0,0.0]

    ## Go through each of the axes of symmetry generated by gen_axes and locate the c-symm axis
    for axis in pdb_axes:
        ## Check each of the axes against a symmetric avg center axis created by CA atoms
        anglelog = angleLog( [[center]+[[x for x in axis]]], [[center]+[pdb_center[:-1]]] )
        ## Anglelog of 1 or 0 corresponds to 0degrees or 180 degrees.
        if anglelog > 0.80 or anglelog < 0.5:
            ## align_axis is the symmetric axis that corresponds to the c-symm axis 
            align_axis = [x for x in axis]
    ## This creates a vector orthogonal to the symmetric axis that points toward the atom type provided of the scaffold
    axis_obj = Transform([pdb_center[:-1]] + [pdb_vector])
    ## This translates the above orthogonal vector to the origin so that we may rotate the cofactor
    axis_obj.translate([0.0,0.0,0.0], pdb_center[:-1])
    ## make a copy of the cofactor object first by returning it to pdb format
    cofactor_pre_transform = cofactor_obj.return_to_pdb()
    variant = Transform(cofactor_pre_transform)
    ## rotate the cofactor until the c-symmetric axis of the scaffold 'align_axis' and 
    ## the c-symmetric axis of the cofacotr 'cof_rot[3]' are aligned 
    variant.symmetric_axis_align(pdb_obj, cof_rot[3], align_axis)
    ## return cofactor object back to pdb format
    cofactor = variant.return_to_pdb()
    ## collect the avg atom center of the cof_atom_type given
    cofactor_center, cofactor_radii, cofactor_vector = avg_atom_center(cofactor, cof_atom_type, subunits, tolerance)
    ## translate cofactor variant's cof_atom_center to origin before transformations
    variant.translate([0.0,0.0,0.0], cofactor_center[0])
    ## rotate the cofactor about the c-symmetric axis 'align_axis' until the cof_atom_type given atoms of 
    ## scaffold 'axis_obj.xyz_total[1]' and cofactor 'cofactor_vector[0]' are aligned
    variant.align_object_to_vector(cofactor_vector[0], axis_obj.xyz_total[1], align_axis)

   
    ## Create several variants rotated about the c-symmetric axis keeping a close distance 
    ## between the cof_atom_type and the pdb_atom_type
    plane_theta = math.acos( (2*(pdb_atomtype_radius**2) - 1) / (2*(pdb_atomtype_radius**2)) )
    rotations = [x*(float(plane_theta)/(5.0)) for x in range(sample_num+1)] + \
                [-1*x*(float(plane_theta)/(5.0)) for x in range(sample_num+1)]

    pre_rot = variant.return_to_pdb()
    rotated_groups = []
    rotation_params = []
    ## go through the rotations and make rotated groups
    for angle in rotations:
        rot_variant = Transform(pre_rot)
        ## rotate 'angle' degrees about the vector created from first input and second input
        rot_variant.rotate_object([0.0,0.0,0.0], align_axis, angle)
        ## translate our rotated cofactor out to the pdb-residue center of cof_atom_type
        rot_variant.translate(pdb_center[:-1], [0.0,0.0,0.0])
    
        cofactor1 = rot_variant.return_to_pdb()
        ## separate the cofactor into 'subunits' number of branches
        branches = symmetric_chain_separate(cofactor1, subunits)
        ## pull side chain of chain A from branch
        side_chain = branches[0][4:]

    
        ## branch_res_ID stores the residue number pulled from the first atom in branch
        branch_res_ID = int(side_chain[0][22:26])

        ## perform a clash check of the side_chain with the scaffold at a given 'clash_factor' threshold
        clash = clash_check(side_chain, pdb_obj.return_to_pdb(), clash_factor, \
                            [ [int(pdb_center[-1]) -1, 'A'], \
                              [int(pdb_center[-1]), 'A'], \
                              [int(pdb_center[-1]) +1, 'A'] ])
        ## if no clashes found, append rotation angle and cofactor to rotated_params and rotated_groups
        if clash == False:
            rotated_groups.append([cofactor1, branch_res_ID])
            rotation_params.append(angle)
        else:
            print 'Clash check failed'
        
        ## we must also assume that we can invert the c-symmetric axis of the cofactor and align the same way
        invert_variant = Transform(pre_rot)
        invert_variant.rotate_object([0.0,0.0,0.0], pdb_vector, math.radians(180))
        invert_variant.rotate_object([0.0,0.0,0.0], align_axis, angle)
        invert_variant.translate(pdb_center[:-1], [0.0,0.0,0.0])

        cofactor2 = invert_variant.return_to_pdb()
        branches2 = symmetric_chain_separate(cofactor2, subunits)
        side_chain2 = branches2[0][4:]
        ## branch_res_ID stores the residue number pulled from the first atom in branch
        branch_res_ID2 = int(side_chain2[0][22:26])
        clash2 = clash_check(side_chain2, pdb_obj.return_to_pdb(), clash_factor, \
                            [ [int(pdb_center[-1]) -1, 'A'], \
                              [int(pdb_center[-1]), 'A'], \
                              [int(pdb_center[-1]) +1, 'A'] ])
        if clash2 == False:
            rotated_groups.append([cofactor2, branch_res_ID2])
            rotation_params.append(angle)
        else:
            print 'Clash check failed'
            

    '''if not rotated_groups:
        print 'No rotated groups found'
        sys.exit()
    else:
        print rotated_groups
        sys.exit()'''
    return rotated_groups, rotation_params

## Attempts to rotate cofactor to see if LASA matches scaff atomtype
## Outputs rigid body rotated cofactor
## **IN PROGRESS** channge all CB atoms to atom_type eventually
def Dsymm_fitter(cofactor_obj, pdb_obj, subunits, residue_in_range, clash_factor=2.8):
    
    rotated_cofs = []
    ## Create a variant of the cof_obj so that changes are not saved.
    variant = cofactor_obj

    ## Prealign the cofactor axis maximum (0) of symmetry to the same axis of the scaffold
    variant.symmetric_axis_align(pdb_obj, 0)
    
    ## Prealign the cofactor axis orthoganl (1) of symmetry to the same axis of the scaffold
    variant.symmetric_axis_align(pdb_obj, 1)
    
    ## Regenerate the new axes of symmetry of the rotated cofactor
    rotate_axes = variant.generate_axes()

    ## Store this postion as the first of 6 D2 cofactor positions
    cofactor1 = variant.return_to_pdb()
    rotated_cofs.append(cofactor1)

    ## Rotate about the third axis, eig minimum by 90 degrees
    variant.rotate_object([0.0,0.0,0.0], rotate_axes[2], math.radians(90))

    ## Store this position as the second of 6 D2 cofactor positions
    cofactor2 = variant.return_to_pdb()
    rotated_cofs.append(cofactor2)

    ## For each axis of symmetry
    for axis in xrange(2):
        
        ## Create a variant of cofactor1 so that changes are not saved.
        variant1 = Transform(cofactor1)
        
        ## Rotate first of two variants 90 degrees about each axis of symmetry
        variant1.rotate_object([0.0,0.0,0.0], rotate_axes[axis], math.radians(90))
        
        ## Store this position as the third and fifth of 6 D2 cofactor positions
        rotated_cofs.append(variant1.return_to_pdb())
        
        ## Create a variant of cofactor2 so that changes are not saved.
        variant2 = Transform(cofactor2)
        
        ## Rotate second of two variants 90 degrees about each axis of symmetry
        variant2.rotate_object([0.0,0.0,0.0], rotate_axes[axis], math.radians(90))
        
        ## Store this position as the fourth and sixth of 6 D2 cofactor positions
        rotated_cofs.append(variant2.return_to_pdb())
    
    ## Find the CB atom coordinates of residue_in_range on chain A (the master chain)
    scaff_chainA_cb_coord = atom_coordinate(pdb_obj.return_to_pdb(), 'CB', residue_in_range, 'A')[:-1]
    
    rotated_groups = []
    rotation_params = []
    ## Iterate over all 6 of the rotated_cofs
    for cof in range(len(rotated_cofs)):
        
        ## Break up rotated cofs into branches equal to the number of subunits
        branches = symmetric_chain_separate(rotated_cofs[cof], subunits)

        ## cofactor_bbless_branch_set = a list of all branches without backbone atoms
        cofactor_bbless_branch_set = []

        ## cb_distances = a group of magnitudes between each cofactor CB and the chain A scaffold CB
        cb_distances = []

        ## cof_res_IDs = a list of residue numbers corresponding to the cb_distances
        cof_res_IDs  = []

        ## Iterate over each branch within the cofactor
        for branch in branches:
            
            ## branch_res_ID stores the residue number pulled from the first atom in branch
            branch_res_ID = int(branch[0][22:26])
            
            ## Store resID in cof_res_IDs
            cof_res_IDs.append(branch_res_ID)
            if type(branch_res_ID) != int:
                print "not grabbing residue for branch bug here"
                sys.exit()
                
            ## Obtain cofactor coordinate for CB of each branch    
            cof_cb = atom_coordinate(branch, 'CB', branch_res_ID)

            ## Calculate and store magnitude difference between scaffold and cofactor CB
            dist = math.sqrt((scaff_chainA_cb_coord[0] - cof_cb[0]) **2 +\
                             (scaff_chainA_cb_coord[1] - cof_cb[1]) **2 +\
                             (scaff_chainA_cb_coord[2] - cof_cb[2]) **2)
            cb_distances.append(dist)

            ## Store a branch without backbone atoms in cofactor_bbless_branch_set for clash check
            bbless_branch = branch[4:]
            cofactor_bbless_branch_set.append(bbless_branch)

        ## Find minimum magnitde difference between CB of scaffold and cofactor
        min_dist = min(cb_distances)
        print 'my min dist is:  ' + str(min_dist)
        print cb_distances

        ## Determine the index of the lowest cb_dist and use index to obtain the residue
        ## of the cofactor that matches chain A of the scaffold
        for ind, cb_dist in enumerate(cb_distances):
            if cb_dist == min_dist:
                
                ## Save the cof_res_ID corresponding to minimum CB difference
                corresponding_cof_res = cof_res_IDs[ind]

                ## Save the back boneless branch corresponding to minimum CB difference
                corresponding_cof_bbless_branch = cofactor_bbless_branch_set[ind]
        
        ## Perform a clash check of the back boneless branch with chain A
        ## exlcuding + or - 1 residue postion of the matching residue
        ## Return True if clash, false otherwise
        clash = clash_check(corresponding_cof_bbless_branch, pdb_obj.return_to_pdb(), clash_factor, \
                            [ [int(residue_in_range) -1, 'A'], \
                              [int(residue_in_range), 'A'], \
                              [int(residue_in_range) +1, 'A'] ])
        
        print 'clash is: ' + str(clash)
        ## If the cofactor branch does not clash with chain A continue
        if clash == False:
            
            ## Store the rotated cofactors that pass a clash check as 1st member
            ## Store the corresponding cofactor residue number needed for chi sampling
            rotated_groups.append([rotated_cofs[cof], corresponding_cof_res])

            ## Store which of the 6 rotational cofactor rotations produced a non-clash
            rotation_params.append('D2_%s' % cof)
    
    ## rotated_groups = [[rotated_cofs[cof], corresponding_cof_res],[]]
    ## rotation_params = [D2_#, D2_#, ...]
    return rotated_groups, rotation_params

################################################################################
#
#
################################################################################
def ion_sampler(tolerance, \
                pdb, \
                pdb_filename, \
                cofactor, \
                subunits, \
                symm_type, \
                gen_match, \
                gen_scaff, \
                rmsd_threshold, \
                avgAlog_threshold, \
                sample_num, \
                degree_factor, \
                cof_atom_type = '', \
                scaff_atom_type = '', \
                cofactor_rotation = (), \
                chi_file = [], \
                automate = False):

    ## If the pdb file is input with a path and '/' separate and grab last item (the file)
    pdb_filename = pdb_filename.split('/')[-1]
    
    ## If the last four characters are '.pdb' save filename without '.pdb'
    if pdb_filename[-4:] == '.pdb':
        pdb_filename = pdb_filename[:-4]
    print pdb_filename

    #center variable contains geometric center for the atom of choice as well as the residue indices generating center
    #the radii is just the radius from geometric center to the atom_type chosen 
    #if match found is set True otherwise it will be false (atom_type too far for complex)
    
    ## Create a transform object from our scaffold
    pdb_obj = Transform(pdb)
    
    ## Translate objects geometric center to a zero origin
    pdb_obj.translate([0.0,0.0,0.0])
    
    ## Return new coordinates back to pdb format
    pdb = pdb_obj.return_to_pdb()
    
    ## Copy above three steps for the input cofactor
    cofactor_obj = Transform(cofactor)
    cofactor_obj.translate([0.0,0.0,0.0])
    cofactor = cofactor_obj.return_to_pdb()

    ## cofactor_center          = list the geometric center of the given atom_type
    ## cofactor_radii           = list the radial distance from the atom_type geometric center to atom_type coordinate (corresponds with center)
    ## atomtype_vector          = list the vector from atom_type geometric center to atom_type coordinate (corresponds with center and radii)
    ## The fifth arguement is not supplied for the cofactor setting base = 0.0 and accepting any cof_atom_type occurances
    cofactor_center, cofactor_radii, atomtype_vector = avg_atom_center(cofactor, cof_atom_type, subunits, tolerance)
    ## pdb_atom_centers         = list of all geometric centers of the given atom_type for scaffold
    ## pdb_atom_radii           = list of all radial distances from the atom_type geometric center to atom_type coordinate (corresponds with center)
    ## pdb_atomtype_vectors     = list of all vectors from atom_type geometric center to atom_type coordinate (corresponds with center and radii)
    ## The fifth arguement is the determined cofactor radii. Will be used to remove residue sets that are not within range of cof_radius + or - tolerance
    pdb_atom_centers, pdb_atom_radii, pdb_atomtype_vectors = avg_atom_center(pdb, scaff_atom_type, subunits, tolerance, cofactor_radii[0])
    ## if pdb_atom_centers contains an empty list then the scaffold contains no matches
    if not pdb_atom_centers:
        print 'no residues close enough'
        sys.exit()
    
    ## For any potential homomeric residue sets that lie within the tolerance allowed distance from the cofactor LASA
    ## first check the symmetry type, perform rigid body rotations, and finally submit for samplinge with chi_dist_file
    '''for i,x in enumerate(pdb_atom_centers):
        if int(x[-1]) == 249:
            pdb_atom_radii = [pdb_atom_radii[i]]
            pdb_atomtype_vectors = [pdb_atomtype_vectors[i]]
            pdb_atom_centers = [x]
    print pdb_atom_radii
    print pdb_atomtype_vectors
    print pdb_atom_centers'''

    match_found = False

    for index, res_set_center in enumerate(pdb_atom_centers):    
        
        ## This check is for any Cn symmetry type
        if symm_type == 'c':
            
            ## Rotated_groups   = cofacotr objects in position for dihedral sampling
            ## Rotations        = the degree rotation within the plane of symmetry
            ## Perform the rigid body rotations and translations on the cofactor with respect to acceptable residue match sets
            ## pdb_atom_radii and pdb_atomtype_vectors correspond to pdb_atom_centers, index represents similar members
            ## A 2.8A clash check is performed on the non-mobile cofactor atoms and the surrounding scaffold backbones
            ## C symmetry cases force the first residue and chain A to align
            ## **IN PROGRESS** add options for clash check distance
            rotated_groups, rotations = Csymm_fitter(cofactor_obj, pdb_obj, cofactor_rotation, res_set_center, pdb_atom_radii[index], pdb_atomtype_vectors[index], cof_atom_type, subunits, tolerance, sample_num)
            scaff_chain = "A"
        
        ## This check is for any Dn symmetry type
        elif symm_type == 'd':

            ## Rotated_groups   = cofacotr objects in position for dihedral sampling paired scaffold chain ID **IN PROGRESS**
            ## Rotations        = the degree rotation within the plane of symmetry
            ## Perform the rigid body rotations and translations on the cofactor with respect to acceptable residue match sets
            ## res_set_center[-1] = the scaffold residue number of the 'in-range' residue for a given center.
            ## A 2.8A clash check is performed on the non-mobile cofactor atoms and the surrounding scaffold backbones
            ## **IN PROGRESS** add options for clash check distance
            ##                 also, outputs chain for matching but should output correct branch residue not chain
            rotated_groups, rotations = Dsymm_fitter(cofactor_obj, pdb_obj, subunits, res_set_center[-1])
            
        ## After performing rigid body rotations and translations check to see if a rotated cofactor exists, if not re-loop
        if not rotated_groups:
            continue
        
        ## Iterate through each accepted rotated_cofactor and submit to sample_chi
        for rot_cof_ind, rotated_cofactor in enumerate(rotated_groups):
            
            ## rotated_cofactor = a nonclashing rigid rotation in 0th term and the paired chain identity 1st term
            '''if symm_type == 'c':
                rotated_cofactor[1] = []'''

            ## tolerance           = the value that scales acceptance. Higher value is more acceptance.
            ## rotated_cofactor[0] = the non-clashing cofactor with performed rigid body rotation/translations.
            ## res_set_center[-1]  = the scaffold residue number of the 'in-range' residue for a given center.
            ## rotated_cofactor[1] = the cofactor residue number of the residue matching chain A.
            match_ensemble, master_cof = sample_chi(tolerance, rotated_cofactor[0], pdb_obj.return_to_pdb(), scaff_atom_type, subunits, chi_file, res_set_center[-1], rotated_cofactor[1], sample_num, degree_factor)
       
            if not match_ensemble:
                continue
            else:
                csvfile = open('SyPRIS_%s.csv' % pdb_filename, 'a')
                writer = csv.writer(csvfile)
                with open('tolPDBs.log', 'a') as myPMS:
                    myPMS.write('%s \n' % pdb_filename)
                for ind_cof, matched_cof in enumerate(match_ensemble):
                    print float(matched_cof[0])
                    print rmsd_threshold
                    print float(matched_cof[1])
                    if float(matched_cof[0]) < rmsd_threshold and \
                       float(matched_cof[1]) < avgAlog_threshold:
                        writer.writerow([pdb_filename]+[res_set_center[-1]]+[rotations[rot_cof_ind]]+["%.11f" % float(x) for x in matched_cof])
                        if gen_match == True:
                            match_found = True
                            with open('%s_%.11f_%.11f.pdb' % (pdb_filename, matched_cof[0], matched_cof[1]), 'w') as mycofFile:
                                mycofFile.writelines(master_cof[ind_cof])
                csvfile.close()
        #if match_count == 0:
            #os.remove('SyPRIS_%s.csv' % pdb_filename)


    if match_found == True and gen_scaff == True:
        write_file(pdb, '%s.scaff' % pdb_filename)


    
    ###IN PROGRESS###
    ##This currently does not work at all
    if automate == True:
        os.system('/home/wah49/Testing_grounds/mutateCofRes_inprog3.py --pdb-path "/lab/symm_scaffolds/d2/INPUT_files/" --cofactor "/home/wah49/SyPRIS_hits/SF4tet/" --SyPRIS-file "/home/wah49/SyPRIS_hits/SF4tet/SyPRIS_%s2.csv" --resType "C" "CYZ" --chi-dist "/home/wah49/chi_dist_files/SF4.chi" --chinum 2 --cst-path "/home/wah49/coordCSTfiles/SF4tet/" -o "/home/wah49/SyPRIS_outPDBs/"' % pdb_filename)
        
        pdbs, residues, rigid_rots, rmsds, anglelogs, chis = read_SyPRIS_file(("/home/wah49/SyPRIS_hits/SF4tet/SyPRIS_%s2.csv" % pdb_filename).strip('\n\r'))
        for ind, match in enumerate(pdbs):
            os.remove(match + '_' + rmsds[ind] + '_' + anglelogs[ind] + '.pdb')
                    
    #sys.exit()

def main(argv):
    
    parser = optparse.OptionParser(usage="\n\nGenerate csv of symmetric cofactor matches.")
    
    parser.add_option('--path', dest = 'path_filename',
        help = 'The input pdb')
    
    parser.add_option('--cofactor', dest = 'cofactor_filename',
        help = 'The input cofactor')
    
    parser.add_option('--s-units', type="int", nargs=1, dest = 'subunits',
        default = 2,
        help = 'Number of symmetric subunits ( default = 2 )')

    parser.add_option('--s-type', type='str', nargs=1, dest = 'symm_type',
        default = 'c', 
        help="Set symmetry type C,D,H, etc. ( default = 'c' )")

    parser.add_option('--tol', type='float', nargs=1, dest = 'tolerance',
        default = 1.0, 
        help="Set tolerance magnitude for matching. ( tol > 1 is low tol < 1 is high )\n \
        Directly relates to LASA atom distance. If tol = 1.0 is 1A allowed. \n \
        Default = 1.0")

    parser.add_option('--sample-num', type='int', nargs=1, dest = 'sample_num',
        default = 0, 
        help="Set a sampling deviation. ( sample_num = 0 is closest match only )\n \
        Sample_num > 0 increases the rotational sampling on either side of the best match torsion angle. \n \
        Larger sampling number yields greater variability. \n \
        Default = 0")

    parser.add_option('--degree-factor', type='float', nargs=1, dest = 'degree_factor',
        default = 1.0, 
        help="Set a sampling degree. ( sample_num = 1.0 sample_num deviations are 1x )\n \
        Increasing the degree factor multiples the sample number deviation. \n \
        Set degree_option > 0.0. \n \
        Default = 1.0")

    parser.add_option('--gen-match',
        help="If set true will generate a pdb file with matched cofactor",
        default = False, 
        action = "store_true", dest='gen_match')

    parser.add_option('--gen-scaff',
        help="If set true will generate the scaffold pdb",
        default = False, 
        action = "store_true", dest='gen_scaff')

    parser.add_option('--automate',
        help="If set true will automatically run post SyPRIS program",
        default = False, 
        action = "store_true", dest='automate')

    parser.add_option('--rmsd', type='float', nargs=1,  dest = 'rmsd_threshold',
        default = '1.0',
        help="Threshold set for maximum rmsd across matched atoms: ( default = 1.0 )")

    parser.add_option('--alog', type='float', nargs=1,  dest = 'avgAlog_threshold',
        default = '0.0',
        help="Threshold set for maximum average anglog: ( default = 0.0 )")
    
    parser.add_option('-j', type='int', nargs=1,  dest = 'num_cpu',
        default = '1',
        help="Number of cpus you wish to run this program: ( default = 1 )")
    
    parser.add_option('--cofact-a', dest = 'cof_atom_type',
        default = '',
        help="From sidechain to N the last static atom: ( 'CB' generally used )")
    
    parser.add_option('--scaff-a', dest = 'scaff_atom_type',
        default = 'CB',
        help="From sidechain to N the last static atom: ( default = 'CB' )")
    
    parser.add_option('--cofactor-rotation', type='int', nargs=4, dest = 'cofactor_rotation',
        default = (0,360,10,2),
        help="For C symmetry, rotation about the symmetric axis chosen (4th int) \
              while sampling the dof between the first and second int at 3rd int intervals. \
              ( default = (0, 360, 10, 2) ) 4th (axis) 2=eig min")
   
    parser.add_option('--chi-file', dest = 'chi_file',
        default = '',
        help="For any symmetry, it is a supplied list of chi ranges and iterations within the\n \
              range to sample for a given dihedral.\n\n Example:\n\nCG CB CA N  30  60  10\n\n")

    (options,args) = parser.parse_args()
    
    start = time.time()
    pdb_list = read_file(str(options.path_filename))
    cofactor = read_file(str(options.cofactor_filename))
    if options.symm_type.lower() == 'c':
        print "Symmetry type listed as:     C%s" % (str(options.subunits))
    elif options.symm_type.lower() == 'd':
        if int(options.subunits) == 4:
            print "Symmetry type listed as:     D2"
            options.cofactor_rotation = ()
        else:
            print "Submitted D symmetry but incorrect number of symmetric subunits"
            sys.exit()

    else:
        print "No symmetry type given."
        sys.exit()

    ffm = run_SyPRIS( pdb_list, options.num_cpu, \
       

        ## 1.  Tolerance impacts code at every distance check. Values of 1.0 will allow 1.0A.
        ##     Anglelog will be affected by a logrithmic factor proportional to distance change.
        ##     Used in:
        ##     avg_atom_center      tolerance added to scaffold residue radii
        ##     sample_chi           tolerance used to check anglelog with:  < math.log10(float(tolerance))/1.5
        ##     sample_chi           tolerance used to check atom_mag with:  < tolerance
        [options.tolerance, \
       

        ## 2.  A list of lines of the input pdb.
        ##     Used in:
        ##     ion_sampler          pdb used as arguemet: Transform(pdb)
        ##     avg_atom_center      pdb needed to determine radii
        #pdb, \
       

        ## 3.  The pdb filename.
        #options.pdb_filename, \
       

        ## 4.  A list of lines of the input cofactor.
        cofactor, \


        ## 5.  Number of subunits of symmetric protein.
        ##     Used in:
        ##     ion_sampler          pdb used as arguemet: Transform(pdb)
        ##     Csymm_fitter         >> avg_atom_center
        ##     Dsymm_fitter         >> symmetric_chain_separate
        ##     avg_atom_center      pdb needed to determine radii
        ##     sample_chi           tolerance used to check atom_mag with:  < tolerance
        int(options.subunits), \


        ## 6.  Symmetry type i.e. 'c', 'd', etc.
        options.symm_type, \
        
        ## 7.  Bool for generating the matching cofactor as a file
        options.gen_match, \

        ## 8.  Bool for generating the matching cofactor as a file
        options.gen_scaff, \
        
        ## 9.  Maximum rmsd allowed for match acceptance
        options.rmsd_threshold, \

        ## 10.  Maximum anglog averaged over all placed atoms for match acceptance
        options.avgAlog_threshold, \

        ## 11.  Number of deviations away from the perfect rotation
        options.sample_num, \
        
        ## 12.  Degree of deviation away from the perfect rotation per sample_num
        options.degree_factor, \

        ## 13. The last static atom (LASA) of the cofactor, generally 'CB'
        ##     Used in:
        ##     Csymm_fitter         >> avg_atom_center
        ##     avg_atom_center      Used to make atom_type_list (var name atom_type)    
        options.cof_atom_type, \
        
        ## 14. The last static atom (LASA) of the pdb scaffold, generally 'CB'
        ##     Used in:
        ##     sample_chi           Used to check if scaffold LASA is the same as placed cofactor atom
        ##                          >>avg_atom_center
        ##                          >>half_angle_sampler
        ##                          >>atom_coordinate
        ##     avg_atom_center      Used to makek atom_type_list (var name atom_type)
        ##     half_angle_sampler   >>atom_coordinate
        options.scaff_atom_type, \

        ## 15. Rotation applied to the cofactor in only C-symmetric cofactors
        ##     Used in:
        ##     Csymm_fitter         axis used for rotation cofactor_rotation[3]
        options.cofactor_rotation, \

        ## 16. The chi distribution file necesary for sampling the dihedral space.
        ##     Used in:
        ##     sample_chi           Iterate through chi-block (lines) and apply rotations 
        read_chi_dist_file(options.chi_file), \

        ## 17. Turn on automation (run post sypris program immeadiately after).
        False ])
    
    end = time.time()
    print end-start
    exit_SyPRIS(ffm)

if __name__ == "__main__":
    main(sys.argv[1:])

"""
Future work and known issues:
1) Cofactor LASA geometric expansion should be by user input
2) Comments missing for several definitions
3) All **IN PROGRESS** lines
4) Post SyPRIS program accesible through -automate option is not working
5) User set maximum hits produced
6) Cloud hit storage (append branch to cofactor as models) per residue matched
7) Plan to setup up SyPRIS as a group of classes for information transfer
8) Add multi-sidechain capabilities and multithreading
"""
