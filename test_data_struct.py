#!/usr/bin/env python
#from rosetta import *
#from toolbox import *
import math, numpy, csv, sys, os, time
#import rosetta.core.conformation.symmetry
#import rosetta.core.pose.symmetry
#import rosetta.core.scoring.symmetry
#import rosetta.protocols.simple_moves.symmetry

def read_file(file_path):
    with open(file_path, 'r') as f:
        my_file = f.readlines()
    return my_file

def read_write_file(file_path):
    with open(file_path, 'r+') as f:
        my_file = f.readlines()
    return my_file

def write_file(input_pdb, path='unamed_pdb.pdb', chain=''):
        
    if not chain:
        with open(path, 'w') as myfile:
            myfile.writelines(input_pdb)
    else:
        pdb = []
        for line in input_pdb:
            if line.startswith("ATOM  ") or line.startswith("HETATM"):
                pdb.append(line)
        if type(chain) == list or type(pose) == tuple:
            pdb = [x for x in pdb if str(x[21]) in chain]
        elif type(chain) == str:
            pdb = [x for x in pdb if str(x[21]) == chain]
        with open(path, 'w') as myfile:
            myfile.writelines(pdb)

def read_chi_dist_file(file_path):
    with open(file_path, 'r') as f:
        my_file = f.readlines()
    blocks = [i.split() for i in my_file]
    final_blocks = []
    block = []
    for item in blocks:
        if item == blocks[-1]:
            if not block:
                continue
            else:
                block.append(item)
                final_blocks.append(block)
                block = []
        elif not item:
            if not block:
                continue
            else:
                final_blocks.append(block)
                block = []
            continue
        else:
            block.append(item)
    return final_blocks

def read_SyPRIS_file(file_path):
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
    

def pdb_to_xyz(file_input): #was xyz_matrix_gen
    xyz1 = []
    #print("")
    for line in file_input:
        if line.startswith("ATOM") or line.startswith("HETATM"):
            x = float(line[30:38].strip())
            y = float(line[38:46].strip())
            z = float(line[46:54].strip())
            key = str(int(line[22:30]))+'_' +\
                          line[21] + '_' +\
                          line[17:20] + '_' +\
                          line[11:17].strip()
            #if line[12:16].strip() == "CA" or line[17:21].strip() == "****":
            #xyz.append([x, y, z])
            xyz1.append([x, y, z])
    nums = xrange(len(xyz1))
    xyz2 = dict(zip(nums, xyz1))
    return xyz1, xyz2

def pose_to_xyz(pose):
    #rosetta.init()
    xyz = []
    total_res = pose.n_residue()
    for res in xrange(total_res):
        for atom in xrange(len(pose.residue(res+1))):
            xyz.append( [pose.residue(res+1)[atom][0], pose.residue(res+1)[atom][1], pose.residue(res+1)[atom][2]] )
    
    return xyz

def xyz_to_pdb(main, rep_coords): #was rep_coords
    atms = []
    pdb = []
    for line in main:
        if line[0:6] in ["ATOM  ", "HETATM"]:
            atms.append(line)
        elif line[0:6] == "REMARK":
            pdb.append(line)

    for index, new_coord in enumerate(rep_coords):
        coord = atms[index]#for num, line in enumerate(atms):
        x_comp = "%.3f" % new_coord[0]
        y_comp = "%.3f" % new_coord[1]
        z_comp = "%.3f" % new_coord[2]

        line_out = ''.join([coord[0:30],('{0:>8}').format(x_comp),('{0:>8}').format(y_comp),('{0:>8}').format(z_comp),coord[54:]])
        pdb.append(line_out)

    return pdb


def create_arbPDB_from_xyz(coords): #was rep_coords
    
    preffix = "ATOM      1  X   UNK Y   1    "
    suffix  = "  1.00  0.00           Z\n"
    out_pdb = []
    for index, coord in enumerate(coords):
        x_comp = "%.3f" % coord[0]
        y_comp = "%.3f" % coord[1]
        z_comp = "%.3f" % coord[2]

        line_out = ''.join([preffix, ('{0:>8}').format(x_comp), \
                                     ('{0:>8}').format(y_comp), \
                                     ('{0:>8}').format(z_comp), \
                            suffix])
        out_pdb.append(line_out)

    return out_pdb

def make_CA_trace(pdb):
    pdb = [ x for x in pdb if x[11:17].strip() == "CA" ]
    return pdb

'''def xyz_to_pose
    unfinished part of code, will get around to it'''


def axis_align(xyz_total1, xyz_total2, axis, symm='D3'):
    #determin protein once
    pose1_obj = Transform(xyz_total1)
    if symm == 'D3':
        pose1_axes = pose1_obj.generate_axes('D3')
    else:
        pose1_axes = pose1_obj.generate_axes()
    print "pose1_axes are:", pose1_axes
    pose2_obj = Transform(xyz_total2)
    if symm == 'D3':
        pose2_axes = pose2_obj.generate_axes('D3')
    else:
        pose2_axes = pose2_obj.generate_axes()
    print "pose2_axes are:", pose2_axes

    ortho_vector = cross(pose2_axes[axis], pose1_axes[axis])
    theta = vector_angle(pose2_axes[axis], pose1_axes[axis])

    rotated_pose2 = []
    for line in xyz_total2:
        #use rotation matrix to find new xyz coords for each atom in ion complex
        #xyz, xyz_radius = normalize(line)
        print "line that we are taking:", line 
        coord_obj = Transform(line)
        new_coord = coord_obj.rotate_vector(ortho_vector, theta)
        print "new_coord is:", new_coord
        rotated_pose2.append(new_coord)

    return rotated_pose2

def generate_axes(xyz, symm = ''):
    xyz_coord = numpy.array(xyz, float)
    xyz_center = numpy.mean(xyz_coord, 0)
    #center each file with geometric center
    delta = xyz_coord - xyz_center

    #return principal axis matrix
    inertia = numpy.dot(delta.transpose(), delta)
    e_values, e_vectors = numpy.linalg.eig(inertia)

    #The axis we are concerned with at this time is the smallest eigan value (3) but writing in 1 and 2 for later
    for i in xrange(len(e_values)):
        if e_values[i] == max(e_values):
            eval1 = e_values[i]
            axis1 = e_vectors[:,i]
        elif e_values[i] == min(e_values):
            eval3 = e_values[i]
            axis3 = e_vectors[:,i]
        else:
            eval2 = e_values[i]
            axis2 = e_vectors[:,i]

    if symm == 'D3':
        #axis1 and axis3 are correct vectors, need to rotate axis1 about axis3 two more times
       
        half_axis = [-1.0 *(axis1[0] + axis2[0]), \
                     -1.0 *(axis1[1] + axis2[1]), \
                     -1.0 *(axis1[2] + axis2[2])]
 

        axis_coords = [axis1, axis2, axis3, half_axis]
        return axis_coords

    axis_coords = [axis1, axis2, axis3]

    return axis_coords

def get_angle(xyz, xyz1, xyz2):
    xyz = [x - xyz1[i] for i,x in enumerate(xyz)]
    xyz2 = [x - xyz1[i] for i,x in enumerate(xyz2)]
    angle = vector_angle(xyz, xyz2)
    return angle

def vector_angle(xyz, xyz2):

    #print "Vecotr_angle vector inuputs:"
    #print xyz
    #print xyz2
    #print ''
    radius_xyz = math.sqrt((xyz[0]**2)+(xyz[1]**2)+(xyz[2]**2))
    radius_xyz2 = math.sqrt((xyz2[0]**2)+(xyz2[1]**2)+(xyz2[2]**2))
    
    dot_over_mag = ((xyz[0]*xyz2[0]) + (xyz[1]*xyz2[1]) + (xyz[2]*xyz2[2])) / (radius_xyz * radius_xyz2)
    #print dot_over_mag
    if dot_over_mag > 0.9999999 and dot_over_mag < 1.0000001: 
        return 0.0
    elif dot_over_mag < -0.9999999 and dot_over_mag > -1.0000001: 
        return 0.0
    else:
        theta_angle = math.acos(dot_over_mag)
    return theta_angle

def get_mag(xyz1, xyz2):
    mag = math.sqrt((xyz1[0] - xyz2[0])**2 +\
                    (xyz1[1] - xyz2[1])**2 +\
                    (xyz1[2] - xyz2[2])**2)
   
    #mag = math.sqrt( sum(n * n for n in rotate_axes[0]) )

    return mag


#for angle_log to work it requires a list of vectors
#example: [ [[1,2,3],[4,5,6]], [[],[]] ]
#where paired atom coords would produce a vector
#from atom1 to atom2 threshold default is 20 degrees, 
#and must be supplied as radians
def angleLog(test, scaffold, threshold=(math.pi/9.0)):
    compare_vectors = zip(test, scaffold)
    vector_ang_sum = 0.00001
    #print 'test is\n'
    #print test
    #print 'scaffold is\n'
    #print scaffold
    #if len(compare_vectors) > 1:
    for num_vectors, vector_pairs in enumerate(compare_vectors):
        vector1 = Transform(vector_pairs[0])
        vector1.translate([0.0,0.0,0.0], vector_pairs[0][0])
        vector2 = Transform(vector_pairs[1])
        vector2.translate([0.0,0.0,0.0], vector_pairs[1][0])
        vector_ang_sum += vector_angle(vector1.xyz_total[1], vector2.xyz_total[1])
    #print vector_ang_sum
    angle_log = math.log10((vector_ang_sum) / (len(compare_vectors) * threshold))
    return angle_log

#if no chain is given it assumes you mean the first name/residue number found
def atom_coordinate(pdb, name, residue=1, chain=[]):
    for ind, line in enumerate(pdb):
        if int(line[22:26]) == int(residue):
            #print 'hi'
            if type(name) == str:
                name = [name]
            if line[12:16].strip(' ') in name:
                #print 'low'
                if not chain:
                    return [ float(line[30:38]),
                             float(line[38:46]),
                             float(line[46:54]), 
                             str(line[21]) ] 
                else:
                    if str(line[21]) == chain:
                        return [ float(line[30:38]),
                                 float(line[38:46]),
                                 float(line[46:54]),
                                 str(line[21]) ]

#if chain not specified then it will include all the atoms of all chains
def list_of_resAtoms(pdb, residue, chain=[]):
    atom_names = []
    for ind, line in enumerate(pdb):
        if line.startswith('ATOM') or line.startswith('HETATM'):
            if int(line[22:26]) == int(residue):
                #print 'hi'
                if not chain:
                    if line[11:17].strip(' ') not in atom_names:
                        atom_names.append(line[11:17].strip(' '))
                else:
                    if str(line[21]) == chain:
                        atom_names.append(line[11:17].strip(' '))
    return atom_names
    
                        
def replace_coordLine(line, coord):
    newline = ''.join([line[0:30], ('{0:>8}').format("%.3f" % coord[0]), \
                                   ('{0:>8}').format("%.3f" % coord[1]), \
                                   ('{0:>8}').format("%.3f" % coord[2]), line[54:]])
    return newline

def symmetrize_pose(pose, path_to_symm):
    pose_symm_data = core.conformation.symmetry.SymmData(pose.n_residue(), pose.num_jump())
    pose_symm_data.read_symmetry_data_from_file("%s" % (path_to_symm))
    core.pose.symmetry.make_symmetric_pose(pose, pose_symm_data)

def symmetric_chain_separate(pdb, subunits, HETATM=False):
    if HETATM == False:
        master = []
        for line in pdb:
            if line[0:6] in ["ATOM  "]:
                master.append(line)
    else:
        master = []
        for line in pdb:
            if line[0:6] in ["ATOM  ", "HETATM"]:
                master.append(line)
    chains = zip(*(iter(master),)*(len(master)/subunits))
    return chains

def gen_residue_list(pdb):
    res_nums = []
    for line in pdb:
        if int(line[22:26]) not in res_nums:
            res_nums.append(int(line[22:26]))
    return res_nums

def block_pdb_by_res(pdb):
    pdb_block = []
    indiv_block = []
    for i, line in enumerate(pdb):
        if i == 0:
            last_res = int(line[22:26])
        current_res = int(line[22:26])
        if current_res == last_res:
            indiv_block.append(line)
        else:
            pdb_block.append(indiv_block)
            indiv_block = []
            indiv_block.append(line)
            last_res = int(line[22:26])
    #append last block
    pdb_block.append(indiv_block)
    return pdb_block

def unblock_pdb_by_res(pdb_block):
    pdb = []
    for block in pdb_block:
        for line in block:
            pdb.append(line)
    return pdb

#must be a monomer to do this
def convert_to_pose_num(pdb):
    #current_res_nums = gen_residue_list(pdb)
    header = [x for x in pdb if x[0:6]  == "REMARK"]
    pdb = [x for x in pdb if x[0:6] == "ATOM  " or x[0:6] == "HETATM"]
    res_separated_blocks = block_pdb_by_res(pdb)
    fixed_blocks = []
    for i_block, block in enumerate(res_separated_blocks):
        new_block = []
        for i_line, line in enumerate(block):
            #force the line at residue position to change to the block index +1
            newline = ''.join([res_separated_blocks[i_block][i_line][:22], ('{0:>4}').format( "%i" % (i_block + 1) ), res_separated_blocks[i_block][i_line][26:]])
            new_block.append(newline)
        fixed_blocks.append(new_block)

    pdb = unblock_pdb_by_res(fixed_blocks)
    if not header:
        return pdb
    else:
        print header
        sys.exit()
        return header + pdb

def combine_all_chains(pdb):
    header = [x for x in pdb if x[0:6]  == "REMARK"]
    pdb = [x for x in pdb if x[0:6] == "ATOM  " or x[0:6] == "HETATM"]
    out_pdb = []
    for line in pdb:
        newline = ''.join([line[:21], 'A', line[22:]])
        out_pdb.append(newline)

    if not header:
        return out_pdb
    else:
        return header + out_pdb
    
#This definition should be used to check the geometry of new loops
#This definition works best determining stretches of atleast 4 residues
#seq_range takes list or string arguements
#It is not necessary but adding a few residues before or after is welcome
def geometric_peptide_angle_check(pdb, seq_range=''):

    pdb = [x for x in pdb if x.startswith("ATOM  ") or x.startswith("HETATM")]

    if type(seq_range) == str:
        seq_range = seq_range.split('-')
    if type(seq_range) == list or type(seq_range) == tuple:
        if len(seq_range) == 2:
            if int(seq_range[0]) < int(seq_range[1]):
                frag = [x for x in pdb if int(x[22:26]) in xrange(int(seq_range[0]), int(seq_range[1]))]
            else:
                frag = [x for x in pdb if int(x[22:26]) in xrange(int(seq_range[1]), int(seq_range[0]))]
        elif len(seq_range) >= 4:
            frag = [x for x in pdb if int(x[22:26]) in seq_range]
        else:
            print "Did not provide a range of residues or a list of residues greater than 3 residues"
            sys.exit()
    else:
        print "seq_range must be either str(#-#) or list[#,#]"
        sys.exit()

    #Need to check the information for the given fragment
    #ten pass_checks per residue

    pdb_res_blocks = block_pdb_by_res(frag)

    N_coords = [atom_coordinate(x, 'N', int(x[0][22:26]))[:-1] for x in pdb_res_blocks]
    CA_coords = [atom_coordinate(x, 'CA', int(x[0][22:26]))[:-1] for x in pdb_res_blocks]
    C_coords = [atom_coordinate(x, 'C', int(x[0][22:26]))[:-1] for x in pdb_res_blocks]

    geo_fail_count = 0
    
    for i in xrange(len(pdb_res_blocks)):
        #C-term angle A
        try: 
            ang = math.degrees(get_angle(CA_coords[i], C_coords[i], CA_coords[i+1]))
            if float(ang) < 141.5 or float(ang) > 150.5:  #141.5 150.5
                geo_fail_count += 1
        except IndexError: 
            continue

        try: 
            ang = math.degrees(get_angle(CA_coords[i], C_coords[i], CA_coords[i+2]))
            if float(ang) < 115.0 or float(ang) > 155.0: #110.0 160.0
                geo_fail_count += 1
        except IndexError: 
            continue

        try: 
            ang = math.degrees(get_angle(CA_coords[i], C_coords[i], CA_coords[i+3]))
            if float(ang) < 90.0 or float(ang) > 160.0: #80.0  170.0
                geo_fail_count += 1
        except IndexError: 
            continue

        #N-term angle A
        try: 
            ang = math.degrees(get_angle(CA_coords[i], N_coords[i+1], CA_coords[i+1]))
            if float(ang) < 152.0 or float(ang) > 159.0: #152.0  159.0
                geo_fail_count += 1
        except IndexError: 
            continue
        
        try: 
            ang = math.degrees(get_angle(CA_coords[i], N_coords[i+2], CA_coords[i+2]))
            if float(ang) < 125.0 or float(ang) > 150.0: #120.0  155.0
                geo_fail_count += 1
        except IndexError: 
            continue

        try: 
            ang = math.degrees(get_angle(CA_coords[i], N_coords[i+3], CA_coords[i+3]))
            if float(ang) < 100.0 or float(ang) > 160.0: #90.0  170.0
                geo_fail_count += 1
        except IndexError: 
            continue

        #C-term angle B
        try: 
            ang = math.degrees(get_angle(C_coords[i], CA_coords[i+1], CA_coords[i+2]))
            if float(ang) < 84.0 or float(ang) > 155.0: #75.0  160.0
                geo_fail_count += 1
        except IndexError:
            continue

        try:
            ang = math.degrees(get_angle(C_coords[i], CA_coords[i+2], CA_coords[i+3]))
            if float(ang) < 62.0 or float(ang) > 160.0: #60.0  163.0
                geo_fail_count += 1
        except IndexError:
            continue

        #N-term angle B
        try:
            ang = math.degrees(get_angle(CA_coords[i], CA_coords[i+1], N_coords[i+2]))
            if float(ang) < 80.0 or float(ang) > 155.0: #75.0  158.0
                geo_fail_count += 1
        except IndexError:
            continue

        try:
            ang = math.degrees(get_angle(CA_coords[i], CA_coords[i+1], N_coords[i+3]))
            if float(ang) < 60.0 or float(ang) > 160.0: # 57.0  162.0
                geo_fail_count += 1
        except IndexError:
            continue
    
    if geo_fail_count == 0:
        return False
    else:
        print "%s failures were predicted" % geo_fail_count
        return True

#exclude res should be a list composed of two arguements
#args are: [ [residue #, chain letter], [], [] ]
def clash_check(test, scaffold, mag_check=2.8, exclude_res=[]):
    back_bone = []
    test_obj = Transform(test)
    for line in scaffold:
        if [int(line[22:26]), str(line[21])] not in \
           [[int(x[0]),str(x[1])] for x in exclude_res]:
            if line[12:16].strip(' ') in ['CA', 'N', 'C', 'O']:
                back_bone.append(line)
    back_bone_obj = Transform(back_bone)
    s1 = time.clock()
    #print test_obj.xyz_total['3_A_THR_N']
    #print back_bone_obj.xyz_total['3_A_ASP_N']

    test_obj.get_geo_center()
    for name,array in test_obj.xyz_total.iteritems():
        test_obj.get_geo_center()
        s2 = time.clock()
        test_obj.translate([0.,0.,0.])
        #print test_obj.xyz_total[0][0]
        #print test_obj.xyz_total[1][1]
        #print test_obj.xyz_total[2][2]
        #print back_bone_obj.xyz_total[0]
        test_obj.rotate_object([0.,0.,0.],[1.,1.,1.],math.radians(30))
        test_obj.translate([0.,0.,0.])
        e2 = time.clock()
        print "1st", e2-s2
        #print test_obj.xyz_total[0]
        #print back_bone_obj.xyz_total[0]
        s3 = time.clock()
        test_obj.rotate_object([0.,0.,0.],[1.,1.,1.],math.radians(30))
        test_obj.translate([0.,0.,0.])
        e3 = time.clock()
        print "2nd", e3-s3
        s4 = time.clock()
        test_obj.rotate_object([0.,0.,0.],[1.,1.,1.],math.radians(30))
        test_obj.translate([0.,0.,0.])
        test_obj.write_to_file('Ummmmmm')
        e4 = time.clock()
        print "3rd", e4-s4
        for line2 in back_bone_obj.xyz_total:
            mag = sum((line2 - line)**2)
            #print line
            #print line2
            #print mag
            #mag =  (diff[0]**2)+(diff[1]**2)+(diff[2]**2)
            if mag < mag_check**2:
                return True #if clash se true
            else:
                continue
        if i == 1:
            e1 = time.clock()
            print "It took",e1-s1,"time to perform one line"
            print mag
            sys.exit()
    #no clashes found
    return False

#The proteins must be symmetric
def rough_symmetric_clash_check(test, test_symm, scaffold, scaffold_symm, mag_check=10.0):
    test_chains = symmetric_chain_separate(test, test_symm)
    scaffold_chains = symmetric_chain_separate(scaffold, scaffold_symm)
    test_points = []
    scaffold_points = []
    for chain in test_chains:
        chain_obj = Transform(chain)
        chain_center = chain_obj.get_geo_center() 
        test_points.append(chain_center)
    for chain in scaffold_chains:
        chain_obj = Transform(chain)
        chain_center = chain_obj.get_geo_center() 
        scaffold_points.append(chain_center)
   
    for line in test_points:
        for line2 in scaffold_points:
            diff = [ line[0] - line2[0], \
                     line[1] - line2[1], \
                     line[2] - line2[2] ]
            mag =  (diff[0]**2)+(diff[1]**2)+(diff[2]**2)
            if mag < mag_check**2:
                return True #if clash se true
            else:
                continue
    #no clashes found
    return False
 
def normalize(xyz, tolerance =0.00001):
    mag2 = sum(n * n for n in xyz)
    if abs(mag2 - 1.0) > tolerance:
        mag = math.sqrt(mag2)
        xyz = tuple(n / mag for n in xyz)
    return xyz

def cross(xyz, xyz2):

    xyz3 = [((xyz[1] * xyz2[2]) - (xyz[2] * xyz2[1])), \
            ((xyz[2] * xyz2[0]) - (xyz[0] * xyz2[2])), \
            ((xyz[0] * xyz2[1]) - (xyz[1] * xyz2[0]))]

    return xyz3

def q_mult(q1, q2):
    w1, x1, y1, z1 = q1
    w2, x2, y2, z2 = q2
    w = (w1 * w2) - (x1 * x2) - (y1 * y2) - (z1 * z2)
    x = (w1 * x2) + (x1 * w2) + (y1 * z2) - (z1 * y2)
    y = (w1 * y2) + (y1 * w2) + (z1 * x2) - (x1 * z2)
    z = (w1 * z2) + (z1 * w2) + (x1 * y2) - (y1 * x2)
    return w, x, y, z

def q_conjugate(q):
    q = normalize(q)
    w, x, y, z = q
    return (w, -x, -y, -z)

def q_xyz_mult(q1, xyz):
    v = normalize(xyz)
    p = (0.0, ) + v
    q2 = q_conjugate(q1)
    return q_mult(q_mult(q1, p), q2)[1:]

def axisangle_to_q(v, theta):
    v = normalize(v)
    x, y, z, = v
    theta /= 2
    w = math.cos(theta)
    x = x * math.sin(theta)
    y = y * math.sin(theta)
    z = z * math.sin(theta)
    return w, x, y, z


class Transform(object):    

    def __init__(self, pose):

        if type(pose) == list or type(pose) == tuple: 
            
            s1 = time.clock()
            self.pdb = pose
            if type(pose[0]) == str:
                self.xyz_total, self.xyz_total_dic = pdb_to_xyz(pose)
                e1 = time.clock()
                print "Took", e1-s1,"time to initialize xyz_total"
            else:
                self.xyz_total = numpy.array(pose)
            
            #self.xyz_atom_arrays = []
            #self.geo_center = numpy.zeros(1)
            #self.atom_mags = numpy.zeros(1)
        #In Progress
        '''else:
            self.pose = pose
            self.xyz_total = pose_to_xyz(pose)'''
    
    def update_vars(self):
        self.get_geo_center()
        self.get_atom_arrays()
        self.get_atom_mags()
        return
    def get_atom_arrays(self):
        self.xyz_atom_arrays = list(self.xyz_total)
        return self.xyz_atom_arrays

    def get_geo_center(self):
        arr = 0.0
        dic = 0.0
        for i in xrange(10000):
            s1 = time.clock()
            x = self.xyz_total[1000][0]
            e1 = time.clock()
            arr += (e1-s1)
            s2 = time.clock()
            x = self.xyz_total_dic[1000][0]
            e2 = time.clock()
            dic += (e2-s2)
        print arr
        print dic
        sys.exit()
        for line in self.xyz_total:
            val += line[0]
        '''for i,k in self.xyz_total_dic.iteritems():
            val += k[0]'''
        print val
        '''array = numpy.array([x for k,x in self.xyz_total.iteritems()], float)
        #array = numpy.array([self.xyz_total[k] for k in self.xyz_total.keys()], float)
        self.geo_center = numpy.mean(array, 0)'''
        e = time.clock()
        print "Took", e-s,"time to initialize geo_center"
        sys.exit()
        #sys.exit()
        return self.geo_center
    
    def get_atom_mags(self):
        if self.geo_center.all() == 0:
            self.get_geo_center()
        if not self.xyz_atom_arrays:
            self.get_atom_arrays()
        
        a_mags = []    
        for atom in self.xyz_atom_arrays:
            atom -= self.geo_center
            value = math.sqrt(sum(atom**2))
            a_mags.append(value)

        atom_mags = numpy.array(a_mags, "float128")[numpy.newaxis] #[math.sqrt(sum((atom - self.geo_center)**2)) for atom in self.xyz_atom_array])
        self.atom_mags = atom_mags.T
        return self.atom_mags
    
    #rotates object around a given bond designated by two atoms.
    #for object rotations at the origin supply atom1 as [0.0, 0.0, 0.0]
    def rotate_object(self, atom1, atom2, theta):
        s= time.clock()
        axis = [x-y for x,y in zip(atom2,atom1)]
        self.translate([0.0, 0.0, 0.0], atom1)
        rotated_vector = axisangle_to_q(axis, theta)
        if self.atom_mags.all() == 0:
            self.get_atom_mags()
        
        rotated_normals = numpy.array([q_xyz_mult(rotated_vector, x) for x in self.xyz_atom_arrays], "float128")
        self.xyz_total = self.atom_mags * rotated_normals
        #self.xyz_total = numpy.array([ y * numpy.array(q_xyz_mult(rotated_vector, x)) for x,y in zip(self.xyz_total, self.atom_mags)])
        
        #print self.xyz_total
        '''for i, atom in enumerate(self.xyz_total):
            
            if atom.all() != 0.0:
                
                mag =  math.sqrt(sum(atom**2))
                Rxyz = [x*mag for x in q_xyz_mult(rotated_vector, atom)]
                rotated_object.append(Rxyz)

        self.xyz_total = numpy.array(rotated_object)'''
        self.translate(atom1, [0.0, 0.0, 0.0])
        
        print self.xyz_total[0]
        e= time.clock()
        print "rotate", e-s, "seconds" 
        return self.xyz_total


    def translate(self, final, initial = []):
        s = time.clock() 
        if initial:
            translation = numpy.array(final,"float128") - numpy.array(initial, "float128")
        
        else:
            #if self.geo_center.all() == 0:
            self.get_geo_center()
            #print self.geo_center
            #print numpy.array(final, "float128")
            translation = numpy.array(final, "float128") - self.geo_center

        for k, v in self.xyz_total.iteritems():
            #v = translation + numpy.array(v)
            v = [v[0] + translation[0],\
                 v[1] + translation[1],\
                 v[2] + translation[2] ]
        e = time.clock()
        print "translate", e-s, "seconds"
        sys.exit() 
        return


    def generate_axes(self, symm = ''):
        #center each file with geometric center
        if self.geo_center.all() == 0:
            self.get_geo_center
        delta = self.xyz_total - self.geo_center

        #return principal axis matrix
        inertia = numpy.dot(delta.transpose(), delta)
        e_values, e_vectors = numpy.linalg.eig(inertia)

        #The axis we are concerned with at this time is the smallest eigan value (3) but writing in 1 and 2 for later
        for i in xrange(len(e_values)):
            if e_values[i] == max(e_values):
                eval1 = e_values[i]
                axis1 = e_vectors[:,i]
            elif e_values[i] == min(e_values):
                eval3 = e_values[i]
                axis3 = e_vectors[:,i]
            else:
                eval2 = e_values[i]
                axis2 = e_vectors[:,i]
        
        axis_coords = [axis1, axis2, axis3]
        
        if symm == 'D3':
            #axis1 and axis3 are correct vectors, need to rotate axis1 about axis3 two more times
           
            half_axis = [-1.0 *(axis1[0] + axis2[0]), \
                         -1.0 *(axis1[1] + axis2[1]), \
                         -1.0 *(axis1[2] + axis2[2])]
     

            axis_coords.append(half_axis)

        return axis_coords

    #only keeping this for old code (plan to remove soon)

    #def get_back_bone_

    def clash_check(self, test, mag_check=2.8, exclude_res=[]):
        back_bone = []
        test_obj = Transform(test)
        for line in scaffold:
            if [int(line[22:26]), str(line[21])] not in \
               [[int(x[0]),str(x[1])] for x in exclude_res]:
                if line[12:16].strip(' ') in ['CA', 'N', 'C', 'O']:
                    back_bone.append(line)
        back_bone_obj = Transform(back_bone)
        for line in test_obj.xyz_total:
            for line2 in back_bone_obj.xyz_total:
                diff = [ line[0] - line2[0], \
                         line[1] - line2[1], \
                         line[2] - line2[2] ]
                mag =  (diff[0]**2)+(diff[1]**2)+(diff[2]**2)
                if mag < mag_check**2:
                    return True #if clash se true
                else:
                    continue
        #no clashes found
        return False

    '''def get_terminal_residues(self):
        pdb = self.return_to_pdb()
        CA_xyz = []
        for ind, line in enumerate(pdb):
            if line[11:17].strip() == "CA":
                CA_xyz.append(self.xyz[ind])'''

    #obj1 will align to obj2
    def symmetric_axis_align(self, obj2, axis, axis2 = [], symm1='', symm2=''):
        #determin protein once
        if symm1 == 'D3':
            obj1_axes = self.generate_axes('D3')
        else:
            obj1_axes = self.generate_axes()
        print "obj1_axes are:", obj1_axes
        if not axis2:
            if symm2 == 'D3':
                obj2_axes = obj2.generate_axes('D3')
            else:
                obj2_axes = obj2.generate_axes()
            print "obj2_axes are:", obj2_axes

            ortho_vector = cross(obj1_axes[axis], obj2_axes[axis])
            theta = vector_angle(obj1_axes[axis], obj2_axes[axis])
        elif type(axis2) == int:
            if symm2 == 'D3':
                obj2_axes = obj2.generate_axes('D3')
            else:
                obj2_axes = obj2.generate_axes()
            print "obj2_axes are:", obj2_axes

            ortho_vector = cross(obj1_axes[axis], obj2_axes[axis2])
            theta = vector_angle(obj1_axes[axis], obj2_axes[axis2])
        
        else:
            print "obj2_axes are:", axis2
            ortho_vector = cross(obj1_axes[axis], axis2)
            theta = vector_angle(obj1_axes[axis], axis2)
            

        self.rotate_object([0.0,0.0,0.0], ortho_vector, theta)
       
    def align_object_to_vector(self, axis1, axis2, ortho_vector=[]):
        #determin protein once
        if ortho_vector:
            theta = vector_angle(axis1, axis2)
            self.rotate_object([0.0,0.0,0.0], ortho_vector, theta)
            
        else:
            ortho_vector = cross(axis1, axis2)
            theta = vector_angle(axis1, axis2)
            self.rotate_object([0.0,0.0,0.0], ortho_vector, theta)

    def return_to_pdb(self):
        return xyz_to_pdb(self.pdb, self.xyz_total)

    def write_to_file(self, path, chain=''):
        if not chain:
            with open(path, 'w') as myfile:
                myfile.writelines(self.return_to_pdb())
        else:
            pdb = self.return_to_pdb()
            pdb = [x for x in pdb if str(x[21]) == chain]
            with open(path, 'w') as myfile:
                myfile.writelines(pdb)
            
    
    #if any('MODEL' in s for s in old_pdb): might still be used
    #in place of the try except
    def append_to_file(self, path, model=False):
        if model == False:
            with open(path, 'a') as myfile:
                myfile.writelines(self.return_to_pdb())
        else:    
            try:
                old_pdb = read_file(path)

                try:
                    indices = [i for i, s in enumerate(old_pdb) if 'MODEL' in s] 
                    last_model_index = indices[-1]
                    last_model = int(old_pdb[last_model_index].split()[1])
                    print "%s contains MODELS, appending %s to file." % (path, last_model + 1)
                    output = ['MODEL  %s\n' % (last_model +1)] + self.return_to_pdb() + ['ENDMDL\n']
                except ValueError:
                    print "Value Error Raised, path submitted does not contain MODELS"
                    output = ['MODEL  1\n'] + old_pdb + ['ENDMDL\n'] + \
                             ['MODEL  2\n'] + self.return_to_pdb() + ['ENDMDL\n']

                    with open(path, 'w') as myfile:
                        myfile.writelines(output)
                    return

            except IOError:
                print "%s does not exist. Writing as new file." % path
                output = ['MODEL  1\n'] + self.return_to_pdb() + ['ENDMDL\n']
                with open(path, 'w') as myfile:
                    myfile.writelines(output)
                return
                
            with open(path, 'a') as myfile:
                myfile.writelines(output)

class Transform_vector(Transform):
    
    def __init__(self, pose):
        
        if len(pose) == 3 and type(pose[0]) != list:
            
            Transform.__init__(self, pose)
            self.xyz = self.xyz_total
            self.mag = math.sqrt(sum(self.xyz**2))

    #Rotate a single point about an axis
    #this is assuming axis is coming from origin but can be given an origin
    #for larger objects it is best to use rotate_object
    def rotate_vector(self, axis, theta, initial=[]):
        
        rotated_vector = axisangle_to_q(axis, theta)
        Rxyz = [x*self.mag for x in q_xyz_mult(rotated_vector, self.xyz)]

        return Rxyz  
