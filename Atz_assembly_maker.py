#!/usr/bin/env python
from transform import *
from AssemB_pkgs import *
import os, sys, math, numpy, time, csv, optparse, collections

## Attempts to rotate cofactor to see if LASA matches scaff atomtype
## Outputs rigid body rotated cofactor
#IN PROGRESS** channge all CB atoms to atom_type eventually
def get_dn_c2_axes(axes_dict):
    clash_chains = set()
    trans_axes = {}
    for key in axes_dict.keys():
        if len(key) > 1 and key[1] == '_':
            chains = key.split('_')
            if 'A' in chains or '1' in chains:
                clash_chains.add(chains[0])
                clash_chains.add(chains[1])
                val = axes_dict[key]
                trans_axes[str(key)] = val
    return clash_chains, trans_axes

def assem_maker(D3, \
                d3_filename, \
                D2, \
                d2_filename):

    Two_comp_obj = D3_D2_AssemMaker(d3_filename, d2_filename)
    Two_comp_obj.generate_c2_axes_positional_ensemble(360.0, 60.0, 90.0, 150.0, 5.0)
    for name, my_obj in collections.OrderedDict(sorted(Two_comp_obj.ensemble.items())).iteritems():
        if name.split('_')[-1].startswith('trans'):
            #my_obj.write_to_file('%s.pdb' % name)
            my_obj.append_to_file('sagars_video.pdb', True)
        #if len(name.split('_')) > 3:
            #my_obj.append_to_file('sagars_video.pdb', True)
        #write_arbPDB_from_xyz([x for i,x in my_obj.get_symmetric_axes().iteritems()], '%s.axis' % name)
    sys.exit()
    ## If the pdb file is input with a path and '/' separate and grab last item (the file)
    d3_filename = d3_filename.split('/')[-1]
    
    ## If the last four characters are '.pdb' save filename without '.pdb'
    if d3_filename[-4:] == '.pdb':
        d3_filename = d3_filename[:-4]
    print d3_filename

    ## Create a transform object from our scaffold
    d3_obj = Transform(D3)
    
    ## Translate objects geometric center to a zero origin
    d3_obj.translate([0.0,0.0,0.0])
    d3_axes = d3_obj.get_symmetric_axes()
    d3_clash_chains, d3_trans_axes = get_dn_c2_axes(d3_axes)

    d3_backbone = TransformCopy(d3_obj)
    d3_backbone.assign(d3_obj)
    d3_backbone.get_xyz_by_chain(list(d3_clash_chains), True)
    d3_backbone.get_xyz_by_atomtype(['N','CA','C','O','OXT'], True)
    
    ## Copy above three steps for the input cofactor
    d2_obj = Transform(D2)
    d2_obj.translate([0.0,0.0,0.0])
    d2_axes = d2_obj.get_symmetric_axes()
    print d2_axes
    sys.exit()
    d2_clash_chains, d2_trans_axes = get_dn_c2_axes(d2_axes)
    
    def major_align(d2_obj, d2_axis, d3_axis):
        d2_obj_variant1 = TransformCopy(d2_obj)
        d2_obj_variant1.assign(d2_obj)
        d2_obj_variant1.align_object_to_vector(min_and_max_d2[0], list(d3_axes['major1']))

    #d2_obj.write_to_file('d2_centered.pdb')
    #print d2_axes
    #sys.exit()
    #d2_obj.major_symmetric_axis_align(d3_obj)
    '''
    d3_obj.write_to_file('check1_d3.pdb')
    d2_obj.write_to_file('check2_d2.pdb')
    sys.exit()
    '''
    min_and_max_d2 = [x for i,x in d2_trans_axes.iteritems()]
    print min_and_max_d2
    meh = create_arbPDB_from_xyz(min_and_max_d2)
    write_file(meh, 'hihi.pdb')
    #sys.exit()
    min_and_max_d2_names = [i for i,x in d2_trans_axes.iteritems()]

    d2_variant_set = {}
    ## Prealign the d2 axis maximum (0) of symmetry to the d3 c2 axis of symmetry
    for d3c2_axis_name, d3c2_axis_coords in d3_trans_axes.iteritems():
        d2_obj_variant1 = TransformCopy(d2_obj)
        d2_obj_variant1.assign(d2_obj)
        d2_obj_variant1.align_object_to_vector(min_and_max_d2[0], list(d3_axes['major1']))
        d2_obj_variant1.align_object_to_vector(min_and_max_d2[1], d3c2_axis_coords, list(d3_axes['major1']))
        d2_variant_set['%s_with_%s_final' % (d3c2_axis_name, min_and_max_d2_names[0])] = d2_obj_variant1
        
        d2_obj_variant2 = TransformCopy(d2_obj)
        d2_obj_variant2.assign(d2_obj)
        d2_obj_variant2.align_object_to_vector(min_and_max_d2[1], list(d3_axes['major1']), list(min_and_max_d2[0]))
        d2_obj_variant2.align_object_to_vector(min_and_max_d2[0], d3c2_axis_coords, list(d3_axes['major1']))
        #d2_variant_set.append(d2_obj_variant)
        d2_variant_set['%s_with_%s' % (d3c2_axis_name, min_and_max_d2_names[1])] = d2_obj_variant2
    
    print len(d2_variant_set.keys())
    for name, d2_var in d2_variant_set.iteritems():
        d2_var.write_to_file('%s.pdb' % name)
    sys.exit()
    
    ## Regenerate the new axes of symmetry of the rotated cofactor
    rotate_axes = d2_obj.generate_axes()


    rotational_axis = rotate_axes[0]
    mag = get_mag(rotate_axes[0], [0.0,0.0,0.0])
    print mag
    #sys.exit()

    translational_axis = tuple(n / mag for n in rotate_axes[0])
    
    for ang in [0.0, 54.75, 90.0, 125.25, 180, 234.75, 270, 305.25]: #[0, 10, 40, 340, 350]: #[120,130,140,150,170,180,190,200]: # #xrange(0, 360, 10):
        model = TransformCopy(d2_obj)
        model.assign(d2_obj)
        model.rotate_object([0.0,0.0,0.0], rotational_axis, math.radians(ang))
        for x in xrange(100,170): #[135]: #[126,129,132,135]: #xrange(108, 150, 3): #range(117, 159, 3): positive ones
            pos_trans_model = TransformCopy(model)
            pos_trans_model.assign(model)
            neg_trans_model = TransformCopy(model)
            neg_trans_model.assign(model)

            neg_trans_model.rotate_object([0.0,0.0,0.0], rotate_axes[2], math.radians(180))
            
            #print "trans axis:", translational_axis
            print x

            pos_trans_model.translate([x * pt for pt in translational_axis], pos_trans_model.get_geo_center())
            neg_trans_model.translate([-1.0 * x * nt for nt in translational_axis], neg_trans_model.get_geo_center())

            if pos_trans_model.clash_check(d3_backbone) == False:
                pos_trans_model.write_to_file('Pos_ang_%s_dist_%s.pdb' % (str(ang).replace('.','_'), x), 'A')
                
            if neg_trans_model.clash_check(d3_backbone) == False:
                neg_trans_model.write_to_file('Neg_ang_%s_dist_%s.pdb' % (str(ang).replace('.','_'), x), 'A')
                
    
    return

################################################################################
#
#
################################################################################

def main(argv):
    
    parser = optparse.OptionParser(usage="\n\nGenerate csv of symmetric cofactor matches.")
    
    parser.add_option('--d3', dest = 'd3_filename',
        help = 'The input pdb')
    
    parser.add_option('--d2', dest = 'd2_filename',
        help = 'The input cofactor')
    
    parser.add_option('--s-units', type="int", nargs=1, dest = 'subunits',
        default = 2,
        help = 'Number of symmetric subunits ( default = 2 )')

    (options,args) = parser.parse_args()
    
    
    D3 = read_file(str(options.d3_filename))
    D2 = read_file(str(options.d2_filename))

    assem_maker(D3, options.d3_filename, D2, options.d2_filename)

    
if __name__ == "__main__":
    main(sys.argv[1:])

