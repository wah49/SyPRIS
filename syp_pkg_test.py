#!/bin/usr/env python

import math, sys, os, time, fork_manager
from SyPRIS_pkgs import *


def syp_pkg_tester(cof_file, scaff_file):
    cof = CsymmCofactor(cof_file, scaff_file)
    cof.std_mult = 2.0
    #cof.use_input_rotamer = True
    cof.origin_align(flip=True)
    #cof.lock_rotamer_addition = True
    cof.generate_inverse_rot_ensemble(std_periodicity=2, lasa='CB', my_rotlib='./rotlib.txt')
    '''for t,y in cof.bb_rot_ensemble.iteritems():
        y.write_to_file('%s.pdb' % t)'''
    cof.inverse_rotamer_sampler(primary_score=2.0,secondary_score=0.0)
    #print cof.match_data_set


    if cof.match_data_set:
        name_scaff_list = scaff_file.split('/')
        scaff_file = name_scaff_list[-1]
        out_lines = []
        for k,v in cof.match_cofactors.iteritems():
            exclude_res_list = [x[:-3] for x in k.split('_')[::2]]
            included_res = set()
            scaffy = copy_transform_object(cof.scaffold)
            scaffy.get_xyz_by_atomtype(['C','N','CA','O'], update=True)
            for name in scaffy.xyz_names:
                name_list = name.split('_')
                res_name = name_list[0]
                if res_name not in exclude_res_list:
                    included_res.add(res_name)
            if v.clash_check(test_scaff=scaffy, include_resnum=list(included_res)) == False:
                out_line = [cof_file.split('_cof.pdb')[0][5:], k]
                sep_name = k.split('_')
                res_to_check = []
                if len(sep_name) == 4:
                    res_to_check.append('_'.join(sep_name[:2]))
                    res_to_check.append('_'.join(sep_name[2:]))
                else:
                    res_to_check.append(k)

                for k2 in res_to_check:
                    alog = cof.match_data_set[k][k2]['alog']
                    rmsd = cof.match_data_set[k][k2]['rmsd']
                    score = cof.match_data_set[k][k2]['score']
                    tors = []
                    for t in cof.match_data_set[k][k2]['torsions'][:3]:
                        if t == 0.0:
                            tors.append('X')
                        else:
                            tors.append(t)
                    scaff_res = k2.split('_')[0][:-3]
                    scaff_chis = get_first_three_chis_of_resnum(copy_transform_object(cof.scaffold), scaff_res)

                    out_line += scaff_chis[:] + tors[:] + [round(alog,5), round(rmsd,5), round(score,5)]
                    v.write_to_file('%s_%s_alog_%s_rmsd_%s_score_%s.pdb' % ( scaff_file[:-4], k, round(alog,10), round(rmsd, 10), round(score,10) ) )
                out_lines.append(out_line)
                #cof.scaffold.write_to_file('%s_syp_centered.pdb' % scaff_file[:-4])
        if out_lines:
            with open('%s.csv' % out_lines[0][0], 'a') as csvfile:
                writer = csv.writer(csvfile)
                for line in out_lines:
                    writer.writerow( line )

        del cof
        return True
    else:
        del cof
        return False



        '''cof_lines = read_file(cof_file)
        res_nums = gen_residue_list(cof_lines)
        scaff = Transform(scaff_file)
        for res_num in res_nums: 
            copy_scaff = copy_transform_object(scaff)
            copy_scaff.get_xyz_by_resnum(res_num, update=True)'''


    #print cof.match_data_set
    #print '\n'
    #print cof.match_cofactors


def main():
    args = sys.argv
    scaff_file_list = args[1]
    s = time.clock()
    counter = 0
    pass_pdbs = []
    with open(scaff_file_list, 'r') as mf:
        scaff_list = mf.readlines()
    
    for line in scaff_list:
        line = line.rstrip('\r\n')
        print line
        split_line = line.split('_')
        if split_line[0] == '!':
            continue
        if split_line[1] == 'v1' or split_line[1] == 'v2':
            cof = 'cofs/' + '_'.join(split_line[:2]) + '_cof.pdb'
        else:
            cof = 'cofs/' + split_line[0] + '_cof.pdb'

        hit = syp_pkg_tester(cof, line)
        if hit:
            counter += 1
            pass_pdbs.append(line)
    
    print 'Took', time.clock() - s, 'seconds. %s of 78 worked.' % counter
    for line in pass_pdbs:
        print line
if __name__ == "__main__":
    main()
