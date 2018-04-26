#!/bin/usr/env python

import math, sys, os, time, fork_manager
from SyPRIS_pkgs import *

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


def run_SyPRIS(scaff_list, num_cpu):

    print "Preparing to run SyPRIS on", len(scaff_list), "files"
    
    ffm = FileForkManager()
    fm = fork_manager.ForkManager(num_cpu)
    fm.error_callback = ffm.handle_failed_run
    fm.success_callback = ffm.handle_successful_run

    for line in scaff_list:
        line = line.rstrip('\r\n')
        print line
        split_line = line.split('_')
        if split_line[0] == '!':
            continue
        if split_line[1] == 'v1' or split_line[1] == 'v2':
            cof = './symmed_cofs/' + ('_'.join(split_line[:2])).lower() + '_cof_symmed.pdb'
        else:
            cof = './symmed_cofs/' + (split_line[0]).lower() + '_cof_symmed.pdb'

        pid = fm.mfork()
        if pid == 0:
            hit = syp_pkg_tester(cof, line)
    
            sys.exit(0)
        else:
            ffm.file_for_job[pid] = line
    fm.wait_for_remaining_jobs()
    return ffm

def exit_SyPRIS( ffm ):
    if ffm.all_files_ran:
        sys.exit(0)
    else:
        for file_name in ffm.files_that_failed:
            print "File", file_name, "did not finish"
        sys.exit(1)


def syp_pkg_tester(cof_file, scaff_file):
    cof = CsymmCofactor(cof_file, scaff_file)
    cof.std_mult = 2.0
    #cof.use_input_rotamer = True
    cof.origin_align(flip=True)
    #cof.lock_rotamer_addition = True
    cof.generate_inverse_rot_ensemble(std_periodicity=2, lasa='CB', my_rotlib='./rotlib_cont.txt')
    cof.inverse_rotamer_sampler(primary_score=1.0,secondary_score=1.0)
    #print cof.match_data_set

    if cof.match_data_set:
        name_scaff_list = scaff_file.split('/')
        scaff_file = name_scaff_list[-1]
        out_lines = []
        for k,v in cof.match_cofactors.iteritems():
            pre_exclude_res_list = [x[:-3] for x in k.split('_')[::2]]
            res_nums = [int(x) for x in pre_exclude_res_list]
            exclude_res_list = set()
            for res_number in res_nums:
                exclude_res_list.add(str(res_number))
                exclude_res_list.add(str(res_number+1))
                exclude_res_list.add(str(res_number-1))
            included_res = set()
            scaffy = copy_transform_object(cof.scaffold)
            scaffy.get_xyz_by_atomtype(['C','N','CA','O'], update=True)
            for name in scaffy.xyz_names:
                name_list = name.split('_')
                res_name = name_list[0]
                if res_name not in exclude_res_list:
                    included_res.add(res_name)
            scaffy.get_xyz_by_resnum(included_res, True)
            if v.clash_check(test_scaff=scaffy, include_resnum=list(included_res)) == False:
                out_line = [cof_file.split('_')[0], k, scaff_file]
                sep_name = k.split('_')
                res_to_check = []
                if len(sep_name) == 4:
                    res_to_check.append('_'.join(sep_name[:2]))
                    res_to_check.append('_'.join(sep_name[2:]))
                else:
                    res_to_check.append(k) 
                #for k2,v2 in cof.match_data_set[k].iteritems():
                for k2 in res_to_check:
                    rmsa = cof.match_data_set[k][k2]['rmsa']
                    rmsd = cof.match_data_set[k][k2]['rmsd']
                    score = cof.match_data_set[k][k2]['score']
                    tors = []
                    for t in cof.match_data_set[k][k2]['torsions'][:3]:
                        if t == 0.0:
                            tors.append('X')
                        else:
                            tors.append(t)
                    scaff_resi = k2.split('_')[0][:-3]
                    cof_resi = ''
                    for char in k2.split('_')[1]:
                        try:
                            int(char)
                            cof_resi += char
                        except ValueError:
                            break
                    if scaff_resi == cof_resi:
                        scaff_chis = get_first_three_chis_of_resnum(copy_transform_object(cof.scaffold), scaff_resi)
                        #write to csv instead of file

                        out_line += scaff_chis[:] + tors[:] + [round(rmsa,5), round(rmsd,5), round(score,5)]
                    else:
                        out_line += ['NA','NA','NA'] + tors[:] + [round(rmsa,5), round(rmsd,5), round(score,5)]
                    v.write_to_file('sub_%s_%s_rmsa_%s_rmsd_%s_score_%s.pdb' % ( scaff_file[:-4], k, round(rmsa,10), round(rmsd, 10), round(score, 10) ) )
                #cof.scaffold.write_to_file('%s_syp_centered.pdb' % scaff_file[:-4])
                out_lines.append(out_line)
        if out_lines:
            with open('sub_%s.csv' % out_lines[0][0], 'a') as csvfile:
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
    sys.exit()


def main():
    args = sys.argv
    scaff_file_list = args[1]
    jobs = args[2]

    with open(scaff_file_list, 'r') as mf:
        scaff_list = mf.readlines()
    ffm = run_SyPRIS( scaff_list, jobs )
if __name__ == "__main__":
    s = time.clock()
    main()
    print 'Took', time.clock() - s, 'seconds.'
