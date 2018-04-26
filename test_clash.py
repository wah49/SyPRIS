from transform import *
import sys, time

total_time = 0.
obj1 = Transform(read_file(sys.argv[1]))
obj2 = Transform(read_file(sys.argv[2]))
#obj3 = TransformCopy(obj2)
#obj3.assign(obj2)
#obj3.get_xyz_by_atomtype(['N','CA','C','O'], True)
#obj4 = TransformCopy(obj1)
#obj4.assign(obj1)
#obj4.get_xyz_by_atomtype(['N','CA','C','O'], True)
#print len(obj2.xyz_total)
#print obj2.xyz_total[:7]
#print obj2.xyz_names[:7]
#
#print len(obj3.xyz_total)
#print obj3.xyz_total[:7]
#print obj3.xyz_names[:7]
axes1 = obj1.get_symmetric_axes()
axes2 = obj2.get_symmetric_axes()
out_1 = create_arbPDB_from_xyz([x for i,x in axes1.iteritems()])
out_2 = create_arbPDB_from_xyz([x for i,x in axes2.iteritems()])
write_file(out_1, 'lala.pdb')
write_file(out_2, 'toto.pdb')
sys.exit()
for i in xrange(100):
    start = time.clock()
    print "Starting %sth round" % i
    value = obj1.clash_check(obj2)
    #value = clash_check(obj4, obj3)
    print value
    #clash_check(read_file(sys.argv[1]),read_file(sys.argv[2]),3.0)

    end = time.clock()
    print "%sth round took: %f seconds" % (i, (end-start))
    total_time += (end - start)

print "Avg time was", total_time/100.0, "seconds to complete this clash check"
