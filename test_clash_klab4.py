from testy import *
import sys, time

total_time = 0.
obj1 = Transform(read_file(sys.argv[1]))
obj2 = Transform(read_file(sys.argv[2]))
for i in xrange(3):
    start = time.clock()
    print "Starting %sth round" % i
    for i2 in xrange(5):
        obj1.translate([2.,2.,2.])
        obj1.rotate_object([0.,0.,0.],[1.,1.,1.],math.radians(30))
        obj2.translate([4.,4.,4.])
        obj2.rotate_object([0.,0.,0.],[1.,1.,1.],math.radians(30))

    #clash_check(read_file(sys.argv[1]),read_file(sys.argv[2]),3.0)

    end = time.clock()
    print "%sth round took: %f seconds" % (i, (end-start))
    total_time += (end - start)

print "Avg time was", total_time/3.0, "seconds to complete this clash check"
