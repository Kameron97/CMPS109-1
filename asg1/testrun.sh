#!/bin/sh -x
# $Id: testrun.sh,v 1.1 2017-08-21 17:29:41-07 - - $


./ydc <qa/input/test0-stk.ydc >test0-stk.ydc.output 2>&1
./ydc <qa/input/test1-stk.ydc >test1-stk.ydc.output 2>&1
./ydc <qa/input/test2-add.ydc >test2-add.ydc.output 2>&1
./ydc <qa/input/test3-sub.ydc >test3-sub.ydc.output 2>&1
./ydc <qa/input/test4-mul.ydc >test4-mul.ydc.output 2>&1
./ydc <qa/input/test5-div.ydc >test5-div.ydc.output 2>&1
./ydc <qa/input/test6-exp.ydc >test6-exp.ydc.output 2>&1
./ydc <qa/input/testx-leaks.ydc >testx-leaks.ydc.output 2>&1

mv *.output myoutput/
 diff -qr /afs/cats.ucsc.edu/users/w/apatnam/CMPS109/asg1/myoutput/ /afs/cats.ucsc.edu/users/w/apatnam/CMPS109/asg1/qa/output/
#diff /afs/cats.ucsc.edu/users/w/apatnam/CMPS109/asg1/myoutput/test0-stk.ydc.output /afs/cats.ucsc.edu/users/w/apatnam/CMPS109/asg1/qa/output/test0-stk.ydc.output
#diff /afs/cats.ucsc.edu/users/w/apatnam/CMPS109/asg1/myoutput/test1-stk.ydc.output /afs/cats.ucsc.edu/users/w/apatnam/CMPS109/asg1/qa/output/test1-stk.ydc.output
#diff /afs/cats.ucsc.edu/users/w/apatnam/CMPS109/asg1/myoutput/test2-add.ydc.output /afs/cats.ucsc.edu/users/w/apatnam/CMPS109/asg1/qa/output/test2-add.ydc.output
#diff /afs/cats.ucsc.edu/users/w/apatnam/CMPS109/asg1/myoutput/test3-sub.ydc.output /afs/cats.ucsc.edu/users/w/apatnam/CMPS109/asg1/qa/output/test3-sub.ydc.output
#diff /afs/cats.ucsc.edu/users/w/apatnam/CMPS109/asg1/myoutput/test4-mul.ydc.output /afs/cats.ucsc.edu/users/w/apatnam/CMPS109/asg1/qa/output/test4-mul.ydc.output
#diff /afs/cats.ucsc.edu/users/w/apatnam/CMPS109/asg1/myoutput/test5-div.ydc.output /afs/cats.ucsc.edu/users/w/apatnam/CMPS109/asg1/qa/output/test5-div.ydc.output
#diff /afs/cats.ucsc.edu/users/w/apatnam/CMPS109/asg1/myoutput/test6-exp.ydc.output /afs/cats.ucsc.edu/users/w/apatnam/CMPS109/asg1/qa/output/test6-exp.ydc.output

