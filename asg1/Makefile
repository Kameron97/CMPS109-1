# Aravind Patnam
# apatnam@ucsc.edu
# ASG1 - DC Calc
# January 22nd, 2018


# $Id: Makefile,v 1.11 2014-07-02 20:01:17-07 - - $

MKFILE      = Makefile
DEPFILE     = ${MKFILE}.dep
NOINCL      = ci clean spotless
NEEDINCL    = ${filter ${NOINCL}, ${MAKECMDGOALS}}
GMAKE       = ${MAKE} --no-print-directory

COMPILECPP  = g++ -g -O0 -Wall -Wextra -std=gnu++11
MAKEDEPCPP  = g++ -MM

CPPHEADER = bigint.h scanner.h debug.h util.h iterstack.h ubigint.h
CPPSOURCE=bigint.cpp scanner.cpp debug.cpp util.cpp main.cpp ubigint.cpp
EXECBIN     = ydc
OBJECTS     = ${CPPSOURCE:.cpp=.o}
OTHERS      = ${MKFILE} README
ALLSOURCES  = ${CPPHEADER} ${CPPSOURCE} ${OTHERS}

all : ${EXECBIN}
	- checksource ${ALLSOURCES}

${EXECBIN} : ${OBJECTS}
	${COMPILECPP} -o $@ ${OBJECTS}

%.o : %.cpp
	${COMPILECPP} -c $<

ci : ${ALLSOURCES}
	- checksource ${ALLSOURCES}
	cid + ${ALLSOURCES}

clean :
	- rm ${OBJECTS} ${DEPFILE} 

spotless : clean
	- rm ${EXECBIN}

dep : ${CPPSOURCE} ${CPPHEADER}
	@ echo "# ${DEPFILE} created `LC_TIME=C date`" >${DEPFILE}
	${MAKEDEPCPP} ${CPPSOURCE} >>${DEPFILE}

submit :
	submit cmps109-wm.w15 asg2 ${ALLSOURCES}
	rm -rf latest_submit
	mkdir latest_submit
	cp ${ALLSOURCES} .score/* latest_submit

${DEPFILE} :
	@ touch ${DEPFILE}
	${GMAKE} dep

again :
	${GMAKE} spotless dep ci all lis

ifeq (${NEEDINCL}, )
include ${DEPFILE}
endif

