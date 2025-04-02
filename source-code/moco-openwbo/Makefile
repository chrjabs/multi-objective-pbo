# Define which solver to use as backend, this can be a name of a file in the
# solvers directory.
SOLVER     ?= glucose4.1
default: build
# The following values should be defined in the included file:
# VERSION    = core or simp 
# SOLVERNAME = name of the SAT solver
# SOLVERDIR  = subdirectory of the SAT solver
# NSPACE     = namespace of the SAT solver
#

# settings pertaining to the cplex libraries. You can sever the
# dependency by calling make like so:
#  make CPLEX_BOUNDS_DEPENDENCY=False 
ifeq ($(cplex_bounds_dependency),True)
include cplex.make
CFLAGS += -DCPLEX_BOUNDS_DEPENDENCY
else
CPLEX_O=$(realpath ./solvers/$(SOLVER))/../../bounds/BoundsCalculatorCPLEX.o
# create empty object value of each type
$(CPLEX_O:.o=.o) $(CPLEX_O:.o=.or) $(CPLEX_O:.o=.od):
	touch $@
endif

# dependencies
build: build_deps build_maxpre s

deps := $(wildcard ./deps/*)

build_deps: 
	@for dep in $(deps);\
		 do \
			 make -C $$dep build;  done
clean_deps:
	@for dep in $(deps);\
		 do \
			 make -C $$dep clean;  done
	
build_maxpre:
	$(MAKE) -C maxpre lib
	
clean_maxpre:
	$(MAKE) -C maxpre clean
	
s: build_maxpre
p: build_maxpre
d: build_maxpre
r: build_maxpre
rs: build_maxpre

CFLAGS += -DPREPRO -I maxpre/src
LFLAGS += -lmaxpre -L maxpre/src/lib

touch_cplex:
	touch cplex.make
clean: touch_cplex
clean_total: clean_deps clean clean_maxpre
# solver
include $(PWD)/solvers/$(SOLVER).mk
# THE REMAINING OF THE MAKEFILE SHOULD BE LEFT UNCHANGED
EXEC       = open-wbo
DEPDIR     += mtl utils core
DEPDIR     +=  ../../encodings ../../algorithms ../../graph ../../classifier ../../bounds
MROOT      ?= $(PWD)/solvers/$(SOLVERDIR)
LFLAGS     += -lgmpxx -lgmp
CFLAGS     += -Wall -Wno-parentheses -Wno-class-memaccess -std=c++17 -DNSPACE=$(NSPACE) -DSOLVERNAME=$(SOLVERNAME) -DVERSION=$(VERSION)
ifeq ($(SANITIZER),asan)
CFLAGS     += -fsanitize=address
LFLAGS     += -fsanitize=address
LFLAGS     += -fuse-ld=gold
endif
ifeq ($(SANITIZER),undef)
CFLAGS     += -fsanitize=undefined -fsanitize-undefined-trap-on-error
LFLAGS     += -fsanitize=undefined
LFLAGS     += -fuse-ld=gold
endif
ifeq ($(VERSION),simp)
DEPDIR     += simp
CFLAGS     += -DSIMP=1 
ifeq ($(SOLVERDIR),glucored)
LFLAGS     += -pthread
CFLAGS     += -DGLUCORED
DEPDIR     += reducer glucored
endif
endif

# Some solvers do not have a template.mk file any more
# E.g.: Minisat or Riss
ifeq ($(SOLVERDIR),$(filter $(SOLVERDIR),minisat riss))
include $(PWD)/mtl/template.mk
else
include $(MROOT)/mtl/template.mk
endif
