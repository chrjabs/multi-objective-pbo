SYSTEM     = x86-64_linux
LIBFORMAT  = static_pic
CPLEXDIR      =/opt/ibm/ILOG/CPLEX_Studio201/cplex
CONCERTDIR    =/opt/ibm/ILOG/CPLEX_Studio201/concert
# ---------------------------------------------------------------------
# Compiler options
# ---------------------------------------------------------------------

CCOPT = -m64 -fPIC -fno-strict-aliasing -fexceptions

# ---------------------------------------------------------------------
# Link options and libraries
# ---------------------------------------------------------------------

CPLEXBINDIR   =$(CPLEXDIR)/bin/$(BINDIST)
CPLEXLIBDIR   =$(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CONCERTLIBDIR =$(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)

CPLEXBINDIR	=$(CPLEXDIR)/bin/$(SYSTEM)
CPLEXLIB	=cplex
CCLNDIRS	=-L$(CPLEXLIBDIR) -L$(CONCERTLIBDIR)
LFLAGS += $(CCLNDIRS) -lconcert -lilocplex -l$(CPLEXLIB) -lm -lpthread -ldl

# For dynamic linking
ifeq ($(dynamic), True)
$(info dynamic linking activated)
CPLEXLIB :=$(CPLEXLIB)2010
CCLNDIRS :=$(CCLNDIRS) -L$(CPLEXBINDIR)
else
$(info static linking, by default)
endif

CONCERTINCDIR = $(CONCERTDIR)/include
CPLEXINCDIR   = $(CPLEXDIR)/include

CFLAGS += $(CCOPT) -I$(CPLEXINCDIR) -I$(CONCERTINCDIR) 
$(firstword $(MAKEFILE_LIST)): cplex.make
	@if [ -z $(CPLEXDIR) ] || [ -z $(CPLEXDIR) ]; then \
		echo "ERROR: empty CPLEX configuration variables"; \
	exit 1; \
	fi

	@if [ ! -d $(CPLEXDIR) ]; then \
		echo "ERROR: configure CPLEXDIR: $(CPLEXDIR) does not exist"; \
	fi
	@if [ ! -d $(CONCERTDIR) ]; then \
		echo "ERROR: configure CONCERTDIR, edit @^: $(CONCERTDIR) does not exist"; \
	fi
	@if [ ! -z $(CPLEXDIR) ] && [ -d $(CPLEXDIR) ] && [ ! -z $(CPLEXDIR) ] && [ -d $(CONCERTDIR) ]; then \
		touch $@; \
	else \
	echo "Please make sure IBM ILOG CPLEX Optimization Studio is installed";\
	echo "Aborting now.";\
	exit 1;\
	fi
