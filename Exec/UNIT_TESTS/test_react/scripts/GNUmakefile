NDEBUG := t
MPI    :=
OMP    :=
MKVERBOSE := t 
COMP := gfortran

# include the main Makefile stuff
include $(AMREX_HOME)/Tools/F_mk/GMakedefs.mak

# to make python libraries, we need to do -fPIC (in gfortran)
F90FLAGS += -fPIC
FFLAGS += -fPIC
CFLAGS += -fPIC

# core AMReX directories
AMREX_CORE := Src/F_BaseLib

Fmpack += $(foreach dir, $(AMREX_CORE), $(AMREX_HOME)/$(dir)/GPackage.mak)
Fmlocs += $(foreach dir, $(AMREX_CORE), $(AMREX_HOME)/$(dir))

# include the necessary GPackage.mak files that define this setup
include $(Fmpack)

# vpath defines the directories to search for the source files
VPATH_LOCATIONS += $(Fmlocs)

programs += fsnapshot

all: python_module

python_module: $(objects)
	f2py --fcompiler=gfortran --f90flags="-J t/Linux.gfortran/m/" -c fsnapshot.f90 -m fsnapshot $(objects)

include $(AMREX_HOME)/Tools/F_mk/GMakerules.mak

clean::
	$(RM) fsnapshot.so


