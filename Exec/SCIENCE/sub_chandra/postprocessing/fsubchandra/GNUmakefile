NDEBUG := t
MPI    :=
OMP    :=

MKVERBOSE := t 

COMP := gfortran

# some routines need an eos/network (i.e. to compute thermodynamic 
# quantities.  If that is the case, set NEED_EOS_NETWORK := t
NEED_EOS_NETWORK := t

# include linear algebra packages
NEED_BLAS := t
NEED_LINPACK := t

# define the location of the MAESTRO root directory
MAESTRO_TOP_DIR := $(MAESTRO_HOME)

# include the main Makefile stuff
include $(BOXLIB_HOME)/Tools/F_mk/GMakedefs.mak

# core BoxLib directories
BOXLIB_CORE := Src/F_BaseLib

# other packages needed for data_processing
Fmdirs := 

# directories containing files that are 'include'-d via Fortran
Fmincludes :=

ifdef NEED_EOS_NETWORK
  EOS_TOP_DIR := $(MAESTRO_TOP_DIR)/Microphysics/EOS
  EOS_DIR := helmeos

  NETWORK_TOP_DIR := $(ASTRODEV_DIR)/networks
  NETWORK_DIR := triple_alpha_plus_cago


  MICROPHYS_CORE := $(MAESTRO_TOP_DIR)/Microphysics/EOS

  MICROPHYS_CORE += $(EOS_TOP_DIR)/$(EOS_DIR) \
                    $(NETWORK_TOP_DIR)/$(NETWORK_DIR)

  ifeq ($(findstring helmeos, $(EOS_DIR)), helmeos)
    Fmincludes += Microphysics/EOS/helmeos
    EOS_PATH := $(MAESTRO_TOP_DIR)/Microphysics/EOS/$(strip $(EOS_DIR))
    ##ALL: table
	 #programs += table
  endif

  f90sources += probin.f90

  include $(NETWORK_TOP_DIR)/$(strip $(NETWORK_DIR))/NETWORK_REQUIRES

  ifdef NEED_VODE
    Fmdirs += Util/VODE
  endif

  ifdef NEED_BLAS
    Fmdirs += Util/BLAS
  endif
  
  ifdef NEED_LINPACK
    Fmdirs += Util/LINPACK
  endif
endif

#Maestro stuff
Fmpack := $(foreach dir, $(Fmdirs), $(MAESTRO_TOP_DIR)/$(dir)/GPackage.mak)
Fmlocs := $(foreach dir, $(Fmdirs), $(MAESTRO_TOP_DIR)/$(dir))

#BoxLib stuff
Fmpack += $(foreach dir, $(BOXLIB_CORE), $(BOXLIB_HOME)/$(dir)/GPackage.mak)
Fmlocs += $(foreach dir, $(BOXLIB_CORE), $(BOXLIB_HOME)/$(dir))

#Microphysics
Fmpack += $(foreach dir, $(MICROPHYS_CORE), $(dir)/GPackage.mak)
Fmlocs += $(foreach dir, $(MICROPHYS_CORE), $(dir))

# any include files
Fmincs := $(foreach dir, $(Fmincludes), $(MAESTRO_TOP_DIR)/$(dir))

# include the necessary GPackage.mak files that define this setup
include $(Fmpack)

#Add module
f90sources += fsubchandra_mod.f90

# vpath defines the directories to search for the source files
VPATH_LOCATIONS += $(Fmlocs)

# list of directories to put in the Fortran include path
FINCLUDE_LOCATIONS += $(Fmincs)

programs += fsubchandra

all: $(pnames) table

# probin stuff
PROBIN_PARAMETER_DIRS :=
EXTERN_PARAMETER_DIRS += $(MICROPHYS_CORE)
EXTERN_PARAMETERS := $(shell $(BOXLIB_HOME)/Tools/F_scripts/findparams.py $(EXTERN_PARAMETER_DIRS))

PROBIN_TEMPLATE := $(MAESTRO_TOP_DIR)/Util/parameters/dummy.probin.template

probin.f90: $(PROBIN_PARAMETERS) $(EXTERN_PARAMETERS) $(PROBIN_TEMPLATE)
	@echo " "
	$(BOXLIB_HOME)/Tools/F_scripts/write_probin.py \
           -t $(PROBIN_TEMPLATE) -o probin.f90 -n probin \
           --pa "$(PROBIN_PARAMETERS)" --pb "$(EXTERN_PARAMETERS)"

include $(BOXLIB_HOME)/Tools/F_mk/GMakerules.mak

%.$(suf).exe:%.f90 $(objects)
ifdef MKVERBOSE
	$(LINK.f90) -o $@ $< $(objects) $(libraries)
else	
	@echo "Linking $@ ... "
	@$(LINK.f90) -o $@ $< $(objects) $(libraries)
endif

table:
	@if [ ! -f helm_table.dat ]; then echo ${bold}Linking helm_table.dat${normal}; ln -s $(EOS_PATH)/helm_table.dat .;  fi