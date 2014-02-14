## This makefile must be executed with gmake (gnu make).

#set to g77 or absoft to choose a compiler under linux
LINUX_COMPILER=g77
#LINUX_COMPILER=absoft

simcdir = .
#CERN_ROOT = /site/cernlib/pc_linux/99
#CERN_ROOT = /cern/pro
#export ABSOFT=/apps/absoft/PRO/usr/absoft
#CERN_ROOT=/u/site/cernlib/x86_64_rhel5/2005

RM        = rm -f 
SHELL     = /bin/sh

SPEC      = ./hms/

SHARE	  = ./shared/

my_objs = $(SPEC)mc_hms.o     $(SPEC)mc_hms_hut.o  $(SPEC)mc_hms_recon.o \
	   project.o           rotate_haxis.o      check_dipole.o        \
	   transp.o            musc.o              musc_ext.o            \
	   stringlib.o         ranecu.o            loren.o               \
	   locforunt.o         gauss1.o            mt19937.o             \
	   trg_track.o         hms_track.o           \
	   qfs_new13_sub.o      init.o                \
	   brem.o              radc.o              hcf2r.o	\
           hallc2h.o enerloss_new.o engine_eloss.o christy_rss.o resmod.o ressf.o \
           qfs_deut.o call_pb_ext_subroutine.o externals_subroutine.o \
           fitemc_2006.o  quasiy8_2006.o  F1F2IN06.o F1F209_rss.o cer_effcorr.o dc_effcorr.o \

my_dobjs = $(SPEC)mc_hms.do     $(SPEC)mc_hms_hut.do  $(SPEC)mc_hms_recon.do \
	   project.do           rotate_haxis.do      check_dipole.do        \
	   transp.do            musc.do              musc_ext.do            \
	   stringlib.do         ranecu.do            loren.do               \
	   locforunt.do         gauss1.do            mt19937.do             \
	   trg_track.do         hms_track.do           \
	   qfs_new13_sub.do      init.do                \
	   brem.do              radc.do              hcf2r.do	\
           hallc2h.do enerloss_new.do engine_eloss.do christy_rss.do  resmod.do ressf.do \
           qfs_deut.do call_pb_ext_subroutine.do externals_subroutine.do \


MYOS := $(subst -,,$(shell uname))
CERNLIBS = -lgeant$(GEANTVER) -lpawlib -lgraflib -lgrafX11 -lpacklib -lmathlib

ifeq ($(MYOS),HPUX)
  ifneq (,$(findstring 09,$(shell uname -r)))
    HPUXVERSION := 09
  else
    HPUXVERSION := 10
  endif
  LIBROOT = $(Csoft)/../$(MYOS)$(HPUXVERSION)/lib
else
  LIBROOT = $(Csoft)/../$(MYOS)/lib
endif

ifeq ($(MYOS),HPUX)
  CERN_ROOT = /site/cernlib/hp700_ux90/96a
  FFLAGS=+U77 +ppu -C +e +es +FPVZOU -O +Onolimit -R8
  LDFLAGS=-Wl,-a archive
  OTHERLIBS = \
	-Wl,-L$(CERN_ROOT)/lib -lpacklib $(CERNLIBS) \
	-Wl,-L/usr/lib/X11R5 -lX11 -lm
endif


ifeq ($(MYOS),ULTRIX)
  FFLAGS=-check_bounds
  LDFLAGS=
  OTHERLIBS = -L$(CERN_ROOT)/lib -lpacklib $(CERNLIBS)
endif

ifeq ($(MYOS),SunOS)
  CERN_ROOT = /site/cernlib/sun4_solaris2/97a
  FFLAGS=-g -e  -I$(Csoft)/SRC/INCLUDE 
  ifeq ($(OSTYPE),SunOS4)
    OTHERLIBS = -L$(CERN_ROOT)/lib $(CERNLIBS) -lnsl -lX11
  else
    OTHERLIBS = -L$(CERN_ROOT)/lib $(CERNLIBS) -lnsl -lsocket -lX11
  endif
endif

ifeq ($(MYOS),AIX)
  F77=f77
  FFLAGS=-g -qfixed=132 -qextname -O -I$(Csoft)/SRC/INCLUDE
  OTHERLIBS = -L$(CERN_ROOT)/lib -lpacklib $(CERNLIBS) -lX11
endif

ifeq ($(MYOS),OSF1)
  F77=f77
  CERN_ROOT = /disk1/lib/cern/new
  LIBROOT = $(Csoft)/OSF1/lib
  FFLAGS= -r8 -extend_source -Wl,-taso -I -warn argument_checking \
        -warn declarations -warn truncated_source -warn unused
  LDFLAGS= 
  OTHERLIBS = -Wl,-L$(CERN_ROOT)/lib \
        -lpacklib $(CERNLIBS) -Wl,-L/usr/lib/X11R5 -lX11 -lm 
endif

ifeq ($(MYOS),Linux)
  ifeq ($(LINUX_COMPILER),absoft)
    FABSFLAGS=-V -W -f -s -N1 -B108 -B100 -N90 -N22 -N2 -N113
    INCLUDES=-I$(Csoft)/SRC/INCLUDE
    EXTRAFLAGS=-DABSOFTFORTRAN
    FFLAGS=-O $(INCLUDES) $(FABSFLAGS) $(EXTRAFLAGS)
    FFLAG1=$(FFLAGS) -c
    OTHERLIBS = -L$(CERN_ROOT)/lib $(CERNLIBS) -lV77 -lU77 -lg2c -lc -lm -lnsl
    FC  := $(ABSOFT)/bin/f77
    F77 := $(ABSOFT)/bin/f77
  endif
  ifeq ($(LINUX_COMPILER),g77)
    FC=gfortran
    F77=gfortran
    OTHERLIBS=-L$(CERN_ROOT)/lib -lpacklib -lmathlib -lnss_nis -lnsl
    OTHERLIBS=-L/usr/lib64/cernlib/2006/lib -Wl,-static -lmathlib -lpacklib -lkernlib -Wl,-dy -llapack -lm -lnsl -lcrypt -ldl
##  FFLAGSO=-O1 -malign-double -fno-automatic -I. -ffixed-line-length-none -fno-second-underscore
##  FFLAGSO=-O1 -malign-double -fno-automatic -I. -ffixed-line-length-none -fno-second-underscore -Wall -Wsurprising
##    FFLAGSO=-O1 -malign-double -fno-automatic -I. -ffixed-line-length-none -fno-second-underscore -Wimplicit -Wall -Wsurprising
##    FFLAGSD=-g -malign-double -fno-automatic -I. -ffixed-line-length-none -fno-second-underscore
    FFLAGSO=-O1  -I. -ffixed-line-length-none   -Wall -Wsurprising
    FFLAGSD=-g  -fno-automatic -I. -ffixed-line-length-none 
  endif
endif

%.o: %.f
	$(F77) $(FFLAGSO) -c $< -o $@
%.do: %.f
	$(F77) $(FFLAGSD) -c $< -o $*.do

mc_hms_single: $(my_objs) Makefile mc_hms_single.o
	$(F77) -o $@ $(FFLAGSO) mc_hms_single.o $(my_objs) $(OTHERLIBS)

debug: $(my_dobjs) Makefile mc_hms_single.do
	$(F77) -o mc_hms_single_dbg $(FFLAGSD) mc_hms_single.do $(my_dobjs) $(OTHERLIBS)

clean:
	find . -name '*.o' -exec rm {} \;
	find . -name '*.do' -exec rm {} \;
	rm -f mc_hms_single mc_hms_single_dbg

#the rule below updates a file, Makefile.dep, that lists which include files
#each fortran file uses.  If an include file is changed, all the fortran files
#that use it will automatically be recompiled.   p.mckee dec03
Makefile.dep : *.inc gen_constants.par
	@ echo "Updating dependencies on include files"
	@ echo -e "#specify include files used by each .f file\n" > Makefile.dep
	@ for file in `ls *.f | sed 's/\.f//'` ; do\
            inc_files=`cat $$file.f | sed -n 's/^[[:blank:]]*[Ii][Nn][Cc][Ll][Uu][Dd][Ee][ ]*['\''"]\(.*\)['\''"]/\1/p' | sort | uniq | xargs echo -n` ;\
            if [ -n "$$inc_files" ]; then \
              echo "$$file.o : $$file.f $$inc_files" >> Makefile.dep;\
              echo "$$file.do : $$file.f $$inc_files" >> Makefile.dep;\
            fi;\
          done

include Makefile.dep
