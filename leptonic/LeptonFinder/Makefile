ARCH          = linux

CXX           =
ObjSuf        = o
SrcSuf        = cxx
ExeSuf        = exe
DllSuf        = so
OutPutOpt     = -o 

ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)
ROOTGLIBS    := $(shell root-config --glibs)

#FASTJETFLAGS := $(shell /afs/cern.ch/user/s/sredford/Software/fastjet-install/bin/fastjet-config --cxxflags)
#FASTJETLIBS  := $(shell /afs/cern.ch/user/s/sredford/Software/fastjet-install/bin/fastjet-config --libs)

ifeq ($(ARCH),linux)
# Linux with egcs, gcc 2.9x, gcc 3.x (>= RedHat 5.2)
CXX           = g++
CXXFLAGS      = -Wall -fPIC -ggdb
LD            = g++ 
LDFLAGS       = -O
SOFLAGS       = -shared
endif



ifeq ($(CXX),)
x(error $(ARCH) invalid architecture)
endif

CXXFLAGS     += $(ROOTCFLAGS) -I$(LCIO)/include
LIBS          = $(ROOTLIBS) $(SYSLIBS) -L$(LCIO)/lib -llcio -lsio
GLIBS         = $(ROOTGLIBS) $(SYSLIBS)

CXXFLAGS     += $(FASTJETFLAGS)
LIBS         += $(FASTJETLIBS)

#------------------------------------------------------------------------------

MAINEVENTO    = tth_analysis.$(ObjSuf) main.$(ObjSuf) 
MAINEVENTS    = tth_analysis.$(SrcSuf) main.$(SrcSuf) 
MAINEXE       = tth.$(ExeSuf)

OBJS          = $(MAINEVENTO)

PROGRAMS      = $(MAINEXE)

.SUFFIXES: .$(SrcSuf) .$(ObjSuf) .$(DllSuf)
.PHONY:   

all:            $(PROGRAMS)

$(MAINEXE):      $(MAINEVENTO)
		$(LD) $(LDFLAGS) $^ $(LIBS) $(OutPutOpt)$@
		@echo "$@ done"



clean:
		@rm -f $(OBJS) core

distclean:      clean
		@rm -f $(PROGRAMS) *Dict.* *.def *.exp \
		   *.root *.ps *.so .def so_locations
		@rm -rf cxx_repository

.SUFFIXES: .$(SrcSuf)

###
.$(SrcSuf).$(ObjSuf):
	$(CXX) $(CXXFLAGS) -c $<