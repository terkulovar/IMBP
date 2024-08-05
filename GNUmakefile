# $Id: GNUmakefile,v 1.2 2003/01/23 15:31:39 maire Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------
#CXX := g++
name := crystal
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = .
endif



.PHONY: all
all: lib bin

include $(G4INSTALL)/config/binmake.gmk
#include /home/terkulov/geant4/counter/config/binmake.gmk

#ifdef ROOTSYS
#  CXXFLAGS += -I$(shell root-config --incdir)
#  CXXFLAGS += -std=c++17
  CPPFLAGS += -I$(shell /home/terkulov/FairSoft/install/bin/root-config --incdir)
#  CPPFLAGS += -std=c++17
#  CPPFLAGS += -I/home/terkulov/litrani/Litraniff
#  CPPFLAGS += -I/home/terkulov/geant4/CRY/cry_v1.7/src
  LDLIBS   += $(shell root-config --libs)
#  LDLIBS   += $(shell /home/terkulov/FairSoft/install/bin/root-config --libs)
#  LDLIBS   += -L/home/terkulov/litrani/Litraniff -lLitrani
#  LDLIBS   += -L/home/terkulov/geant4/CRY/cry_v1.7/lib -lCRY
#endif
G4BINDIR := .

tar:
	tar zcvf crystal.tgz *.cc src/*.cc include/*.hh GNUmakefile *.mac *.in visTutor/* data/* *.loop

visclean:
	rm -f g4*.prim g4*.eps g4*.wrl
	rm -f .DAWN_*
