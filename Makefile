CXX=g++
RM=rm -rfv
#OPTLVL=-g
OPTLVL=-O2 -march=native

#get rootflags from root-config --cflags

ROOTFLAGS=-pthread -I$(ROOTSYS)/include

CPPFLAGS=$(OPTLVL) -Wall -Wextra -pedantic-errors -fmax-errors=8 $(ROOTFLAGS)
LDFLAGS=$(OPTLVL) $(shell root-config --ldflags)
LDLIBS=$(shell root-config --libs) -lgsl -lgslcblas

SRCS=met.cc

OBJS=$(subst .cc,.o,$(SRCS))


ws: met.o
	$(CXX) $(LDFLAGS) -o met met.o $(LDLIBS)

mrproper: 
	$(RM) *.pcm *.d *.so *.o met

