
IDIR= -I ../../../lapack-3.7.0/include/
CXX=g++
CPPFLAGS=-std=c++14 $(IDIR) -llapacke -llapack -lgfortran -lblas -ltmglib

ODIR=.
LDIR= -L ../../../lapack-3.7.0/
LDLIBS= ${LDIR}

_DEPS= krylov.h tools.h test.h problems.h tools_Test.cpp krylov_Test.cpp 
DEPS=$(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ=krylov.o tools.o problems.o
OBJ=$(patsubst %,$(ODIR)/%,$(_OBJ))

_OBJ_T=tools_Test.o krylov_Test.o
OBJ_T=$(patsubst %,$(ODIR)/%,$(_OBJ_T))

$(ODIR)/%.o: %.c $(DEPS)
	$(CXX) -c -o $@ $< $(CPPFLAGS)

all: main test

main: main.cpp $(OBJ)
	g++ -o $@ $^ ${LDLIBS} $(CPPFLAGS)

test: test.cpp $(OBJ) $(OBJ_T)
	g++ -o $@ $^ ${LDLIBS} $(CPPFLAGS) 
	

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~


