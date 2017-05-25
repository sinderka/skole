
IDIR= -I ../../../lapack/LAPACKE/include/ -I ../../../lapack/CBLAS/include/ 
CXX=g++
CPPFLAGS=-std=c++14 $(IDIR) -llapack -lcblas

ODIR=.
LDIR= -L ../../../lapack/CBLAS/src/ -L ../../../lapack/LAPACKE/src/ -L ../../../lapack/SRC/
LDLIBS= ${LDIR}

_DEPS= io.h krylov.h random_Arrayi.h tools.h test.h tools_Test.cpp krylov_Test.cpp
DEPS=$(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ=io.o krylov.o random_Array.o tools.o tools_Test.o krylov_Test.o
OBJ=$(patsubst %,$(ODIR)/%,$(_OBJ))

$(ODIR)/%.o: %.c $(DEPS)
	$(CXX) -c -o $@ $< $(CPPFLAGS)

all: main test

main: main.cpp $(OBJ)
	g++ -o $@ $^ ${LDLIBS} $(CPPFLAGS)

test: test.cpp $(OBJ)
	g++ -o $@ $^ ${LDLIBS} $(CPPFLAGS) 
	

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~


