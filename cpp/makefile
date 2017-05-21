
IDIR= cblas/include/
CXX=g++
CPPFLAGS=-std=c++14 -I $(IDIR)

ODIR=.
LDIR= cblas/lib/
LDLIBS=-L ${LDIR} -lopenblas

_DEPS= io.h krylov.h random_Arrayi.h tools.h test.h tools_Test.cpp krylov_Test.cpp
DEPS=$(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ=io.o krylov.o random_Array.o tools.o tools_Test.o krylov_Test.o
OBJ=$(patsubst %,$(ODIR)/%,$(_OBJ))

$(ODIR)/%.o: %.c $(DEPS)
	$(CXX) -c -o $@ $< $(CPPFLAGS)

all: main test

main: main.cpp $(OBJ)
	g++ -o $@ $^ $(CPPFLAGS) ${LDLIBS}

test: test.cpp $(OBJ)
	g++ -o $@ $^ $(CPPFLAGS) ${LDLIBS}
	

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~




#$(program_NAME): $(program_OBJS)
#    $(LINK.cc) $(program_OBJS) -o $(program_NAME)
